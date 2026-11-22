#' Initialize a UBiQTree Epistemic SHAP Explainer
#'
#' @param model A trained tree ensemble (`ranger`, `randomForest`, or `xgb.Booster`).
#' @param X_train Training feature `data.frame` or `matrix`.
#' @param y_train Training target vector (factor/character for classification, numeric for regression).
#' @param beta Temperature parameter for softmax tree weighting (default: 5.0).
#' @param seed Optional integer random seed.
#'
#' @return An object of class `ubiq_explainer`.
#' @export
ubiq_explainer <- function(model, X_train, y_train, beta = 5.0, seed = NULL) {
  X_df <- as.data.frame(X_train)
  is_classif <- is.factor(y_train) || is.character(y_train)
  task <- if (is_classif) "classification" else "regression"
  classes <- if (is_classif) levels(as.factor(y_train)) else NULL
  
  unified <- treeshap::unify(model, X_df)
  
  # Handle multi-output unified object (multi-class classification)
  base_tree_df <- if (inherits(unified, "model_unified_multioutput")) {
    unified[[1]]$model
  } else {
    unified$model
  }
  tree_ids <- sort(unique(base_tree_df$Tree))
  n_trees <- length(tree_ids)
  
  weights <- compute_tree_weights(model, n_trees, X_df, y_train, task, beta)
  
  structure(
    list(
      model = model,
      unified = unified,
      task = task,
      classes = classes,
      tree_ids = tree_ids,
      n_trees = n_trees,
      feature_names = colnames(X_df),
      weights = weights,
      beta = beta,
      seed = seed
    ),
    class = "ubiq_explainer"
  )
}

#' @keywords internal
compute_tree_weights <- function(model, n_trees, X_df, y_train, task, beta) {
  if (n_trees <= 1) return(1.0)
  
  if (inherits(model, "ranger")) {
    preds <- stats::predict(model, data = X_df, predict.all = TRUE)$predictions
  } else if (inherits(model, "randomForest")) {
    preds <- stats::predict(model, newdata = X_df, predict.all = TRUE)$individual
  } else {
    return(rep(1 / n_trees, n_trees))
  }
  
  scores <- vapply(seq_len(n_trees), function(t_idx) {
    if (task == "classification") {
      if (length(dim(preds)) == 3L) {
        # probability = TRUE ranger model: (n_obs, n_classes, n_trees)
        prob_mat <- preds[, , t_idx, drop = FALSE]
        pred_idx <- max.col(matrix(prob_mat, nrow = nrow(X_df)))
        tree_pred <- dimnames(preds)[[2]][pred_idx]
      } else {
        tree_pred <- preds[, t_idx]
        if (is.numeric(tree_pred) && !is.null(model$forest$levels)) {
          tree_pred <- model$forest$levels[tree_pred]
        }
      }
      mean(as.character(tree_pred) == as.character(y_train))
    } else {
      tree_pred <- preds[, t_idx]
      ss_res <- sum((y_train - as.numeric(tree_pred))^2)
      ss_tot <- sum((y_train - mean(y_train))^2)
      if (ss_tot == 0) 0 else 1 - (ss_res / ss_tot)
    }
  }, numeric(1))
  
  max_score <- max(scores)
  exp_vals <- exp(beta * (scores - max_score))
  exp_vals / sum(exp_vals)
}

#' Compute SHAP Values with Epistemic Uncertainty
#'
#' @param explainer An object of class `ubiq_explainer`.
#' @param x A single-row `data.frame` representing the instance to explain.
#' @param n_samples Number of Dirichlet hypothesis samples (default: 500).
#' @param alpha Dirichlet concentration parameter (default: 1.0).
#' @param class_idx Integer index of the class to explain for multi-class models (default: 1).
#'
#' @return An object of class `ubiq_explanation`.
#' @export
explain_uncertainty <- function(explainer, x, n_samples = 500, alpha = 1.0, class_idx = 1L) {
  stopifnot(inherits(explainer, "ubiq_explainer"))
  if (!is.null(explainer$seed)) set.seed(explainer$seed)
  
  x_df <- as.data.frame(x)[1, , drop = FALSE]
  n_trees <- explainer$n_trees
  n_features <- length(explainer$feature_names)
  
  # Select class-specific unified model if multi-class
  target_unified <- if (inherits(explainer$unified, "model_unified_multioutput")) {
    explainer$unified[[class_idx]]
  } else {
    explainer$unified
  }
  class_name <- if (!is.null(explainer$classes)) explainer$classes[class_idx] else NULL
  
  # Precompute per-tree SHAP matrix (n_trees x n_features) once
  tree_shaps <- vapply(explainer$tree_ids, function(t_id) {
    sub_unified <- target_unified
    sub_df <- sub_unified$model[sub_unified$model$Tree == t_id, , drop = FALSE]
    # Reset Tree ID to 0L so treeshap's C++ backend treats it as a single-tree model
    sub_df$Tree <- 0L
    sub_unified$model <- sub_df
    res <- treeshap::treeshap(sub_unified, x_df, verbose = FALSE)
    as.numeric(res$shaps[1, ])
  }, numeric(n_features))
  
  tree_shaps <- t(tree_shaps)
  
  # Sample Dirichlet hypothesis space
  phi_dist <- matrix(0, nrow = n_samples, ncol = n_features)
  colnames(phi_dist) <- explainer$feature_names
  
  for (s in seq_len(n_samples)) {
    gamma_draws <- stats::rgamma(n_trees, shape = pmax(alpha * explainer$weights, 1e-12), rate = 1)
    dirichlet_weights <- gamma_draws / sum(gamma_draws)
    tree_counts <- as.numeric(stats::rmultinom(1, size = n_trees, prob = dirichlet_weights))
    phi_dist[s, ] <- (tree_counts %*% tree_shaps) / n_trees
  }
  
  mean_shap <- colMeans(phi_dist)
  std_shap <- apply(phi_dist, 2, stats::sd)
  ci_95 <- apply(phi_dist, 2, stats::quantile, probs = c(0.025, 0.975))
  entropy_vals <- compute_entropy(phi_dist)
  sign_stab <- compute_sign_stability(phi_dist, mean_shap)
  
  summary_df <- data.frame(
    feature = explainer$feature_names,
    mean = mean_shap,
    std = std_shap,
    ci_low = ci_95[1, ],
    ci_high = ci_95[2, ],
    entropy = entropy_vals,
    sign_stability = sign_stab,
    row.names = NULL
  )
  
  structure(
    list(
      summary = summary_df,
      samples = phi_dist,
      task = explainer$task,
      class_name = class_name
    ),
    class = "ubiq_explanation"
  )
}

#' @keywords internal
compute_entropy <- function(phi_dist, bins = 10) {
  apply(phi_dist, 2, function(col_vals) {
    if (diff(range(col_vals)) == 0) return(0)
    h <- graphics::hist(col_vals, breaks = bins, plot = FALSE)
    probs <- h$density * diff(h$breaks)
    probs <- probs[probs > 0]
    -sum(probs * log(probs))
  })
}

#' @keywords internal
compute_sign_stability <- function(phi_dist, mean_shap) {
  mean_sign <- sign(mean_shap)
  vapply(seq_along(mean_shap), function(j) {
    mean(sign(phi_dist[, j]) == mean_sign[j])
  }, numeric(1))
}