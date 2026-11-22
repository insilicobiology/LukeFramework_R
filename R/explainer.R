#' Initialize a UBiQTree Epistemic SHAP Explainer
#'
#' Prepares a trained tree ensemble (`ranger`, `randomForest`, or `xgb.Booster`)
#' to quantify epistemic uncertainty in SHAP feature attributions.
#'
#' @param model A trained tree ensemble (`ranger`, `randomForest`, or `xgb.Booster`).
#' @param X_train Training feature `data.frame` or `matrix` (predictors only, without the target column).
#' @param y_train Training target vector (factor/character for classification, numeric for regression).
#' @param beta Softmax temperature parameter for weighting trees by training accuracy (default: 5.0).
#' @param seed Optional integer random seed for reproducible results.
#'
#' @return An object of class `ubiq_explainer`.
#' @export
ubiq_explainer <- function(model, X_train, y_train, beta = 5.0, seed = NULL) {
  if (missing(model) || missing(X_train) || missing(y_train)) {
    stop("Please provide `model`, `X_train` (features only), and `y_train` (outcome vector).", call. = FALSE)
  }

  X_df <- as.data.frame(X_train)
  if (nrow(X_df) != length(y_train)) {
    stop(sprintf(
      "`X_train` has %d rows, but `y_train` has %d values. Please ensure they match.",
      nrow(X_df), length(y_train)
    ), call. = FALSE)
  }

  if (inherits(model, "ranger") && is.null(model$forest)) {
    stop("Your `ranger` model was trained without saving trees. Please re-train with `write.forest = TRUE` (the default).", call. = FALSE)
  }

  is_classif <- is.factor(y_train) || is.character(y_train)
  task <- if (is_classif) "classification" else "regression"
  classes <- if (is_classif) levels(as.factor(y_train)) else NULL

  if (inherits(model, "ranger") && is_classif && length(classes) > 2L) {
    unified <- unify_ranger_multiclass(model, X_df, classes)
  } else {
    unified <- treeshap::unify(model, X_df)
  }

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

#' @export
print.ubiq_explainer <- function(x, ...) {
  cat("UBiQTree Epistemic SHAP Explainer\n")
  cat("  Task:     ", x$task, "\n", sep = "")
  if (!is.null(x$classes)) {
    cat("  Classes:  ", paste(paste0(seq_along(x$classes), "=", x$classes), collapse = ", "), "\n", sep = "")
  }
  cat("  Trees:    ", x$n_trees, "\n", sep = "")
  cat("  Features: ", length(x$feature_names), "\n", sep = "")
  invisible(x)
}

#' @keywords internal
unify_ranger_multiclass <- function(model, X_df, classes) {
  n_trees <- model$num.trees
  tree_infos <- lapply(seq_len(n_trees), function(t_idx) {
    ti <- ranger::treeInfo(model, tree = t_idx)
    ti$Tree <- t_idx - 1L
    ti$Node <- ti$nodeID
    ti$Feature <- ifelse(ti$terminal, NA_character_, as.character(ti$splitvarName))
    ti$Split <- ti$splitval
    ti$Yes <- ti$leftChild
    ti$No <- ti$rightChild
    ti$Cover <- 1.0
    ti
  })

  unified_list <- lapply(seq_along(classes), function(c_idx) {
    pred_col <- paste0("pred.", classes[c_idx])
    if (!pred_col %in% names(tree_infos[[1]])) {
      pred_col <- paste0("pred.", c_idx)
    }
    models_c <- lapply(tree_infos, function(ti) {
      ti$Prediction <- if (pred_col %in% names(ti)) {
        ti[[pred_col]]
      } else {
        ifelse(ti$terminal, as.numeric(as.character(ti$prediction) == classes[c_idx]), NA_real_)
      }
      ti[, c("Tree", "Node", "Feature", "Split", "Yes", "No", "Prediction", "Cover")]
    })
    list(model = do.call(rbind, models_c), data = X_df)
  })

  class(unified_list) <- "model_unified_multioutput"
  unified_list
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

#' Exact Pure-R TreeSHAP for a Single Unified Tree Data Frame
#' @keywords internal
single_tree_shap <- function(tree_df, x_row, feature_names) {
  n_features <- length(feature_names)
  phi <- numeric(n_features)
  feat_map <- stats::setNames(seq_along(feature_names), feature_names)

  get_row <- function(node_id) {
    match(node_id, tree_df$Node)
  }

  extend <- function(m, p_z, p_o) {
    l <- length(m)
    m <- c(m, 0)
    for (idx in seq(l + 1L, 1L, by = -1L)) {
      val <- if (idx > 1L) m[idx - 1L] * p_z * (idx - 1L) / l else 0
      if (idx <= l) val <- val + m[idx] * p_o * (l - idx + 1L) / l
      m[idx] <- val
    }
    m
  }

  walk <- function(r, path_feats, path_pz, path_po) {
    if (is.na(r)) return(invisible(NULL))
    feat <- as.character(tree_df$Feature[r])

    if (is.na(feat) || feat == "Leaf") {
      leaf_val <- tree_df$Prediction[r]
      if (is.na(leaf_val) || length(path_feats) == 0L) return(invisible(NULL))

      uniq_f <- unique(path_feats)
      for (uf in uniq_f) {
        pos <- which(path_feats == uf)
        pz_f <- prod(path_pz[pos])
        po_f <- prod(path_po[pos])

        other_f <- setdiff(uniq_f, uf)
        m <- c(1)
        for (of in other_f) {
          opos <- which(path_feats == of)
          m <- extend(m, prod(path_pz[opos]), prod(path_po[opos]))
        }
        M <- length(m)
        shap_w <- vapply(seq_len(M) - 1L, function(s) {
          1 / (M * choose(M - 1L, s))
        }, numeric(1))
        phi[uf] <<- phi[uf] + leaf_val * (pz_f - po_f) * sum(m * shap_w)
      }
      return(invisible(NULL))
    }

    f_idx <- feat_map[[feat]]
    split_val <- tree_df$Split[r]
    yes_r <- get_row(tree_df$Yes[r])
    no_r  <- get_row(tree_df$No[r])

    cov_p <- tree_df$Cover[r]
    cov_y <- if (!is.na(yes_r)) tree_df$Cover[yes_r] else NA_real_
    cov_n <- if (!is.na(no_r))  tree_df$Cover[no_r]  else NA_real_

    r_y <- if (!is.na(cov_p) && cov_p > 0 && !is.na(cov_y)) cov_y / cov_p else 0.5
    r_n <- if (!is.na(cov_p) && cov_p > 0 && !is.na(cov_n)) cov_n / cov_p else 0.5

    x_val <- as.numeric(x_row[[feat]])
    go_yes <- !is.na(x_val) && x_val <= split_val

    walk(yes_r, c(path_feats, f_idx), c(path_pz, if (go_yes) 1 else 0), c(path_po, r_y))
    walk(no_r,  c(path_feats, f_idx), c(path_pz, if (!go_yes) 1 else 0), c(path_po, r_n))
  }

  root_row <- get_row(min(tree_df$Node))
  walk(root_row, integer(0), numeric(0), numeric(0))
  phi
}

#' Compute SHAP Values with Epistemic Uncertainty for One Sample
#'
#' @param explainer An object created by [ubiq_explainer()].
#' @param x A single-row `data.frame` representing one patient or biological sample.
#' @param n_samples Number of Dirichlet hypothesis samples (default: 500).
#' @param alpha Dirichlet concentration parameter (default: 1.0; smaller values focus on top-performing trees).
#' @param class_idx Integer index or character name of the target class to explain (default: 1).
#'
#' @return An object of class `ubiq_explanation`.
#' @export
explain_uncertainty <- function(explainer, x, n_samples = 500, alpha = 1.0, class_idx = 1L) {
  if (!inherits(explainer, "ubiq_explainer")) {
    stop("`explainer` must be created using `ubiq_explainer()`.", call. = FALSE)
  }
  if (!is.null(explainer$seed)) set.seed(explainer$seed)

  x_df <- as.data.frame(x)
  if (nrow(x_df) > 1L) {
    warning("`x` has multiple rows; explaining the first row only.", call. = FALSE)
  }
  x_df <- x_df[1, , drop = FALSE]

  missing_cols <- setdiff(explainer$feature_names, colnames(x_df))
  if (length(missing_cols) > 0L) {
    stop("Missing feature columns in `x`: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  x_df <- x_df[, explainer$feature_names, drop = FALSE]

  # Allow character class name in addition to integer index
  if (is.character(class_idx)) {
    if (is.null(explainer$classes) || !class_idx %in% explainer$classes) {
      stop("Unknown class name. Available classes: ", paste(explainer$classes, collapse = ", "), call. = FALSE)
    }
    class_idx <- match(class_idx, explainer$classes)
  }

  n_trees <- explainer$n_trees
  n_features <- length(explainer$feature_names)

  target_unified <- if (inherits(explainer$unified, "model_unified_multioutput")) {
    explainer$unified[[class_idx]]
  } else {
    explainer$unified
  }
  class_name <- if (!is.null(explainer$classes)) explainer$classes[class_idx] else NULL

  tree_shaps <- vapply(explainer$tree_ids, function(t_id) {
    sub_df <- target_unified$model[target_unified$model$Tree == t_id, , drop = FALSE]
    single_tree_shap(sub_df, x_df, explainer$feature_names)
  }, numeric(n_features))

  tree_shaps <- t(tree_shaps)

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
  ci_95 <- apply(phi_dist, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  entropy_vals <- compute_entropy(phi_dist)
  sign_stab <- compute_sign_stability(phi_dist, mean_shap)

  status <- ifelse(
    sign_stab < 0.70, "Unstable Direction",
    ifelse(sign_stab < 0.90 | std_shap > 0.15, "Review Recommended", "Reliable")
  )

  summary_df <- data.frame(
    feature = explainer$feature_names,
    mean = mean_shap,
    std = std_shap,
    ci_low = ci_95[1, ],
    ci_high = ci_95[2, ],
    entropy = entropy_vals,
    sign_stability = sign_stab,
    status = status,
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

#' @export
print.ubiq_explanation <- function(x, top_n = 10L, ...) {
  cat("UBiQTree Epistemic SHAP Explanation")
  if (!is.null(x$class_name)) cat(" (Class: ", x$class_name, ")", sep = "")
  cat("\n\n")

  ord <- order(-abs(x$summary$mean))
  df <- utils::head(x$summary[ord, , drop = FALSE], top_n)
  df$mean <- round(df$mean, 3)
  df$std <- round(df$std, 3)
  df$ci_low <- round(df$ci_low, 3)
  df$ci_high <- round(df$ci_high, 3)
  df$entropy <- round(df$entropy, 3)
  df$sign_stability <- sprintf("%.1f%%", df$sign_stability * 100)

  print(df, row.names = FALSE)
  invisible(x)
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
