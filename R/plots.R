#' Plot SHAP Values with Epistemic Uncertainty Bars & Violins
#'
#' @param explanation An object of class `ubiq_explanation`.
#' @param title Optional plot title.
#' @return A `ggplot` object.
#' @export
plot_uncertainty_bars <- function(explanation, title = "SHAP with Epistemic Uncertainty") {
  if (!is.null(explanation$class_name)) {
    title <- paste0(title, " - Class: ", explanation$class_name)
  }

  df <- explanation$summary |>
    dplyr::mutate(
      abs_mean = abs(.data$mean),
      feature = stats::reorder(.data$feature, .data$mean)
    )

  samples_long <- as.data.frame(explanation$samples) |>
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "feature", values_to = "shap") |>
    dplyr::mutate(feature = factor(.data$feature, levels = levels(df$feature)))

  ggplot2::ggplot() +
    ggplot2::geom_violin(
      data = samples_long,
      ggplot2::aes(x = .data$shap, y = .data$feature),
      fill = "dodgerblue", alpha = 0.25, color = NA, scale = "width"
    ) +
    ggplot2::geom_col(
      data = df,
      ggplot2::aes(x = .data$mean, y = .data$feature, fill = .data$abs_mean),
      alpha = 0.85, color = "black", width = 0.6
    ) +
    ggplot2::geom_errorbar(
      data = df,
      ggplot2::aes(xmin = .data$ci_low, xmax = .data$ci_high, y = .data$feature),
      orientation = "y",
      width = 0.25, linetype = "dashed", color = "black"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.6) +
    ggplot2::scale_fill_viridis_c(name = "Abs SHAP\nMagnitude") +
    ggplot2::labs(title = title, x = "SHAP Value (Impact on Prediction)", y = NULL) +
    ggplot2::theme_minimal(base_size = 13)
}

#' Plot Kernel Density Estimate for a Single Feature
#'
#' @param explanation An object of class `ubiq_explanation`.
#' @param feature Name or integer index of the feature to plot.
#' @param title Optional plot title.
#' @return A `ggplot` object.
#' @export
plot_uncertainty_distribution <- function(explanation, feature, title = "SHAP Value Distribution") {
  if (!is.null(explanation$class_name)) {
    title <- paste0(title, " - Class: ", explanation$class_name)
  }

  feat_name <- if (is.numeric(feature)) explanation$summary$feature[feature] else feature
  row <- explanation$summary[explanation$summary$feature == feat_name, ]
  vals <- data.frame(shap = explanation$samples[, feat_name])

  unc_label <- sprintf(
    "Epistemic Uncertainty: %s\n(Std = %.3f)",
    ifelse(row$std > 0.15, "High", "Low"), row$std
  )
  stab_label <- sprintf(
    "%s\n(Sign Stability = %.1f%%)",
    ifelse(row$sign_stability > 0.9, "Direction Certain",
           ifelse(row$sign_stability < 0.7, "Direction Uncertain", "Direction Mostly Stable")),
    row$sign_stability * 100
  )

  ggplot2::ggplot(vals, ggplot2::aes(x = .data$shap)) +
    ggplot2::annotate("rect", xmin = row$ci_low, xmax = row$ci_high,
                      ymin = -Inf, ymax = Inf, fill = "skyblue", alpha = 0.25) +
    ggplot2::geom_density(color = "navy", linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = row$mean, color = "firebrick", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
    ggplot2::annotate("label", x = -Inf, y = Inf, label = paste(unc_label, stab_label, sep = "\n\n"),
                      hjust = -0.05, vjust = 1.1, size = 3.8) +
    ggplot2::labs(
      title = paste0(title, "\n", feat_name),
      x = paste("SHAP Value for", feat_name),
      y = "Probability Density"
    ) +
    ggplot2::theme_minimal(base_size = 13)
}

#' Plot 3-Panel Uncertainty Metric Comparison (Std, Entropy, Sign Stability)
#'
#' @param explanation An object of class `ubiq_explanation`.
#' @param title Overall figure title.
#' @return A `patchwork` plot object.
#' @export
plot_uncertainty_comparison <- function(explanation, title = "Uncertainty Metric Comparison") {
  if (!is.null(explanation$class_name)) {
    title <- paste0(title, " - Class: ", explanation$class_name)
  }

  df <- explanation$summary |>
    dplyr::mutate(feature = stats::reorder(.data$feature, .data$mean))

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$std, y = .data$feature, fill = .data$feature)) +
    ggplot2::geom_col(show.legend = FALSE, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = mean(df$std), color = "red", linetype = "dashed") +
    ggplot2::scale_fill_viridis_d(option = "plasma") +
    ggplot2::labs(title = "Standard Deviation of SHAP Values", x = "Magnitude of Uncertainty", y = NULL) +
    ggplot2::theme_minimal()

  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$entropy, y = .data$feature, fill = .data$feature)) +
    ggplot2::geom_col(show.legend = FALSE, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = mean(df$entropy), color = "red", linetype = "dashed") +
    ggplot2::scale_fill_viridis_d(option = "plasma") +
    ggplot2::labs(title = "Explanation Entropy (Information Uncertainty)", x = "Entropy Value", y = NULL) +
    ggplot2::theme_minimal()

  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = .data$sign_stability, y = .data$feature, fill = .data$feature)) +
    ggplot2::geom_col(show.legend = FALSE, alpha = 0.85) +
    ggplot2::geom_vline(xintercept = c(0.9, 0.7, 0.4), color = c("green", "orange", "red"), linetype = "dashed") +
    ggplot2::scale_fill_viridis_d(option = "plasma") +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(title = "Sign Stability (Direction Consistency)", x = "Probability of Consistent Direction", y = NULL) +
    ggplot2::theme_minimal()

  patchwork::wrap_plots(p1, p2, p3, ncol = 1) +
    patchwork::plot_annotation(title = title)
}
