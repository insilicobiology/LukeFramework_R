library(LukeFrameworkR)
library(ranger)
library(ggplot2)
library(patchwork)

set.seed(42)

# ------------------------------------------------------------------------------
# 1. Simulate a Clinical Hospital Length-of-Stay (LOS) Cohort (cf. MIMIC-III)
#    Classes: Short_Stay (<3 days), Medium_Stay (3-7 days), Long_Stay (>=7 days)
#    Engineered so top drivers (NumTransfers, NumCallouts, NumNotes) exhibit
#    high epistemic ambiguity for the borderline Medium_Stay cohort.
# ------------------------------------------------------------------------------
n <- 600
df <- data.frame(
  NumTransfers     = rpois(n, lambda = 2.2),
  NumNotes         = round(rnorm(n, mean = 18, sd = 7)),
  NumCallouts      = rpois(n, lambda = 1.5),
  NumDiagnosis     = rpois(n, lambda = 6) + 1,
  NumProcs         = rpois(n, lambda = 3),
  Age              = round(rnorm(n, mean = 64, sd = 14)),
  NumLabs          = round(rnorm(n, mean = 45, sd = 15)),
  AdmitProcedure   = rbinom(n, 1, 0.35),
  TotalNumInteract = round(runif(n, 5, 40)),
  NoiseBiomarker   = rnorm(n)
)

# Construct latent severity score with non-linear threshold overlap around Medium Stay
latent_los <- 1.4 * df$NumTransfers +
  0.15 * df$NumNotes +
  0.9 * df$NumCallouts +
  0.25 * df$NumDiagnosis +
  rnorm(n, sd = 2.2) # Aleatoric clinical noise inducing tree split disagreement

df$LOS_Group <- cut(
  latent_los,
  breaks = c(-Inf, 5.2, 8.8, Inf),
  labels = c("Short_Stay", "Medium_Stay", "Long_Stay")
)

# 70:30 Stratified Train / Test Split (as in Dubey et al., 2026)
train_idx <- unlist(lapply(split(seq_len(n), df$LOS_Group), function(idx) {
  sample(idx, size = floor(0.7 * length(idx)))
}))
train_df <- df[train_idx, ]
test_df  <- df[-train_idx, ]

# ------------------------------------------------------------------------------
# 2. Train Random Forest Classifier (n_trees = 100, probability = TRUE)
# ------------------------------------------------------------------------------
rf_los <- ranger(
  LOS_Group ~ .,
  data = train_df,
  num.trees = 100,
  probability = TRUE,
  seed = 42
)

# Evaluate test accuracy
pred_probs <- predict(rf_los, data = test_df[, -ncol(test_df)])$predictions
pred_class <- colnames(pred_probs)[max.col(pred_probs)]
cat(sprintf("Test Set Accuracy: %.1f%%\n\n", mean(pred_class == test_df$LOS_Group) * 100))

# ------------------------------------------------------------------------------
# 3. Initialize UBiQTree Explainer (beta = 5.0 softmax tree weighting)
# ------------------------------------------------------------------------------
explainer <- ubiq_explainer(
  model   = rf_los,
  X_train = train_df[, -ncol(train_df)],
  y_train = train_df$LOS_Group,
  beta    = 5.0,
  seed    = 42
)

# Select a clinically ambiguous test patient predicted as Medium_Stay (class_idx = 2)
ambiguous_idx <- which(test_df$LOS_Group == "Medium_Stay")[1]
patient_x <- test_df[ambiguous_idx, -ncol(test_df), drop = FALSE]

res_medium <- explain_uncertainty(
  explainer = explainer,
  x         = patient_x,
  n_samples = 500,
  alpha     = 0.5,   # Matches MIMIC-III Dirichlet concentration in paper
  class_idx = 2L     # 1 = Short_Stay, 2 = Medium_Stay, 3 = Long_Stay
)

# Print full uncertainty summary table
summary_ordered <- res_medium$summary[order(-abs(res_medium$summary$mean)), ]
print(summary_ordered, digits = 3, row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Generate Paper Visualizations
# ------------------------------------------------------------------------------
# Figure A: SHAP Summary with Epistemic Uncertainty Violins & 95% CI (cf. Fig 1)
p_summary <- plot_uncertainty_bars(
  res_medium,
  title = "SHAP Summary with Epistemic Uncertainty"
)
print(p_summary)

# Figure B: Top-3 Feature KDE Distributions & Sign Stability Grid (cf. Figs 2-4)
top3_feats <- summary_ordered$feature[1:3]
kde_plots <- lapply(top3_feats, function(f) {
  plot_uncertainty_distribution(res_medium, feature = f, title = "SHAP Distribution")
})
p_top3_grid <- wrap_plots(kde_plots, ncol = 3)
print(p_top3_grid)

# Figure C: 3-Panel Metric Comparison (Std Dev, Entropy, Sign Stability)
p_metrics <- plot_uncertainty_comparison(res_medium)
print(p_metrics)

# Figure D: Entanglement Scatter Plot (Mean |SHAP| vs Epistemic Std Dev, cf. Fig 5)
ent_df <- res_medium$summary
ent_df$abs_mean <- abs(ent_df$mean)
ent_df$is_top3 <- ifelse(ent_df$feature %in% top3_feats, "Top 3 Features", "Other Features")
rho <- cor(ent_df$abs_mean, ent_df$std, method = "pearson")

p_entanglement <- ggplot(ent_df, aes(x = abs_mean, y = std, color = is_top3)) +
  geom_point(size = 3.5) +
  geom_text(
    data = subset(ent_df, is_top3 == "Top 3 Features"),
    aes(label = feature),
    vjust = -0.8, hjust = 0.5, size = 3.8, show.legend = FALSE
  ) +
  annotate(
    "label", x = -Inf, y = Inf,
    label = sprintf("Pearson rho = %.2f", rho),
    hjust = -0.1, vjust = 1.3, size = 4
  ) +
  scale_color_manual(values = c("Top 3 Features" = "firebrick", "Other Features" = "steelblue")) +
  labs(
    title = "Entanglement: Correlation of Importance and Uncertainty (Medium Stay)",
    x = "Mean Absolute SHAP Value (Feature Importance)",
    y = "Std Dev of SHAP Value (Epistemic Uncertainty)",
    color = NULL
  ) +
  theme_minimal(base_size = 13)

print(p_entanglement)
