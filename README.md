# LukeFrameworkR (`Luke`)

### Trustworthy AI explanations for biologists

**Luke** (*Leverage UBiQTree, Keep Explanations*) is an R package that adds **error bars and reliability checks** to the explanations of tree-based machine learning models. It implements the **UBiQTree** framework ([Dubey et al., *Patterns*, 2026](https://www.cell.com/patterns/fulltext/S2666-3899(25)00302-2); original Python code: [dubeyakshat07/UBiQTree](https://github.com/dubeyakshat07/UBiQTree)).

It works with models from `ranger`, `randomForest`, and `xgboost`, and is designed to help life-science and clinical researchers separate **true biological drivers** from **random noise**.

---

## Contents

- [Why use LukeFrameworkR?](#why-use-lukeframeworkr)
- [Installation](#installation)
- [5-minute quick start](#5-minute-quick-start)
- [How to interpret the output](#how-to-interpret-the-output)
- [Key functions and parameters](#key-functions-and-parameters)
- [Best practices for biological and clinical studies](#best-practices-for-biological-and-clinical-studies)
- [Full walkthrough](#full-walkthrough)
- [Citation and license](#citation-and-license)

---

## Why use LukeFrameworkR?

Standard SHAP values give you **one number** per gene, biomarker, or clinical measurement, with no indication of how confident that number is. In small or noisy cohorts, a measurement that is pure noise can look just as important as a real driver.

> **The lab analogy.** Think of a Random Forest as a panel of 100 reviewers, each of whom saw a slightly different subset of your patients. Standard SHAP simply averages their votes. LukeFrameworkR repeatedly re-weights the reviewers according to how well they perform, then checks whether they still agree.

**What you get**, alongside the usual average SHAP value:

- a **95% credible interval** (error bar)
- the **spread (SD)** across model versions
- the **entropy** of the SHAP distribution
- **sign stability**: how often model versions agree on the direction of effect (up- vs. down-regulation)
- an automatic **reliability status**: `Reliable`, `Review Recommended`, or `Unstable Direction`

---

## Installation

Install the development version from GitHub using `remotes` (or `devtools`):

```r
# install.packages("remotes")
remotes::install_github("insilicobiology/LukeFramework_R")
```

**Dependencies.** `LukeFrameworkR` automatically installs `dplyr`, `ggplot2`, `patchwork`, `ranger`, `rlang`, `tidyr`, and `treeshap`. To explain `randomForest` or `xgboost` models, make sure those packages are installed too.

---

## 5-minute quick start

This minimal workflow trains a Random Forest (`ranger`) to predict a clinical or biological outcome, then explains a single sample. Replace `train_df` and `test_df` with your own data: one column per predictor (gene, biomarker, measurement) plus a factor column named `Outcome`.

### 1. Train your model

For classification, set `probability = TRUE` so the explainer can access class probabilities.

```r
library(LukeFrameworkR)
library(ranger)

set.seed(42)
rf_model <- ranger(
  Outcome ~ .,
  data        = train_df,
  num.trees   = 100,
  probability = TRUE,
  seed        = 42
)
```

### 2. Create the UBiQTree explainer

Provide the fitted model, the **training predictors** (`X_train`), and the **training outcome** (`y_train`). This step scores each tree by its predictive performance.

```r
feature_cols <- setdiff(names(train_df), "Outcome")

explainer <- ubiq_explainer(
  model   = rf_model,
  X_train = train_df[, feature_cols],
  y_train = train_df$Outcome,
  beta    = 5.0,   # how strongly accurate trees are favored
  seed    = 42
)
```

### 3. Explain one patient or sample, with uncertainty

Select a single row of predictors and name the outcome class you want to explain.

```r
patient_x <- test_df[1, feature_cols, drop = FALSE]

explanation <- explain_uncertainty(
  explainer = explainer,
  x         = patient_x,
  n_samples = 500,     # number of alternative tree weightings
  alpha     = 0.5,     # Dirichlet concentration (see parameters below)
  class_idx = "Case"   # class name or integer index
)

# Top 10 features with 95% intervals and reliability flags
print(explanation)
```

### 4. Visualize the results

```r
# Bar chart with uncertainty violins and 95% credible intervals
plot_uncertainty_bars(explanation)

# Full SHAP distribution for one gene or biomarker
plot_uncertainty_distribution(explanation, feature = "Biomarker_A")

# 3-panel comparison of SD, entropy, and sign stability
plot_uncertainty_comparison(explanation)
```

---

## How to interpret the output

`print(explanation)` and `explanation$summary` contain the following columns:

| Column | Plain-language meaning | What to look for |
| :--- | :--- | :--- |
| `mean` | Average SHAP push on the model's prediction (the predicted probability for classification, the predicted value for regression). | Positive values push *toward* the target class (or a higher value); negative values push *away* from it. |
| `std` | Standard deviation across the re-weighted model versions. | Smaller means the model versions agree more. |
| `ci_low` / `ci_high` | 95% credible interval (error bar around `mean`). | If the interval crosses `0`, the direction of effect is uncertain. |
| `sign_stability` | Percentage of model versions that agree on the direction (+ or −). | Above 90% is consistent; near 50% is a coin flip. |
| `entropy` | How spread out the SHAP distribution is. | Lower entropy means a sharper, more concentrated signal. |
| `status` | Automatic triage flag: `Reliable`, `Review Recommended`, or `Unstable Direction`. | Prioritize `Reliable` features; treat `Unstable Direction` with caution. |

---

## Key functions and parameters

### `ubiq_explainer(model, X_train, y_train, beta = 5.0, seed = NULL)`

| Argument | Description |
| :--- | :--- |
| `model` | A trained `ranger`, `randomForest`, or `xgb.Booster` object. |
| `X_train` | Training predictors (data frame). |
| `y_train` | Training outcome vector. |
| `beta` | Softmax temperature (default `5.0`). Higher values give more weight to the best-performing trees. |
| `seed` | Random seed for reproducibility. |

### `explain_uncertainty(explainer, x, n_samples = 500, alpha = 1.0, class_idx = 1L)`

| Argument | Description |
| :--- | :--- |
| `explainer` | Object returned by `ubiq_explainer()`. |
| `x` | A single-row data frame of predictors. |
| `n_samples` | Number of alternative tree weightings (Dirichlet draws) to generate. Default `500`. |
| `alpha` | Dirichlet concentration parameter. Default `1.0`. Smaller values (e.g. `0.5`) produce more uneven weightings, where a few trees dominate each draw. |
| `class_idx` | Name or integer index of the class to explain. Ignored for regression. |

### Plotting functions

| Function | What it shows |
| :--- | :--- |
| `plot_uncertainty_bars(explanation, top_n = 15, title = ...)` | Mean SHAP bars overlaid with distribution violins and dashed 95% error bars. |
| `plot_uncertainty_distribution(explanation, feature, title = ...)` | Full density curve for one feature, with the shaded 95% region and a zero line. |
| `plot_uncertainty_comparison(explanation, top_n = 15, title = ...)` | Three panels comparing standard deviation, entropy, and sign stability (reference cutoffs at 90%, 70%, and 40%). |

---

## Best practices for biological and clinical studies

1. **Include a negative control.** Add a column of random noise to your dataset, the in-silico equivalent of a vehicle or scrambled control. If it shows high `sign_stability` or a narrow interval, your results are likely overconfident: be more conservative about which features you trust, and re-examine your model.
2. **Check direction before magnitude.** A biomarker whose direction flips across model versions (`sign_stability` below 70%) should not support a biological claim, even if its average SHAP value looks large.
3. **Evaluate across a cohort.** Single-patient explanations can be noisy. Loop `explain_uncertainty()` over many test samples and report how often each feature is flagged `Reliable`.
4. **Association is not causation.** SHAP with uncertainty tells you what the *model* reliably depends on, not necessarily the underlying biological mechanism.

---

## Full walkthrough

For a complete tutorial, including a simulated clinical length-of-stay (LOS) cohort with known true drivers, decoy variables, and a template for aggregating results across many patients, open the vignette in R:

```r
vignette("clinical_los_experiment", package = "LukeFrameworkR")
```

---

## Citation and license

If you use `LukeFrameworkR` or the UBiQTree method in your research, please cite the original paper:

> Dubey A, Anžel A, İlgen B, Hattab G. **UbiQTree: Uncertainty quantification in XAI with tree ensembles.** *Patterns* (2026). DOI: [10.1016/j.patter.2025.101454](https://doi.org/10.1016/j.patter.2025.101454)

- **Original Python implementation:** [github.com/dubeyakshat07/UBiQTree](https://github.com/dubeyakshat07/UBiQTree)
- **License:** MIT (see `LICENSE.md`)
