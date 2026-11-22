# LukeFrameworkR (`Luke`): Trustworthy AI Explanations for Biologists



**Luke** (*Leverage UBiQTree, Keep Explanations*) is a biologist-friendly R package implementing the **UBiQTree** framework ([Dubey et al., *Patterns*, 2026](https://www.cell.com/patterns/fulltext/S2666-3899(25)00302-2); original Python code at [dubeyakshat07/UBiQTree](https://github.com/dubeyakshat07/UBiQTree)). 



It adds **error bars and reliability checks** to SHAP feature explanations from tree-based machine learning models (`ranger`, `randomForest`, and `xgboost`), helping life sciences and clinical researchers separate true biological drivers from random noise.



---



## Why Use LukeFrameworkR?



Standard SHAP values give **one single number** per gene, biomarker, or clinical measurement—with no measure of confidence. In noisy or small biological cohorts, a pure-noise measurement can easily look just as important as a real driver.



* **The Lab Analogy:** Think of a Random Forest as a panel of 100 reviewers who each saw a slightly different subset of patients. Standard SHAP averages their votes. **LukeFrameworkR** repeatedly re-weights those reviewers based on how well they perform and checks whether they still agree.

* **What You Get:** Alongside the average SHAP push, you get a **95% credible interval (error bar)**, **epistemic spread (SD)**, **entropy**, **sign stability** (how often model versions agree on up vs. down regulation), and an automated reliability status (`Reliable`, `Review Recommended`, or `Unstable Direction`).



---



## Installation



Install the development version directly from GitHub using `remotes` (or `devtools`):



```r

# install.packages("remotes")

remotes::install_github("insilicobiology/LukeFramework_R")

```



### Dependencies

`LukeFrameworkR` automatically imports `dplyr`, `ggplot2`, `patchwork`, `ranger`, `rlang`, `tidyr`, and `treeshap`. To explain `randomForest` or `xgboost` models, ensure those packages are also installed.



---



## 5-Minute Quick Start



Below is a minimal workflow using a Random Forest (`ranger`) to predict a clinical or biological outcome and explain a single sample.



### 1. Train your model

Make sure `probability = TRUE` is set for classification tasks so the explainer can access class probabilities.



```r

library(LukeFrameworkR)

library(ranger)



# Example data: predictors + factor outcome column

set.seed(42)

rf_model <- ranger(

  Outcome ~ .,

  data        = train_df,

  num.trees   = 100,

  probability = TRUE,

  seed        = 42

)

```



### 2. Initialize the UBiQTree Explainer

Pass your fitted model along with the **training predictors** (`X_train`) and **training outcome vector** (`y_train`). This step weights each tree by its predictive performance.



```r

feature_cols <- setdiff(names(train_df), "Outcome")



explainer <- ubiq_explainer(

  model   = rf_model,

  X_train = train_df[, feature_cols],

  y_train = train_df$Outcome,

  beta    = 5.0,   # Controls how strongly accurate trees are favored

  seed    = 42

)

```



### 3. Explain a Patient or Sample with Uncertainty

Select a single row of predictors (`patient_x`) and specify which outcome class (`class_idx`) you want to explain.



```r

patient_x <- test_df[1, feature_cols, drop = FALSE]



explanation <- explain_uncertainty(

  explainer = explainer,

  x         = patient_x,

  n_samples = 500,          # Number of alternative tree weightings

  alpha     = 0.5,          # Dirichlet concentration parameter

  class_idx = "Case"        # Class name or integer index

)



# Print the top 10 features with 95% intervals and reliability flags

print(explanation)

```



### 4. Visualize the Results



```r

# 1. Bar chart with uncertainty violins and 95% credible intervals

plot_uncertainty_bars(explanation)



# 2. Posterior SHAP density curve for a specific biomarker/gene

plot_uncertainty_distribution(explanation, feature = "Biomarker_A")



# 3. 3-panel diagnostic comparing SD, Entropy, and Sign Stability

plot_uncertainty_comparison(explanation)

```



---



## How to Interpret the Output



When you run `print(explanation)` or inspect `explanation$summary`, you will see the following columns:



| Column | Plain-Language Meaning | What to Look For |

| :--- | :--- | :--- |

| **`mean`** | Average SHAP push on the predicted probability. | Positive pushes *toward* the target class; negative pushes *away*. |

| **`std`** | Epistemic standard deviation across re-weighted ensembles. | Smaller values mean higher agreement across model versions. |

| **`ci_low` / `ci_high`** | 95% credible interval (error bar around `mean`). | If the range crosses `0`, the direction of effect is uncertain. |

| **`sign_stability`** | Percentage of model versions agreeing on direction (`+` vs `-`). | `> 90%` is consistent; near `50%` is a coin flip. |

| **`entropy`** | How spread out or flat the SHAP distribution is. | Lower entropy indicates a sharp, concentrated signal. |

| **`status`** | Automated triage flag (`Reliable`, `Review Recommended`, `Unstable Direction`). | Prioritize `Reliable` features; treat `Unstable Direction` with caution. |



---



## Key Functions & Parameters



* **`ubiq_explainer(model, X_train, y_train, beta = 5.0, seed = NULL)`**

  * `model`: A trained `ranger`, `randomForest`, or `xgb.Booster` object.

  * `beta`: Softmax temperature parameter (default `5.0`). Higher values give more weight to top-performing trees.

* **`explain_uncertainty(explainer, x, n_samples = 500, alpha = 1.0, class_idx = 1L)`**

  * `x`: A single-row `data.frame` of predictors.

  * `n_samples`: Number of Dirichlet hypothesis draws (default `500`).

  * `alpha`: Dirichlet concentration parameter (default `1.0`; smaller values like `0.5` concentrate weight on the strongest trees).

  * `class_idx`: Name or integer index of the target class to explain (ignored for regression).

* **`plot_uncertainty_bars(explanation, top_n = 15, title = ...)`**

  * Plots mean SHAP bars overlaid with posterior violins and dashed 95% error bars.

* **`plot_uncertainty_distribution(explanation, feature, title = ...)`**

  * Plots the full kernel density distribution for one feature alongside its 95% shaded region and zero-line.

* **`plot_uncertainty_comparison(explanation, top_n = 15, title = ...)`**

  * Generates a 3-panel comparison of **Standard Deviation**, **Entropy**, and **Sign Stability** (with reference cutoffs at 90%, 70%, and 40%).



---



## Best Practices for Biological & Clinical Studies



1. **Include a Negative Control:** Add a column of random noise (like a scrambled siRNA or vehicle control) to your dataset. If it achieves high `sign_stability` or narrow intervals, tighten your thresholds or inspect your model.

2. **Check Direction Before Magnitude:** A biomarker whose direction flips across model versions (`sign_stability < 70%`) should not be used to support a biological claim, even if its average SHAP looks large.

3. **Evaluate Across a Cohort:** Single-patient explanations can be noisy. Loop `explain_uncertainty()` over multiple test samples and report how often each feature is flagged as `Reliable` across the cohort.

4. **Remember: Association != Causation:** SHAP with uncertainty tells you what the *model* reliably depends on—not necessarily the underlying causal mechanism.



---



## Full Walkthrough Vignette



For a complete, step-by-step tutorial—including a simulated clinical Length-of-Stay (LOS) cohort with known true drivers, decoy variables, and a multi-patient aggregation template—open the package vignette in R:



```r

vignette("clinical_los_experiment", package = "LukeFrameworkR")

```



---



## Citation & License



If you use `LukeFrameworkR` or the UBiQTree methodology in your research, please cite the original paper:



> Dubey A, Anžel A, İlgen B, Hattab G. **UbiQTree: Uncertainty quantification in XAI with tree ensembles.** *Patterns* (2026). DOI: [10.1016/j.patter.2025.101454](https://doi.org/10.1016/j.patter.2025.101454)



* **Original Python Implementation:** [github.com/dubeyakshat07/UBiQTree](https://github.com/dubeyakshat07/UBiQTree)

* **License:** MIT License (see `LICENSE.md`)
