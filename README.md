# LukeFramework_R


**Luke** (Leverage UbiQtree Keep Explanations) is a biologist-friendly wrapper 
around [UBiQTree](https://github.com/dubeyakshat07/UBiQTree). It takes advanced machine learning 
uncertainty and explainable AI (XAI) from tree ensembles and turns them into clean,
accessible tools for life sciences researchers.


Package: UBiQTreeR
Type: Package
Title: Epistemic Uncertainty Quantification in XAI for Tree Ensembles
Version: 0.1.0
Description: Pure R port of UBiQTree. Quantifies epistemic uncertainty in 
    TreeSHAP explanations using softmax performance weighting and Dirichlet 
    Process sub-ensemble sampling.
License: CC BY-NC-SA 4.0
Imports:
    R6,
    treeshap,
    ranger,
    ggplot2,
    patchwork,
    stats,
    graphics
Suggests:
    randomForest,
    xgboost,
    testthat
Encoding: UTF-8
RoxygenNote: 7.3.1