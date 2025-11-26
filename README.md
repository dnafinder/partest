[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/partest)

🌐 Overview
This repository contains the MATLAB function `partest.m`, which performs an extensive diagnostic, predictive, similarity, accuracy, association, and entropy-based analysis of a binary diagnostic test given a 2×2 confusion matrix.

⭐ Features
- Confusion matrix display and basic counts (TP, FN, FP, TN, totals)
- Imbalance test between groups (Z test)
- Sensitivity, Specificity, FNR, FPR with Newcombe–Wilson confidence intervals
- AUC, Distance Index, and Gini Index with confidence intervals
- Youden’s J (Informedness) and Number Needed to Diagnose (NND)
- Likelihood ratios (PLR, NLR) with confidence intervals and qualitative interpretation
- Predictive values (PPV, NPV), FDR, FOR, prevalence threshold, Lift and Information Score
- Markedness (MK) and Number Needed to Predict (NNP)
- Similarity and overlap indices (Bray–Curtis, Jaccard, Otsuka–Ochiai, Overlap, Braun–Blanquet)
- Accuracy metrics: ACC, bACC, Optimized Precision, F1, adjusted F, G-measure, adjusted G
- Matthews Correlation Coefficient (MCC) and normalized Matthews index (with MK/J formulation if prevalence is provided)
- Skill scores: Gilbert Skill Score, True Skill Statistic, Heidke Skill Score (Cohen’s kappa)
- Confusion Entropy (CEN) and modified CEN with derived indices
- Association measures: Test Bias, Error Odds Ratio, Diagnostic Odds Ratio (DOR), Yule coefficient, tetrachoric coefficient, Cramer’s V, discriminant power, Bayesian credibility assessment
- Two summary plots (stacked “PARTEST GRAPH” and “ROSEPLOT PARTEST GRAPH”)

🛠️ Installation
1. Download or clone the repository:
   https://github.com/dnafinder/partest
2. Add the folder to your MATLAB path:
   addpath('path_to_folder');

▶️ Usage
Call the function with a 2×2 confusion matrix:
   partest(X)

Optionally specify the confidence level (alpha):
   partest(X, alpha)

Optionally specify both alpha and prevalence:
   partest(X, alpha, pr)

🔣 Inputs
X:
   2×2 numeric matrix of non-negative integers with layout

```text

      Affected   Healthy
    ------------------------
T+       TP        FN
T-       FP        TN

```

Example matrix:
   x = [A B;
        C D];


ALPHA (optional):
   Scalar in (0,1), significance level for confidence intervals.
   Default: 0.05 (95% confidence).

PR (optional):
   Scalar in [0,1], disease prevalence.
   • If provided, it is treated as the true prevalence.
   • If omitted or set to NaN, prevalence is estimated from the confusion matrix as:

     pr = (TP + FN) / (TP + FN + FP + TN)

📤 Outputs
If called without output:
   Results are printed to the Command Window and two figures are generated:

   • Tables for:
     - basic parameters and imbalance test
     - performance metrics (Sensitivity, Specificity, AUC, etc.)
     - predictive metrics (PLR, NLR, PPV, NPV, FDR, FOR, MK, NNP, etc.)
     - similarity indices (Bray–Curtis, Jaccard, Overlap, Otsuka–Ochiai, etc.)
     - accuracy and skill scores (ACC, bACC, F1, MCC, TSS, HSS, etc.)
     - entropy measures (CEN, modified CEN)
     - association and credibility measures (TB, Error OR, DOR, Yule, Cramer’s V, COR, etc.)

   • Plots:
     - PARTEST GRAPH: stacked areas of TP, FP, FN, TN vs subject proportion
     - ROSEPLOT PARTEST GRAPH: rose/coxcomb plot of error and success components

If called with output:
   The function is currently designed as a display/visualization tool only and does not return a structure.
   If you need programmatic access to the metrics, you can adapt the code so that all computed indices are stored in a struct (e.g. stats) before being displayed.

📘 Confusion matrix structure
The expected confusion matrix is:

   X = [ TP  FN
         FP  TN ];

where:
   TP = True Positives (affected subjects with positive test)
   FN = False Negatives (affected subjects with negative test)
   FP = False Positives (healthy subjects with positive test)
   TN = True Negatives (healthy subjects with negative test)

From this matrix, `partest` derives:
   • Column sums: Condition Positive (P), Condition Negative (N)
   • Row sums: Test Outcome Positive (TOP), Test Outcome Negative (TON)
   • Total population (POP)

📘 Metrics and interpretation (summary)

Basic performance:
   • Sensitivity (TPR) = TP / (TP + FN)
   • Specificity (TNR) = TN / (FP + TN)
   • False Negative Rate (FNR) = 1 − TPR
   • False Positive Rate (FPR) = 1 − TNR
   • Confidence intervals for proportions via Newcombe–Wilson score method

AUC and related indices:
   • AUC (approximate) = (TPR + TNR) / 2
   • Distance Index d_ind = sqrt(FNR^2 + FPR^2)
   • Gini Index G = 2·AUC − 1

Imbalance & diagnostic metrics:
   • Imbalance ratio (IR) between groups
   • Youden’s J (Informedness): J = TPR + TNR − 1
   • Number Needed to Diagnose: NND = 1 / J

Predictive values & likelihood ratios:
   • Positive Likelihood Ratio: PLR = TPR / FPR
   • Negative Likelihood Ratio: NLR = FNR / TNR
   • Positive Predictive Value (PPV) and Negative Predictive Value (NPV) via Bayes’ theorem
   • False Discovery Rate (FDR) and False Omission Rate (FOR)
   • Prevalence Threshold (PT), Lift Score (LS), Information Score (IS)
   • Markedness: MK = PPV + NPV − 1
   • Number Needed to Predict: NNP = 1 / MK

Similarity indices:
   • Bray–Curtis Dissimilarity (BCD)
   • Similarity Index (sInd) derived from distance index
   • Jaccard Index (JI)
   • Overlap Coefficient (OC)
   • Braun–Blanquet similarity (BB)
   • Otsuka–Ochiai Coefficient (OOC), equivalent to cosine similarity

Accuracy & skill scores:
   • Random Accuracy (RACC) and Unbiased Random Accuracy (RACCU)
   • Accuracy (ACC), Balanced Accuracy (bACC), Optimized Precision (OP)
   • F1 score (harmonic mean of Sensitivity and Precision)
   • Adjusted F-score (AGF)
   • G-measure (geometric mean of Sensitivity and Precision) and adjusted G-measure
   • Matthews Correlation Coefficient (MCC) and normalized Matthews index:
       - if PR is provided, MCC is derived from Youden’s J and Markedness MK
       - otherwise, MCC is computed from the 2×2 determinant formula
   • Gilbert Skill Score (GSS)
   • True Skill Statistic (TSS)
   • Heidke Skill Score (HSS, Cohen’s kappa)

Derived indices for each accuracy-like metric:
   • Misclassification Rate (MCR) = 1 − ACC
   • Identification Index (II) = 2·ACC − 1
   • Number Needed to Screen (NNS) = 1 / II (when defined)
   • Number Needed to Misdiagnose (NNMD) = 1 / MCR
   • Qualitative interpretation of accuracy strength (negligible, weak, moderate, strong, very strong)

Entropy measures:
   • Confusion Matrix Entropy (CEN)
   • Modified CEN
   • For each: confidence intervals, MCR, II, NNS, NNMD (when meaningful)

Association measures:
   • Test Bias (TB) and interpretation (perfect, over-/under-estimation)
   • Error Odds Ratio (error distribution between groups)
   • Diagnostic Odds Ratio (DOR) with confidence interval and qualitative assessment
   • Yule coefficient (normalized DOR) with strength classification
   • Tetrachoric Coefficient
   • Cramer’s V with association strength classification (weak, moderate, strong)
   • Discriminant Power
   • Bayesian Credibility Assessment via Critical Odds Ratio (COR) and “credible / not credible” text

📝 Notes
• The function is non-interactive: there are no dialog boxes; prevalence is either provided as argument or estimated from the table.  
• No specialized toolbox is required: geometric and harmonic means are computed via explicit formulas.  
• Zeros in the 2×2 table are not automatically corrected; if needed, continuity corrections must be applied by the user before calling `partest`.  
• The main purpose is an in-depth interpretative report. For a compact numeric API, you can wrap or modify the function to return a `stats` structure.

📚 Citation
If you use this function in scientific work, you may cite it as:

Cardillo G. (2025). PARTEST: performance analysis of a binary diagnostic test on a 2×2 confusion matrix.  
GitHub: https://github.com/dnafinder/partest

👤 Author
Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

⚖️ License
This project is released under the GNU GPL-3.0 license.
