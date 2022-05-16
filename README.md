# partest
This function calculate the performance, based on Bayes theorem, of a
clinical test.

 Syntax: 	PARTEST(X,ALPHA)
      
     Input:
           X is the following 2x2 matrix.
           ALPHA - significance level for confidence intervals (default = 0.05).
     Outputs:
           - Prevalence
           - Sensibility and false negative rate
           - Specificity and false positive rate
           - Youden's Index and Number Needed to Diagnose (NDD)
           - False discovery and discovery rates
           - Prevalence threshold
           - Positive and Negative Likelihood Ratio
           - Positive and Negative predictivity
           - False discovery and omission rate
           - Predictive Summary Index (PSI) and Number Needed to Predict (NNP)
           - Test Accuracy, Identification Index and Number Needed to Screen (NNS)
           - Misclassification Rate and Number Needed to Misdiagnose (NNMD)
           - F-Measure; G- Measure; Matthews Index
           - Threat score (TS) or critical success index (CSI)
           - Gilbert Skill Score; True Skill Statistic; Heidke Skill Score
           - Test bias
           - Error odds ratio and Cramer's V
           - Risk Ratio
           - Diagnostic odds ratio
           - Yule's Coefficient
           - Tetrachoric Coefficient
           - Discriminant Power
 
      Example: 
 
       x=[731 270;78 1500]
           Calling on Matlab the function: partest(x)
           Answer is:
 
 DIAGNOSTIC TEST PERFORMANCE PARAMETERS
 ----------------------------------------------------------------------------------------------------
 Estimed Prevalence: 0.3137 (0.2959 - 0.3321)
  
 Sensitivity, Recall, Hit rate or True Positive rate (TPR): 0.9036 (0.8914 - 0.9146)
 Miss rate or False Negative rate (FNR): 0.0964 (0.0854 - 0.1086)
  
 Specificity, Selectivity or True Negative Rate (TNR): 0.8475 (0.8329 - 0.8610)
 Fall-out or False Positive rate (FPR): 0.1525 (0.1390 - 0.1671)
  
 Informedness - Youden's Index - J: 0.7510
 Number Needed to Diagnose (NND): 1.33
 Around 134 patients need to be tested to correctly detect 100 positive tests for the presence of disease
  
 Prevalence Threshold (PT): 0.2912
  
 Positive Likelihood Ratio (PLR): 5.9235 (5.6953 - 6.1609)
 Moderate increase in possibility of disease presence
  
 Negative Likelihood Ratio (NLR): 0.1138 (0.1094 - 0.1183)
 Moderate increase in possibility of disease absence
  
 Precision or Positive Predictive Value (PPV): 0.7303 (0.7126 - 0.7472)
 False discovery rate (FDR): 0.2697 (0.2528 - 0.2874)
  
 Negative Predictive Value (NPV): 0.9506 (0.9413 - 0.9585)
 False omission rate (FOR): 0.0494 (0.0415 - 0.0587)
  
 Markedness - Predictive Summary Index (PSI): 0.6808
 Number Needed to Predict (NNP): 1.47
 Around 147 patients need to be screened to correctly predict 100 diagnosis
  
 Accuracy or Potency (ACC): 0.8651 (0.8511 - 0.8779)
 Identification Index (II): 0.7301
 Number Needed to Screen (NNS): 1.37
 Around 137 patients need to be screened to prevent 100 deaths by disease.
  
 Mis-classification Rate (MCR): 0.1349 (0.1221 - 0.1489)
 Number Needed to Mis-diagnose (NNMD): 7.41
 Around 742 patients need to be screened to have 100 mis-diagnosis
  
 Balanced Accuracy or Potency (BA): 0.8755 (0.8620 - 0.8879)
 Balanced Identification Index (BII): 0.7510
 Balanced Number Needed to Screen (NNS): 1.33
 Around 134 patients need to be screened to prevent 100 deaths by disease.
  
 Balanced Mis-classification Rate (BMCR): 0.1349 (0.1221 - 0.1489)
 Number Needed to Mis-diagnose (BNNMD): 7.41
 Around 742 patients need to be screened to have 100 mis-diagnosis
  
 F1-measure (Sørensen–Dice index): 0.8077 (0.7919 - 0.8227)
 G-measure (Fowlkes–Mallows index): 0.8123 (0.7966 - 0.8271)
 Matthews index: 0.7151 - Normalized Index 0.8575 (0.8433 - 0.8707)
  
 Threat score (TS) or critical success index (CSI): 0.6775 (0.6590 - 0.6954)
 Gilbert Skill Score: 0.5451 (0.5256 - 0.5644)
 True Skill Statistic (Hanssen-Kuipper skill score; Pierce's skill score): 0.7510 (0.7338 - 0.7675)
 Heidke Skill Score (Cohen's Kappa): 0.7056 (0.6768 - 0.7344)
  
 Test bias: 1.2373 (0.9474 - 1.6160) - Test overestimates the phenomenon
 Error odds ratio: 1.6869 (1.2916 - 2.2032)
 Diagnostic odds ratio: 52.0655 (39.8649 - 68.0002)
 Normalized Diagnostic odds ratio: 0.9623
 Cramer's V: 0.7141
 Strong positive association (risk factor)
 
 Bayesian Credibility Assessment
 Critical Diagnostic Odds Ratio: 1.0182
 DOR>COR. Test is credible at the 95%
  
 Risk Ratio: 11.3119<14.7739<19.2955
 Absolute risk reduction: 68.1%
 Relative risk reduction: 93.2%
  
 Yule's Coefficient: 0.9623
 Tetrachoric Coefficient: 0.9278
 Discriminant Power: 2.1791
    
 
           Created by Giuseppe Cardillo
           giuseppe.cardillo.75@gmail.com
 
 To cite this file, this would be an appropriate format:
 Cardillo G. (2006). Clinical test performance: the performance of a
 clinical test based on the Bayes theorem. 
 http://www.mathworks.com/matlabcentral/fileexchange/12705
