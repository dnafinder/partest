# partest
This function calculate the performance, based on Bayes theorem, of a clinical test.<br/>
Syntax: 	PARTEST(X,ALPHA)
     
    Input:
          X is the following 2x2 matrix. [true positive false positive; false negative true negative];
          ALPHA - significance level for confidence intervals (default = 0.05).

    Outputs:
          - Prevalence
          - Sensibility
          - Specificity
          - False positive and negative proportions
          - False discovery and discovery rates
          - Youden's Index and Number Needed to Diagnose (NDD)
          - Positive predictivity
          - Positive Likelihood Ratio
          - Negative predictivity
          - Negative Likelihood Ratio
          - Predictive Summary Index (PSI) and Number Needed to Screen (NNS)
          - Test Accuracy
          - Mis-classification Rate
          - F-Measure
          - Test bias
          - Error odds ratio
          - Diagnostic odds ratio
          - Discriminant Power

     Example: 

      x=[731 270;78 1500]
          Calling on Matlab the function: partest(x)
          Answer is:

DIAGNOSTIC TEST PERFORMANCE PARAMETERS
----------------------------------------------------------------------------------------------------
Prevalence: 31.4% (29.6% - 33.2%)
 
Sensitivity (probability that test is positive on unhealthy subject): 90.4% (89.1% - 91.5%)
False negative proportion: 9.6% (8.5% - 10.9%)
False discovery rate: 27.0% (25.3% - 28.7%)
 
Specificity (probability that test is negative on healthy subject): 84.7% (83.3% - 86.1%)
False positive proportion: 15.3% (13.9% - 16.7%)
False omission rate: 4.9% (4.2% - 5.9%)
 
Youden's Index (a perfect test would have a Youden's index of +1): 0.7510
Number Needed to Diagnose (NND): 1.33
Around 14 persons need to be tested to return 10 positive tests for the presence of disease
 
Precision or Predictivity of positive test
(probability that a subject is unhealthy when test is positive): 73.0% (71.3% - 74.7%)
Positive Likelihood Ratio: 5.9 (5.7 - 6.2)
Moderate increase in possibility of disease presence
 
Predictivity of negative test
(probability that a subject is healthy when test is negative): 95.1% (94.1% - 95.8%)
Negative Likelihood Ratio: 0.1138 (0.1094 - 0.1183)
Moderate increase in possibility of disease absence
 
Predictive Summary Index: 0.6808
Number Needed to Screen (NNS): 1.47
Around 15 persons need to be screened to avoid 10 events (i.e. death) for the presence of disease
 
Accuracy or Potency: 86.5% (85.1% - 87.8%)
Mis-classification Rate: 13.5% (12.2% - 14.9%)
F-measure: 80.8% (79.2% - 82.3%)
 
Test bias: 1.2373 (0.9474 - 1.6160)
Test overestimates the phenomenon
Error odds ratio: 1.6869 (1.2916 - 2.2032)
Diagnostic odds ratio: 52.0655 (39.8649 - 68.0002)1.0968
Discriminant Power: 2.2
     A test with a discriminant value of 1 is not effective in discriminating between affected and unaffected individuals.
     A test with a discriminant value of 3 is effective in discriminating between affected and unaffected individuals.
   

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2006). Clinical test performance: the performance of a
clinical test based on the Bayes theorem. 
http://www.mathworks.com/matlabcentral/fileexchange/12705
