# partest
This function calculate the performance, based on Bayes theorem, of a
clinical test.

 Syntax: 	PARTEST(X,ALPHA)
      
     Input:
           X is the following 2x2 matrix.
           ALPHA - significance level for confidence intervals (default = 0.05).
     Outputs:
           All that you could compute onto a 2x2 matrix
 
      Example: 
 
       x=[731 270;78 1500]
           Calling on Matlab the function: partest(x)
           Answer is:
 
 CONFUSION MATRIX
                     Affected    Healthy
                     ________    _______

    Positive_test      731         270  
    Negative_test       78        1500  

----------------------------------------------------------------------------------------------------
BASIC PARAMETERS
                                 Parameter
                                 _________

    True_Positive_TP                731   
    False_Negative_FN                78   
    False_Positive_FP               270   
    True Negative_TN               1500   
    Condition_Positive_P            809   
    Condition_Negative_N           1770   
    Test_Outcome_Positive_TOP      1001   
    Test_Outcome_Negative_TON      1578   
    Population_POP                 2579   

----------------------------------------------------------------------------------------------------
IMBALANCE TEST
 
Imbalance ratio (IR): 2.1879
Z = 18.9233 p-value = 0.0000 - Groups are unbalanced
----------------------------------------------------------------------------------------------------
PERFORMANCE PARAMETERS
                                Parameter      LB        UB        Comment  
                                _________    ______    ______    ___________

    Sensitivity_TPR              0.9036      0.8914    0.9146    'None'     
    False_Negative_Rate_FNR      0.0964      0.0854    0.1086    'None'     
    Specificity_TNR              0.8475      0.8329     0.861    'None'     
    False_Positive_Rate_FPR      0.1525       0.139    0.1671    'None'     
    Area_Under_the_Curve_AUC     0.8755      0.8626    0.8885    'Good test'
    Distance Index_dInd          0.1805         NaN       NaN    'None'     
    Gini_Index_GI                 0.751      0.7251     0.777    'None'     

                                     Parameter
                                     _________

    Automatic_Manual                     192  
    Youden_Index_J                     0.751  
    Number_needed_to_diagnose_NND     1.3315  

Around 134 patients need to be tested to correctly detect 100 positive tests for the presence of disease
 
----------------------------------------------------------------------------------------------------
PREDICTIVE PARAMETERS
 
Estimed Prevalence: 0.3137 (0.2959 - 0.3321)
 
                                     Parameter      LB        UB                             Comment                        
                                     _________    ______    ______    ______________________________________________________

    Positive_Likelihood_Ratio_PLR     5.9235      5.6953    6.1609    'Moderate increase in possibility of disease presence'
    Negative_Likelihood_Ratio_NLR     0.1138      0.1094    0.1183    'Moderate increase in possibility of disease absence' 
    Positive_Prediction_Value_PPV     0.7303      0.7126    0.7472    'None'                                                
    False_Discovery_Rate_FDR          0.2697      0.2528    0.2874    'None'                                                
    Negative_Predictive_Value_NPV     0.9506      0.9413    0.9585    'None'                                                
    False_Omission_Rate_FOR           0.0494      0.0415    0.0587    'None'                                                

                                    Parameter
                                    _________

    Prevalence_Threshold_PT          0.2912  
    Lift_Score_LS                     2.328  
    Information_Score_IS             1.2191  
    Markedness_MK                    0.6808  
    Number_needed_to_predict_NNP     1.4688  

Around 147 patients need to be screened to correctly predict 100 diagnosis
 
----------------------------------------------------------------------------------------------------
SIMILARITY PARAMETERS
                                     Parameter      LB        UB  
                                     _________    ______    ______

    Bray_Curtis_Dissimilarity_BCD     0.0372      0.0304    0.0454
    Similarity_Index_sInd             0.8724      0.8588    0.8849
    Jaccard_Index_JI                  0.6775       0.659    0.6954
    Overlap_Coefficient_OC            0.9036      0.8914    0.9146
    Braun_Blanquet_similarity_BB      0.7303      0.7126    0.7472
    Otsuka_Ochiai_Coefficient_OOC     0.8123      0.7966    0.8271

----------------------------------------------------------------------------------------------------
ACCURACY PARAMETERS
                                      Parameter      LB        UB  
                                      _________    ______    ______

    Random_Accuracy_RACC               0.1218      0.1095    0.1351
    Random_Accuracy_Unbiased_RACCU     0.1231      0.1108    0.1366

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    Accuracy_ACC                          0.8651      0.8511    0.8779    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR            0.1349      0.1221    0.1489    'None'                                                                     
    Identification_Index_II               0.7301      0.7023    0.7558    'None'                                                                     
    Number_needed_to_screen_NNS           1.3696      1.3231     1.424    'Around 137 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     7.4109      6.7181    8.1896    'Around 742 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    Balanced_Accuracy_bACC                0.8755       0.862    0.8879    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR            0.1245      0.1121     0.138    'None'                                                                     
    Identification_Index_II                0.751       0.724    0.7758    'None'                                                                     
    Number_needed_to_screen_NNS           1.3315       1.289    1.3812    'Around 134 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     8.0335      7.2482    8.9199    'Around 804 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    Optimized_Precision_OP                 0.833      0.8179    0.8471    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR             0.167      0.1529    0.1821    'None'                                                                     
    Identification_Index_II                0.666      0.6358    0.6942    'None'                                                                     
    Number_needed_to_screen_NNS           1.5015      1.4405    1.5727    'Around 151 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     5.9884      5.4926      6.54    'Around 599 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    F1_Measure                            0.8077      0.7919    0.8227    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR            0.1923      0.1773    0.2081    'None'                                                                     
    Identification_Index_II               0.6155      0.5837    0.6454    'None'                                                                     
    Number_needed_to_screen_NNS           1.6248      1.5495    1.7131    'Around 163 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     5.2011       4.805    5.6392    'Around 521 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    Adjusted_F_Score                      0.8947      0.8821    0.9062    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR            0.1053      0.0938    0.1179    'None'                                                                     
    Identification_Index_II               0.7894      0.7642    0.8124    'None'                                                                     
    Number_needed_to_screen_NNS           1.2667       1.231    1.3086    'Around 127 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     9.4982      8.4814    10.657    'Around 950 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    G_Measure                             0.8123      0.7966    0.8271    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR            0.1877      0.1729    0.2034    'None'                                                                     
    Identification_Index_II               0.6246      0.5932    0.6542    'None'                                                                     
    Number_needed_to_screen_NNS           1.6009      1.5285    1.6858    'Around 161 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     5.3282      4.9164    5.7839    'Around 533 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    Adjusted_G_Measure                    0.8266      0.8113    0.8409    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR            0.1734      0.1591    0.1887    'None'                                                                     
    Identification_Index_II               0.6532      0.6227    0.6819    'None'                                                                     
    Number_needed_to_screen_NNS           1.5308      1.4665     1.606    'Around 154 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     5.7677      5.3005    6.2865    'Around 577 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    Normalized_Matthews_Index_nMCI        0.8575      0.8433    0.8707    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR            0.1425      0.1293    0.1567    'None'                                                                     
    Identification_Index_II               0.7151      0.6866    0.7414    'None'                                                                     
    Number_needed_to_screen_NNS           1.3984      1.3488    1.4564    'Around 140 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     7.0195      6.3829     7.733    'Around 702 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                   
                                         _________    ______    ______    ____________________________________________________________________________

    Gilbert Skill Score                   0.5451      0.5256    0.5644    'Accuracy is moderate'                                                      
    Misclassification_Rate_MCR            0.4549      0.4356    0.4744    'None'                                                                      
    Identification_Index_II               0.0902      0.0513    0.1289    'None'                                                                      
    Number_needed_to_screen_NNS           11.087      7.7605    19.507    'Around 1109 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     2.1983      2.1081    2.2958    'Around 220 patients need to be screened to have 100 mis-diagnosis'         

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    True Skill Statistic                   0.751      0.7338    0.7675    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR             0.249      0.2325    0.2662    'None'                                                                     
    Identification_Index_II               0.5021      0.4676    0.5351    'None'                                                                     
    Number_needed_to_screen_NNS           1.9917      1.8688    2.1387    'Around 200 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     4.0167      3.7566    4.3018    'Around 402 patients need to be screened to have 100 mis-diagnosis'        

                                         Parameter      LB        UB                                        Comment                                  
                                         _________    ______    ______    ___________________________________________________________________________

    Heidke Skill Score                    0.7056      0.6768    0.7344    'Accuracy is strong'                                                       
    Misclassification_Rate_MCR            0.2944       0.277    0.3125    'None'                                                                     
    Identification_Index_II               0.4112      0.3536    0.4687    'None'                                                                     
    Number_needed_to_screen_NNS           2.4321      2.1335    2.8278    'Around 244 patients need to be screened to prevent 100 deaths by disease.'
    Number_needed_to_misdiagnose_NNMD     3.3965      3.2001    3.6107    'Around 340 patients need to be screened to have 100 mis-diagnosis'        

----------------------------------------------------------------------------------------------------
ENTROPY MEASURE
                                         Parameter      LB        UB                                    Comment                              
                                         _________    ______    ______    ___________________________________________________________________

    Confusion_Matrix_Entropy_CEN          0.4844       0.465    0.5039    'Accuracy is weak'                                                 
    Misclassification_Rate_MCR            0.5156      0.4961     0.535    'None'                                                             
    Identification_Index_II                  NaN         NaN       NaN    'None'                                                             
    Number_needed_to_screen_NNS              NaN         NaN       NaN    'Not computable'                                                   
    Number_needed_to_misdiagnose_NNMD     1.9396      1.8691    2.0158    'Around 194 patients need to be screened to have 100 mis-diagnosis'

                                         Parameter      LB        UB                                    Comment                              
                                         _________    ______    ______    ___________________________________________________________________

    Modified_CEN                          0.4801      0.4606    0.4996    'Accuracy is weak'                                                 
    Misclassification_Rate_MCR            0.5199      0.5004    0.5394    'None'                                                             
    Identification_Index_II                  NaN         NaN       NaN    'None'                                                             
    Number_needed_to_screen_NNS              NaN         NaN       NaN    'Not computable'                                                   
    Number_needed_to_misdiagnose_NNMD     1.9233       1.854    1.9982    'Around 193 patients need to be screened to have 100 mis-diagnosis'

----------------------------------------------------------------------------------------------------
ASSOCIATION MEASURE
                                 Parameter      LB         UB                              Comment                            Critical_Bayesian_OR    Bayesian_Credibility
                                 _________    _______    _______    ______________________________________________________    ____________________    ____________________

    Test_Bias_TB                   1.2373      0.9474      1.616    'There is not test bias'                                            NaN            'Not needed'       
    Error_OR                       1.6869      1.2916     2.2032    'Test wrongly classifies in positives'                           1.1718            'Test is credible' 
    Diagnostic_Odds_Ratio_DOR      52.066      39.865         68    'Test discriminates between diseased and not diseased'           1.0182            'Test is credible' 
    Yule_Coefficient              0.96231     0.95106    0.97101    'Strong'                                                            NaN            'Not needed'       
    Tetrachoric_Coefficient       0.92777     0.91054    0.94501    'None'                                                              NaN            'Not needed'       
    Cramer_V                      0.71409         NaN        NaN    'Strong association'                                                NaN            'Not needed'       
    Discriminant_Power             2.1791         NaN        NaN    'Fair'  
    
  
 To cite this file, this would be an appropriate format:
 Cardillo G. (2006). Clinical test performance: the performance of a
 clinical test based on the Bayes theorem. 
 http://www.mathworks.com/matlabcentral/fileexchange/12705
