function partest(x,varargin)
%This function calculate the performance, based on Bayes theorem, of a
%clinical test.
%
% Syntax: 	PARTEST(X,ALPHA)
%      
%     Input:
%           X is the following 2x2 matrix.
%           ALPHA - significance level for confidence intervals (default = 0.05).
% 
% ....................Affected(D+)..Healthy(D-)
%                    _______________________
% Positive Test(T+) |   True    |  False    |
%                   | positives | positives |
%                   |___________|___________|
%                   |  False    |   True    |
% Negative Test(T-) | negatives | negatives |
%                   |___________|___________|
% 
%     Outputs:
%           - all that you could compute onto a 2x2 matrix
% 
%      Example: 
% 
%       x=[731 270;78 1500]
%           Calling on Matlab the function: partest(x)
%           Answer is:
% 
% CONFUSION MATRIX
%                      Affected    Healthy
%                      ________    _______
% 
%     Positive_test      731         270  
%     Negative_test       78        1500  
% 
% ----------------------------------------------------------------------------------------------------
% BASIC PARAMETERS
%                                  Parameter
%                                  _________
% 
%     True_Positive_TP                731   
%     False_Negative_FN                78   
%     False_Positive_FP               270   
%     True Negative_TN               1500   
%     Condition_Positive_P            809   
%     Condition_Negative_N           1770   
%     Test_Outcome_Positive_TOP      1001   
%     Test_Outcome_Negative_TON      1578   
%     Population_POP                 2579   
% 
% ----------------------------------------------------------------------------------------------------
% IMBALANCE TEST
%  
% Imbalance ratio (IR): 2.1879
% Z = 18.9233 p-value = 0.0000 - Groups are unbalanced
% ----------------------------------------------------------------------------------------------------
% PERFORMANCE PARAMETERS
%                                 Parameter      LB        UB        Comment  
%                                 _________    ______    ______    ___________
% 
%     Sensitivity_TPR              0.9036      0.8914    0.9146    'None'     
%     False_Negative_Rate_FNR      0.0964      0.0854    0.1086    'None'     
%     Specificity_TNR              0.8475      0.8329     0.861    'None'     
%     False_Positive_Rate_FPR      0.1525       0.139    0.1671    'None'     
%     Area_Under_the_Curve_AUC     0.8755      0.8626    0.8885    'Good test'
%     Distance Index_dInd          0.1805         NaN       NaN    'None'     
%     Gini_Index_GI                 0.751      0.7251     0.777    'None'     
% 
%                                      Parameter
%                                      _________
% 
%     Automatic_Manual                     192  
%     Youden_Index_J                     0.751  
%     Number_needed_to_diagnose_NND     1.3315  
% 
% Around 134 patients need to be tested to correctly detect 100 positive tests for the presence of disease
%  
% ----------------------------------------------------------------------------------------------------
% PREDICTIVE PARAMETERS
%  
% Estimed Prevalence: 0.3137 (0.2959 - 0.3321)
%  
%                                      Parameter      LB        UB                             Comment                        
%                                      _________    ______    ______    ______________________________________________________
% 
%     Positive_Likelihood_Ratio_PLR     5.9235      5.6953    6.1609    'Moderate increase in possibility of disease presence'
%     Negative_Likelihood_Ratio_NLR     0.1138      0.1094    0.1183    'Moderate increase in possibility of disease absence' 
%     Positive_Prediction_Value_PPV     0.7303      0.7126    0.7472    'None'                                                
%     False_Discovery_Rate_FDR          0.2697      0.2528    0.2874    'None'                                                
%     Negative_Predictive_Value_NPV     0.9506      0.9413    0.9585    'None'                                                
%     False_Omission_Rate_FOR           0.0494      0.0415    0.0587    'None'                                                
% 
%                                     Parameter
%                                     _________
% 
%     Prevalence_Threshold_PT          0.2912  
%     Lift_Score_LS                     2.328  
%     Information_Score_IS             1.2191  
%     Markedness_MK                    0.6808  
%     Number_needed_to_predict_NNP     1.4688  
% 
% Around 147 patients need to be screened to correctly predict 100 diagnosis
%  
% ----------------------------------------------------------------------------------------------------
% SIMILARITY PARAMETERS
%                                      Parameter      LB        UB  
%                                      _________    ______    ______
% 
%     Bray_Curtis_Dissimilarity_BCD     0.0372      0.0304    0.0454
%     Similarity_Index_sInd             0.8724      0.8588    0.8849
%     Jaccard_Index_JI                  0.6775       0.659    0.6954
%     Overlap_Coefficient_OC            0.9036      0.8914    0.9146
%     Braun_Blanquet_similarity_BB      0.7303      0.7126    0.7472
%     Otsuka_Ochiai_Coefficient_OOC     0.8123      0.7966    0.8271
% 
% ----------------------------------------------------------------------------------------------------
% ACCURACY PARAMETERS
%                                       Parameter      LB        UB  
%                                       _________    ______    ______
% 
%     Random_Accuracy_RACC               0.1218      0.1095    0.1351
%     Random_Accuracy_Unbiased_RACCU     0.1231      0.1108    0.1366
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     Accuracy_ACC                          0.8651      0.8511    0.8779    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR            0.1349      0.1221    0.1489    'None'                                                                     
%     Identification_Index_II               0.7301      0.7023    0.7558    'None'                                                                     
%     Number_needed_to_screen_NNS           1.3696      1.3231     1.424    'Around 137 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     7.4109      6.7181    8.1896    'Around 742 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     Balanced_Accuracy_bACC                0.8755       0.862    0.8879    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR            0.1245      0.1121     0.138    'None'                                                                     
%     Identification_Index_II                0.751       0.724    0.7758    'None'                                                                     
%     Number_needed_to_screen_NNS           1.3315       1.289    1.3812    'Around 134 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     8.0335      7.2482    8.9199    'Around 804 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     Optimized_Precision_OP                 0.833      0.8179    0.8471    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR             0.167      0.1529    0.1821    'None'                                                                     
%     Identification_Index_II                0.666      0.6358    0.6942    'None'                                                                     
%     Number_needed_to_screen_NNS           1.5015      1.4405    1.5727    'Around 151 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     5.9884      5.4926      6.54    'Around 599 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     F1_Measure                            0.8077      0.7919    0.8227    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR            0.1923      0.1773    0.2081    'None'                                                                     
%     Identification_Index_II               0.6155      0.5837    0.6454    'None'                                                                     
%     Number_needed_to_screen_NNS           1.6248      1.5495    1.7131    'Around 163 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     5.2011       4.805    5.6392    'Around 521 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     Adjusted_F_Score                      0.8947      0.8821    0.9062    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR            0.1053      0.0938    0.1179    'None'                                                                     
%     Identification_Index_II               0.7894      0.7642    0.8124    'None'                                                                     
%     Number_needed_to_screen_NNS           1.2667       1.231    1.3086    'Around 127 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     9.4982      8.4814    10.657    'Around 950 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     G_Measure                             0.8123      0.7966    0.8271    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR            0.1877      0.1729    0.2034    'None'                                                                     
%     Identification_Index_II               0.6246      0.5932    0.6542    'None'                                                                     
%     Number_needed_to_screen_NNS           1.6009      1.5285    1.6858    'Around 161 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     5.3282      4.9164    5.7839    'Around 533 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     Adjusted_G_Measure                    0.8266      0.8113    0.8409    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR            0.1734      0.1591    0.1887    'None'                                                                     
%     Identification_Index_II               0.6532      0.6227    0.6819    'None'                                                                     
%     Number_needed_to_screen_NNS           1.5308      1.4665     1.606    'Around 154 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     5.7677      5.3005    6.2865    'Around 577 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     Normalized_Matthews_Index_nMCI        0.8575      0.8433    0.8707    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR            0.1425      0.1293    0.1567    'None'                                                                     
%     Identification_Index_II               0.7151      0.6866    0.7414    'None'                                                                     
%     Number_needed_to_screen_NNS           1.3984      1.3488    1.4564    'Around 140 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     7.0195      6.3829     7.733    'Around 702 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                   
%                                          _________    ______    ______    ____________________________________________________________________________
% 
%     Gilbert Skill Score                   0.5451      0.5256    0.5644    'Accuracy is moderate'                                                      
%     Misclassification_Rate_MCR            0.4549      0.4356    0.4744    'None'                                                                      
%     Identification_Index_II               0.0902      0.0513    0.1289    'None'                                                                      
%     Number_needed_to_screen_NNS           11.087      7.7605    19.507    'Around 1109 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     2.1983      2.1081    2.2958    'Around 220 patients need to be screened to have 100 mis-diagnosis'         
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     True Skill Statistic                   0.751      0.7338    0.7675    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR             0.249      0.2325    0.2662    'None'                                                                     
%     Identification_Index_II               0.5021      0.4676    0.5351    'None'                                                                     
%     Number_needed_to_screen_NNS           1.9917      1.8688    2.1387    'Around 200 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     4.0167      3.7566    4.3018    'Around 402 patients need to be screened to have 100 mis-diagnosis'        
% 
%                                          Parameter      LB        UB                                        Comment                                  
%                                          _________    ______    ______    ___________________________________________________________________________
% 
%     Heidke Skill Score                    0.7056      0.6768    0.7344    'Accuracy is strong'                                                       
%     Misclassification_Rate_MCR            0.2944       0.277    0.3125    'None'                                                                     
%     Identification_Index_II               0.4112      0.3536    0.4687    'None'                                                                     
%     Number_needed_to_screen_NNS           2.4321      2.1335    2.8278    'Around 244 patients need to be screened to prevent 100 deaths by disease.'
%     Number_needed_to_misdiagnose_NNMD     3.3965      3.2001    3.6107    'Around 340 patients need to be screened to have 100 mis-diagnosis'        
% 
% ----------------------------------------------------------------------------------------------------
% ENTROPY MEASURE
%                                          Parameter      LB        UB                                    Comment                              
%                                          _________    ______    ______    ___________________________________________________________________
% 
%     Confusion_Matrix_Entropy_CEN          0.4844       0.465    0.5039    'Accuracy is weak'                                                 
%     Misclassification_Rate_MCR            0.5156      0.4961     0.535    'None'                                                             
%     Identification_Index_II                  NaN         NaN       NaN    'None'                                                             
%     Number_needed_to_screen_NNS              NaN         NaN       NaN    'Not computable'                                                   
%     Number_needed_to_misdiagnose_NNMD     1.9396      1.8691    2.0158    'Around 194 patients need to be screened to have 100 mis-diagnosis'
% 
%                                          Parameter      LB        UB                                    Comment                              
%                                          _________    ______    ______    ___________________________________________________________________
% 
%     Modified_CEN                          0.4801      0.4606    0.4996    'Accuracy is weak'                                                 
%     Misclassification_Rate_MCR            0.5199      0.5004    0.5394    'None'                                                             
%     Identification_Index_II                  NaN         NaN       NaN    'None'                                                             
%     Number_needed_to_screen_NNS              NaN         NaN       NaN    'Not computable'                                                   
%     Number_needed_to_misdiagnose_NNMD     1.9233       1.854    1.9982    'Around 193 patients need to be screened to have 100 mis-diagnosis'
% 
% ----------------------------------------------------------------------------------------------------
% ASSOCIATION MEASURE
%                                  Parameter      LB         UB                              Comment                            Critical_Bayesian_OR    Bayesian_Credibility
%                                  _________    _______    _______    ______________________________________________________    ____________________    ____________________
% 
%     Test_Bias_TB                   1.2373      0.9474      1.616    'There is not test bias'                                            NaN            'Not needed'       
%     Error_OR                       1.6869      1.2916     2.2032    'Test wrongly classifies in positives'                           1.1718            'Test is credible' 
%     Diagnostic_Odds_Ratio_DOR      52.066      39.865         68    'Test discriminates between diseased and not diseased'           1.0182            'Test is credible' 
%     Yule_Coefficient              0.96231     0.95106    0.97101    'Strong'                                                            NaN            'Not needed'       
%     Tetrachoric_Coefficient       0.92777     0.91054    0.94501    'None'                                                              NaN            'Not needed'       
%     Cramer_V                      0.71409         NaN        NaN    'Strong association'                                                NaN            'Not needed'       
%     Discriminant_Power             2.1791         NaN        NaN    'Fair'  
%    
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo.75@gmail.com
% 
% To cite this file, this would be an appropriate format:
% Cardillo G. (2006). Clinical test performance: the performance of a
% clinical test based on the Bayes theorem. 
% http://www.mathworks.com/matlabcentral/fileexchange/12705

%Input Error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'real','finite','integer','nonnegative','nonnan','size',[2 2]}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
parse(p,x,varargin{:});
alpha=p.Results.alpha;
clear p varargin
z=-realsqrt(2)*erfcinv(2-alpha); 
tr=repmat('-',1,100);
cs=sum(x); %columns sum
rs=sum(x,2); %rows sums
N=sum(x(:)); %numbers of elements
dx=det(x); %matrix determinant

clc
disp('CONFUSION MATRIX')
disp(array2table(x,'RowNames',{'Positive_test';'Negative_test'},'VariableNames',{'Affected';'Healthy'}))
disp(tr)
disp('BASIC PARAMETERS')
disp(array2table([x(:);cs(:);rs(:);N],'RowNames',{'True_Positive_TP';'False_Negative_FN';...
    'False_Positive_FP';'True Negative_TN';'Condition_Positive_P';'Condition_Negative_N';...
    'Test_Outcome_Positive_TOP';'Test_Outcome_Negative_TON';...
    'Population_POP'},'VariableNames',{'Parameter'}))
disp(tr)
disp('IMBALANCE TEST')
disp(' ')
IR=cs(1)/cs(2); if IR<1; IR=1/IR; end
fprintf('Imbalance ratio (IR): %0.4f\n',IR)
zval=diff(cs)/sqrt(N);
pvalue=2 * normcdf(-abs(zval),0,1);
if pvalue<0.05
    fprintf('Z = %0.4f p-value = %0.4f - Groups are unbalanced\n',zval,pvalue)
else
     fprintf('Z = %0.4f p-value = %0.4f - Groups are balanced\n',zval,pvalue)
end
clear zval pvalue
disp(tr)

disp('PERFORMANCE PARAMETERS')
d=diag(x); %true positives and true negatives
%Sensitivity and Specificity
%The Sensitivity is the probability that the test is positive on sick subjects: P(T+|D+) 
%The Specificity is the probability that the test is negative on healthy subjects: P(T-|D-) 
%In Matlab both parameters are obtained with only one instruction:
SS=d./cs'; %Sensitivity and Specificity 
%Of course the false proportion is the complement to 1
fp=1-SS; %false proportions
% 95% confidence intervals
Seci=newcombe(SS(1)); Spci=newcombe(SS(2));
fnci=newcombe(fp(1)); fpci=newcombe(fp(2));

%Area under the ROC curve
%The area under the curve (often referred to as simply the AUC) is equal to
%the probability that a classifier will rank a randomly chosen positive
%instance higher than a randomly chosen negative one (assuming 'positive'
%ranks higher than 'negative'). Thus, AUC corresponds to the arithmetic
%mean of sensitivity and specificity values of each class    
AUC=mean(SS);
%Hanley and McNeil formula
q0=AUC*(1-AUC);
q1=2*AUC^2/(1+AUC)-AUC^2;
q2=AUC/(2-AUC)-AUC^2;
SE=realsqrt([1 (cs-1)]*[q0;q1;q2]/prod(cs));
AUCci=AUC+[-1 1].*z*SE;
clear q0 q1 q2 SE

if AUC==1
    str='Perfect test';
elseif AUC>=0.90 && AUC<1
    str='Excellent test';
elseif AUC>=0.80 && AUC<0.90
    str='Good test';
elseif AUC>=0.70 && AUC<0.80
    str='Fair test';
elseif AUC>=0.60 && AUC<0.70
    str='Poor test';
elseif AUC>=0.50 && AUC<0.60
    str='Fail test';
else
    str='Failed test - less than chance';
end

%Distance Index
%Euclidean distance of a ROC point from the top left corner of the ROC
%space, which can take values between 0 (perfect classification) and
%sqrt(2)
dind=realsqrt(sum(fp.^2));

%Gini Index
%A chance-standardized variant of the AUC is given by Gini coefficient,
%taking values between 0 (no difference between the score distributions of
%the two classes) and 1 (complete separation between the two
%distributions). Gini coefficient is widespread use metric in imbalanced
%data learning    
GI=2*AUC-1;

Parameter=round([SS(1);fp(1);SS(2);fp(2);AUC;dind;GI],4);
LB=round([Seci(1);fnci(1);Spci(1);fpci(1);AUCci(1);NaN;2*AUCci(1)-1],4);
UB=round([Seci(2);fnci(2);Spci(2);fpci(2);AUCci(2);NaN;2*AUCci(2)-1],4);
Comment={'None';'None';'None';'None';str;'None';'None'};

disp(table(Parameter,LB,UB,Comment,...
    'RowNames',{'Sensitivity_TPR';'False_Negative_Rate_FNR';...
    'Specificity_TNR';'False_Positive_Rate_FPR';'Area_Under_the_Curve_AUC';...
    'Distance Index_dInd';'Gini_Index_GI'}...
    ))

clear Spci fnci fpci d AUC AUCci sind GI str parameter LB UB Comment

%AM (Automatic/Manual)
%Difference between automatic and manual classification i.e., the
%difference between positive outcomes and of positive samples. 
AM=rs(1)-cs(1);

%Youden's Index
%Youden's J Index (also called Informedness) is a single statistic that
%captures the performance of a diagnostic test. The use of such a single index
%is "not generally to be recommended". It is equal to the risk difference for a
%dichotomous test and it defined as: J = Sensitivity + Specificity - 1. 
%A perfect test has J=1.
%The informedness of a prediction method as captured by a contingency
%matrix is defined as the probability that the prediction method will make
%a correct decision as opposed to guessing and is calculated using the
%bookmaker algorithm   
J=sum(SS)-1; %Youden's index

%The inverse of J has been defined as the “number needed to diagnose”
%(NND), that is, the number of patients who need to be examined in order
%to correctly detect one person with the disease of interest in a study
%population of persons with and without the known disease    
NND=1/J; %Number needed to Diagnose (NND)

disp(array2table(round([AM;J;NND],4),...
    'RowNames',{'Automatic_Manual';'Youden_Index_J';'Number_needed_to_diagnose_NND'},...
    'VariableNames',{'Parameter'}))
fprintf('Around %i patients need to be tested to correctly detect 100 positive tests for the presence of disease\n',ceil(NND*100))
disp(' ')

clear NND PT

disp(tr)
disp('PREDICTIVE PARAMETERS')
disp(' ')
ButtonName = questdlg('Do you want to input the true prevalence?', 'Prevalence Question', 'Yes', 'No', 'No');
if strcmp(ButtonName,'Yes')
    ButtonName2 = questdlg('Do you want to input the true prevalence as:', 'Prevalence Question', 'Ratio', 'Probability', 'Ratio');
    switch ButtonName2
        case 'Ratio'
            prompt={'Enter the Numerator or the prevalence ratio:','Enter the denominator or the prevalence ratio:'};
            name='Input for Ratio prevalence';
            Ratio=str2double(inputdlg(prompt,name));
            pr=Ratio(1)/Ratio(2);
            clear Ratio
        case 'Probability'
            prompt={'Enter the prevalence probability comprised between 0 and 1:'};
            name='Input for prevalence';
            pr=str2double(inputdlg(prompt,name));
    end
    clear prompt name
    if pr<1e-4
        fprintf('Reported Prevalence: %0.4e\n',pr)
    else
        fprintf('Reported Prevalence: %0.4f\n',pr)
    end
    clear prompt name Ratio 
else
    %Prevalence
    %the prevalence of disease is estimated all D+/N
    pr=(cs(1)/N);
    %95% confidence interval critical value for Prevalence
    prci=newcombe(pr);
    if pr<1e-4
        fprintf('Estimed Prevalence: %0.4e (%0.4e - %0.4e)\n',pr,prci(1),prci(2))
    else
        fprintf('Estimed Prevalence: %0.4f (%0.4f - %0.4f)\n',pr,prci(1),prci(2))
    end
    clear prci
end
disp(' ')
clear ButtonName2

%Positive and Negative Likelihood Ratio
%When we decide to order a diagnostic test, we want to know which test (or
%tests) will best help us rule-in or rule-out disease in our patient.  In the
%language of clinical epidemiology, we take our initial assessment of the
%likelihood of disease ("pre-test probability"), do a test to help us shift our
%suspicion one way or the other, and then determine a final assessment of the
%likelihood of disease ("post-test probability"). 
%Likelihood ratios tell us how much we should shift our suspicion for a
%particular test result. Because tests can be positive or negative, there are at
%least two likelihood ratios for each test.  The "positive likelihood ratio"
%(LR+) tells us how much to increase the probability of disease if the test is
%positive, while the "negative likelihood ratio" (LR-) tells us how much to
%decrease it if the test is negative.
%You can also define the LR+ and LR- in terms of sensitivity and specificity:
%LR+ = sensitivity / (1-specificity)
%LR- = (1-sensitivity) / specificity
plr=SS(1)/fp(2); %Positive Likelihood Ratio
nlr=fp(1)/SS(2); %Negative Likelihood Ratio
cv=realsqrt(sum(1./diag(x))-sum(1./sum(x,2)))*z;
plrci=exp(log(plr)+[-1 1].*cv);
nlrci=exp(log(nlr)+[-1 1].*cv);
plrtxt=dlr(plr);
nlrtxt=dlr(nlr);

%Positive and Negative predictivity
%Positive predictivity is the probability that a subject is sick when test is positive: P(D+|T+)
%Negative predictivity is the probability that a subject is healthy when test is negative: P(D-|T-)
%Positive predictivity=Precision
%apply Bayes Theorem
po=pr/(1-pr);%prior odds
PNp=[1/(1+1/(plr*po)) 1/(1+nlr*po)]; %Positive and Negative predictivity
%The false discovery rate is the probability that the disease is absent when the test is positive: P(D-|T+) 
%The false omission rate is the probability that the disease is present when the test is negative: P(D+|T-)
fd=1-PNp; %false discovery and omission rate
% 95% confidence interval critical value for Positive and Negative predictivity
Ppci=newcombe(PNp(1));
Npci=newcombe(PNp(2));
fdci=newcombe(fd(1)); 
foci=newcombe(fd(2));

Parameter=round([plr;nlr;PNp(1);fd(1);PNp(2);fd(2)],4);
LB=round([plrci(1);nlrci(1);Ppci(1);fdci(1);Npci(1);foci(1)],4);
UB=round([plrci(2);nlrci(2);Ppci(2);fdci(2);Npci(2);foci(2)],4);
Comment={plrtxt;nlrtxt;'None';'None';'None';'None'};

disp(table(Parameter,LB,UB,Comment,...
    'RowNames',{'Positive_Likelihood_Ratio_PLR';'Negative_Likelihood_Ratio_NLR';...
    'Positive_Prediction_Value_PPV';'False_Discovery_Rate_FDR';...
    'Negative_Predictive_Value_NPV';'False_Omission_Rate_FOR'}))

%Below a disease prevalence of PT, the PPV for a screening test with
%these sensitivities and specificities drops significantly and is therefore
%more unreliable.
PT=(sqrt(SS(1)*fp(2))+SS(2)-1)/J;

%LIFT SCORE
%In the context of classification, lift compares model predictions to
%randomly generated predictions. Lift is often used in marketing research
%combined with gain and lift charts as a visual aid   
LS=PNp(1)/pr;

%Information Score
%The amount of information needed to correctly classify an example into
%class C, whose prior probability is p(C), is defined as −log2(p(C)) 
IS=-log2(cs(1)/N)+log2(PNp(1));

%Markedness
%In statistics and psychology, the social science concept of markedness is
%quantified as a measure of how much one variable is marked as a predictor
%or possible cause of another and is also known as △P in simple two-choice
%cases    
MK=sum(PNp)-1;
NNP=1/MK; %Number needed to screen
disp(array2table(round([PT;LS;IS;MK;NNP],4),...
    'RowNames',{'Prevalence_Threshold_PT';'Lift_Score_LS';...
    'Information_Score_IS';'Markedness_MK';'Number_needed_to_predict_NNP'},...
    'VariableNames',{'Parameter'}))
fprintf('Around %i patients need to be screened to correctly predict 100 diagnosis\n',ceil(NNP*100))
disp(' ')

clear plrci nlrci cv  fd po pr Npci fdci foci NNP LS IS

disp(tr)
disp('SIMILARITY PARAMETERS')

%BCD (Bray-Curtis dissimilarity)
%In ecology and biology, the Bray–Curtis dissimilarity, named after J.
%Roger Bray and John T. Curtis, is a statistic used to quantify the
%compositional dissimilarity between two different sites, based on counts
%at each site  
BCD=abs(AM)/(2*N); BCDci=newcombe(BCD);

%Similarity index
%sInd is comprised between 0 (no correct classifications) and 1 (perfect classification)
sind=1-dind/sqrt(2); sindci=newcombe(sind);

%Jaccard Index
%The Jaccard index, also known as Intersection over Union and the Jaccard
%similarity coefficient (originally coined coefficient de communauté by
%Paul Jaccard), is a statistic used for comparing the similarity and
%diversity of sample sets.
%It is also known as critical success index (CSI) or the threat score (TS),
%and it is a verification measure of categorical forecast performance equal
%to the total number of correct event forecasts (hits) divided by the total
%number of storm forecasts plus the number of misses (hits + false alarms +
%misses). The CSI is not affected by the number of non-event forecasts that
%verify (correct rejections). However, the CSI is a biased score that is
%dependent upon the frequency of the event. For an unbiased version of the
%CSI, see  the Gilbert skill score (GS). 
JI=x(1)/(rs(1)+cs(1)-x(1)); JIci=newcombe(JI);

%The overlap coefficient, or Szymkiewicz–Simpson coefficient, is a
%similarity measure that measures the overlap between two finite sets. It
%is defined as the size of the intersection divided by the smaller of the
%size of the two sets   
[OC,I]=max([PNp(1) SS(1)]);
if I==1
    OCci=Ppci;
else
    OCci=Seci;
end

%The Braun-Blanquet coefficient is a similarity measure that is mostly used
%in botany. It is defined as the size of the intersection divided by the
%larger of the size of the two sets   
[BB,I]=min([PNp(1) SS(1)]);
if I==1
    BBci=Ppci;
else
    BBci=Seci;
end

%In biology, there is a similarity index, known as the Otsuka-Ochiai
%coefficient named after Yanosuke Otsuka and Akira Ochiai, also known as
%the Ochiai-Barkman or Ochiai coefficient. If sets are represented as bit
%vectors, the Otsuka-Ochiai coefficient can be seen to be the same as the
%cosine similarity    
OOC=x(1)/geomean([rs(1) cs(1)]); OOCci=newcombe(OOC);

Parameter=round([BCD;sind;JI;OC;BB;OOC],4);
LB=round([BCDci(1);sindci(1);JIci(1);OCci(1);BBci(1);OOCci(1)],4);
UB=round([BCDci(2);sindci(2);JIci(2);OCci(2);BBci(2);OOCci(2)],4);

disp(table(Parameter,LB,UB,...
    'RowNames',{'Bray_Curtis_Dissimilarity_BCD';'Similarity_Index_sInd';...
    'Jaccard_Index_JI';'Overlap_Coefficient_OC';'Braun_Blanquet_similarity_BB';...
    'Otsuka_Ochiai_Coefficient_OOC'}))

clear AM BCD* dind* sind* JI* OC* BB* Ppci Seci

disp(tr)
disp('ACCURACY PARAMETERS')

%Random Accuracy
%The expected accuracy from a strategy of randomly guessing categories
%according to reference and response distributions  
RACC=rs(1)*cs(1)/N^2; RACCci=newcombe(RACC);

%Random Accuracy Unbiased
%The expected accuracy from a strategy of randomly guessing categories
%according to the average of the reference and response distributions 
RACCU=((rs(1)+cs(1))/(2*N))^2; RACCUci=newcombe(RACCU);

Parameter=round([RACC;RACCU;],4);
LB=round([RACCci(1);RACCUci(1);],4);
UB=round([RACCci(2);RACCUci(2);],4);
clear RACC*

disp(table(Parameter,LB,UB,'RowNames',{'Random_Accuracy_RACC';'Random_Accuracy_Unbiased_RACCU';}))

%Accuracy and Mis-classification rate
%Diagnostic accuracy (or Power) is defined as the proportion of all tests
%that give a correct result. The Mis-classification rate is its complement to 1. 

acc=trace(x)/N; %accuracy
acc=[acc;mean(SS)]; %balanced accuracy
acc=[acc;acc(1)-abs(diff(SS))/abs(sum(SS))]; %optimized precision

%In statistics, the F1 score (also F-score or F-measure or Sørensen–Dice coefficient) is a measure of a
%test's accuracy. It considers both the Precision (positive predictivity) 
%and the Sensitivity of the test to compute the score.
%The F1 score can be interpreted as the harmonic mean of the Precision and
%Sensitivity, where an F1 score reaches its best value at 1 and worst score at 0. 
%Similarly, G-Measure (Fowlkes–Mallows index) is the geometric mean of Precision and Sensitivity
acc=[acc;harmmean([SS(1) PNp(1)])]; %F-measure

%The F-measures used only three of the four elements of the confusion
%matrix and hence two classifiers with different TNR values may have the
%same F-score. Therefore, the AGF metric is introduced to use all elements
%of the confusion matrix and provide more weights to samples which are
%correctly classified in the minority class    
acc=[acc;geomean([5*PNp(1)*SS(1)/(4*PNp(1)+SS(1)) (1+0.5^2)*PNp(2)*SS(2)/(0.5^2*PNp(2)+SS(2))])];

%An adjusted version of the geometric mean of specificity and sensitivity 
acc=[acc;geomean([SS(1) PNp(1)])]; %G-measure
Nn=cs(2)/N; acc=[acc;(acc(end)+SS(2)*Nn)/(1+Nn)]; clear Nn 

%Matthews Index
%The Matthews correlation coefficient is used in machine learning as a
%measure of the quality of binary (two-class) classifications, introduced
%by biochemist Brian W. Matthews in 1975. It takes into account true and
%false positives and negatives and is generally regarded as a balanced
%measure that can be used even if the classes are of very different sizes.
%The MCC is, in essence, a correlation coefficient between the observed and
%predicted binary classifications; it returns a value between −1 and +1. A
%coefficient of +1 represents a perfect prediction, 0 no better than random
%prediction and −1 indicates total disagreement between prediction and
%observation         
if strcmp(ButtonName,'Yes')
    MMeasure=geomean([J PSI]);
else
    MMeasure=dx/sqrt(prod(rs)*prod(cs));
end
acc=[acc;(MMeasure+1)/2]; %normalized Matthews
clear MMeasure nMMeasure J PSI PNp

%Gilbert skill score. See Jaccard index
ar=rs(1)*cs(1)/N;
acc=[acc;(x(1)-ar)/(N-x(4)-ar)]; 

clear ar GS

%true skill statistic
%popular combination of sensitivity and false positive rate; measures the
%ability to separate yes and no cases; range from -1 to 1; 
%perfect score=1; no skill level=0;
acc=[acc;dx/prod(cs)];
clear dx TSS 

%Heidke skill score
% It is another form of Cohen's kappa
%measures fractional improvements over random chance; range from -inf to 1;
%perfect score=1; no skill level=0;
f=diag(ones(1,2));
po=sum(sum(x./N.*f));
pe=sum(sum(rs*cs.*f./N^2));
HSS=(po-pe)/(1-pe); %Cohen's kappa
acc=[acc;HSS];
sek=sqrt((po*(1-po))/(N*(1-pe)^2)); %kappa standard error for confidence interval
clear f po pe 

acctxt={'Accuracy_ACC';'Balanced_Accuracy_bACC';'Optimized_Precision_OP';...
    'F1_Measure';'Adjusted_F_Score';'G_Measure';'Adjusted_G_Measure';'Normalized_Matthews_Index_nMCI';...
    'Gilbert Skill Score';'True Skill Statistic';'Heidke Skill Score'};

L=length(acc);
for I=1:L
    if I<L
        accci=newcombe(acc(I)); %Newcombe confidence interval
    else
        accci=HSS+([-1 1].*(z*sek));
        clear HSS sek
    end
    mcr=1-acc(I); mcrci=newcombe(mcr); %Mis-classification rate
    II=2*acc(I)-1; IIci=2.*accci-1; %Identification index
    NNS=1/II; NNSci=fliplr(1./IIci); %Number needed to screen 
    NNMD=1/mcr; NNMDci=fliplr(1./mcrci); %Number needed to misdiagnose
    if acc(I)<0.3
        txt2=' negligible';
    elseif acc(I)>=0.3 && acc(I)<0.5
        txt2=' weak';
    elseif acc(I)>=0.5 && acc(I)<0.7
         txt2=' moderate';
    elseif acc(I)>=0.7 && acc(I)<0.9
         txt2=' strong';
    else
        txt2=' very strong';
    end
    txt=strcat('Accuracy is',txt2);
    Parameter=round([acc(I);mcr;II;NNS;NNMD;],4);
    LB=round([accci(1);mcrci(1);IIci(1);NNSci(1);NNMDci(1);],4);
    UB=round([accci(2);mcrci(2);IIci(2);NNSci(2);NNMDci(2);],4);
    Comment={txt;'None';'None';...
        sprintf('Around %i patients need to be screened to prevent 100 deaths by disease.',ceil(NNS*100));...
        sprintf('Around %i patients need to be screened to have 100 mis-diagnosis',ceil(NNMD*100));};

    disp(table(Parameter,LB,UB,Comment,...
    'RowNames',{acctxt{I};'Misclassification_Rate_MCR';...
    'Identification_Index_II';'Number_needed_to_screen_NNS';...
    'Number_needed_to_misdiagnose_NNMD';}))
end

clear txt* acc* mcr* II NNS NNMD L acc Parameter LB UB Comment

disp(tr)
disp('ENTROPY MEASURE')
%Confusion Entropy
%CEN based upon the concept of entropy for evaluating classifier
%performances. By exploiting the misclassification information of confusion
%matrices, the measure evaluates the confusion level of the class
%distribution of misclassified samples. Both theoretical analysis and
%statistical results show that the proposed measure is more discriminating
%than accuracy and RCI while it remains relatively consistent with the two
%measures. Moreover, it is more capable of measuring how the samples of
%different classes have been separated from each other. Hence the proposed
%measure is more precise than the two measures and can substitute for them
%to evaluate classifiers in classification applications
AF=x(2)+x(3); 
loga=x(2)*log2(x(2))+x(3)*log2(x(3));
trx=trace(rot90(x));
acc(1)=(AF*log2(N^2-diff(diag(x))^2)/2-loga)/N;
acc(2)=2*(trx*log2((N-x(4))*(N-x(1)))-2*loga)/(3*N+trx);
clear AF loga trx

acctxt={'Confusion_Matrix_Entropy_CEN';'Modified_CEN'};

L=length(acc);
for I=1:L
    accci=newcombe(acc(I)); %Newcombe confidence interval
    mcr=1-acc(I); mcrci=newcombe(mcr); %Mis-classification rate
    if acc(I)>0.5
        II=2*acc(I)-1; IIci=2.*accci-1; %Identification index
        NNS=1/II; NNSci=fliplr(1./IIci); %Number needed to screen
    else
        II=NaN; IIci=NaN(1,2); NNS=NaN; NNSci=NaN(1,2);
    end
    NNMD=1/mcr; NNMDci=fliplr(1./mcrci); %Number needed to misdiagnose
    if acc(I)<0.3
        txt2=' negligible';
    elseif acc(I)>=0.3 && acc(I)<0.5
        txt2=' weak';
    elseif acc(I)>=0.5 && acc(I)<0.7
         txt2=' moderate';
    elseif acc(I)>=0.7 && acc(I)<0.9
         txt2=' strong';
    else
        txt2=' very strong';
    end
    txt=strcat('Accuracy is',txt2);
    Parameter=round([acc(I);mcr;II;NNS;NNMD;],4);
    LB=round([accci(1);mcrci(1);IIci(1);NNSci(1);NNMDci(1);],4);
    UB=round([accci(2);mcrci(2);IIci(2);NNSci(2);NNMDci(2);],4);
    if ~isnan(II)
        txt3=sprintf('Around %i patients need to be screened to prevent 100 deaths by disease.',ceil(NNS*100));
    else
        txt3='Not computable';
    end
    Comment={txt;'None';'None';txt3;...
        sprintf('Around %i patients need to be screened to have 100 mis-diagnosis',ceil(NNMD*100));};

    disp(table(Parameter,LB,UB,Comment,...
    'RowNames',{acctxt{I};'Misclassification_Rate_MCR';...
    'Identification_Index_II';'Number_needed_to_screen_NNS';...
    'Number_needed_to_misdiagnose_NNMD';}))
end
clear txt* acc* mcr* II NNS NNMD L acc Parameter LB UB Comment

disp(tr)
disp('ASSOCIATION MEASURE')

orse=realsqrt(sum(1./x(:))); %standard error of log(OR)
cv=([-1 1].*(z*orse)); 
Comment=cell(3,1);
txt=cell(3,1);

%Test Bias (TB)
%A test which shows provable and systematic differences in the results of people
%based on group membership. For example, a test might be considered biased if
%members of one particular gender or race consistently and systematic have
%statistically different results from the rest of the testing population. 
%It is defined as (T+)/(D+)=(TP+FP)/(TP+FN)
%A perfect test has a TB=1;
%If TB<1 the test underestimates the disease because there are more affected than positive test
%If TB>1 the test overestimates the disease because there are more positive test than affected
Parameter(1)=round(rs(1)/cs(1),4); %Test Bias
orci=exp(reallog(Parameter(1))+cv); %OR confidence interval
LB(1)=round(orci(1),4); UB(1)=round(orci(2),4);
flag=0;
if Parameter(1)>1
    if orci(1)<1 && orci(2)>1
        Comment{1}='There is not test bias';
        flag=0;
    else
        Comment{1}='Test overestimates the phenomenon';
        flag=1;
    end
elseif Parameter(1)<1
    if orci(1)<1 && orci(2)>1
        Comment{1}='There is not test bias';
        flag=0;
    else
        Comment{1}='Test underestimates the phenomenon';
        flag=1;
    end
else
    Comment{1}='Test is perfect';
    flag=0;
end
if flag==1
     [cor(1),txt{1}]=bca(Parameter(1),orci);
else
    cor(1)=NaN; txt{1}='Not needed';
end

%Error Odds Ratio. 
%Indicates if the probability of being wrongly classified is highest in the
%diseased or in the non-diseased group. If the error odds is higher than one the
%probability is highest in the diseased group (and the specificity of the test
%is better than the sensitivity), if the value is lower than one the probability
%of an incorrect classification is highest in the non-diseased group (and the
%sensitivity of the test is better than the specificity).     
%It is defined as (Sensitivity/(1-Sensitivity))/(Specificity/(1-Specificity));
Parameter(2)=(SS(1)/fp(1))/(SS(2)/fp(2)); %Error odds ratio
orci=exp(reallog(Parameter(2))+cv); %OR confidence interval
LB(2)=round(orci(1),4); UB(2)=round(orci(2),4);
flag=0;
if Parameter(2)>1
    if orci(1)<1 && orci(2)>1
        Comment{2}='There is not misclassification';
        flag=0;
    else
        Comment{2}='Test wrongly classifies in positives';
        flag=1;
    end
elseif Parameter(2)<1
    if orci(1)<1 && orci(2)>1
        Comment{2}='There is not misclassification';
        flag=0;
    else
        Comment{2}='Test wrongly classifies in negatives';
        flag=1;
    end
else
    Comment{2}='Test is perfect';
    flag=0;
end

if flag==1
     [cor(2),txt{2}]=bca(Parameter(2),orci);
else
    cor(2)=NaN; txt{2}='Not needed';
end

%Diagnostic Odds Ratio. 
%Diagnostic odds ratio is defined as how much more likely will the test
%make a correct diagnosis than an incorrect diagnosis in patients with the
%disease (Scott et al. 2008).  
%Often used as a measure of the discriminative power of the test. Has the value
%one if the test does not discriminate between diseased and not diseased. Very
%high values above one means that a test discriminates well. Values lower than
%one mean that there is something wrong in the application of the test.   

DOR=plr/nlr; Parameter(3)=DOR; %Diagnostic odds ratio 
orci=exp(reallog(DOR)+cv); %OR confidence interval
LB(3)=round(orci(1),4); UB(3)=round(orci(2),4);

flag=0;
if DOR>1
    if DOR<=3 || (orci(1)<1 && orci(2)>1)
        Comment{3}='Test does not discriminate between diseased and not diseased';
    else
        Comment{3}='Test discriminates between diseased and not diseased';
        flag=1;
    end
elseif DOR<1
    if orci(1)<1 && orci(2)>1
        Comment{3}='Test does not discriminate between diseased and not diseased';
    else
        Comment{3}='There is something wrong in the application of the test';
    end
else
    Comment{3}='Test does not discriminate between diseased and not diseased';
end
if flag==1
     [cor(3),txt{3}]=bca(Parameter(3),orci);
else
    cor(3)=NaN; txt{3}='Not needed';
end
clear plr nlr 

%DOR can be normalized into Yule's coefficient.
%If Yule's coefficient=0.5 the test works by chance. 
%This means that a DOR<3 the test in not effective into discriminate healthy from not healthy.
QI=(DOR-1)/(DOR+1); %Yule coefficient
Parameter(4)=QI; LB(4)=(orci(1)-1)/(orci(1)+1); UB(4)=(orci(2)-1)/(orci(2)+1); 
cor(4)=NaN; txt{4}='Not needed';
if QI<=0.25
    Comment{4}='Negligible'; 
elseif QI>0.25 && QI<=0.50
    Comment{4}='Weak';
elseif QI>0.5 && QI<=0.75
    Comment{4}='Moderate';
else
    Comment{4}='Strong';
end
clear QI

TC=cos(pi/(1+sqrt(DOR))); %Tetrachoric Coefficient
ac=DOR^(pi/4); 
SE=realsqrt((pi*ac/(2*(1+ac)^2))^2)*orse;
tcci=TC+[-1 +1].*(SE*z);
Parameter(5)=TC; LB(5)=tcci(1); UB(5)=tcci(2); 
Comment{5}='None'; cor(5)=NaN; txt{5}='Not needed';
clear TC ac SE orse tcci

Phi=(det(x)-N/2)/realsqrt(prod(cs)*prod(rs));
phi_hat=max(0,Phi^2-1/(N-1));
V=realsqrt(phi_hat*(N-1)/(N-2)); %Cramer V
clear Phi phi_hat k_hat
Parameter(6)=V; LB(6)=NaN; UB(6)=NaN;
if V<=0.3
    Comment{6}='Weak association';
elseif (V>0.3 && V<=0.7)
    Comment{6}='Moderate association';
else
    Comment{6}='Strong association';
end
cor(6)=NaN; txt{6}='Not needed';

%Discriminant power
%The discriminant power for a test, also termed the test effectiveness, is a
%measure of how well a test distinguishes between affected and unaffected
%persons. It is the sum of logs of Sensivity and Specificity over own false
%proportion, scaled by the standard deviation of the logistic normal
%distribution curve (square root of 3 divided by π). Test effectiveness is
%interpreted as the standardized distance between the means for both populations.     
dpwr=(realsqrt(3)/pi)*log(DOR); %Discriminant power
Parameter(7)=dpwr; LB(7)=NaN; UB(7)=NaN;
cor(7)=NaN; txt{7}='Not needed';
if dpwr<=1
    Comment{7}='Poor'; 
elseif dpwr>1 && dpwr<=2
    Comment{7}='Limited';
elseif dpwr>2 && dpwr<=3
    Comment{7}='Fair';
else
    Comment{7}='Good';
end
    

disp(table(Parameter',LB',UB',Comment,cor',txt,...
    'RowNames',{'Test_Bias_TB';'Error_OR';'Diagnostic_Odds_Ratio_DOR';...
    'Yule_Coefficient';'Tetrachoric_Coefficient';'Cramer_V';'Discriminant_Power'},...
    'VariableNames',{'Parameter','LB','UB','Comment',...
    'Critical_Bayesian_OR','Bayesian_Credibility'}))

clear flag Parameter LB UB Comment cor txt V dpwr alpha DOR orci flag cv z

%Display graphs
xg=cs./N;
subplot(121)
hold on
H=zeros(1,4);
H(1)=fill(xg(1)+[0 xg(2) xg(2) 0],SS(2)+[0 0 fp(2) fp(2)],'r');
H(2)=fill([0 xg(1) xg(1) 0],fp(1)+[0 0 SS(1) SS(1)],'g');
H(3)=fill([0 xg(1) xg(1) 0],[0 0 fp(1) fp(1)],'y');
H(4)=fill(xg(1)+[0 xg(2) xg(2) 0],[0 0 SS(2) SS(2)],'b');
hold off
axis square
title('PARTEST GRAPH')
xlabel('Subjects proportion')
ylabel('Parameters proportion')
legend(H,'False Positive','True Positive (Sensibility)','False Negative','True Negative (Specificity)','Location','NorthOutside')

%The rose plot is a variation of the common pie chart. For both, we have k data
%points where each point denotes a frequency or a count. Pie charts and rose
%plots both use the area of segments of a circle to convey amounts. 
%The pie chart uses a common radius and varies the central angle according to
%the data. That is, the angle is proportional to the frequency. So if the i-th
%point has count X and the total count is N, the i-th angle is 360*(X/N). 
%For the rose plot, the angle is constant (i.e, divide 360 by the number of
%groups, k) and it is the square root of the radius that is proportional to the
%data. 
%According to Wainer (Wainer (1997), Visual Revelations: Graphical Tales of Fate
%and Deception from Napolean Bonaporte to Ross Perot, Copernicus, Chapter 11),
%the use of a common angle is the strength of the rose plot since it allows us
%to easily compare a sequence of rose plots (i.e, the corresponding segments in
%different rose plots are always in the same relative position). 
%In particular, this makes rose plots an effective technique for displaying the
%data in contingency tables.  
%As an interesting historical note, Wainer points out that rose plots were used
%by Florence Nightingale (she referred to them as coxcombs).  
subplot(122)
axis square
color={'r','g','y','b'};
rp=[fp(2) SS(1) fp(1) SS(2)];
k=length(rp); H=zeros(1,k);
ang=(0:1:k)./k.*2.*pi;
axis equal
hold on
for I=1:k
    theta=[ang(I) linspace(ang(I),ang(I+1),500) ang(I+1)];
    rho=[0 repmat(realsqrt(rp(I)),1,500) 0];
    [xg,yg]=pol2cart(theta,rho);
    H(I)=patch(xg,yg,color{I});
end
hold off
title('ROSEPLOT PARTEST GRAPH')
legend(H,'False Positive','True Positive (Sensibility)','False Negative','True Negative (Specificity)','Location','NorthOutside')

function txt=dlr(lr) %Likelihood dialog
if lr==1
    txt='Test is not suggestive of the presence/absence of disease';
    return
end

if lr>10 || lr<0.1
    p1='Large (often conclusive)';
elseif (lr>5 && lr<=10) || (lr>0.1 && lr<=0.2)
    p1='Moderate';
elseif (lr>2 && lr<=5) || (lr>0.2 && lr<=0.5)
    p1='Low';
elseif (lr>1 && lr<=2) || (lr>0.5 && lr<=1)
    p1='Poor';
end

p2=' increase in possibility of disease';

if lr>1
    p3=' presence';
elseif lr<1
    p3=' absence';
end
txt=strcat(p1,p2,p3);
end 

function ci=newcombe(p)
a=2*N*p+z^2;
b=z*realsqrt(z^2-2-1/N+4*p*(N*(1-p)+1));
c=2*(N+z^2);
%Of course, the critical interval lower bound cannot be less than 0 and the
%upper bound cant be greater than 1 and so:
ci(1)=max([0 (a-b-1)/c]);
ci(2)=min([1 (a+b+1)/c]);
end

function [cor,txt]=bca(or,orci)
lorci=reallog(orci);
if or<1
    sk=-1;
else
    sk=1;
end
cor=exp(sk*diff(lorci)^2/(4*realsqrt(prod(lorci)))); %Critical odds ratio (COR)
if or>cor
    txt='Test is credible';
else
    txt='Test is not credible';
end
end

end
