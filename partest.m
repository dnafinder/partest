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
%           - Prevalence
%           - Imbalance ratio
%           - Sensibility and false negative rate
%           - Specificity and false positive rate
%           - Youden's Index and Number Needed to Diagnose (NDD)
%           - False discovery and discovery rates
%           - Prevalence threshold
%           - Positive and Negative Likelihood Ratio
%           - Positive and Negative predictivity
%           - False discovery and omission rate
%           - Predictive Summary Index (PSI) and Number Needed to Predict (NNP)
%           - Test Accuracy, Identification Index and Number Needed to Screen (NNS)
%           - Misclassification Rate and Number Needed to Misdiagnose (NNMD)
%           - F-Measure; G- Measure; Matthews Index
%           - Threat score (TS) or critical success index (CSI)
%           - Gilbert Skill Score; True Skill Statistic; Heidke Skill Score
%           - Test bias
%           - Error odds ratio 
%           - Diagnostic odds ratio
%           - Yule's Coefficient
%           - Tetrachoric Coefficient
%           - Cramer's V
%           - Discriminant Power
% 
%      Example: 
% 
%       x=[731 270;78 1500]
%           Calling on Matlab the function: partest(x)
%           Answer is:
% 
% DIAGNOSTIC TEST PERFORMANCE PARAMETERS
% ----------------------------------------------------------------------------------------------------
% Estimed Prevalence: 0.3137 (0.2959 - 0.3321)
%
% Imbalance ratio (IR): 2.1879
% Z = -18.9233 p-value = 0.0000 - Groups are unbalanced
%  
% Sensitivity, Recall, Hit rate or True Positive rate (TPR): 0.9036 (0.8914 - 0.9146)
% Miss rate or False Negative rate (FNR): 0.0964 (0.0854 - 0.1086)
%  
% Specificity, Selectivity or True Negative Rate (TNR): 0.8475 (0.8329 - 0.8610)
% Fall-out or False Positive rate (FPR): 0.1525 (0.1390 - 0.1671)
%  
% Informedness - Youden's Index - J: 0.7510
% Number Needed to Diagnose (NND): 1.33
% Around 134 patients need to be tested to correctly detect 100 positive tests for the presence of disease
%  
% Prevalence Threshold (PT): 0.2912
%  
% Positive Likelihood Ratio (PLR): 5.9235 (5.6953 - 6.1609)
% Moderate increase in possibility of disease presence
%  
% Negative Likelihood Ratio (NLR): 0.1138 (0.1094 - 0.1183)
% Moderate increase in possibility of disease absence
%  
% Precision or Positive Predictive Value (PPV): 0.7303 (0.7126 - 0.7472)
% False discovery rate (FDR): 0.2697 (0.2528 - 0.2874)
%  
% Negative Predictive Value (NPV): 0.9506 (0.9413 - 0.9585)
% False omission rate (FOR): 0.0494 (0.0415 - 0.0587)
%  
% Markedness - Predictive Summary Index (PSI): 0.6808
% Number Needed to Predict (NNP): 1.47
% Around 147 patients need to be screened to correctly predict 100 diagnosis
%  
% Accuracy or Potency (ACC): 0.8651 (0.8511 - 0.8779)
% Identification Index (II): 0.7301
% Number Needed to Screen (NNS): 1.37
% Around 137 patients need to be screened to prevent 100 deaths by disease.
%  
% Mis-classification Rate (MCR): 0.1349 (0.1221 - 0.1489)
% Number Needed to Mis-diagnose (NNMD): 7.41
% Around 742 patients need to be screened to have 100 mis-diagnosis
%  
% Balanced Accuracy or Potency (BA): 0.8755 (0.8620 - 0.8879)
% Balanced Identification Index (BII): 0.7510
% Balanced Number Needed to Screen (NNS): 1.33
% Around 134 patients need to be screened to prevent 100 deaths by disease.
%  
% Balanced Mis-classification Rate (BMCR): 0.1245 (0.1121 - 0.1380)
% Number Needed to Mis-diagnose (BNNMD): 8.03
% Around 804 patients need to be screened to have 100 mis-diagnosis
%  
% F1-measure (Sørensen–Dice index): 0.8077 (0.7919 - 0.8227)
% G-measure (Fowlkes–Mallows index): 0.8123 (0.7966 - 0.8271)
% Matthews index: 0.7151 - Normalized Index 0.8575 (0.8433 - 0.8707)
%  
% Threat score (TS) or critical success index (CSI): 0.6775 (0.6590 - 0.6954)
% Gilbert Skill Score: 0.5451 (0.5256 - 0.5644)
% True Skill Statistic (Hanssen-Kuipper skill score; Pierce's skill score): 0.7510 (0.7338 - 0.7675)
% Heidke Skill Score (Cohen's Kappa): 0.7056 (0.6768 - 0.7344)
%  
% Test bias: 1.2373 (0.9474 - 1.6160) - There is not test bias
%  
% Error odds ratio: 1.6869 (1.2916 - 2.2032) - Test wrongly classifies in positives
% Bayesian Credibility Assessment
% Critical Error Odds Ratio (EOR): 1.1718
% Error Odds Ratio (EOR)>COR. Test is credible at the 95%
%  
% Diagnostic odds ratio: 52.0655 (39.8649 - 68.0002) - Test discriminates between diseased and not diseased
% Normalized Diagnostic odds ratio (Yule's coefficient): 0.9623 (0.9511 - 0.9710)
% Tetrachoric Coefficient: 0.9278 (0.9105 - 0.9450)
% Cramer's V: 0.7141 - Strong positive association (risk factor)
%  
% Bayesian Credibility Assessment
% Critical Diagostic Odds Ratio (DOR): 1.0182
% Diagostic Odds Ratio (DOR)>COR. Test is credible at the 95%
% 
% Discriminant Power: 2.1791
% 
%  
%    
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
cs=sum(x); %columns sum
rs=sum(x,2); %rows sums
N=sum(x(:)); %numbers of elements

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
    disp('DIAGNOSTIC TEST PERFORMANCE PARAMETERS')
    disp(repmat('-',1,100))
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
    disp('DIAGNOSTIC TEST PERFORMANCE PARAMETERS')
    disp(repmat('-',1,100))
    if pr<1e-4
        fprintf('Estimed Prevalence: %0.4e (%0.4e - %0.4e)\n',pr,prci(1),prci(2))
    else
        fprintf('Estimed Prevalence: %0.4f (%0.4f - %0.4f)\n',pr,prci(1),prci(2))
    end
    clear prci
end
disp(' ')
clear ButtonName2

IR=cs(1)/cs(2); if IR<1; IR=1/IR; end
fprintf('Imbalance ratio (IR): %0.4f\n',IR)
zval=(cs(1)/N-0.5)/sqrt(0.25/N);
pvalue=2 * normcdf(-abs(zval),0,1);
if pvalue<0.05
    fprintf('Z = %0.4f p-value = %0.4f - Groups are unbalanced\n',zval,pvalue)
else
     fprintf('Z = %0.4f p-value = %0.4f - Groups are balanced\n',zval,pvalue)
end
clear zval pvalue
disp(' ')

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

fprintf('Sensitivity, Recall, Hit rate or True Positive rate (TPR): %0.4f (%0.4f - %0.4f)\n',SS(1),Seci(1),Seci(2))
fprintf('Miss rate or False Negative rate (FNR): %0.4f (%0.4f - %0.4f)\n',fp(1),fnci(1),fnci(2))
disp(' ')
fprintf('Specificity, Selectivity or True Negative Rate (TNR): %0.4f (%0.4f - %0.4f)\n',SS(2),Spci(1),Spci(2))
fprintf('Fall-out or False Positive rate (FPR): %0.4f (%0.4f - %0.4f)\n',fp(2),fpci(1),fpci(2))
disp(' ')
clear Seci Spci fnci fpci d

%Youden's Index
%Youden's J statistics (also called Youden's index) is a single statistic that
%captures the performance of a diagnostic test. The use of such a single index
%is "not generally to be recommended". It is equal to the risk difference for a
%dichotomous test and it defined as: J = Sensitivity + Specificity - 1. 
%A perfect test has J=1. 
J=sum(SS)-1; %Youden's index
fprintf('Informedness - Youden''s Index - J: %0.4f\n', J)

%The inverse of J has been defined as the “number needed to diagnose”
%(NND), that is, the number of patients who need to be examined in order
%to correctly detect one person with the disease of interest in a study
%population of persons with and without the known disease    
NND=1/J; %Number needed to Diagnose (NND)
fprintf('Number Needed to Diagnose (NND): %0.2f\n',NND);
fprintf('Around %i patients need to be tested to correctly detect 100 positive tests for the presence of disease\n',ceil(NND*100)) 
disp(' ')
clear NND

%Below a disease prevalence of PT, the PPV for a screening test with
%these sensitivities and specificities drops significantly and is therefore
%more unreliable.
PT=(sqrt(SS(1)*fp(2))+SS(2)-1)/J;
fprintf('Prevalence Threshold (PT): %0.4f\n', PT)
disp(' ')
clear PT

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
fprintf('Positive Likelihood Ratio (PLR): %0.4f (%0.4f - %0.4f)\n',plr,plrci(1),plrci(2))
dlr(plr)
disp(' ')
if nlr<1e-4 || nlrci(1)<1e-4
    fprintf('Negative Likelihood Ratio (NLR): %0.4e (%0.4e - %0.4e)\n',nlr,nlrci(1),nlrci(2))
else
    fprintf('Negative Likelihood Ratio (NLR): %0.4f (%0.4f - %0.4f)\n',nlr,nlrci(1),nlrci(2))
end
dlr(nlr)
disp(' ')
clear plrci nlrci cv

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
fprintf('Precision or Positive Predictive Value (PPV): %0.4f (%0.4f - %0.4f)\n', PNp(1),Ppci(1),Ppci(2))
fprintf('False discovery rate (FDR): %0.4f (%0.4f - %0.4f)\n',fd(1),fdci(1),fdci(2))
disp(' ')
fprintf('Negative Predictive Value (NPV): %0.4f (%0.4f - %0.4f)\n', PNp(2),Npci(1),Npci(2))
fprintf('False omission rate (FOR): %0.4f (%0.4f - %0.4f)\n',fd(2),foci(1),foci(2))
disp(' ')
clear fd po pr Ppci Npci fdci foci

%Predictive Summary Index (similar to Youden's Index)
PSI=sum(PNp)-1;
NNP=1/PSI; %Number needed to screen
fprintf('Markedness - Predictive Summary Index (PSI): %0.4f\n', PSI)
fprintf('Number Needed to Predict (NNP): %0.2f\n',NNP);
fprintf('Around %i patients need to be screened to correctly predict 100 diagnosis\n',ceil(NNP*100)) 
disp(' ')
clear NNP

%Accuracy and Mis-classification rate
%Diagnostic accuracy (or Power) is defined as the proportion of all tests
%that give a correct result. The Mis-classification rate is its complement to 1. 
acc=trace(x)/N; %Accuracy
accci=newcombe(acc);
mcr=1-acc; %Mis-classification rate
mcrci=newcombe(mcr);
fprintf('Accuracy or Potency (ACC): %0.4f (%0.4f - %0.4f)\n',acc,accci(1),accci(2))
II=2*acc-1; NNS=1/II;
fprintf('Identification Index (II): %0.4f\n', II)
fprintf('Number Needed to Screen (NNS): %0.2f\n',NNS);
fprintf('Around %i patients need to be screened to prevent 100 deaths by disease.\n',ceil(NNS*100)) 
disp(' ')
fprintf('Mis-classification Rate (MCR): %0.4f (%0.4f - %0.4f)\n',mcr,mcrci(1),mcrci(2))
fprintf('Number Needed to Mis-diagnose (NNMD): %0.2f\n',1/mcr);
fprintf('Around %i patients need to be screened to have 100 mis-diagnosis\n',ceil(100/mcr)) 
disp(' ')
clear acc accci mcr mcrci II NNS

bacc=mean(SS); %Balanced Accuracy
baccci=newcombe(bacc);
bmcr=1-bacc; %Balanced Mis-classification rate
bmcrci=newcombe(bmcr);
fprintf('Balanced Accuracy or Potency (BA): %0.4f (%0.4f - %0.4f)\n',bacc,baccci(1),baccci(2))
bII=2*bacc-1; bNNS=1/bII;
fprintf('Balanced Identification Index (BII): %0.4f\n', bII)
fprintf('Balanced Number Needed to Screen (NNS): %0.2f\n',bNNS);
fprintf('Around %i patients need to be screened to prevent 100 deaths by disease.\n',ceil(bNNS*100)) 
disp(' ')
fprintf('Balanced Mis-classification Rate (BMCR): %0.4f (%0.4f - %0.4f)\n',bmcr,bmcrci(1),bmcrci(2))
fprintf('Number Needed to Mis-diagnose (BNNMD): %0.2f\n',1/bmcr);
fprintf('Around %i patients need to be screened to have 100 mis-diagnosis\n',ceil(100/bmcr)) 
disp(' ')
clear bacc baccci bmcr bmcrci bII bNNS

%In statistics, the F1 score (also F-score or F-measure or Sørensen–Dice coefficient) is a measure of a
%test's accuracy. It considers both the Precision (positive predictivity) 
%and the Sensitivity of the test to compute the score.
%The F1 score can be interpreted as the harmonic mean of the Precision and
%Sensitivity, where an F1 score reaches its best value at 1 and worst score at 0. 
%Similarly, G-Measure (Fowlkes–Mallows index) is the geometric mean of Precision and Sensitivity
FMeasure=harmmean([SS(1) PNp(1)]); %F-measure
FMci=newcombe(FMeasure);
GMeasure=geomean([SS(1) PNp(1)]); %G-measure
GMci=newcombe(GMeasure);
if strcmp(ButtonName,'Yes')
    MMeasure=geomean([J PSI]); %Matthews
else
    MMeasure=det(x)/sqrt(prod(rs)*prod(cs));
end
nMMeasure=(MMeasure+1)/2; %normalized Matthews
MMci=newcombe(nMMeasure);

fprintf('F1-measure (Sørensen–Dice index): %0.4f (%0.4f - %0.4f)\n',FMeasure,FMci(1),FMci(2))
fprintf('G-measure (Fowlkes–Mallows index): %0.4f (%0.4f - %0.4f)\n',GMeasure,GMci(1),GMci(2))
fprintf('Matthews index: %0.4f - Normalized Index %0.4f (%0.4f - %0.4f)\n',MMeasure,nMMeasure,MMci(1),MMci(2))
disp(' ')
clear FMeasure FMci GMeasure GMci MMeasure nMMeasure MMci J PSI PNp
% Also called the threat score (TS), is a verification measure of categorical forecast performance
% equal to the total number of correct event forecasts (hits) divided by the total number of storm
% forecasts plus the number of misses (hits + false alarms + misses). The CSI is not affected by the
% number of non-event forecasts that verify (correct rejections). However, the CSI is a biased
% score that is dependent upon the frequency of the event. For an unbiased version of the CSI, see
% the Gilbert skill score (GS). With respect to the 2x2 verification problem example outlined in the
% definition of contingency table, CSI= (TP)/(TP+FN+FP).
TS=x(1)/(N-x(4));
TSci=newcombe(TS);
ar=rs(1)*cs(1)/N;
GS=(x(1)-ar)/(N-x(4)-ar); 
GSci=newcombe(GS);
fprintf('Threat score (TS) or critical success index (CSI): %0.4f (%0.4f - %0.4f)\n', TS,TSci(1),TSci(2))
fprintf('Gilbert Skill Score: %0.4f (%0.4f - %0.4f)\n', GS,GSci(1),GSci(2))
clear TS TSci ar GS GSci
%true skill statistic
%popular combination of sensitivity and false positive rate; measures the
%ability to separate yes and no cases; range from -1 to 1; 
%perfect score=1; no skill level=0;
dx=det(x); TSS=dx/prod(cs); TSSci=newcombe(TSS);
fprintf('True Skill Statistic (Hanssen-Kuipper skill score; Pierce''s skill score): %0.4f (%0.4f - %0.4f)\n',TSS,TSSci(1),TSSci(2))
clear dx TSS TSSci

%Heidke skill score
% It is another form of Cohen's kappa
%measures fractional improvements over random chance; range from -inf to 1;
%perfect score=1; no skill level=0;
f=diag(ones(1,2));
po=sum(sum(x./N.*f));
pe=sum(sum(rs*cs.*f./N^2));
HSS=(po-pe)/(1-pe); %Cohen's kappa
sek=sqrt((po*(1-po))/(N*(1-pe)^2)); %kappa standard error for confidence interval
HSSci=HSS+([-1 1].*(z*sek)); %k confidence interval
fprintf('Heidke Skill Score (Cohen''s Kappa): %0.4f (%0.4f - %0.4f)\n',HSS,HSSci(1),HSSci(2))
disp(' ')
clear f po pe HSS sek HSSci

orse=realsqrt(sum(1./x(:))); %standard error of log(OR)
cv=([-1 1].*(z*orse)); 

%Test Bias (TB)
%A test which shows provable and systematic differences in the results of people
%based on group membership. For example, a test might be considered biased if
%members of one particular gender or race consistently and systematic have
%statistically different results from the rest of the testing population. 
%It is defined as (T+)/(D+)=(TP+FP)/(TP+FN)
%A perfect test has a TB=1;
%If TB<1 the test underestimates the disease because there are more affected than positive test
%If TB>1 the test overestimates the disease because there are more positive test than affected
TB=rs(1)/cs(1); %Test Bias
orci=exp(reallog(TB)+cv); %OR confidence interval
fprintf('Test bias: %0.4f (%0.4f - %0.4f) - ',TB,orci(1),orci(2))
if TB>1
    if orci(1)<1 && orci(2)>1
        disp('There is not test bias')
        flag=0;
    else
        disp('Test overestimates the phenomenon')
        flag=1;
    end
elseif TB<1
    if orci(1)<1 && orci(2)>1
        disp('There is not test bias')
        flag=0;
    else
        disp('Test underestimates the phenomenon')
        flag=1;
    end
else
    disp('Test is perfect')
    flag=0;
end
if flag==1
    bca(TB,orci,'Test Bias (TB)')
end
disp(' ')
clear TB orci flag

%Error Odds Ratio. 
%Indicates if the probability of being wrongly classified is highest in the
%diseased or in the non-diseased group. If the error odds is higher than one the
%probability is highest in the diseased group (and the specificity of the test
%is better than the sensitivity), if the value is lower than one the probability
%of an incorrect classification is highest in the non-diseased group (and the
%sensitivity of the test is better than the specificity).     
%It is defined as (Sensitivity/(1-Sensitivity))/(Specificity/(1-Specificity));
EOR=(SS(1)/fp(1))/(SS(2)/fp(2)); %Error odds ratio
orci=exp(reallog(EOR)+cv); %OR confidence interval
fprintf('Error odds ratio: %0.4f (%0.4f - %0.4f) - ',EOR,orci(1),orci(2))
if EOR>1
    if orci(1)<1 && orci(2)>1
        disp('There is not misclassification')
        flag=0;
    else
        disp('Test wrongly classifies in positives')
        flag=1;
    end
elseif EOR<1
    if orci(1)<1 && orci(2)>1
        disp('There is not misclassification')
        flag=0;
    else
        disp('Test wrongly classifies in negatives')
        flag=1;
    end
else
    disp('Test is perfect')
    flag=0;
end
if flag==1
    bca(EOR,orci,'Error Odds Ratio (EOR)')
end
clear EOR orci flag

%Diagnostic Odds Ratio. 
%Diagnostic odds ratio is defined as how much more likely will the test
%make a correct diagnosis than an incorrect diagnosis in patients with the
%disease (Scott et al. 2008).  
%Often used as a measure of the discriminative power of the test. Has the value
%one if the test does not discriminate between diseased and not diseased. Very
%high values above one means that a test discriminates well. Values lower than
%one mean that there is something wrong in the application of the test.   
%DOR can be normalized into Yule's coefficient.
%If Yule's coefficient=0.5 the test works by chance. 
%This means that a DOR<3 the test in not effective into discriminate healthy from not healthy.

DOR=plr/nlr; %Diagnostic odds ratio
orci=exp(reallog(DOR)+cv); %OR confidence interval
fprintf('Diagnostic odds ratio: %0.4f (%0.4f - %0.4f) - ',DOR,orci(1),orci(2))
flag=0;
if DOR>1
    if DOR<=3 || (orci(1)<1 && orci(2)>1)
        disp('Test does not discriminate between diseased and not diseased')
    else
        disp('Test discriminates between diseased and not diseased')
        flag=1;
    end
elseif DOR<1
    if orci(1)<1 && orci(2)>1
        disp('Test does not discriminate between diseased and not diseased')
    else
        disp('There is something wrong in the application of the test')
    end
else
    disp('Test does not discriminate between diseased and not diseased')
end
clear plr nlr 
fprintf('Normalized Diagnostic odds ratio (Yule''s coefficient): %0.4f (%0.4f - %0.4f)\n',(DOR-1)/(DOR+1),(orci(1)-1)/(orci(1)+1),(orci(2)-1)/(orci(2)+1))
TC=cos(pi/(1+sqrt(DOR)));
ac=DOR^(pi/4); 
SE=realsqrt((pi*ac/(2*(1+ac)^2))^2)*orse;
tcci=TC+[-1 +1].*(SE*z);
fprintf('Tetrachoric Coefficient: %0.4f (%0.4f - %0.4f)\n',TC,tcci(1),tcci(2))
clear TC ac SE orse tcci
Phi=(det(x)-N/2)/realsqrt(prod(cs)*prod(rs));
phi_hat=max(0,Phi^2-1/(N-1));
k_hat=2-1/(N-1);
V=sqrt(phi_hat/(k_hat-1));
clear Phi phi_hat k_hat
fprintf('Cramer''s V: %0.4f - ',V)
switch sign(V)
    case -1
        txt2='negative association (protective factor)';
    case 1
        txt2='positive association (risk factor)';
end
V=abs(V);
if V<=0.3
    txt1='Weak ';
elseif (V>0.3 && V<=0.7)
    txt1='Moderate ';
else
    txt1='Strong ';
end
disp([txt1 txt2])
disp(' ')
clear V txt1 txt2
if flag==1
    bca(DOR,orci,'Diagostic Odds Ratio (DOR)')
end

%Discriminant power
%The discriminant power for a test, also termed the test effectiveness, is a
%measure of how well a test distinguishes between affected and unaffected
%persons. It is the sum of logs of Sensivity and Specificity over own false
%proportion, scaled by the standard deviation of the logistic normal
%distribution curve (square root of 3 divided by π). Test effectiveness is
%interpreted as the standardized distance between the means for both populations.     
dpwr=(realsqrt(3)/pi)*log(DOR); %Discriminant power
fprintf('Discriminant Power: %0.4f\n',dpwr)
clear dpwr
disp(' ')
clear alpha DOR orci flag cv z

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

function dlr(lr) %Likelihood dialog
if lr==1
    disp('Test is not suggestive of the presence/absence of disease')
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

p2=' increase in possibility of disease ';

if lr>1
    p3='presence';
elseif lr<1
    p3='absence';
end
disp([p1 p2 p3])
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

function bca(or,orci,txt)
disp('Bayesian Credibility Assessment')
lorci=reallog(orci);
cor=exp(-diff(lorci)^2/(4*realsqrt(prod(lorci)))); %Critical odds ratio (COR)
if or<1
    fprintf('Critical %s: %0.4f\n',txt,cor)
    if or<cor
        fprintf('%s<COR. Test is credible at the %d%%\n',txt,(1-alpha)*100)
    else
        fprintf('%s>=COR. Test isn''t credible at the %d%%\n',txt,(1-alpha)*100)
    end
else
    cor=1/cor; %correct cor
    fprintf('Critical %s: %0.4f\n',txt,cor)
    if or>cor
        fprintf('%s>COR. Test is credible at the %d%%\n',txt,(1-alpha)*100)
    else
        fprintf('%s<=COR. Test isn''t credible at the %d%%\n',txt,(1-alpha)*100)
    end
end
disp(' ')
end

end
