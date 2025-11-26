function partest(x,varargin)
%PARTEST  Performance analysis of a binary diagnostic test (2x2 confusion matrix).
%
%   PARTEST(X) computes a wide set of diagnostic, predictive, similarity,
%   accuracy, association and entropy-based indices from a 2x2 confusion
%   matrix X:
%
%       X = [ TP  FN
%             FP  TN ]
%
%   where:
%       TP = True Positive
%       FN = False Negative
%       FP = False Positive
%       TN = True Negative
%
%   PARTEST(X,ALPHA) uses the specified significance level ALPHA (0<ALPHA<1)
%   for the computation of (1-ALPHA)*100% confidence intervals (default 0.05
%   for 95% confidence intervals).
%
%   PARTEST(X,ALPHA,PR) additionally uses the specified disease prevalence PR
%   (0<=PR<=1). If PR is not provided (or is NaN), the prevalence is estimated
%   from the confusion matrix as P(D+)=P/(P+N).
%
%   The function:
%     - prints the confusion matrix and basic counts;
%     - performs an imbalance test;
%     - computes Sensitivity, Specificity, their complements, AUC and CI,
%       distance and Gini indices;
%     - computes Youden's J, Number Needed to Diagnose, predictive values,
%       likelihood ratios and several related measures;
%     - computes similarity indices (Jaccard, Otsuka-Ochiai, etc.);
%     - computes accuracy-related indices (F1, G-measure, adjusted F, TSS,
%       Matthews index, HSS, etc.);
%     - computes entropy-based indices (CEN and modified CEN);
%     - computes association measures (Test Bias, Error OR, DOR, Yule,
%       tetrachoric coefficient, Cramer's V, discriminant power);
%     - shows two diagnostic plots.
%
%   INPUT:
%       X      - 2x2 matrix of nonnegative integers (TP,FN;FP,TN).
%       ALPHA  - optional scalar in (0,1), default 0.05.
%       PR     - optional scalar in [0,1]; if omitted or NaN, prevalence is
%                estimated from the confusion matrix.
%
%   OUTPUT:
%       No variables are returned; results are displayed in the Command
%       Window and two figures are produced.
%
%   Author:  Giuseppe Cardillo
%   Email:   giuseppe.cardillo.75@gmail.com
%   GitHub:  https://github.com/dnafinder/partest
%   General: https://github.com/dnafinder
%
%   Version: 2.0.1
%   Last update: 2025-11-26
%
%   This code is distributed under a BSD-style license.
%

% Input Error handling
narginchk(1,3);

p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},...
    {'real','finite','integer','nonnegative','nonnan','size',[2 2]}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},...
    {'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'pr',NaN, @(x) validateattributes(x,{'numeric'},...
    {'scalar','real','finite','nonnan','>=',0,'<=',1}));
parse(p,x,varargin{:});
alpha = p.Results.alpha;
pr    = p.Results.pr;
userProvidedPr = ~isnan(pr);
clear p varargin

% Critical value
z  = -realsqrt(2)*erfcinv(2-alpha); 
tr = repmat('-',1,100);

% Basic counts
cs = sum(x);      % columns sum
rs = sum(x,2);    % rows sums
N  = sum(x(:));   % numbers of elements
dx = det(x);      % matrix determinant

clc
disp('CONFUSION MATRIX')
disp(array2table(x,...
    'RowNames',{'Positive_test';'Negative_test'},...
    'VariableNames',{'Affected';'Healthy'}))
disp(tr)

disp('BASIC PARAMETERS')
disp(array2table([x(:);cs(:);rs(:);N],...
    'RowNames',{'True_Positive_TP';'False_Negative_FN';...
                'False_Positive_FP';'True Negative_TN';...
                'Condition_Positive_P';'Condition_Negative_N';...
                'Test_Outcome_Positive_TOP';'Test_Outcome_Negative_TON';...
                'Population_POP'},...
    'VariableNames',{'Parameter'}))
disp(tr)

%% IMBALANCE TEST
disp('IMBALANCE TEST')
disp(' ')
IR = cs(1)/cs(2);
if IR<1
    IR = 1/IR;
end
fprintf('Imbalance ratio (IR): %0.4f\n',IR)

zval   = diff(cs)/sqrt(N);
pvalue = 2 * normcdf(-abs(zval),0,1);
if pvalue<0.05
    fprintf('Z = %0.4f p-value = %0.4f - Groups are unbalanced\n',zval,pvalue)
else
    fprintf('Z = %0.4f p-value = %0.4f - Groups are balanced\n',zval,pvalue)
end
clear zval pvalue
disp(tr)

%% PERFORMANCE PARAMETERS
disp('PERFORMANCE PARAMETERS')
d = diag(x); % true positives and true negatives

% Sensitivity and Specificity
SS = d./cs';  % [Sensitivity; Specificity]
% False proportions
fp = 1-SS;

% 95% confidence intervals
Seci = newcombe(SS(1)); 
Spci = newcombe(SS(2));
fnci = newcombe(fp(1)); 
fpci = newcombe(fp(2));

% Area under the ROC curve
AUC = mean(SS);

% Hanley and McNeil formula
q0 = AUC*(1-AUC);
q1 = 2*AUC^2/(1+AUC)-AUC^2;
q2 = AUC/(2-AUC)-AUC^2;
SE = realsqrt([1 (cs-1)]*[q0;q1;q2]/prod(cs));
AUCci = AUC+[-1 1].*z*SE;
clear q0 q1 q2 SE

if AUC==1
    str = 'Perfect test';
elseif AUC>=0.90 && AUC<1
    str = 'Excellent test';
elseif AUC>=0.80 && AUC<0.90
    str = 'Good test';
elseif AUC>=0.70 && AUC<0.80
    str = 'Fair test';
elseif AUC>=0.60 && AUC<0.70
    str = 'Poor test';
elseif AUC>=0.50 && AUC<0.60
    str = 'Fail test';
else
    str = 'Failed test - less than chance';
end

% Distance Index (Euclidean distance from top-left corner in ROC space)
dind = realsqrt(sum(fp.^2));

% Gini Index
GI = 2*AUC-1;

Parameter = round([SS(1);fp(1);SS(2);fp(2);AUC;dind;GI],4);
LB        = round([Seci(1);fnci(1);Spci(1);fpci(1);AUCci(1);NaN;2*AUCci(1)-1],4);
UB        = round([Seci(2);fnci(2);Spci(2);fpci(2);AUCci(2);NaN;2*AUCci(2)-1],4);
Comment   = {'None';'None';'None';'None';str;'None';'None'};

disp(table(Parameter,LB,UB,Comment,...
    'RowNames',{'Sensitivity_TPR';'False_Negative_Rate_FNR';...
                'Specificity_TNR';'False_Positive_Rate_FPR';...
                'Area_Under_the_Curve_AUC';...
                'Distance_Index_dInd';'Gini_Index_GI'}))

clear Spci fnci fpci d AUC AUCci GI str Parameter LB UB Comment

%% AM (Automatic/Manual)
AM = rs(1)-cs(1);

% Youden's Index (Informedness)
J   = sum(SS)-1;      % Youden's index
NND = 1/J;            % Number needed to Diagnose

disp(array2table(round([AM;J;NND],4),...
    'RowNames',{'Automatic_Manual';'Youden_Index_J';...
                'Number_needed_to_diagnose_NND'},...
    'VariableNames',{'Parameter'}))
fprintf(['Around %i patients need to be tested to correctly detect 100 ',...
         'positive tests for the presence of disease\n'],ceil(NND*100))
disp(' ')

clear NND

disp(tr)
disp('PREDICTIVE PARAMETERS')
disp(' ')

%% Prevalence from argument or estimated

if userProvidedPr
    % Prevalence provided by the user
    if pr<1e-4
        fprintf('Reported Prevalence: %0.4e\n',pr)
    else
        fprintf('Reported Prevalence: %0.4f\n',pr)
    end
else
    % Prevalence estimated from confusion matrix
    pr   = cs(1)/N;
    prci = newcombe(pr);
    if pr<1e-4
        fprintf('Estimed Prevalence: %0.4e (%0.4e - %0.4e)\n',pr,prci(1),prci(2))
    else
        fprintf('Estimed Prevalence: %0.4f (%0.4f - %0.4f)\n',pr,prci(1),prci(2))
    end
    clear prci
end
disp(' ')

%% Likelihood ratios, predictive values and related measures

% Positive and Negative Likelihood Ratio
plr = SS(1)/fp(2); % Positive Likelihood Ratio
nlr = fp(1)/SS(2); % Negative Likelihood Ratio

cv    = realsqrt(sum(1./diag(x))-sum(1./sum(x,2)))*z;
plrci = exp(log(plr)+[-1 1].*cv);
nlrci = exp(log(nlr)+[-1 1].*cv);
plrtxt = dlr(plr);
nlrtxt = dlr(nlr);

% Positive and Negative predictivity (Bayes theorem)
po   = pr/(1-pr); % prior odds
PNp  = [1/(1+1/(plr*po))  1/(1+nlr*po)]; % PPV, NPV
fd   = 1-PNp;  % false discovery and omission rate
Ppci = newcombe(PNp(1));
Npci = newcombe(PNp(2));
fdci = newcombe(fd(1)); 
foci = newcombe(fd(2));

Parameter = round([plr;nlr;PNp(1);fd(1);PNp(2);fd(2)],4);
LB        = round([plrci(1);nlrci(1);Ppci(1);fdci(1);Npci(1);foci(1)],4);
UB        = round([plrci(2);nlrci(2);Ppci(2);fdci(2);Npci(2);foci(2)],4);
Comment   = {plrtxt;nlrtxt;'None';'None';'None';'None'};

disp(table(Parameter,LB,UB,Comment,...
    'RowNames',{'Positive_Likelihood_Ratio_PLR';...
                'Negative_Likelihood_Ratio_NLR';...
                'Positive_Prediction_Value_PPV';...
                'False_Discovery_Rate_FDR';...
                'Negative_Predictive_Value_NPV';...
                'False_Omission_Rate_FOR'}))

% Prevalence threshold
PT = (sqrt(SS(1)*fp(2))+SS(2)-1)/J;

% LIFT SCORE
LS = PNp(1)/pr;

% Information Score
IS = -log2(cs(1)/N)+log2(PNp(1));

% Markedness
MK  = sum(PNp)-1;
NNP = 1/MK; % Number needed to predict

disp(array2table(round([PT;LS;IS;MK;NNP],4),...
    'RowNames',{'Prevalence_Threshold_PT';'Lift_Score_LS';...
                'Information_Score_IS';'Markedness_MK';...
                'Number_needed_to_predict_NNP'},...
    'VariableNames',{'Parameter'}))
fprintf('Around %i patients need to be screened to correctly predict 100 diagnosis\n',...
    ceil(NNP*100))
disp(' ')

clear plrci nlrci cv fd po Npci fdci foci NNP LS IS

disp(tr)
disp('SIMILARITY PARAMETERS')

%% Similarity / dissimilarity indices

% Bray-Curtis dissimilarity
BCD   = abs(AM)/(2*N); 
BCDci = newcombe(BCD);

% Similarity index
sind   = 1-dind/sqrt(2); 
sindci = newcombe(sind);

% Jaccard Index
JI   = x(1)/(rs(1)+cs(1)-x(1)); 
JIci = newcombe(JI);

% Overlap coefficient (Szymkiewicz–Simpson)
[OC,I] = max([PNp(1) SS(1)]);
if I==1
    OCci = Ppci;
else
    OCci = Seci;
end

% Braun-Blanquet coefficient
[BB,I] = min([PNp(1) SS(1)]);
if I==1
    BBci = Ppci;
else
    BBci = Seci;
end

% Otsuka-Ochiai coefficient
OOCden = realsqrt(rs(1)*cs(1));
OOC    = x(1)/OOCden; 
OOCci  = newcombe(OOC);

Parameter = round([BCD;sind;JI;OC;BB;OOC],4);
LB        = round([BCDci(1);sindci(1);JIci(1);OCci(1);BBci(1);OOCci(1)],4);
UB        = round([BCDci(2);sindci(2);JIci(2);OCci(2);BBci(2);OOCci(2)],4);

disp(table(Parameter,LB,UB,...
    'RowNames',{'Bray_Curtis_Dissimilarity_BCD';'Similarity_Index_sInd';...
                'Jaccard_Index_JI';'Overlap_Coefficient_OC';...
                'Braun_Blanquet_similarity_BB';...
                'Otsuka_Ochiai_Coefficient_OOC'}))

clear AM BCD* dind* sind* JI* OC* BB* Ppci Seci OOCden

disp(tr)
disp('ACCURACY PARAMETERS')

%% Random accuracy measures

RACC   = rs(1)*cs(1)/N^2; 
RACCci = newcombe(RACC);

RACCU   = ((rs(1)+cs(1))/(2*N))^2; 
RACCUci = newcombe(RACCU);

Parameter = round([RACC;RACCU],4);
LB        = round([RACCci(1);RACCUci(1)],4);
UB        = round([RACCci(2);RACCUci(2)],4);
clear RACC*

disp(table(Parameter,LB,UB,...
    'RowNames',{'Random_Accuracy_RACC';...
                'Random_Accuracy_Unbiased_RACCU'}))

%% Accuracy, F-measures and skill scores

acc = trace(x)/N;       % accuracy
acc = [acc;mean(SS)];   % balanced accuracy
acc = [acc;acc(1)-abs(diff(SS))/abs(sum(SS))]; % optimized precision

% F1 Measure (harmonic mean of Sensitivity and Precision)
if (SS(1)+PNp(1))>0
    F1 = 2*SS(1)*PNp(1)/(SS(1)+PNp(1));
else
    F1 = 0;
end
acc = [acc;F1];

% Adjusted F-score (AGF metric): geometric mean of two components
AF1 = 5*PNp(1)*SS(1)/(4*PNp(1)+SS(1));
AF2 = (1+0.5^2)*PNp(2)*SS(2)/(0.5^2*PNp(2)+SS(2));
AGF = realsqrt(AF1*AF2);
acc = [acc;AGF];

% G-measure (geometric mean of Sensitivity and Precision)
Gm = realsqrt(SS(1)*PNp(1));
acc = [acc;Gm];

% Adjusted G-measure
Nn  = cs(2)/N; 
acc = [acc;(acc(end)+SS(2)*Nn)/(1+Nn)];
clear Nn AF1 AF2 AGF Gm F1

% Matthews Index (MCC) and normalized Matthews
if userProvidedPr
    MMeasure = realsqrt(J*MK); % MK e J sono entrambe probabilità / indici [0,1]
else
    MMeasure = dx/sqrt(prod(rs)*prod(cs));
end
acc = [acc;(MMeasure+1)/2]; % normalized Matthews

clear MMeasure J MK PNp

% Gilbert skill score
ar  = rs(1)*cs(1)/N;
acc = [acc;(x(1)-ar)/(N-x(4)-ar)];
clear ar

% True skill statistic
acc = [acc;dx/prod(cs)];
clear dx

% Heidke skill score (Cohen's kappa)
f  = diag(ones(1,2));
po = sum(sum(x./N.*f));
pe = sum(sum(rs*cs.*f./N^2));
HSS = (po-pe)/(1-pe);
acc = [acc;HSS];
sek = sqrt((po*(1-po))/(N*(1-pe)^2));
clear f po pe 

acctxt = {'Accuracy_ACC';'Balanced_Accuracy_bACC';'Optimized_Precision_OP';...
          'F1_Measure';'Adjusted_F_Score';'G_Measure';'Adjusted_G_Measure';...
          'Normalized_Matthews_Index_nMCI';'Gilbert_Skill_Score';...
          'True_Skill_Statistic';'Heidke_Skill_Score'};

L = length(acc);
for I = 1:L
    if I < L
        accci = newcombe(acc(I));
    else
        accci = HSS+([-1 1].*(z*sek));
        clear HSS sek
    end
    mcr   = 1-acc(I); 
    mcrci = newcombe(mcr);
    II    = 2*acc(I)-1; 
    IIci  = 2.*accci-1;
    NNS   = 1/II; 
    NNSci = fliplr(1./IIci);
    NNMD  = 1/mcr; 
    NNMDci= fliplr(1./mcrci);

    if     acc(I)<0.3
        txt2 = ' negligible';
    elseif all([acc(I)>=0.3 acc(I)<0.5])
        txt2 = ' weak';
    elseif all([acc(I)>=0.5 acc(I)<0.7])
        txt2 = ' moderate';
    elseif all([acc(I)>=0.7 acc(I)<0.9])
        txt2 = ' strong';
    else
        txt2 = ' very strong';
    end
    txt = strcat('Accuracy is',txt2);

    Parameter = round([acc(I);mcr;II;NNS;NNMD],4);
    LB        = round([accci(1);mcrci(1);IIci(1);NNSci(1);NNMDci(1)],4);
    UB        = round([accci(2);mcrci(2);IIci(2);NNSci(2);NNMDci(2)],4);
    Comment   = {txt;'None';'None';...
        sprintf(['Around %i patients need to be screened to prevent 100 ',...
                 'deaths by disease.'],ceil(NNS*100));...
        sprintf(['Around %i patients need to be screened to have 100 ',...
                 'mis-diagnosis'],ceil(NNMD*100))};

    disp(table(Parameter,LB,UB,Comment,...
        'RowNames',{acctxt{I};'Misclassification_Rate_MCR';...
                    'Identification_Index_II';...
                    'Number_needed_to_screen_NNS';...
                    'Number_needed_to_misdiagnose_NNMD';}))
end

clear txt* acc* mcr* II NNS NNMD L acc Parameter LB UB Comment

disp(tr)
disp('ENTROPY MEASURE')

%% Confusion Entropy (CEN and modified CEN)

AF   = x(2)+x(3); 
loga = x(2)*log2(x(2))+x(3)*log2(x(3));
trx  = trace(rot90(x));
acc  = zeros(2,1);
acc(1) = (AF*log2(N^2-diff(diag(x))^2)/2-loga)/N;
acc(2) = 2*(trx*log2((N-x(4))*(N-x(1)))-2*loga)/(3*N+trx);
clear AF loga trx

acctxt = {'Confusion_Matrix_Entropy_CEN';'Modified_CEN'};

L = length(acc);
for I = 1:L
    accci = newcombe(acc(I));
    mcr   = 1-acc(I); 
    mcrci = newcombe(mcr);
    if acc(I)>0.5
        II   = 2*acc(I)-1; 
        IIci = 2.*accci-1;
        NNS  = 1/II; 
        NNSci= fliplr(1./IIci);
    else
        II   = NaN; 
        IIci = NaN(1,2); 
        NNS  = NaN; 
        NNSci= NaN(1,2);
    end
    NNMD   = 1/mcr; 
    NNMDci = fliplr(1./mcrci);

    if     acc(I)<0.3
        txt2 = ' negligible';
    elseif acc(I)>=0.3 && acc(I)<0.5
        txt2 = ' weak';
    elseif acc(I)>=0.5 && acc(I)<0.7
        txt2 = ' moderate';
    elseif acc(I)>=0.7 && acc(I)<0.9
        txt2 = ' strong';
    else
        txt2 = ' very strong';
    end
    txt = strcat('Accuracy is',txt2);

    Parameter = round([acc(I);mcr;II;NNS;NNMD],4);
    LB        = round([accci(1);mcrci(1);IIci(1);NNSci(1);NNMDci(1)],4);
    UB        = round([accci(2);mcrci(2);IIci(2);NNSci(2);NNMDci(2)],4);
    if ~isnan(II)
        txt3 = sprintf(['Around %i patients need to be screened to ',...
                        'prevent 100 deaths by disease.'],ceil(NNS*100));
    else
        txt3 = 'Not computable';
    end
    Comment = {txt;'None';'None';txt3;...
        sprintf(['Around %i patients need to be screened to have 100 ',...
                 'mis-diagnosis'],ceil(NNMD*100))};

    disp(table(Parameter,LB,UB,Comment,...
        'RowNames',{acctxt{I};'Misclassification_Rate_MCR';...
                    'Identification_Index_II';...
                    'Number_needed_to_screen_NNS';...
                    'Number_needed_to_misdiagnose_NNMD';}))
end
clear txt* acc* mcr* II NNS NNMD L acc Parameter LB UB Comment

disp(tr)
disp('ASSOCIATION MEASURE')

%% Association measures

orse = realsqrt(sum(1./x(:))); % standard error of log(OR)
cv   = ([-1 1].*(z*orse)); 
Comment = cell(3,1);
txt     = cell(3,1);

% Test Bias (TB)
Parameter(1) = round(rs(1)/cs(1),4); % Test Bias
orci         = exp(reallog(Parameter(1))+cv);
LB(1)        = round(orci(1),4); 
UB(1)        = round(orci(2),4);
flag         = 0;
if Parameter(1)>1
    if orci(1)<1 && orci(2)>1
        Comment{1} = 'There is not test bias';
        flag = 0;
    else
        Comment{1} = 'Test overestimates the phenomenon';
        flag = 1;
    end
elseif Parameter(1)<1
    if orci(1)<1 && orci(2)>1
        Comment{1} = 'There is not test bias';
        flag = 0;
    else
        Comment{1} = 'Test underestimates the phenomenon';
        flag = 1;
    end
else
    Comment{1} = 'Test is perfect';
    flag = 0;
end
if flag==1
    [cor(1),txt{1}] = bca(Parameter(1),orci);
else
    cor(1) = NaN; 
    txt{1} = 'Not needed';
end

% Error Odds Ratio
Parameter(2) = (SS(1)/fp(1))/(SS(2)/fp(2)); 
orci         = exp(reallog(Parameter(2))+cv);
LB(2)        = round(orci(1),4); 
UB(2)        = round(orci(2),4);
flag         = 0;
if Parameter(2)>1
    if orci(1)<1 && orci(2)>1
        Comment{2} = 'There is not misclassification';
        flag = 0;
    else
        Comment{2} = 'Test wrongly classifies in positives';
        flag = 1;
    end
elseif Parameter(2)<1
    if orci(1)<1 && orci(2)>1
        Comment{2} = 'There is not misclassification';
        flag = 0;
    else
        Comment{2} = 'Test wrongly classifies in negatives';
        flag = 1;
    end
else
    Comment{2} = 'Test is perfect';
    flag = 0;
end

if flag==1
    [cor(2),txt{2}] = bca(Parameter(2),orci);
else
    cor(2) = NaN; 
    txt{2} = 'Not needed';
end

% Diagnostic Odds Ratio
DOR          = plr/nlr; 
Parameter(3) = DOR;
orci         = exp(reallog(DOR)+cv);
LB(3)        = round(orci(1),4); 
UB(3)        = round(orci(2),4);

flag = 0;
if DOR>1
    if DOR<=3 || (orci(1)<1 && orci(2)>1)
        Comment{3} = 'Test does not discriminate between diseased and not diseased';
    else
        Comment{3} = 'Test discriminates between diseased and not diseased';
        flag = 1;
    end
elseif DOR<1
    if orci(1)<1 && orci(2)>1
        Comment{3} = 'Test does not discriminate between diseased and not diseased';
    else
        Comment{3} = 'There is something wrong in the application of the test';
    end
else
    Comment{3} = 'Test does not discriminate between diseased and not diseased';
end
if flag==1
    [cor(3),txt{3}] = bca(Parameter(3),orci);
else
    cor(3) = NaN; 
    txt{3} = 'Not needed';
end
clear plr nlr 

% Yule coefficient (normalized DOR)
QI        = (DOR-1)/(DOR+1); 
Parameter(4) = QI; 
LB(4)    = (orci(1)-1)/(orci(1)+1); 
UB(4)    = (orci(2)-1)/(orci(2)+1); 
cor(4)   = NaN; 
txt{4}   = 'Not needed';
if QI<=0.25
    Comment{4} = 'Negligible'; 
elseif QI>0.25 && QI<=0.50
    Comment{4} = 'Weak';
elseif QI>0.5 && QI<=0.75
    Comment{4} = 'Moderate';
else
    Comment{4} = 'Strong';
end
clear QI

% Tetrachoric Coefficient
TC  = cos(pi/(1+sqrt(DOR)));
ac  = DOR^(pi/4); 
SE  = realsqrt((pi*ac/(2*(1+ac)^2))^2)*orse;
tcci = TC+[-1 +1].*(SE*z);
Parameter(5) = TC; 
LB(5)        = tcci(1); 
UB(5)        = tcci(2); 
Comment{5}   = 'None'; 
cor(5)       = NaN; 
txt{5}       = 'Not needed';
clear TC ac SE orse tcci

Phi     = (det(x)-N/2)/realsqrt(prod(cs)*prod(rs));
phi_hat = max(0,Phi^2-1/(N-1));
V       = realsqrt(phi_hat*(N-1)/(N-2)); % Cramer V
clear Phi phi_hat
Parameter(6) = V; 
LB(6)        = NaN; 
UB(6)        = NaN;
if V<=0.3
    Comment{6} = 'Weak association';
elseif (V>0.3 && V<=0.7)
    Comment{6} = 'Moderate association';
else
    Comment{6} = 'Strong association';
end
cor(6) = NaN; 
txt{6} = 'Not needed';

% Discriminant power
dpwr        = (realsqrt(3)/pi)*log(DOR); 
Parameter(7)= dpwr; 
LB(7)       = NaN; 
UB(7)       = NaN;
cor(7)      = NaN; 
txt{7}      = 'Not needed';
if dpwr<=1
    Comment{7} = 'Poor'; 
elseif dpwr>1 && dpwr<=2
    Comment{7} = 'Limited';
elseif dpwr>2 && dpwr<=3
    Comment{7} = 'Fair';
else
    Comment{7} = 'Good';
end

disp(table(Parameter',LB',UB',Comment,cor',txt,...
    'RowNames',{'Test_Bias_TB';'Error_OR';'Diagnostic_Odds_Ratio_DOR';...
                'Yule_Coefficient';'Tetrachoric_Coefficient';'Cramer_V';...
                'Discriminant_Power'},...
    'VariableNames',{'Parameter','LB','UB','Comment',...
                     'Critical_Bayesian_OR','Bayesian_Credibility'}))

clear flag Parameter LB UB Comment cor txt V dpwr alpha DOR orci cv ...
      userProvidedPr

%% Display graphs
xg = cs./N;

subplot(121)
hold on
H    = zeros(1,4);
H(1) = fill(xg(1)+[0 xg(2) xg(2) 0],SS(2)+[0 0 fp(2) fp(2)],'r');
H(2) = fill([0 xg(1) xg(1) 0],fp(1)+[0 0 SS(1) SS(1)],'g');
H(3) = fill([0 xg(1) xg(1) 0],[0 0 fp(1) fp(1)],'y');
H(4) = fill(xg(1)+[0 xg(2) xg(2) 0],[0 0 SS(2) SS(2)],'b');
hold off
axis square
title('PARTEST GRAPH')
xlabel('Subjects proportion')
ylabel('Parameters proportion')
legend(H,'False Positive','True Positive (Sensibility)',...
          'False Negative','True Negative (Specificity)',...
          'Location','NorthOutside')

subplot(122)
axis square
color = {'r','g','y','b'};
rp    = [fp(2) SS(1) fp(1) SS(2)];
k     = length(rp); 
H     = zeros(1,k);
ang   = (0:1:k)./k.*2.*pi;
axis equal
hold on
for I = 1:k
    theta = [ang(I) linspace(ang(I),ang(I+1),500) ang(I+1)];
    rho   = [0 repmat(realsqrt(rp(I)),1,500) 0];
    [xg,yg] = pol2cart(theta,rho);
    H(I) = patch(xg,yg,color{I});
end
hold off
title('ROSEPLOT PARTEST GRAPH')
legend(H,'False Positive','True Positive (Sensibility)',...
          'False Negative','True Negative (Specificity)',...
          'Location','NorthOutside')

%% Nested utilities

    function txt = dlr(lr) %#ok<NESTED>
        % DLR  Likelihood ratio qualitative interpretation.
        if lr==1
            txt = 'Test is not suggestive of the presence/absence of disease';
            return
        end

        if lr>10 || lr<0.1
            p1 = 'Large (often conclusive)';
        elseif (lr>5 && lr<=10) || (lr>0.1 && lr<=0.2)
            p1 = 'Moderate';
        elseif (lr>2 && lr<=5) || (lr>0.2 && lr<=0.5)
            p1 = 'Low';
        elseif (lr>1 && lr<=2) || (lr>0.5 && lr<=1)
            p1 = 'Poor';
        else
            p1 = 'Unspecified';
        end

        p2 = ' increase in possibility of disease';
        if lr>1
            p3 = ' presence';
        else
            p3 = ' absence';
        end
        txt = strcat(p1,p2,p3);
    end 

    function ci = newcombe(p) 
        % NEWCOMBE  Newcombe-Wilson score interval for a proportion.
        a = 2*N*p+z^2;
        b = z*realsqrt(z^2-2-1/N+4*p*(N*(1-p)+1));
        c = 2*(N+z^2);
        ci(1) = max([0 (a-b-1)/c]);
        ci(2) = min([1 (a+b+1)/c]);
    end

    function [cor,txt] = bca(or,orci) 
        % BCA  Bayesian credibility assessment for odds ratio.
        lorci = reallog(orci);
        if or<1
            sk = -1;
        else
            sk = 1;
        end
        cor = exp(sk*diff(lorci)^2/(4*realsqrt(prod(lorci)))); % Critical odds ratio
        if or>cor
            txt = 'Test is credible';
        else
            txt = 'Test is not credible';
        end
    end

end
