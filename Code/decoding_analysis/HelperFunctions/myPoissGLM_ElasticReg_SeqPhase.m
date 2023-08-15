function [fitStats, coeff, mdl] = myPoissGLM_ElasticReg_SeqPhase(X, y, y_phase, trainVsTest, trainCVfold, alpha, lambdas)
% X can be a numeric array or a cell array with the first dimension
% representing trial number in a session..
% if X is a cell array, then so must y and y_phase
% y_phase indicates phase/time of each observation within its trial
% trainVsTest can be either a scalar or a vector of indices with length = size(X,1)
% trainCV can be either a scalar or a vector of indices with length = sum(trainVsTest~=1)

fitStats = struct;

if nargin < 4
    trainVsTest = 4; % reciprocal of the fraction of test data
end

if nargin < 5
    trainCVfold = 5;
end

if nargin < 6
    alpha = 0.9;
end

if nargin < 7
    lambdas = fliplr(logspace(-6, -1, 20)); % default grid of lambdas
%     lambdas = fliplr(logspace(-5, -1, 20)); % default grid of lambdas
end

addpath('./glmnet_matlab/');

rng(1); % to ensure the same cross-val folds every run

if numel(trainVsTest) == 1
    trainVsTest = crossvalind('kfold', size(X,1), trainVsTest);
end


if isnumeric(X) % X is obs * features
    X_train = X(trainVsTest~=1,:);
    y_train = y(trainVsTest~=1);
    
    X_test = X(trainVsTest==1,:);
    y_test = y(trainVsTest==1);

else % cell array of trials
   X_train = cell2mat(X(trainVsTest~=1,:));
   y_train = cell2mat(y(trainVsTest~=1));
   
   X_test = cell2mat(X(trainVsTest==1,:));
   y_test = cell2mat(y(trainVsTest==1));
   
   y_phase_test = cell2mat(y_phase(trainVsTest==1));   
   
   % trial identity of test set observations
   tr_test = cell2mat(cellfun(@(y,tr) ones(size(y))*tr, y(trainVsTest==1), num2cell(1:sum(trainVsTest==1),1)', 'unif', 0));
end

fitStats.ntrain = numel(y_train);
fitStats.ntest = numel(y_test);

if numel(trainCVfold)==1
    trainCVfold = crossvalind('kfold', sum(trainVsTest~=1), trainCVfold); % for regularization path
end

if ~isnumeric(X) % trial-wise cross-validation
    trainCVfold = cellfun(@(x, cvi) cvi * ones(size(x,1),1), X(trainVsTest~=1,1), num2cell(trainCVfold), 'unif', 0);
    trainCVfold = vertcat(trainCVfold{:});
end

% determine regularization path (lambda param) by cross-validation within training data
if ~isempty(lambdas)
    if numel(lambdas)>1
        options = glmnetSet(struct('lambda', lambdas, 'nlambda', numel(lambdas), 'alpha', alpha, 'standardize', 0));
    else
        options = glmnetSet(struct('lambda', [lambdas,1], 'nlambda', 2, 'alpha', alpha, 'standardize', 0));
    end
else
    options = glmnetSet(struct('alpha', alpha, 'standardize', 0)); % let glmnet pick its own vector of lambdas
end

mdl = cvglmnet(X_train, y_train, 'poisson', options, [], [], trainCVfold, 1);

[fitStats.pR2_train, fitStats.dev_train] = poissR2(y_train, cvglmnetPredict(mdl, X_train, 'lambda_min', 'response'));
y_pred = cvglmnetPredict(mdl, X_test, 'lambda_min', 'response');
[fitStats.pR2_test, fitStats.dev_test] = poissR2(y_test, y_pred(:));

% rng(1);
% % R2 expected if model captured all variability
% y_poiss = arrayfun(@(n) poissrnd(y_pred),1:100,'unif',0);
% [R2_mdlnoise, dev_mdlnoise] = poissR2(cat(2,y_poiss{:}), y_pred(:)); 
% R2_mdlnoise = [mean(R2_mdlnoise), std(R2_mdlnoise)];

coeff = cvglmnetCoef(mdl, 'lambda_min');

%% goodness of fit as a function of the phase

phase_bins = -.15:.1:0.95;
fitStats.phases = (phase_bins(1:end-1)+phase_bins(2:end))/2;
[~,~,ph_bn] = histcounts(y_phase_test, phase_bins);
fitStats.n_phase = arrayfun(@(b) sum(ph_bn==b), 1:length(phase_bins)-1);
fitStats.sp_phase = arrayfun(@(b) mean(y_test(ph_bn==b)), 1:length(phase_bins)-1);

pR2_phase = NaN(1,length(phase_bins)-1);
dev_phase = NaN(1,length(phase_bins)-1);
for b = 1 : length(phase_bins)-1
    [pR2_phase(b), dev_phase(b)] = poissR2(y_test(ph_bn==b), y_pred(ph_bn==b));
end
pR2_phase(isinf(pR2_phase)) = NaN;
dev_phase(isinf(dev_phase)) = NaN;
fitStats.pR2_phase = pR2_phase;
fitStats.dev_phase = dev_phase;


%% does difference in phase explain residual differences in encoding?

tr_dist = pdist(tr_test);
y_dist = pdist(sqrt(y_test)); % also applied simple variance stabilizing transform 
y_pred_dist = pdist(sqrt(y_pred)); % also applied simple variance stabilizing transform
d_phase = pdist(y_phase_test);

% ignore same-trial differences
y_dist = y_dist(tr_dist~=0)'; 
y_pred_dist = y_pred_dist(tr_dist~=0)'; 
d_phase = d_phase(tr_dist~=0)'; 

fitStats.n_dist = numel(y_dist);

% fit regression models with and without phase differences

CV_dist = round(rand(size(y_dist)))+1;

warning('off', 'stats:regress:RankDefDesignMat');
fitStats.coeff_dist = regress(y_dist(CV_dist==1), cat(2, ones(sum(CV_dist==1),1), y_pred_dist(CV_dist==1)));
fitStats.coeff_dist_dphase = regress(y_dist(CV_dist==1), cat(2, ones(sum(CV_dist==1),1), y_pred_dist(CV_dist==1), d_phase(CV_dist==1)));
warning('on', 'stats:regress:RankDefDesignMat');

y_dist_pred = fitStats.coeff_dist(1) + fitStats.coeff_dist(2)*y_pred_dist(CV_dist==2);
fitStats.R2_dist = 1 - sum((y_dist(CV_dist==2) - y_dist_pred).^2) / sum((y_dist(CV_dist==2) - mean(y_dist(CV_dist==2))).^2);

y_dist_dphase_pred = fitStats.coeff_dist_dphase(1) + fitStats.coeff_dist_dphase(2)*y_pred_dist(CV_dist==2) + fitStats.coeff_dist_dphase(3)*d_phase(CV_dist==2);
fitStats.R2_dist_dphase = 1 - sum((y_dist(CV_dist==2) - y_dist_dphase_pred).^2) / sum((y_dist(CV_dist==2) - mean(y_dist(CV_dist==2))).^2);


d_phase_bins = 0:0.1:1;
fitStats.dphases = (d_phase_bins(1:end-1)+d_phase_bins(2:end))/2;
[~,~,dph_bn] = histcounts(d_phase, d_phase_bins);
fitStats.n_dphases = arrayfun(@(dp) sum(dph_bn==dp), 1:length(d_phase_bins)-1);
fitStats.corr_dist_dphases = arrayfun(@(dp) corr(y_dist(dph_bn==dp), y_pred_dist(dph_bn==dp)), 1:length(d_phase_bins)-1);





