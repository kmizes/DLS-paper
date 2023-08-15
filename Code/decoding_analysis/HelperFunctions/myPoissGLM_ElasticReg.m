function [fitStats, coeff, mdl] = myPoissGLM_ElasticReg(X, y, trainVsTest, trainCVfold, alpha, lambdas)
% X can be a numeric array or a cell array with the first dimension
% representing trial number in a session..
% if X is a cell array, then so must y
% trainVsTest can be either a scalar or a vector of indices with length = size(X,1)
% trainCV can be either a scalar or a vector of indices with length = sum(trainVsTest~=1)

fitStats = struct;

if nargin < 3
    trainVsTest = 4; % reciprocal of the fraction of test data
end

if nargin < 4
    trainCVfold = 5;
end

if nargin < 5
    alpha = 0.9;
end

if nargin < 6
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
[fitStats.pR2_test, fitStats.dev_test] = poissR2(y_test, cvglmnetPredict(mdl, X_test, 'lambda_min', 'response'));

coeff = cvglmnetCoef(mdl, 'lambda_min');