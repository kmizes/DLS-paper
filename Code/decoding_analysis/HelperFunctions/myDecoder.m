function [R2_test, R2_train] = myDecoder(X, y, trainVsTest, trainCVfold, decode_method, R2_sets, varargin)
% X can be a numeric array or a cell array with the first dimension
% representing trial number in a session..
% if X is a cell array, then so must y
% trainVsTest can be either a scalar or a vector of indices with length = size(X,1)
% trainCV can be either a scalar or a vector of indices with length = sum(trainVsTest~=1)
% decode_method: 'GLM', 'SVR', 'NN', 'TreeBoost'

% fitStats = struct;

if nargin < 3
    trainVsTest = 4; % reciprocal of the fraction of test data
end

if nargin < 4
    trainCVfold = 5;
end

if nargin < 5
    decode_method = 'GLM';
end

rng(1); % to ensure the same cross-val folds every run

if numel(trainVsTest) == 1
    trainVsTest = crossvalind('kfold', size(X,1), trainVsTest);
end


if isnumeric(X) % X is obs * features
    X_train = X(trainVsTest~=1,:);
    y_train = y(trainVsTest~=1,:);
    
    X_test = X(trainVsTest==1,:);
    y_test = y(trainVsTest==1,:);

else % cell array of trials
   X_train = cell2mat(X(trainVsTest~=1,:));
   y_train = cell2mat(y(trainVsTest~=1,:));
   
   X_test = cell2mat(X(trainVsTest==1,:));
   y_test = cell2mat(y(trainVsTest==1,:));
   
end

if numel(trainCVfold)==1
    trainCVfold = crossvalind('kfold', sum(trainVsTest~=1), trainCVfold); % for regularization path
end

if ~isnumeric(X) % trial-wise cross-validation
    trainCVfold = cellfun(@(x, cvi) cvi * ones(size(x,1),1), X(trainVsTest~=1,1), num2cell(trainCVfold), 'unif', 0);
    trainCVfold = vertcat(trainCVfold{:});
end


switch decode_method
    case 'GLM'
        if nargin < 7
            alpha = 0.9;
        else
            alpha = varargin{1};
        end

        if nargin < 8
            lambdas = fliplr(logspace(-6, -1, 20)); % default grid of lambdas
        %     lambdas = fliplr(logspace(-5, -1, 20)); % default grid of lambdas
        else
            lambdas = varargin{2};
        end

        addpath('./glmnet_matlab/');

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
        
        if size(y_train,2) == 1
            mdl = cvglmnet(X_train, y_train, 'gaussian', options, [], [], trainCVfold);
        else
            mdl = cvglmnet(X_train, y_train, 'mgaussian', options, [], [], trainCVfold);
        end

        y_train_pred = cvglmnetPredict(mdl, X_train, 'lambda_min', 'response');
        y_pred = cvglmnetPredict(mdl, X_test, 'lambda_min', 'response');
        coeff = cvglmnetCoef(mdl, 'lambda_min');
    
    case 'SVR'
        y_train_pred = NaN(size(y_train));
        y_pred = NaN(size(y_test));        
        for d = 1 : size(y_train,2)
            mdl = fitrsvm(X_train, y_train(:,d), 'BoxConstraint', 3, 'KernelFunction', 'linear', 'Standardize', false);
            y_train_pred(:,d) = predict(mdl, X_train);
            y_pred(:,d) = predict(mdl, X_test);
        end
        
    case 'TreeBoost'
        y_train_pred = NaN(size(y_train));
        y_pred = NaN(size(y_test));        
        for d = 1 : size(y_train,2)
            mdl = fitrensemble(X_train, y_train(:,d), 'Method', 'LSBoost', 'Learners', templateTree('MaxNumSplits',3), 'LearnRate', 0.3, 'NumLearningCycles',300);
            y_train_pred(:,d) = predict(mdl, X_train);
            y_pred(:,d) = predict(mdl, X_test);
        end
        
    case 'NN'
        layers = [
            imageInputLayer([size(X_train,2),1,1], 'Normalization', 'none')
            
            dropoutLayer(0.25)
            fullyConnectedLayer(400)
            reluLayer
            
            dropoutLayer(0.25)
            fullyConnectedLayer(400)
            reluLayer
            
            fullyConnectedLayer(size(y_train,2))
            regressionLayer];
        
        optns = trainingOptions('adam', 'MaxEpochs', 10, 'Verbose', 0);
        
        mdl = trainNetwork(permute(X_train, [2,3,4,1]), y_train, layers, optns);
        y_train_pred = predict(mdl, permute(X_train, [2,3,4,1]));
        y_pred = predict(mdl, permute(X_test, [2,3,4,1]));
end

if nargin < 6
    R2_sets = num2cell(logical(eye(size(y_test,2))),2); % calc R2 on each dim separately
end

R2_train = NaN(size(R2_sets));
R2_test = NaN(size(R2_sets));
for r = 1 : length(R2_sets)    
    R2_train(r) = myR2(y_train(:,R2_sets{r}), y_train_pred(:,R2_sets{r}));
    R2_test(r) = myR2(y_test(:,R2_sets{r}), y_pred(:,R2_sets{r}));
end

end
%% R2 function

function [R2, mse] = myR2(dat, dat_pred)
    dat_null = repmat(nanmean(dat,1),size(dat,1),1);
    R2 = 1 - nansum((dat(:) - dat_pred(:)).^2,1) ./ sum((dat(:)-dat_null(:)).^2,1);
    mse = nanmean((dat(:) - dat_pred(:)).^2,1);
end

