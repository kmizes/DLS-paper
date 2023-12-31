function [fitStats, y_pred_all] = myDecoder_Kfold_R2wts(X, y, trainVsTest, decode_method, R2_sets, R2_wts, varargin)
% X can be a numeric array or a cell array with the first dimension
% representing trial number in a session..
% if X is a cell array, then so must y
% trainVsTest can be either a scalar or a vector of indices with length = size(X,1)
% decode_method: 'GLM', 'SVR', 'NN' or 'TreeBoost'

if nargin < 3
    trainVsTest = 4; % reciprocal of the fraction of test data
end

if nargin < 4
    decode_method = 'NN';
end

rng(1); % to ensure the same cross-val folds every run

if numel(trainVsTest) == 1
    trainVsTest = crossvalind('kfold', size(X,1), trainVsTest);
end

nFolds = max(trainVsTest);
tempStats = struct;

if isnumeric(y)
    y_pred_all = NaN(size(y));
else
    y_pred_all = cell(size(y));
    y_tr = cellfun(@(i,tr) cat(2,ones(numel(i),1)*tr,(1:length(i))'), y(:,1), num2cell(1:size(y,1))', 'unif', 0);
end

for k = 1 : nFolds
    
    if isnumeric(X) % X is obs * features
        X_train = X(trainVsTest~=k,:);
        y_train = y(trainVsTest~=k,:);

        X_test = X(trainVsTest==k,:);
        y_test = y(trainVsTest==k,:);

    else % cell array of trials
       X_train = cell2mat(X(trainVsTest~=k,:));
       y_train = cell2mat(y(trainVsTest~=k,:));

       X_test = cell2mat(X(trainVsTest==k,:));
       y_test = cell2mat(y(trainVsTest==k,:));
       
       y_tr_test = cell2mat(y_tr(trainVsTest==k,:)); 
    end

    switch decode_method
        case 'GLM'
            alpha = 0.9;
            lambdas = fliplr(logspace(-6, -1, 20)); % default grid of lambdas

            if nargin < 7
                trainCVfold = 5;
            else
                trainCVfold = varargin{1};
            end

            if numel(trainCVfold)==1
                trainCVfold = crossvalind('kfold', sum(trainVsTest~=k), trainCVfold); % for regularization path
            end

            if ~isnumeric(X) % trial-wise cross-validation
                trainCVfold = cellfun(@(x, cvi) cvi * ones(size(x,1),1), X(trainVsTest~=k,1), num2cell(trainCVfold), 'unif', 0);
                trainCVfold = vertcat(trainCVfold{:});
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

                dropoutLayer(0.05)
                fullyConnectedLayer(400)
                reluLayer

                dropoutLayer(0.05)
                fullyConnectedLayer(400)
                reluLayer

                fullyConnectedLayer(size(y_train,2))
                regressionLayer];

            %optns = trainingOptions('adam', 'MaxEpochs', 10, 'Verbose', 0, 'ExecutionEnvironment','cpu');
            optns = trainingOptions('sgdm', 'MaxEpochs', 50, 'Verbose', 0, 'ExecutionEnvironment','cpu');
            
            mdl = trainNetwork(permute(X_train, [2,3,4,1]), y_train, layers, optns);
            y_train_pred = predict(mdl, permute(X_train, [2,3,4,1]));
            y_pred = predict(mdl, permute(X_test, [2,3,4,1]));
        case 'NNbig'
            layers = [
                imageInputLayer([size(X_train,2),1,1], 'Normalization', 'none')

                dropoutLayer(0.05)
                fullyConnectedLayer(400)
                reluLayer

                dropoutLayer(0.05)
                fullyConnectedLayer(400)
                reluLayer

                dropoutLayer(0.05)
                fullyConnectedLayer(400)
                reluLayer
                
		fullyConnectedLayer(size(y_train,2))
                regressionLayer];

            %optns = trainingOptions('adam', 'MaxEpochs', 10, 'Verbose', 0, 'ExecutionEnvironment','cpu');
            optns = trainingOptions('sgdm', 'MaxEpochs', 200, 'Verbose', 0, 'ExecutionEnvironment','cpu');
            
            mdl = trainNetwork(permute(X_train, [2,3,4,1]), y_train, layers, optns);
            y_train_pred = predict(mdl, permute(X_train, [2,3,4,1]));
            y_pred = predict(mdl, permute(X_test, [2,3,4,1]));
    end

    if nargin < 5 || size(y_test,2)==1
        R2_sets = num2cell(logical(eye(size(y_test,2))),2); % calc R2 separately for each response dim
    end
    
    if nargin < 6 || size(y_test,2)==1
        R2_wts = cellfun(@(x) ones(size(x)), R2_sets, 'unif', 0);
    end

    for r = 1 : length(R2_sets)
        tempStats(k,r).ntrain = size(y_train,1);
        tempStats(k,r).ntest = size(y_test,1);
        [tempStats(k,r).R2_train, tempStats(k,r).mse_train] = myR2(y_train(:,R2_sets{r}), y_train_pred(:,R2_sets{r}), R2_wts{r});
        [tempStats(k,r).R2_test, tempStats(k,r).mse_test] = myR2(y_test(:,R2_sets{r}), y_pred(:,R2_sets{r}), R2_wts{r});
    end
    
    % collect predicted y
    if isnumeric(y)
        y_pred_all(trainVsTest==k,:) = y_pred;
    else
        for tr = unique(y_tr_test(:,1))'
            y_pred_all(tr,:) = num2cell(y_pred(y_tr_test(:,1)==tr,:),1);
        end
    end
end

fitStats = struct;
flds = fieldnames(tempStats);
for r = 1 : length(R2_sets)
    for f = 1:length(flds)
        fitStats(r).(flds{f}) = mean(vertcat(tempStats(:,r).(flds{f})),1); % mean across folds
    end
end

end

%% R2 function

function [R2, mse] = myR2(dat, dat_pred, wts)
    dat = dat .* repmat(wts, size(dat,1), 1);
    dat_pred = dat_pred .* repmat(wts, size(dat_pred,1), 1);
    
    dat_null = repmat(nanmean(dat,1),size(dat,1),1);
    R2 = 1 - nansum((dat(:) - dat_pred(:)).^2,1) ./ sum((dat(:)-dat_null(:)).^2,1);
    mse = nanmean((dat(:) - dat_pred(:)).^2,1);
end
