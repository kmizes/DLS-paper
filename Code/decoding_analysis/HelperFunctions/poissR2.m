function [poiss_pR2, poiss_dev] = poissR2(y, y_pred, y_null)
% poissR2: calculate pseudoR2 and deviance for poisson regression results based on Cameron and Windmeijer 1997
% size(y,1) should equal number of observations/predictions per set.
% Sets represented by different columns.

if size(y,2)==1  && size(y_pred,2)>1 % >1 set of predictions and 1 set of observations
    yy = repmat(y,1,size(y_pred,2));
    yy_pred = y_pred;
elseif size(y,2)>1 && size(y_pred,2)==1 % >1 of observations and 1 set of predictions
    yy = y;
    yy_pred = repmat(y_pred,1,size(y,2));
elseif isequal(size(y),size(y_pred)) % equal number of observations and predictions
    yy = y;
    yy_pred = y_pred;
else
    error('Mismatch between numbers of observations and predictions');
end
        
if nargin < 3
    y_null = nanmean(y(:));
    yy_null = ones(size(yy))*y_null;
else
    if numel(y_null) == 1 % a single null value has been provided
        yy_null = ones(size(yy))*y_null;
    elseif size(y_null,2)==1 && size(y_null,1)==size(yy,1) % one null value for each observation
        yy_null = repmat(y_null,1, size(yy,2));
    elseif size(y_null,1)==1 && size(y_null,2)==size(yy,2) % one null value for each set
        yy_null = repmat(y_null,size(yy,1),1);
    elseif isequal(size(y_null),size(yy)) % same sized null values and actual data
        yy_null = y_null;
    else
        error('Mismatch between numbers of observations and null predictions');
    end
end


poiss_dev = 2 * sum(yy.*log((yy+eps)./(yy_pred+eps)) - (yy - yy_pred), 1); % deviance of fitted model

null_dev = 2 * sum(yy.*log((yy+eps)./(yy_null+eps)) - (yy - yy_null), 1); % deviance of null model (mean only)

poiss_pR2 = 1 - poiss_dev ./ null_dev;

end

