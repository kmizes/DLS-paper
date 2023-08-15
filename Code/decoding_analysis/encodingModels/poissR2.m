function [poissR2, poissDev] = poissR2(y, y_pred, y_null)
% poissR2: calculate pseudoR2 and deviance for poisson regression results based on Cameron and Windmeijer 1997

y = y(:);
if numel(y_pred) > size(y_pred,1) % there is > 1 set of predictions
    if size(y_pred,1) ~= numel(y)
        y_pred = y_pred';
    end
end
        
if nargin < 3
    y_null = mean(y);
end

y_null = y_null(:);

yy = repmat(y,1,size(y_pred,2));
yy_null = repmat(y_null,1,size(y_pred,2));

poissDev = 2 * sum(yy.*log((yy+eps)./(y_pred+eps)) - (yy-y_pred), 1); % deviance of fitted model

nullDev = 2 * sum(yy.*log((yy+eps)./(yy_null+eps)) - (yy - yy_null), 1); % deviance of null model (mean only)

poissR2 = 1 - poissDev ./ nullDev;

end

