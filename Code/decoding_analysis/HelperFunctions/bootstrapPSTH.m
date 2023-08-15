function [p_psth, z_psth, sig_psth, conf_bnds] = bootstrapPSTH(psth, nboot, alpha)

% psth = ntrials X bin vector of FR/spike counts
% nboot = num of shuffles
% alpha (p-value) = 0.0001 by default (for each psth bin)
% p_psth = nbins sized probability vector that each bin occured by chance
% z_psth = zscores for each bin based on the mean and std of the shuffled
% distribution
% sig_psth = prob vector threshold by alpha value
% conf_bnds = confidence bounds at percentile = alpha and 1-alpha

if nargin < 3
    alpha = 0.0001; % at a prob of 1e-4, the prob of atleast one bin in a 100 bin psth being significant is 0.01
end

% if nboot < 10/alpha
%    disp('nboot needs to be atleast 10 times the reciprocal of alpha.'); 
% end
     
if any(psth(:) > 0)    % check if there is at least 1 spike in the psth

    % Create shuffled versions of the psth
    if nboot*size(psth,1) <= 1e9       
        boot_ind = randi(numel(psth), nboot*size(psth,1), 1);
        f_psth = psth(boot_ind(:));
        f_psth = sort(sum(reshape(f_psth, size(psth,1), nboot)));
        clear boot_ind;
    else % when shuffled data wont all fit in memory
        f_psth = zeros(1,nboot);
        for n = 1 : nboot
            boot_ind = randi(numel(psth), size(psth,1), 1);
            f_psth(n) = sum(psth(boot_ind(:)));
        end
        f_psth = sort(f_psth);
    end

    p_psth = zeros(size(sum(psth)));

    sum_psth = sum(psth);
    mf_psth = mean(f_psth);
    sdf_psth = std(f_psth);

    z_psth = (sum_psth - mf_psth)/sdf_psth;

    for bin = 1 : size(psth,2)
        if sum_psth(bin) >= mf_psth
            score = find(f_psth <= sum_psth(bin),1,'last') / numel(f_psth);
            score = 1 - score;
        else
            score = -find(f_psth > sum_psth(bin),1,'first') / numel(f_psth);
        end
        p_psth(bin) = score;
    end

    % threshold p-values
    sig_psth = p_psth <= alpha;

    conf_bnds = prctile(f_psth, [alpha*100, (1-alpha)*100]);
else
    p_psth = NaN(size(sum(psth)));
    z_psth = NaN(size(sum(psth)));
    sig_psth = NaN(size(sum(psth)));
    conf_bnds = NaN(1,2);
end

