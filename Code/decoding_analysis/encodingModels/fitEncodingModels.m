function encFits = fitEncodingModels(sess, encParams)
%FITENCODINGMODELS

addpath('./HelperFunctions/');

%% unpack params

preTap1 = encParams.preTap1; % 0.2 s
postTap2 = encParams.postTap2; % 0.2 s
minIPI = encParams.minIPI; % 0.3 s
maxIPI = encParams.maxIPI; % 1.2 s
bin = encParams.bin; % 0.025 s
minTrials = encParams.minTrials; % 20 trials
minSpikeFrac = encParams.minSpikeFrac; % 0.25. min ratio of spikes to trials to include units
dt = encParams.dt; % 0 s. delay between neural activity and kinematics for one-to-one regression
t_length = encParams.t_length; % 1. In multiples of encParams.bin. For one-to-one regression.
dt_list = encParams.dts; % (-150:25:150)*1e-3 s. delay between neural activity and kinematics for many-to-one regression
t_lengths = encParams.t_lengths; % 1:2:7. In multiples of encParams.bin. Should be odd numbers.
trainVsTest = encParams.trainVsTest; % 4. reciprocal of fraction of test trials
trainCVfold = encParams.trainCVfold; % 5. num folds for cross-val regularization using training data
timeWarpMode = encParams.timeWarpMode; % 1, 2 or 3. time warp method.

%% global params

samprate = 3e4; % Hz
traj_preTap = 0.5; % s
traj_postTap = 1.5; % s
fr_time = 1/120; % s
acc_time = 1/150; % s

%% define trials

IPIs = (sess.tapTimes.tap2Times(:,5) - sess.tapTimes.tap1Times(:,5)) / samprate;

modes = sess.modes;

trSet = find( IPIs >= minIPI ...
            & IPIs <= maxIPI ...
            & (IPIs+postTap2) <= traj_postTap ...
            & isnan(sess.stimTimes(:,3)) ...
            & (sess.modes > 0));

mode_list = unique(modes(trSet));
mode_list = mode_list(arrayfun(@(m) sum(modes(trSet) == m), mode_list) >= minTrials);

trSet = trSet(ismember(modes(trSet), mode_list));

postTap = max(IPIs(trSet)) + postTap2;
postTap = ceil(postTap/bin) * bin; % round up to bin size
histedges = -preTap1:bin:postTap;
binTimes = histedges(1:end-1) + bin/2;


%% get kinematics

tTraj = -traj_preTap:fr_time:traj_postTap;
tTraj_v = tTraj(1:end-1) + fr_time/2;
tTraj_a = tTraj_v(1:end-1) + fr_time/2;

if contains(sess.ratName, 'Gandhar')
    pawLInd = find(strcmp(sess.traj.Robjs, 'paw2'),1);
    headLInd = find(strcmp(sess.traj.Lobjs, 'nose'),1);
    pawRInd = find(strcmp(sess.traj.Robjs, 'paw'),1);
    headRInd = find(strcmp(sess.traj.Robjs, 'ear'),1);
    traj = double(cat(4, sess.traj.Rtraj(trSet,:,1:2,pawLInd), sess.traj.Rtraj(trSet,:,1:2,pawRInd), ...
            sess.traj.Ltraj(trSet,:,1:2,headLInd), sess.traj.Rtraj(trSet,:,1:2,headRInd)));
        
elseif contains(sess.ratName, 'SW116')
    pawLInd = find(strcmp(sess.traj.Lobjs, 'paw'),1);
    headLInd = find(strcmp(sess.traj.Lobjs, 'eye'),1);
    pawRInd = find(strcmp(sess.traj.Robjs, 'paw'),1);
    headRInd = find(strcmp(sess.traj.Robjs, 'eye'),1);
    traj = double(cat(4, sess.traj.Ltraj(trSet,:,1:2,pawLInd), sess.traj.Rtraj(trSet,:,1:2,pawRInd), ...
            sess.traj.Ltraj(trSet,:,1:2,headLInd), sess.traj.Rtraj(trSet,:,1:2,headRInd)));
    
else
    pawLInd = find(strcmp(sess.traj.Lobjs, 'paw'),1);
    headLInd = find(strcmp(sess.traj.Lobjs, 'ear'),1);
    pawRInd = find(strcmp(sess.traj.Robjs, 'paw'),1);
    headRInd = find(strcmp(sess.traj.Robjs, 'ear'),1);
    traj = double(cat(4, sess.traj.Ltraj(trSet,:,1:2,pawLInd), sess.traj.Rtraj(trSet,:,1:2,pawRInd), ...
            sess.traj.Ltraj(trSet,:,1:2,headLInd), sess.traj.Rtraj(trSet,:,1:2,headRInd)));
end

tap1Fr = find(tTraj >= 0, 1, 'first');
tap2Fr = arrayfun(@(x) find(tTraj > x, 1, 'first'), IPIs(trSet));

startFr = tap1Fr - preTap1/fr_time; 
endFr = tap2Fr + postTap2/fr_time;
    
NanThresh = 10; % do not interpolate on trials with more NaN frames than this
FrNanCnt = cellfun(@(x,y) max(max(sum(isnan(x(1,startFr:y+1,:,:)),2),[],3),[],4), num2cell(traj, [2,3,4]), num2cell(endFr)); % count number of max NaNs over objects in trial

disp([num2str(sum(FrNanCnt > NanThresh)) ' of ' num2str(length(trSet)) ' frames lost due to NaNs in trajectory']);

trSet = trSet(FrNanCnt <= NanThresh);
traj = traj(FrNanCnt <= NanThresh,:,:,:);

% interpolate over remaining NaNs for each object trajectory 
for tr = 1 : size(traj,1)
    for xy = 1 : size(traj, 3)
        for o = 1 : size(traj, 4)
            if any(squeeze(isnan(traj(tr,:,xy,o))))
                nan_pts = squeeze(isnan(traj(tr,:,xy,o)));
                t = 1:size(traj,2);
                traj(tr,nan_pts,xy,o) = interp1(t(~nan_pts), squeeze(traj(tr,~nan_pts,xy,o)), t(nan_pts), 'linear');            
            end
        end
    end
end

% Spline smoothing of trajectories
warning('off', 'SPLINES:CHCKXYWP:NaNs');
traj_s = NaN(size(traj));
for tr = 1 : size(traj,1)
    for xy = 1 : size(traj,3)
        for o = 1 : size(traj,4)
            traj_s(tr,:,xy,o) = csaps(1:size(traj,2), traj(tr,:,xy,o), .1, 1:size(traj,2));
        end
    end
end
warning('on', 'SPLINES:CHCKXYWP:NaNs');

% NaN any trajectories that cross image borders
traj_s(traj_s < 0 | traj_s > 640) = NaN;

% convert to mm scale
if isfield(sess.traj, 'LmmPerPixel')
    traj_s(:,:,:,1:2:end) = traj_s(:,:,:,1:2:end)*sess.traj.LmmPerPixel;
    traj_s(:,:,:,2:2:end) = traj_s(:,:,:,2:2:end)*sess.traj.RmmPerPixel;
end

traj_p = traj_s - repmat(nanmean(traj_s(:,startFr,:,:),1), 1,size(traj_s,2),1,1);

% compute vel. and acc.
traj_v = NaN(size(traj_p));
traj_a = NaN(size(traj_p));
for xy = 1 : size(traj,3)
    for o = 1 : size(traj,4)
        traj_v(:,:,xy,o) = interp1(tTraj_v, diff(traj_p(:,:,xy,o),1,2)', tTraj, 'linear', 'extrap')';
        traj_a(:,:,xy,o) = interp1(tTraj_a, diff(diff(traj_p(:,:,xy,o),1,2),1,2)', tTraj, 'linear', 'extrap')';
    end
end

% resample acc to timescale of trajectory
tAcc = -traj_preTap:acc_time:traj_postTap;

acc = zeros(size(traj,1), length(tTraj), size(sess.acc.filtAcc,3));
for a = 1 : size(sess.acc.filtAcc,3)
    acc(:,:,a) = interp1(tAcc, sess.acc.filtAcc(trSet,:,a)', tTraj, 'linear', 'extrap')';
end
        
% kinematics: vel and acc
kin = cat(3, reshape(traj_v, length(trSet), length(tTraj), []), ...
            reshape(traj_a, length(trSet), length(tTraj), []), ...
            acc);        

% unsigned kinematics (vigor variables)
kinmag = cat(3, squeeze(sqrt(sum(traj_v.^2,3))), ...
            squeeze(sqrt(sum(traj_a.^2,3))), ...
            squeeze(sqrt(sum(acc.^2,3))));
        
% position (state variable)
pos = reshape(traj_p, length(trSet), length(tTraj), []);
        

%% bin spike train

nunits = length(sess.units);
upeth = NaN(length(trSet), length(binTimes), nunits);
for u = 1 : nunits
    upeth(:,:,u) = calcPETH(sess.units(u).spikes, sess.tapTimes.tap1Times(trSet,5), histedges);    
end

for tr = 1:length(trSet)
    upeth(tr,binTimes>(IPIs(trSet(tr))+postTap2),:) = NaN;
end

% get rid of units that do not fire at least 25% as many spikes as there are trials
u_ind = squeeze(nansum(nansum(upeth,1),2)) >= length(trSet)*minSpikeFrac;
upeth = upeth(:,:,u_ind);
u_nums = [sess.units(u_ind).unitnum];

%% Initialize fit struct

encFits = struct('ratName', sess.ratName, 'exptID', sess.exptID, 'defID', sess.defID, ...
    'bin', bin, 'preTap1', preTap1, 'postTap2', postTap2, 'timeWarpMode', timeWarpMode, ...
    'dt', dt, 't_length', t_length, 'dts', dt_list, 't_lengths', t_lengths, ...
    'unitNums', u_nums, 'unitTypes', [sess.units(u_ind).unitType], ...
    'trSet', trSet, 'mode_list', mode_list);

encFits.unitNames = {sess.units(u_ind).unitname};

%% Regression models

if ~isempty(u_nums)
    %% Generate predictor set for phase/time vector co-variates

    targ = nanmedian(IPIs(trSet)); % target IPI
    t_tvec = tTraj(1):fr_time:3; 
    tvec_sd = 0.1; % width (std dev) of a time field (sec)
    tvec_dt = tvec_sd; % spacing between time fields (sec)
    tvec_mus = sort(cat(2, 0:-tvec_dt:-preTap1, tvec_dt:tvec_dt:postTap2+targ));
    n_tvecs = length(tvec_mus);

    tvecs = NaN(length(trSet), length(tTraj), n_tvecs); % three types of time-warping

    switch timeWarpMode
        case 1
            % warp firing fields
            tvecs_template = arrayfun(@(mu) normpdf(t_tvec, mu, tvec_sd)', tvec_mus, 'unif', 0);
            tvecs_template = [tvecs_template{:}];
            for tr = 1 : length(trSet)
                t_warp_IPI = linspace(0, targ, sum(tTraj>=0 & tTraj<=targ)*IPIs(trSet(tr))/targ);    
                t_warp = [tTraj(1):fr_time:-fr_time, t_warp_IPI, targ+cumsum(fr_time*ones(1,length(tTraj)-(tap1Fr+length(t_warp_IPI)-1)))];
                tvecs(tr,:,:) = interp1(t_tvec, tvecs_template, t_warp); % warp FRs (less spikes but same peak FR on shorter trials)
            end

        case 2        
             % warp spike train (higher FR but same number of spikes on shorter trials)
            tvecs_template = arrayfun(@(mu) normpdf(t_tvec, mu, tvec_sd)', tvec_mus, 'unif', 0);
            tvecs_template = [tvecs_template{:}];
            for tr = 1 : length(trSet)
                t_warp_IPI = linspace(0, targ, sum(tTraj>=0 & tTraj<=targ)*IPIs(trSet(tr))/targ);    
                t_warp = [tTraj(1):fr_time:-fr_time, t_warp_IPI, targ+cumsum(fr_time*ones(1,length(tTraj)-(tap1Fr+length(t_warp_IPI)-1)))];
                tvecs_temp = interp1(t_tvec, tvecs_template, t_warp); 
                tvecs(tr,:,:) = tvecs_temp ./ repmat(sum(tvecs_temp,2),1,length(tTraj),1);
            end

        case 3        
            % shift, but not warp firing fields
            for tr = 1 : length(trSet)
                tvec_mus_warp = cat(2, tvec_mus(tvec_mus < 0), linspace(0, IPIs(trSet(tr)), sum(tvec_mus>=0 & tvec_mus<=targ)), tvec_mus(tvec_mus > targ)+(IPIs(trSet(tr))-targ));
                temp = arrayfun(@(mu) normpdf(tTraj, mu, tvec_sd)', tvec_mus_warp, 'unif', 0);
                tvecs(tr,:,:) = [temp{:}];
            end        
    end


    %% Predictor set for phase/time vector co-variates with modes

    % simple on/off mode set
    mode_vec = arrayfun(@(m) mode_list(:)'==m, modes(trSet), 'unif', 0);
    mode_vec = double(vertcat(mode_vec{:})); 

    % mode timevec
    mode_tvecs = zeros(size(tvecs,1),size(tvecs,2),size(tvecs,3),length(mode_list));
    for tr = 1 : length(trSet)
        mode_tvecs(tr,:,:,mode_list==modes(trSet(tr))) = tvecs(tr,:,:);
    end

    mode_tvecs = reshape(mode_tvecs,size(mode_tvecs,1),size(mode_tvecs,2),size(mode_tvecs,3)*length(mode_list));


    %% define train and test indices

    rng(1);

    % generate train versus test ind so as to equally divide modes
    trainVsTestInd = zeros(size(trSet));
    for m = mode_list'
        trainVsTestInd(modes(trSet)==m) = crossvalind('kfold', sum(modes(trSet)==m), trainVsTest);
    end

    trainCVind = zeros(sum(trainVsTestInd~=1),1);
    for m = mode_list'
        trainCVind(modes(trSet(trainVsTestInd~=1))==m) = crossvalind('kfold', sum(modes(trSet(trainVsTestInd~=1))==m), trainCVfold);
    end


    %% prepare datasets for one-to-one regression (dt = 0, t_length = 1)

    num_t = (t_length-1)/2;

    trajPerBin = bin / fr_time; % in units of frame time
    hwind_traj = (((num_t*2+1)*trajPerBin)-1)/2;

    % prepare time vec regressors
    tvec_dt = NaN(length(trSet), length(binTimes), n_tvecs);
    for tv = 1 : n_tvecs
        tvec_dt(:,:,tv) = interp1(tTraj + dt, tvecs(:,:,tv)', binTimes)';           
    end 
    mode_tvec_dt = NaN(length(trSet), length(binTimes), size(mode_tvecs,3));
    for mtv = 1 : size(mode_tvecs,3)
        mode_tvec_dt(:,:,mtv) = interp1(tTraj + dt, mode_tvecs(:,:,mtv)', binTimes)';            
    end

    kin_dt = NaN(length(trSet), length(binTimes), size(kin,3), hwind_traj*2+1);
    for b = 1 : length(binTimes)
        for o = 1 : size(kin,3)
            kin_dt(:,b,o,:) = interp1(tTraj + dt, kin(:,:,o)', binTimes(b) + fr_time*(-hwind_traj:hwind_traj))';
        end
    end 

    pos_dt = NaN(length(trSet), length(binTimes), size(pos,3), hwind_traj*2+1);
    for b = 1 : length(binTimes)
        for o = 1 : size(pos,3)
            pos_dt(:,b,o,:) = interp1(tTraj + dt, pos(:,:,o)', binTimes(b) + fr_time*(-hwind_traj:hwind_traj))';
        end
    end 


    kinmag_dt = NaN(length(trSet), length(binTimes), size(kinmag,3), hwind_traj*2+1);
    for b = 1 : length(binTimes)
        for o = 1 : size(kinmag,3)
            kinmag_dt(:,b,o,:) = interp1(tTraj + dt, kinmag(:,:,o)', binTimes(b) + fr_time*(-hwind_traj:hwind_traj))';
        end
    end        
    kinmag_dt = reshape(kinmag_dt, size(kinmag,1), length(binTimes), numel(kinmag_dt(1,1,:,:)));

    u_dat = cell(length(trSet),size(upeth,3));
    kin_dat = cell(length(trSet), 19); % separate predictor dimensions in different columns of the cell array to make it easier to pick and choose particular sets of predictors
    pos_dat = cell(length(trSet), 1);
    kinmag_dat = cell(length(trSet), 1);
    tvec_dat = cell(length(trSet), 1);
    mode_dat = cell(length(trSet), 1);
    mode_tvec_dat = cell(length(trSet), 1);


    % collect regressors and response data
    NaNcheck = cat(3, reshape(kin_dt,size(kin_dt,1),size(kin_dt,2),numel(kin_dt(1,1,:,:))), upeth);
    for tr = 1 : length(trSet)      
        NaNframes = any(isnan(NaNcheck(tr,:,:)),3);
        u_dat(tr,:) = num2cell(permute(upeth(tr,~NaNframes,:),[2,3,1]), 1);
        kin_dat(tr,:) = num2cell(permute(kin_dt(tr,~NaNframes,:,:), [2,4,3,1]), [1,2]);
        pos_dat{tr} = permute(pos_dt(tr,~NaNframes,:), [2,3,1]);            
        kinmag_dat{tr} = permute(kinmag_dt(tr,~NaNframes,:), [2,3,1]);            
        tvec_dat{tr} = permute(tvec_dt(tr,~NaNframes,:), [2,3,1]);
        mode_tvec_dat{tr} = permute(mode_tvec_dt(tr,~NaNframes,:), [2,3,1]);            
        mode_dat{tr} = repmat(mode_vec(tr,:),sum(~NaNframes),1);
    end  

    %% define kinematic sets

    Lpaw_kinInd = [1:2, 9:10];
    Rpaw_kinInd = [3:4, 11:12];
    head_kinInd = [5:8, 13:16, 17:19];
    vel_kinInd = [1:8];
    acc_kinInd = [9:16, 17:19];


    % define ipsi versus contra paws
    if contains(sess.ratName, {'Gandhar','Hamir','SW158','SW116','SW166'}) % left-side recording
        ipsiPaw_kinInd = Lpaw_kinInd;
        contraPaw_kinInd = Rpaw_kinInd;
    elseif contains(sess.ratName, {'Dhanashri','Gorakh','JaunpuriL','Hindol','SW160','SW163','SW233'}) % right-side recording
        ipsiPaw_kinInd = Rpaw_kinInd;
        contraPaw_kinInd = Lpaw_kinInd;
    end


    %% Get R2 for one-to-one regression models

    encFits.kinR2 = NaN(length(u_nums),1);
    encFits.posR2 = NaN(length(u_nums),1);
    encFits.ipsiPawR2 = NaN(length(u_nums),1);
    encFits.contraPawR2 = NaN(length(u_nums),1);
    encFits.pawsR2 = NaN(length(u_nums),1);
    encFits.headR2 = NaN(length(u_nums),1);
    encFits.velR2 = NaN(length(u_nums),1);
    encFits.accR2 = NaN(length(u_nums),1);
    encFits.tvecR2 = NaN(length(u_nums),1);
    encFits.tvec_modeR2 = NaN(length(u_nums),1);
    encFits.modetvecR2 = NaN(length(u_nums),1);
    encFits.kinmagR2 = NaN(length(u_nums),1);

    tic;
    for u = 1 : length(u_nums)
        disp(['unit ' num2str(u)]);
        encFits.kinR2(u) = myPoissGLM_ElasticReg(kin_dat, u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.posR2(u) = myPoissGLM_ElasticReg(pos_dat, u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.ipsiPawR2(u) = myPoissGLM_ElasticReg(kin_dat(:,ipsiPaw_kinInd), u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.contraPawR2(u) = myPoissGLM_ElasticReg(kin_dat(:,contraPaw_kinInd), u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.pawsR2(u) = myPoissGLM_ElasticReg(kin_dat(:,[ipsiPaw_kinInd,contraPaw_kinInd]), u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.headR2(u) = myPoissGLM_ElasticReg(kin_dat(:,head_kinInd), u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.velR2(u) = myPoissGLM_ElasticReg(kin_dat(:,vel_kinInd), u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.accR2(u) = myPoissGLM_ElasticReg(kin_dat(:,acc_kinInd), u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.tvecR2(u) = myPoissGLM_ElasticReg(tvec_dat, u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.tvec_modeR2(u) = myPoissGLM_ElasticReg(cat(2, tvec_dat, mode_dat), u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.modetvecR2(u) = myPoissGLM_ElasticReg(mode_tvec_dat, u_dat(:,u), trainVsTestInd, trainCVind);
        encFits.kinmagR2(u) = myPoissGLM_ElasticReg(kinmag_dat, u_dat(:,u), trainVsTestInd, trainCVind);
        toc;
    end

    %% prepare datasets for kinematic, time and mode data for multiple dts and lengths of the trajectory

    num_ts = (t_lengths-1)/2; 

    u_dat_n = cell(length(trSet),size(upeth,3),length(num_ts),length(dt_list));
    kin_dat_n = cell(length(trSet),length(num_ts),length(dt_list));
    tvec_dat_n = cell(length(trSet),length(num_ts),length(dt_list));
    mode_dat_n = cell(length(trSet),length(num_ts),length(dt_list));
    mode_tvec_dat_n = cell(length(trSet),length(num_ts),length(dt_list));

    trajPerBin = bin / fr_time; % in units of frame time
    hwind_trajs = (((num_ts*2+1)*trajPerBin)-1)/2;

    for t = 1 : length(dt_list)

        % prepare time vec regressors
        tvec_dt = NaN(length(trSet), length(binTimes), n_tvecs);
        mode_tvec_dt = NaN(length(trSet), length(binTimes), size(mode_tvecs,3));

        for tv = 1 : n_tvecs
            tvec_dt(:,:,tv) = interp1(tTraj + dt_list(t), tvecs(:,:,tv)', binTimes)';           
        end 
        for mtv = 1 : size(mode_tvecs,3)
            mode_tvec_dt(:,:,mtv) = interp1(tTraj + dt_list(t), mode_tvecs(:,:,mtv)', binTimes)';            
        end

        u_dt = upeth;        

        for n = 1 : length(num_ts)                    
            kin_dt_n = NaN(length(trSet), length(binTimes), size(kin,3), hwind_trajs(n)*2+1);
            for b = 1 : length(binTimes)
                for o = 1 : size(kin,3)
                    kin_dt_n(:,b,o,:) = interp1(tTraj + dt_list(t), kin(:,:,o)', binTimes(b) + fr_time*(-hwind_trajs(n):hwind_trajs(n)))';
                end
            end        
            kin_dt_n = reshape(kin_dt_n, size(kin,1), length(binTimes), numel(kin_dt_n(1,1,:,:)));

            NaNcheck = cat(3, kin_dt_n, u_dt);
            for tr = 1 : length(trSet)      
                NaNframes = any(isnan(NaNcheck(tr,:,:)),3);
                u_dat_n(tr,:,n,t) = num2cell(permute(u_dt(tr,~NaNframes,:), [2,3,1]),1);
                kin_dat_n{tr,n,t} = permute(kin_dt_n(tr,~NaNframes,:), [2,3,1]);           
                tvec_dat_n{tr,n,t} = permute(tvec_dt(tr,~NaNframes,:), [2,3,1]);
                mode_tvec_dat_n{tr,n,t} = permute(mode_tvec_dt(tr,~NaNframes,:), [2,3,1]);            
                mode_dat_n{tr,n,t} = repmat(mode_vec(tr,:),sum(~NaNframes),1);
            end  
        end
    end

    %% measure performance of GLM models on different predictor sets (cross-validated pseudo-R2)

    encFits.kin_n_R2s = NaN(length(u_nums),length(num_ts),length(dt_list)); % kinematics
    encFits.tvec_kin_n_R2s = NaN(length(u_nums),length(num_ts),length(dt_list)); % kinematics + time
    encFits.tvec_mode_kin_n_R2s = NaN(length(u_nums),length(num_ts),length(dt_list)); % kinematics + time + simple mode variable
    encFits.modetvec_kin_n_R2s = NaN(length(u_nums),length(num_ts),length(dt_list)); % kinematics + time X mode vector

    tic;
    for t = 1 : length(dt_list)
        disp(['************dt = ' num2str(dt_list(t)) '*****************']);    
        for n = 1 : length(num_ts)
            disp(['n = ' num2str(num_ts(n))]);
            for u = 1 :  length(u_nums)
                if u/5 == fix(u/5); sprintf('Unit #%d: took %.2f sec', u, toc); end
                encFits.kin_n_R2s(u,n,t) = myPoissGLM_ElasticReg(kin_dat_n(:,n,t), u_dat_n(:,u,n,t), trainVsTestInd, trainCVind);
                encFits.tvec_kin_n_R2s(u,n,t) = myPoissGLM_ElasticReg(cat(2, kin_dat_n(:,n,t), tvec_dat_n(:,n,t)), u_dat_n(:,u,n,t), trainVsTestInd, trainCVind);            
                encFits.tvec_mode_kin_n_R2s(u,n,t) = myPoissGLM_ElasticReg(cat(2, kin_dat_n(:,n,t), tvec_dat_n(:,n,t), mode_dat_n(:,n,t)), u_dat_n(:,u,n,t), trainVsTestInd, trainCVind);
                encFits.modetvec_kin_n_R2s(u,n,t) = myPoissGLM_ElasticReg(cat(2, kin_dat_n(:,n,t), mode_tvec_dat_n(:,n,t)), u_dat_n(:,u,n,t), trainVsTestInd, trainCVind);            
            end
            toc;
        end
    end

end

end

