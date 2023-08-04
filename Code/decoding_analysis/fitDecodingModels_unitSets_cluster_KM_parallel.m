

% i've done encoding, onto decoding?



% get unit count
% fpath = 'D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat47\ratBehUnitTrajSess_quality';
% fpath = 'D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat33\ratBehUnitTrajSess_quality';
% fpath = 'D:\Kevin\Sequence_tap\EphysE8_output\Results-EphysE8-Rat45\ratBehUnitTrajSess';
% fnames = dir([fpath '\*.mat']);
% 
% unitCount = [];
% for s = 1:length(fnames)
%     a = fnames(s).name;
%     a = strsplit(a,'.');
%     ss = str2num(a{1});
%     load(fullfile(fpath, fnames(s).name));
%     disp(ss);
%     
%     if ~isempty(ratBehUnitTrajStructSess.spikeTimes)
%     unitCount(ss) = sum(cellfun(@length,ratBehUnitTrajStructSess.spikeTimes) > 10);
%     end
%     
% end

function fitDecodingModels_unitSets_cluster_KM_parallel(sss)

%% paths and files to load
addpath(genpath('/n/holylfs02/LABS/olveczky_lab/Kevin/Analysis/encodingModels'))
addpath('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/Utilities');
addpath('../HelperFunctions');

% addpath(genpath('D:\Kevin\Sequence_tap_analysis\basic_beh_ephys_analysis\encodingModels'))
% 
% cd('D:\Kevin\Sequence_tap_analysis\basic_beh_ephys_analysis\decodingModels');
% % good session, lots of units
%  load('D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat47\ratBehUnitTrajSess_quality\84.mat')
% % lots of units, OT only
%  load('D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat47\ratBehUnitTrajSess_quality\87.mat')
% 
% 
%  load('D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat47\ratBEHstruct.mat')
 
%% loop over sessions


% % 
% files = dir(['D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat47\ratBehUnitTrajSess_quality\*.mat']);
% load('D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat47\unitCount.mat')
% load('D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat47\ratBEHstruct.mat')
% % 
% files = dir('D:\Kevin\Sequence_tap\EphysE8_output\Results-EphysE8-Rat45\ratBehUnitTrajSess\*.mat');
% load('D:\Kevin\Sequence_tap\EphysE8_output\Results-EphysE8-Rat45\unitCount.mat')
% load('D:\Kevin\Sequence_tap\EphysE8_output\Results-EphysE8-Rat45\ratBEHstruct.mat')

%ratname = 'EphysE7-Rat47';
%fpath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat47/ratBehUnitTrajSess_quality/';
%files = dir('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat47/ratBehUnitTrajSess_quality/*.mat');
%load('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat47/unitCount.mat')
%load('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat47/ratBEHstruct.mat')
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat47/decoding_paws_nose_10_20_lag_new/';
%if ~exist(savepath); mkdir(savepath); end
%loadpath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat47/day_matches.mat';
%load(loadpath);


ratname = 'EphysE7-Rat33';
files = dir('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/ratBehUnitTrajSess/*.mat');
load('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/ratBEHstruct.mat')
fpath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/ratBehUnitTrajSess/';
extrinsicpath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/';

%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_10_20_lag_new/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_scaleattempt/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_treeattempt/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_bignnattempt/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_control/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_inversescaleattempt/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_onlymod/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_control_moreunits/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_shuffled/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_uMod_nosesscombined/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_uMod_nosesscombined_SVR/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_nosesscombined_SVR/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_nosesscombined_morespikebins/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_stricterunittoss/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_stricterunittoss_nosesscombined/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_stricterunittoss_nosesscombined_SVR/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_nosesscombined/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_nosesscombined_stricterunittoss/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined_subseq/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined_shuffle/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined_shuffle_allunits/'
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined_allunits/';
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_nosesscombined_allunits_egocentric/';
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_allunits_egocentric/';

%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined_allunits_withphase/';
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined_3d/';
%savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined_3dpaws/';
savepath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/decoding_paws_nose_morespikebins_mediumunittoss_nosesscombined_3dpaws_and_nose/';

if ~exist(savepath); mkdir(savepath); end
loadpath = '/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/day_matches.mat';
load(loadpath);

% dont forget to rerun for rat 35 without seq subset!

doscale = 0;
douMod = 0;
doShuffle = 0;
doMoreSpikeBins = 1;
doStrictUnitToss = 0; %if neuron doesnt fire on 1/4 of all the trials 
useCombinedOne = 0;

doSubSeq = 0;
doMediumUnitToss = 1; % average fr must be 0.25hz across all trials

%% lag 
tuse = [-.5:.1:.5]; % seconds lag
tuse = tuse/(1/40);
unitsPerGrp = [4,8,12,15]; % ??
%unitsPerGrp = 2; % for egocentric only

%unitsPerGrp = [15,12,8,4];
unitsPerGrp = [12]; % run quicker for check
%unitsPerGrp = [20];
%% variables to save
% R2_pos_paws = cell(size(files));
% R2_vel_paws = cell(size(files)); % test x train
% R2_pos_nose = cell(size(files));
% R2_vel_nose = cell(size(files));
% num_units = cell(size(files));% track units using
% R2_phase = cell(size(files));
% R2_taps = cell(size(files));
% R2_vel_paws_lag = cell(size(files,1), length(tuse)); % 11 hard coded...
% R2_paws_lag = cell(size(files,1),length(tuse));
% y_pred_pos = cell(size(files,1),2);  % save y_predictions for a figure?
% y_pred_vel = cell(size(files,1),2);
% R2_pos = cell(size(files));
% R2_vel = cell(size(files));
% 
% R2_pos_paws20 = cell(size(files));
% R2_vel_paws20 = cell(size(files)); % test x train
% R2_pos_nose20 = cell(size(files));
% R2_vel_nose20 = cell(size(files));
% num_units20 = cell(size(files));% track units using
% R2_phase20 = cell(size(files));
% R2_taps20 = cell(size(files));
% R2_vel_paws_lag20 = cell(size(files,1),length(tuse));
% R2_paws_lag20 = cell(size(files,1),length(tuse));
% y_pred_pos20 = cell(size(files,1),2);  % save y_predictions for a figure?
% y_pred_vel20 = cell(size(files,1),2);

R2_ego = cell(size(files,1));

R2_pos_paws = cell(size(files,1),length(unitsPerGrp));
R2_vel_paws = cell(size(files,1),length(unitsPerGrp)); % test x train
R2_pos_nose = cell(size(files,1),length(unitsPerGrp));
R2_vel_nose = cell(size(files,1),length(unitsPerGrp));
num_units = cell(size(files,1),length(unitsPerGrp));% track units using
R2_phase = cell(size(files,1),length(unitsPerGrp));
R2_taps = cell(size(files,1),length(unitsPerGrp));
R2_vel_paws_lag = cell(size(files,1), length(tuse),length(unitsPerGrp)); % 11 hard coded...
R2_paws_lag = cell(size(files,1),length(tuse),length(unitsPerGrp));
y_pred_pos = cell(size(files,1),2,length(unitsPerGrp));  % save y_predictions for a figure?
y_pred_vel = cell(size(files,1),2,length(unitsPerGrp));
R2_pos = cell(size(files,1),length(unitsPerGrp));
R2_vel = cell(size(files,1),length(unitsPerGrp));
R2_vel_lag = cell(size(files,1),length(tuse),length(unitsPerGrp));
R2_lag = cell(size(files,1),length(tuse),length(unitsPerGrp));


%% decoding params

%unitsPerGrp = [1,4,7,11,14];
%unitsPerGrp = 10; % only want this for now...
grpReps = 20;

%unitsPerGrp = [4,8,12,15]; % ??


% parallelize by splitting into sessions?
%%
% %for sss = 1:length(files)
% %     sess = files(sss).name;
% %     sess = strsplit(sess,'.'); sess = str2num(sess{1});
% %     try
% %         load(fullfile(files(sss).folder,files(sss).name));
% %     catch
% %         return;
% %    end
%     % other way loading
%     sess = sss;
%     try
%         load([fpath num2str(sess) '.mat']);
%     catch
%         return
%     end
%     
%     
%     disp(sess)
%     if ratBehUnitTrajStructSess.protocol==8; return; end
%     if unitCount(sess) < 5; return; end
%     if isempty(ratBehUnitTrajStructSess(1).trajMaster); return; end
%     if isempty(ratBehUnitTrajStructSess(1).trajSlave1); return; end
%     if isempty(ratBehUnitTrajStructSess(1).trajSlave2); return; end

%% load sessCat
load(loadpath); % load day matches
% load the combined one
if useCombinedOne
disp('using combined one!')
%load('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat47/day_matches_merge.mat'); % skip OT sessions!
load('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/day_matches_merge.mat');
sessMatch = sessCombined;
end
disp(length(sessMatch))
sessCat = sessMatch{sss}; % this is the dayiter
disp(sessCat);
%unitCat = unitMatch;
% load first session?
load([fpath num2str(sessCat(1)) '.mat'])

%% get unitMatch
unitMatch = [];
for s = sessCat
    %load(['D:\Kevin\Sequence_tap\EphysE7_output\Results-EphysE7-Rat47\ratBehUnitTrajSess_quality\' num2str(s) '.mat'])
   % load(['D:\Kevin\Sequence_tap\EphysE8_output\Results-EphysE8-Rat45\ratBehUnitTrajSess\' num2str(s) '.mat'])
    load([fpath num2str(s) '.mat']);
   % if ratBehUnitTrajStructSess.protocol==8; continue; end
    if isempty(ratBehUnitTrajStructSess.unitLabel); disp(s); end
    if isempty(unitMatch)
        unitMatch = ratBehUnitTrajStructSess.unitLabel;
    else
        unitMatch = intersect(ratBehUnitTrajStructSess.unitLabel, unitMatch);
    end
end

disp(unitMatch);
unitCat = unitMatch;

%% uMod
if douMod
load('/n/holylfs02/LABS/olveczky_lab/Kevin/Video/Tracking/EphysE7-Rat33/uMod_zscore_1.mat');
unitCat = intersect(unitCat, uMod);
disp('only modulated units')
disp(unitCat)

end
 
%% organize data?


offset = .75; % seconds before and after! (was at 2.5, why so large?)
binsize = 25; % milliseconds
offset = .1; % crop closer to taps
% set binisize to sampling rate of camera? (40 hz = 25 millisecond frames)
off = offset * 40; % 40 hz

thresh = 0.95; % probability cutoff of tracking

% asheshs stuff
sm_bin = .04;
bin = .025;


edges = -offset*1000 : binsize : offset*1000;

assert(strcmp(ratBEHstruct(1).name,ratBehUnitTrajStructSess(1).name))

% count for trials
count = 1;

% load structs until has one without empty trajMaster?
val = 1;
while isempty(ratBehUnitTrajStructSess(1).trajMaster)
load([fpath num2str(val) '.mat']);
val = val+1;
end

% take everything
joints1 = fieldnames(ratBehUnitTrajStructSess(1).trajMaster)';
cams1 = cell(1,length(joints1)); cams1(:) = {'Master'};
joints2 = fieldnames(ratBehUnitTrajStructSess(1).trajSlave1)';
cams2 = cell(1,length(joints2)); cams2(:) = {'Slave1'};
joints3 = fieldnames(ratBehUnitTrajStructSess(1).trajSlave2)';
cams3 = cell(1,length(joints3)); cams3(:) = {'Slave2'};
cams = [cams1 cams2 cams3];
joints = [joints1 joints2 joints3];
dimension1 = cell(1,length(joints)); dimension1(:)={'X'};
dimension2 = cell(1,length(joints)); dimension2(:)={'Y'};
dimension = [dimension1 dimension2];



miscount = 0;
count = 1;
s = 1;
% pull out spike times?
uUse = find(cellfun(@length, ratBehUnitTrajStructSess.spikeTimes) > 10);

fr = {};
frAshesh = {};
paws = {};
badtrackfrac = [];
spTimesAll = {}; % use FR or spike times??
spTimesReal = {};
tapTimes = {};
tapTimesVid = {};
unitID = [];
unitType = [];
% better to just make a giant list (independent of this?)
Hit = []; % binary if hit
lever = {}; % nx1 of levers in trial
target = {}; % 3x1 of target sequence
move = {};
protocol = [];
HitVal = {};
WM = [];
sessID = [];
trialID = [];
type = [];gain = -.9:.001:1.1;
trialTimes = {}; % array of time in seconds aligned to first tap

vidNames = {};
vidFrames = {};

for sess = sessCat
    % load day matches 
    disp(['sess: ' num2str(sess) ' of ' num2str(sessCat(end))]);

    try
        load([fpath num2str(sess) '.mat']);
    catch
        continue; %return
    end
    
   % if ratBehUnitTrajStructSess.protocol==8; return; end
 %   if unitCount(sess) < 5; return; end
    if isempty(ratBehUnitTrajStructSess(1).trajMaster); continue; end
    if isempty(ratBehUnitTrajStructSess(1).trajSlave1); continue; end
    if isempty(ratBehUnitTrajStructSess(1).trajSlave2); continue; end
    
    % pull out unit matches
    uid = ismember(ratBehUnitTrajStructSess.unitLabel,unitCat);
    ratBehUnitTrajStructSess.unitLabel(~uid) = [];
    ratBehUnitTrajStructSess.unitType(~uid) = [];
    ratBehUnitTrajStructSess.spikeTimes(~uid) = [];
   % assert(issorted( ratBehUnitTrajStructSess.unitLabel ))
    
    uUse = find(cellfun(@length, ratBehUnitTrajStructSess.spikeTimes) > 10);
    uUse = 1:sum(uid);


for trial = 1:length(ratBehUnitTrajStructSess(s).Hit)

    % skip if not a hit?
    if ratBehUnitTrajStructSess(s).Hit(trial) == 0 ; continue; end

    ports = ratBehUnitTrajStructSess(s).pokeNames{trial};
    if isempty(ports) ; continue; end


    for tap = 1%:length(ports) % for multiple levers?

        % get time range of session to crop spikes by
        tapTime = ratBehUnitTrajStructSess(s).pokeTimes{trial}((tap));
        tapEnd = ratBehUnitTrajStructSess(s).pokeTimes{trial}((tap) + length(ports) - 1);
        allTaps = ratBehUnitTrajStructSess(s).pokeTimes{trial}(((tap)):((tap)+length(ports)-1));
        trange = (tapTime - offset*1000) : (tapEnd + offset*1000);
        edges = (tapTime - offset*1000) : binsize : (tapEnd + offset*1000 + binsize); % get one more bin here


        % calculate scale factor for tapTime?

        % get spikes
        

       %** 20191022
       % temp skip if already tracked

       % get traj info
       % skip if already tracked, for not caring about units
%           temp = [trialID;sessID];
%           if ~isempty(temp); if sum(ismember(temp',[trial,s],'rows'))>=1; continue; end ; end
       %* at end check uniqueness
       %* x = [trialID;sessID];
       %* y = unique(x','rows'); % y should equal x

       %** todo: add in tracking prob!

       skippedflag = 0;
       checkframesstart = []; checkframesstop = [];
        for cid = 1:length(joints)
        if isempty(ratBehUnitTrajStructSess(s).(['pokeFrames' cams{cid}])); continue; end

        framesstart = ratBehUnitTrajStructSess(s).(['pokeFrames' cams{cid}]){trial};
        framesstop = ratBehUnitTrajStructSess(s).(['pokeOutFrames' cams{cid}]){trial};

        checkframesstart(cid) = framesstart(1); checkframesstop(cid) = framesstop(end);

        if length(unique(framesstart))~= length(framesstart); continue ; end
        if any(isnan(framesstart)) | any(isnan(framesstop)); continue ; end

        if ~isfield(ratBehUnitTrajStructSess(s).(['traj' cams{cid}]),(joints{cid})); continue; end
        if ~isfield(ratBehUnitTrajStructSess(s).(['traj' cams{cid}]),(joints{cid})); continue; end

        if isempty(ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid}){1}); continue; end

            % get traj
        if ( any(find(diff(framesstart)<0)) | ((framesstop(end)+off )> length(ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid}){1}))   ) ...
                & length(ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid})) > 1
            if framesstart(1)-off < 1; continue; end % shoudl rarely happen, where rat waits 10 seconds between taps...

            vid1 = ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid}){1};
            vid2 = ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid}){2};
            %** shouldn't have this continue in this loop!!!
            if isempty(vid2); continue; end
            L = length(vid1);
            paws{cid,count} = [vid1((framesstart(1)-off):end,1) ; vid2(1:mod(framesstart(end)+off,L),1)];
            paws{cid+length(joints),count} = [vid1((framesstart(1)-off):end,2) ; vid2(1:mod(framesstart(end)+off,L),2)];

            vid1prob = ratBehUnitTrajStructSess(s).(['traj_prob' cams{cid}])((trial)).(joints{cid}){1};
            vid2prob = ratBehUnitTrajStructSess(s).(['traj_prob' cams{cid}])((trial)).(joints{cid}){2};
            prob = [vid1prob((framesstart(1)-off):end) , vid2prob(1:mod(framesstart(end)+off,L))];


            % replace with nans
            paws{cid,count}(prob < thresh) = nan;
            paws{cid+length(joints),count}(prob<thresh) = nan;

            stitchind = find(diff(framesstart)<0);
            framesstart((stitchind+1):end) = framesstart((stitchind+1):end) + length(vid1);
            stitchind = find(diff(framesstop)<0);
            framesstop((stitchind+1):end) = framesstop((stitchind+1):end) + length(vid1);
            tapTimesVid{cid,count} = [framesstart - framesstart(1)+off];

        % goes into next video but doesn't have trajectory (case due to offset)
        elseif ( any(find(diff(framesstart)<0)) | ((framesstop(end)+off )> length(ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid}){1}))   ) ...
                & length(ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid})) == 1
            1+1; continue;
        % goes into previous video (works??)
        elseif framesstart(1)-off < 1 && trial==1
            continue;
        elseif framesstart(1)-off < 1 && trial>1
            if isempty(ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial-1)).(joints{cid})); continue; end
            vid1 = ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial-1)).(joints{cid}){end};
            vid2 = ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid}){1};
            if isempty(vid1) | isempty(vid2); continue; end
            L = length(vid1);
            paws{cid,count} = [vid1((end-(off-framesstart(1))):end,1) ; vid2(1:mod(framesstart(end)+off,L),1)];
            paws{cid+length(joints),count} = [vid1((end-(off-framesstart(1))):end,2) ; vid2(1:mod(framesstart(end)+off,L),2)];
            % replace with nans bad tracking
            vid1prob = ratBehUnitTrajStructSess(s).(['traj_prob' cams{cid}])((trial-1)).(joints{cid}){end};
            vid2prob = ratBehUnitTrajStructSess(s).(['traj_prob' cams{cid}])((trial)).(joints{cid}){1};
            prob = [vid1prob((end-(off-framesstart(1))):end) , vid2prob(1:mod(framesstart(end)+off,L))];

            paws{cid,count}(prob < thresh) = nan;
            paws{cid+length(joints),count}(prob<thresh) = nan;

            tapTimesVid{cid,count} = [framesstart - framesstart(1)+off];
        else
            vid1 = ratBehUnitTrajStructSess(s).(['traj' cams{cid}])((trial)).(joints{cid}){1};
            paws{cid,count} = vid1((framesstart(1)-off):(framesstart(end)+off),1);
            paws{cid+length(joints),count} = vid1((framesstart(1)-off):(framesstart(end)+off),2);

            vid1prob = ratBehUnitTrajStructSess(s).(['traj_prob' cams{cid}])((trial)).(joints{cid}){1};
            prob = vid1prob((framesstart(1)-off):(framesstart(end)+off));

            paws{cid,count}(prob<thresh) = nan;
            paws{cid+length(joints),count}(prob<thresh) = nan;

            tapTimesVid{cid,count} = [framesstart - framesstart(1)+off];
        end

        if any(prob>1); 
            disp('here')
        end

        skippedflag = skippedflag+1;
        end % end loop over joints

       if ~skippedflag; continue; end

       if any(cellfun(@isempty, paws(:,count))); continue; end


       % check that trajectory length will be the same as sampling
       % spike times!

       % problem is frame variance in trajectories.
       % - should just interpolate to edges length?

       % step 1) fill in nans for bad probs below thresh?
       % - done above
       % - what about trials with many nans at start of sequence?
       % how to interpolate over that?
       % step 1.5) note badly marked points?
       % - if body part not well tracked, do not include in
       % analysis?
       % step 2) interpolate to length of edges

       for cid = 1:length(joints)*2
           % cound fraction bad
           badtrackfrac(cid,count) = sum(isnan(paws{cid,count})) / length(paws{cid,count});
           % fill missing
           paws{cid,count} = fillmissing(paws{cid,count}, 'linear','EndValues','nearest'); %?
           % interpolate to fixed length (change to different data
           % struct since can concatenate?
           paws{cid,count} = interp1(linspace(1,length(edges)-1,length(paws{cid,count})), paws{cid,count}, 1:(length(edges)-1)); % maybe?
       
           % smooth?
           paws{cid,count} =  imgaussfilt(paws{cid,count}, [eps, 250/1000*10]);
       end



       % * old remove
       % lengths can be different depending on exact spacing of
       % sampling. if so, add on extra bin at end of edges
%                if length(paws{1,count}) == length(edges)
%                    edges = [edges edges(end)+binsize];
%                else
%                    
%                end

        if count==50; 
            disp('here')
        end
       % get unit info
       for u = 1:length(uUse)
           spikeTimes = ratBehUnitTrajStructSess(1).spikeTimes{uUse(u)};
           spikes = spikeTimes(spikeTimes > trange(1) & spikeTimes < trange(end));
           spTimesReal{u,count} = spikes-trange(1)-offset*1000; % is this right?
           spTimesAll{u,count} = histcounts(spikes, edges); % these are different sizes???
           % pad spTimesAll to correct length?
           
           if isfield(ratBehUnitTrajStructSess(s),'unitType')
               unitType(u,count) = ratBehUnitTrajStructSess(s).unitType(uUse(u));
           end
           
           fr{u,count} = getFRfromSpikeCount(spTimesAll{u,count}, mean(diff(edges)));
           %frAshesh{u,count} = imgaussfilt(spTimesAll{u,count}, [eps, sm_bin/bin]);
           frAshesh{u,count} = imgaussfilt(spTimesAll{u,count},[eps, 250/1000*10]); 
       end



        tapTimes{count} = find(histcounts(allTaps, edges));
        unitID(count) = u;
        
        protocol(count) = ratBehUnitTrajStructSess(s).protocol;
        Hit(count)=  ratBehUnitTrajStructSess(s).Hit(trial);
        %HitVal(end+1) = ratBehUnitTrajStructSess(s).HitLeverVals{trial}(ind(tap));
        HitVal{count} = ratBehUnitTrajStructSess(s).HitLeverVals{trial};%(ind(tap):(ind(tap)+1));
        if ratBehUnitTrajStructSess(s).protocol==7
            if ~isempty(ratBehUnitTrajStructSess(s).blocknumRepair)
                WM(count) = (ratBehUnitTrajStructSess(s).blocknumRepair(trial) >= 3);
            else
                WM(count) = (ratBehUnitTrajStructSess(s).blocknum(trial)>=3);
            end
        else
            WM(count) = -1;
        end
        %IPI(end+1) = tapEnd - tapTime;
        %IPI(end+1,:) = diff(allTaps);
        sessID(count) = s;
        trialID(count) = trial;

        lever{count} = ports;
        %startTimeBEH =  extractfield(ratBEHstruct, 'startTime');
        %startTimeBEH = [startTimeBEH{:}];
        startTimeBEH = [ratBEHstruct(:).startTime];
        startTimeBTU = ratBehUnitTrajStructSess(s).startTime;
        if isfield(ratBehUnitTrajStructSess(s),'targetNames')
            target{count} = ratBehUnitTrajStructSess(s).targetNames{trial};
        else
            target{count} = ratBEHstruct(find(startTimeBTU==startTimeBEH)).targetNames{trial};
        end

        trialTimes{count} = (edges(1:end-1) - tapTime)/1000; % trial time in seconds

        % * get vid names and frames
        %* this might be pulling the wrong info???

        total_cams = unique(cams);
        for cid = 1:length(total_cams)
            if ~isempty(ratBehUnitTrajStructSess(s).(['VidNames' total_cams{cid}]))
                vidNames{cid,count} = ratBehUnitTrajStructSess(s).(['VidNames' total_cams{cid}]){trial};
                vidFrames{cid,count} = ratBehUnitTrajStructSess(s).(['pokeFrames' total_cams{cid}]){trial};
            else
                vidNames{cid,count} = [];
                vidFrames{cid,count} = [];
            end
        end
        % * 

        count = count+1;






    end % end tap


end % end trial

end % end sessCat

%% remove non firing trials

% fires 0 spikes
a =cellfun(@sum,fr);
badid = find(sum(a,2)==0);



% get rid of units that do not fire at least 25% as many spikes as there are trials
fr_floor = sum(cellfun(@sum,spTimesAll),2);
badid = find(fr_floor < size(spTimesAll,2)*.25);

% get rid of units below 0.25 hz firing rate during trial period?
% - can have needs to do in every trial
% - which feels unreasonable for the different sequences
% - or can have needs to have average over all trials?
if doMediumUnitToss
trialFR = cellfun(@(v) sum(v) / (length(v)*10/1000),spTimesAll);
badid = find(mean(trialFR,2)<0.25);
end
%  - want to toss units that don't fire in > 25 % of all trials
% or could just do 0.25 hz also

if doStrictUnitToss; % if do this, have nothing!!
fr_zero = cellfun(@sum, spTimesAll); 
val = sum(fr_zero==0,2);
badid = find(val > size(spTimesAll,2)*.25); 
end

spTimesReal(badid,:) = [];
spTimesAll(badid,:) = [];
fr(badid,:) = [];
frAshesh(badid,:) = [];
size(uUse)
uUse(badid) = [];


%% skipping session criteria...
if length(uUse)<10; return; end % new min?
if isempty(paws); return; end
if isempty(fr); return; end

%% preprocess for asheshs stuff
%** to do here, setup loop to look at multiple joints at once?

% units are spTimesAll, just bins with spike count...
% dims are trials x unit? x binning size? x time_lags
%  -  unit and bin size might be messed up
% actual input is cell trials x units

X = fr';
X = spTimesAll'; % use a smoothed version of this like ashesh did?
X = frAshesh';
X = cellfun(@transpose,X,'un',0);

% trajectories, cell array, trial x joint
% - each array in cell is the trajectory of that joint

% loop over many joints?

%cam = {'Slave1'}; % also for rat47
%joint = {'pawL'};
% 
%cam = {'Slave2'}; % for rat 47, also rat 45
%joint = {'nose'};

%joint = {'pawL','nose'}; % rat 47
%cam = {'Slave1','Slave2'};

% rat 45
%joint = {'pawL','nose'}; % left, top
%cam = {'Slave1','Slave2'};

% rat 35
%joint = {'pawR','nose'};
%cam = {'Slave1','Slave2'};

%** todo: replace for 3d!

% rat 33
joint = {'pawR','nose'};
cam = {'Master','Slave2'};


id = [];
for c = 1:length(cam)
id(end+1) = find(contains(joints,joint{c}) & contains(cams,cam{c}) );
id(end+1) = [id(end) + length(joints)]; % add in y dimension
end
% add length(joints) to get y-dim
% nose-x, or position 3, is the 


y = paws(id,:);
if size(y,2) ~= size(X,1)
    y(:,end) = []; % hacky fix
    %y(cellfun(@isempty,y)) = [];
end
% get absolute position, normalized?
try % kept crashing below!
for c = 1:size(y,1)
%     hi = max([y{c,:}]); lo = min([y{c,:}]);
%     hi = cellfun(@max,y(c,:)); hi = sort(hi,'descend'); hi = hi(5);
%     lo = cellfun(@min,y(c,:)); lo = sort(lo,'ascend'); lo = lo(5);
%     y(c,:) = cellfun(@(v) (v-lo)./(hi-lo), y(c,:),'un',0);
    % z-score instead
    u = mean([y{c,:}]); sig = std([y{c,:}]);
    y(c,:) = cellfun(@(v) (v-u)./sig, y(c,:), 'un',0);
end
% scale position based on distance to camera
if doscale
if any(contains(cam,'Master')) % looking from right, so scale to positive
    scale = -.5;
    scale = .5; % inverse
    y(1,:) = cellfun(@(v,w) v+scale.*w, y(1,:), y(3,:),'un',0);
    y(2,:) = cellfun(@(v,w) v+scale.*w, y(2,:), y(3,:),'un',0);
else
    scale = .5;
    scale = -.5; % inverse
    y(1,:) = cellfun(@(v,w) v+scale.*w, y(1,:), y(3,:),'un',0);
    y(2,:) = cellfun(@(v,w) v+scale.*w, y(2,:), y(3,:),'un',0);
end
end


y = y';
y = cellfun(@transpose,y,'un',0);
% now get velocity?
y_vel = cell(size(y));
for tr = 1:size(y,1)
    for j = 1:size(y,2);
        y_vel{tr,j} = getVelocityFromTrajectory(y{tr,j});
    end
end
catch
    return
end

%% get traj 3d!
% - could prob put this in an if statement to amke a conditional...



%add trajs

% rat 33
joint = {'nose','nose'};
cam = {'Master','Slave2'};
cam_angle = {'R','T'};


joint = {'pawR','pawR'};
cam = {'Master','Slave1'};
cam_angle = {'R','L'};

joint = {'pawR','pawR','nose','nose'};
cam = {'Master','Slave1','Master','Slave2'};
cam_angle = {'R','L','R','T'};

% pick out id
id = [];
for c = 1:length(cam)
id(end+1) = find(contains(joints,joint{c}) & contains(cams,cam{c}) );
id(end+1) = [id(end) + length(joints)]; % add in y dimension
end

% get positions from side and top
y = paws(id,:);
if size(y,2) ~= size(X,1)
    y(:,end) = []; % hacky fix
end
y = y'; % joint1 x-dim, y-dim , joint2 x-dim y-dim

% toss poorly tracked frames?
badid = any(badtrackfrac(id,:)>.2);
y(badid,:) = [];
X(badid,:) = [];
tapTimes(:,badid) = [];
trialTimes(:,badid) = [];
lever(:,badid) = [];

% load extrinsics
if strcmp(cam_angle{1},'R')
    load(fullfile(extrinsicpath, 'calib_topright','extrinsics.mat'));
    load(fullfile(extrinsicpath, 'calib_topright','camMatrix.mat'));
elseif strcmp(cam_angle{1},'L')
    load(fullfile(extrinsicpath, 'calib_topleft','extrinsics.mat'));
    load(fullfile(extrinsicpath, 'calib_topleft','camMatrix.mat'));
end
if ~strcmp(cam_angle{2},'T')
   disp('loading from sides for tracking paw')
   load(fullfile(extrinsicpath, 'calib','camMatrix.mat'))
end
%camMatrix1 = cameraMatrix(intrinsics{1},rotationMatrix{1},translationVector{1});
%camMatrix2 = cameraMatrix(intrinsics{2},rotationMatrix{2},translationVector{2});
camMatrix1
camMatrix2

% convert everything to 3d
%assert(strcmp(cam_angle{2},'T')); % check
y3d = cell(size(y,1), 3); % for one joint only...
y3d_vel = cell(size(y,1),3);
for trial = 1:size(y,1)
    val = triangulate([y{trial,1}; y{trial,2}]', [y{trial,3}; y{trial,4}]', camMatrix1, camMatrix2);
    y3d{trial,1} = val(:,1);
    y3d{trial,2} = val(:,2);
    y3d{trial,3} = val(:,3);
end
% reorient
y3d = cellfun(@transpose,y3d,'un',0);
y3d = y3d';

% for two joints
if length(joint)>2
% the case where doing both paw and nose!
% 1) load top extrinsics and sie extrinsics
load(fullfile(extrinsicpath,'calib','camMatrix.mat'))
camMatrix1_paw = camMatrix1;
camMatrix2_paw = camMatrix2;
if strcmp(cam_angle{3},'R')
    load(fullfile(extrinsicpath,'calib_topright','extrinsics.mat'))
elseif strcmp(cam_angle{3},'L')
    load(fullfile(extrinsicpath, 'calib_topleft','camMatrix.mat'));
end
camMatrix1_nose = camMatrix1;
camMatrix2_nose = camMatrix2;
% triangulate
disp('ok really traingulating for both paws and nose')
y3d = cell(size(y,1),6); % hard coding at 6 for now...
y3d_vel = cell(size(y,1),6);
for trial = 1:size(y,1)
    val = triangulate([y{trial,1}; y{trial,2}]', [y{trial,3}; y{trial,4}]', camMatrix1_paw, camMatrix2_paw);
    y3d{trial,1} = val(:,1);
    y3d{trial,2} = val(:,2);
    y3d{trial,3} = val(:,3);
    val = triangulate([y{trial,5}; y{trial,6}]', [y{trial,7}; y{trial,8}]', camMatrix1_nose, camMatrix2_nose);
    y3d{trial,4} = val(:,1);
    y3d{trial,5} = val(:,2);
    y3d{trial,6} = val(:,3);
end
% reorient
y3d = cellfun(@transpose,y3d,'un',0);
y3d = y3d';
disp(size(y3d))
end

% scale the velocities and trajs, as did above
for c = 1:size(y3d,1)
    u = mean([y3d{c,:}]); sig = std([y3d{c,:}]);
    y3d(c,:) = cellfun(@(v) (v-u)./sig, y3d(c,:), 'un',0);
end
y3d = y3d';
y3d = cellfun(@transpose,y3d,'un',0);
% now get velocity?
y3d_vel = cell(size(y3d));
for tr = 1:size(y3d,1)
    for j = 1:size(y3d,2)
        y3d_vel{tr,j} = getVelocityFromTrajectory(y3d{tr,j});
    end
end

% assign back
y = y3d;
y_vel = y3d_vel;

%** replace with velocity

%interp1(tTraj_v, diff(traj_p(:,:,xy,o),1,2)', tTraj, 'linear', 'extrap')';


% decoder parameters
%% shuffle

if doShuffle
y = cellfun(@(v) v(randperm(length(v))),y,'un',0);
y_vel = cellfun(@(v) v(randperm(length(v))),y_vel,'un',0);
end

%% add in phase
% want:
% - a phase signal of the entire sequence
% - a phase signal relative to the taps!

% entire sequence
y_phase = trialTimes'; y_phase = cellfun(@transpose, y_phase,'un',0);

% taps?
% - trial x 3 vector for time around each tap
try
y_taps = cell(size(y,1),3);
y_taps(:,1) = cellfun(@(v,w) v-v(w(1)), y_phase, tapTimes', 'un', 0);
y_taps(:,2) = cellfun(@(v,w) v-v(w(2)), y_phase, tapTimes', 'un', 0);
y_taps(:,3) = cellfun(@(v,w) v-v(w(3)), y_phase, tapTimes', 'un', 0);
useTapsFlag = 1;
catch
useTapsFlag = 0;
end

%% add in weird tap thing that I think makes no sense
X_ego = {};
for t = 1:length(tapTimes)
    X_ego{t,1} = zeros(length(X{t,1}),1);
    % label taps
    rr = arrayfun(@(v) v+(-4:4), tapTimes{t}, 'un',0);
    X_ego{t,1}([rr{:}]) = 1;
    % label moves
    X_ego{t,2} = zeros(length(X{t,1}),1);
    if length(tapTimes{t})==2 
        ss = setdiff(1:length(X{t,1}),[rr{:}]);	
        left = contains(lever{t}(1:2),{'RL','RC','CL'});
        left = left*2-1;
        X_ego{t,2}(ss) = left;
    elseif length(tapTimes{t})==3
        ss = setdiff(1:length(X{t,1}),[rr{:}]);
        jj = find(diff(ss)>1);
        
        left = contains(lever{t}(1:2),{'RL','RC','CL'});
        left = left*2-1;
        X_ego{t,2}(ss(1:jj)) = left;
        left = contains(lever{t}(2:3),{'RL','RC','CL'});
        left = left*2-1;
        X_ego{t,2}(ss((jj+1):end)) = left;
    end
end

%% control for number of units!
% - so looking across sessions, use similar # of units
% subselect 10 units? run a bunch of times? what does ashesh do...
u_grps = cell(length(unitsPerGrp),1);
for g = 1 : length(unitsPerGrp)
    if unitsPerGrp(g) > length(uUse); continue; end;
    if nchoosek(length(uUse),unitsPerGrp(g)) <= grpReps
        u_grps{g} = nchoosek(1:length(uUse),unitsPerGrp(g));
    else
        temp = zeros(grpReps,unitsPerGrp(g));
        comb_num = randperm(nchoosek(length(uUse),unitsPerGrp(g)),grpReps);
        for r = 1 : grpReps
           temp(r,:) = nchoosek_m(length(uUse),unitsPerGrp(g),comb_num(r)); 
        end
        u_grps{g} = temp;
    end    
end
u_grps(cellfun(@isempty, u_grps)) = [];

%% lag analysis
% * 20201106
% - might be easier to write my own than use asheshs

% - how to pad?

disp(size(X))
disp(size(y_vel))

tuse = [-.5:.1:.5]; % .5 seconds before to .5 seconds after
tuse = -.5; % to run faster
% convert to frames
tuse = round(tuse /(1/40));

Xlag = cell(length(tuse),1); ylag = cell(length(tuse),1); 
ylag_pos = cell(length(tuse),1);
for iii = 1:length(tuse)
    % remove last t trials of y, remove first t trials of x?
    if tuse(iii) < 0
        Xlag{iii} = cellfun(@(v) v(1:end+tuse(iii),:), X,'un',0);
        ylag{iii} = cellfun(@(v) v(-(tuse(iii)-1):end,:),y_vel,'un',0);
        ylag_pos{iii} = cellfun(@(v) v(-(tuse(iii)-1):end,:),y,'un',0);
    elseif tuse(iii) > 0
        Xlag{iii} = cellfun(@(v) v((tuse(iii)+1):end,:),X,'un',0);
        ylag{iii} = cellfun(@(v) v(1:end-tuse(iii),:),y_vel,'un',0);
        ylag_pos{iii} = cellfun(@(v) v(1:end-tuse(iii),:),y,'un',0);
    else
        Xlag{iii} = X; 
        ylag{iii} = y_vel;
        ylag_pos{iii} = y;
    end
end

%% add in more spikes like ashesh did
% add in more spikes 
% - before and after
if doMoreSpikeBins
X_bins = X;

X_plus = cellfun(@(v) v(3:end), X_bins,'un',0);
X_minus = cellfun(@(v) v(1:end-2), X_bins,'un',0);
X_bins = cellfun(@(v) v(2:end-1), X_bins, 'un',0);
X_bins = [X_minus, X_bins, X_plus];

% othe rvariabels to shorten
X = X_bins;
y = cellfun(@(v) v(2:end-1),y,'un',0);
y_vel = cellfun(@(v) v(2:end-1),y_vel,'un',0);
trialTimes = cellfun(@(v) v(2:end-1), trialTimes, 'un',0);
tapTimes = cellfun(@(v) v-1, tapTimes, 'un',0);
% fix u_grps to actually use the bins!!!!
for g = 1:length(u_grps)
    u_grps{g} = [u_grps{g}, u_grps{g}+length(uUse), u_grps{g}+2*length(uUse)];
end
y_phase = cellfun(@(v) v(2:end-1), y_phase, 'un', 0);

% also for X_ego?
X_bins = X_ego;
X_plus = cellfun(@(v) v(3:end), X_bins,'un',0);
X_minus = cellfun(@(v) v(1:end-2), X_bins,'un',0);
X_bins = cellfun(@(v) v(2:end-1), X_bins, 'un',0);
X_bins = [X_minus, X_bins, X_plus];
X_ego = X_bins;

end

disp(size(X))
disp(size(y_vel))
disp(u_grps)
%%

allmove = {'LCL','LRL','LCR','LRC', 'CLC', 'CLR','CRC','CRL','RLC','RLR','RCR','RCL'};

%% actual decoder
trSet = 1:length(y);
trainVsTest = 2; % 2 to match splits for subseq
dec = 'NN'; % NN
%dec = 'TreeBoost'; % faster to test stuff...
R2_sets = {1}; % this is set of arrays of body part groupings!
R2_wts = {1}; % array of weights for each body part
R2_sets = {[1:2],[3:4],[1:4]}; 
R2_wts = {[1,1],[1,1],[1,1,1,1]};

% for 3d trajs
R2_sets = {[1:3],[1:3],[1:3]};
R2_wts = {[1,1,1],[1,1,1],[1,1,1]};

% for 3d trajs two joints
R2_sets = {[1:3],[4:6],[1:6]};
R2_wts = {[1,1,1],[1,1,1],[1,1,1,1,1,1]};

% - need to edit to add in both paw and nose? also include y component of
% paw too?

% - out here since no units to shuffle over?
trainVsTestInd = crossvalind('kfold', length(trSet), trainVsTest);
[fitStats_ego,y_pred_ego_temp] = myDecoder_Kfold_R2wts...
           (X_ego, y_vel, trainVsTestInd, dec, R2_sets, R2_wts); 
R2_ego{sess}(:,end+1) = [fitStats_ego(3).R2_train, fitStats_ego(3).R2_test];

for uu = 1:length(u_grps)
    for uiter = 1:size(u_grps{uu},1)
        disp(['num u: ' num2str(uu) ', iteration : ' num2str(uiter)])
        %disp('looping through unit groups') 
        trainVsTestInd = crossvalind('kfold', length(trSet), trainVsTest);
       
        if doSubSeq 
	% split train vs test by sequence
        trainVsTestInd = ismember(lever, allmove(randsample(12,6)));
        end



	[fitStats,y_pred] = myDecoder_Kfold_R2wts...
                                (X(:,u_grps{uu}(uiter,:)), y, trainVsTestInd, dec, R2_sets, R2_wts); 
                                %(u_dat(:,:,n,t), kin_sets{k}, trainVsTestInd, dec, R2_sets, R2_wts);
                                
        [fitStats_vel,y_pred_vel_temp] = myDecoder_Kfold_R2wts...
                                (X(:,u_grps{uu}(uiter,:)), y_vel, trainVsTestInd, dec, R2_sets, R2_wts); 
        
        [fitStats_phase,y_pred_phase] = myDecoder_Kfold_R2wts...
                                (X(:,u_grps{uu}(uiter,:)), y_phase, trainVsTestInd, dec, {[1]}, {[1]}); 
        
% 	if useTapsFlag;                    
%         [fitStats_taps,y_pred_taps] = myDecoder_Kfold_R2wts...
%                                 (X(:,u_grps{uu}(uiter,:)), y_taps, trainVsTestInd, dec, {[1,2,3]}, {[1,1,1]});           
%         end
                    
        % store data
        R2_pos_paws{sess,uu}(:,end+1) = [fitStats(1).R2_train, fitStats(1).R2_test];
        R2_vel_paws{sess,uu}(:,end+1) = [fitStats_vel(1).R2_train, fitStats_vel(1).R2_test];
        R2_pos_nose{sess,uu}(:,end+1) = [fitStats(2).R2_train, fitStats(2).R2_test];
        R2_vel_nose{sess,uu}(:,end+1) = [fitStats_vel(2).R2_train, fitStats_vel(2).R2_test];
        num_units{sess,uu}(:,end+1) = uUse(u_grps{uu}(uiter,1:unitsPerGrp));
        R2_phase{sess,uu}(:,end+1) = [fitStats_phase(:).R2_train, fitStats_phase(:).R2_test];
        %if useTapsFlag; R2_taps{sess,uu}(:,end+1) = [fitStats_taps(:).R2_train, fitStats_taps(:).R2_test]; end
        R2_pos{sess,uu}(:,end+1) = [fitStats(3).R2_train, fitStats(3).R2_test];
        R2_vel{sess,uu}(:,end+1) = [fitStats_vel(3).R2_train, fitStats_vel(3).R2_test];
        
         % only take last iteration? no k-fold, just need for a plot
        y_pred_pos{uiter,1,uu} = y_pred;
        y_pred_pos{uiter,2,uu} = y;
        y_pred_vel{uiter,1,uu} = y_pred_vel_temp;
        y_pred_vel{uiter,2,uu} = y_vel;
        
        % store for lag
        for iii = 1:length(tuse)
		disp([tuse(iii)])
            [fitStats_lag,y_pred_lag] = myDecoder_Kfold_R2wts...
                      (Xlag{iii}(:,u_grps{uu}(uiter,1:unitsPerGrp)),ylag{iii},...
                      trainVsTestInd, dec, R2_sets, R2_wts); 
                  
            [fitStats_lag_pos,y_pred_lag_pos] = myDecoder_Kfold_R2wts...
                      (Xlag{iii}(:,u_grps{uu}(uiter,1:unitsPerGrp)),ylag_pos{iii},...
                      trainVsTestInd, dec, R2_sets, R2_wts); 

            R2_vel_paws_lag{sess,iii,uu}(:,end+1) = ...
                [fitStats_lag(1).R2_train, fitStats_lag(1).R2_test,... % paws
                 fitStats_lag(2).R2_train, fitStats_lag(2).R2_test,... % nose
                 fitStats_lag(3).R2_train, fitStats_lag(3).R2_test]; % both
             
            R2_paws_lag{sess,iii,uu}(:,end+1) = ...
                 [fitStats_lag_pos(1).R2_train, fitStats_lag_pos(1).R2_test,... % paws
                 fitStats_lag_pos(2).R2_train, fitStats_lag_pos(2).R2_test,...
                 fitStats_lag(3).R2_train, fitStats_lag(3).R2_test];
        end
        
        
    end % end uiter

   disp('temp save')

 save([savepath num2str(sss) '.mat'],...
     'num_units','R2_pos_paws','R2_vel_paws','R2_pos_nose','R2_vel_nose',...
     'R2_phase','R2_taps','R2_vel_paws_lag','R2_paws_lag',...
     'y_pred_pos','y_pred_vel','R2_pos','R2_vel');
 

end % end units

% ego decoder
% - out here since no units to shuffle over?
[fitStats_ego,y_pred_ego_temp] = myDecoder_Kfold_R2wts...
           (X_ego, y_vel, trainVsTestInd, dec, R2_sets, R2_wts); 
R2_ego{sess}(:,end+1) = [fitStats_ego(3).R2_train, fitStats_ego(3).R2_test];

 
% save([savepath num2str(sess) '.mat'],'num_units',...
%      'R2_pos_paws','R2_vel_paws','R2_pos_nose','R2_vel_nose');
%  
 
% uncell things
disp('saving for 10 units')
disp([savepath num2str(sss) '.mat'])
% R2_pos_paws = R2_pos_paws{sess};
% R2_vel_paws = R2_vel_paws{sess};
% R2_pos_nose = R2_pos_nose{sess};
% R2_vel_nose = R2_vel_nose{sess};
% num_units = num_units{sess};
% R2_phase = R2_phase{sess};
% R2_taps = R2_taps{sess};
% R2_vel_paws_lag = R2_vel_paws_lag(sess,:);
% R2_paws_lag = R2_paws_lag{sess,:};
% disp([savepath num2str(sess) '.mat'])
R2_pos_paws = R2_pos_paws(sess,:);
R2_vel_paws = R2_vel_paws(sess,:);
R2_pos_nose = R2_pos_nose(sess,:);
R2_vel_nose = R2_vel_nose(sess,:);
R2_phase = R2_phase(sess,:);
R2_taps = R2_taps(sess,:);
R2_vel_paws_lag = R2_vel_paws_lag(sess,:,:);
R2_paws_lag = R2_paws_lag(sess,:,:);
R2_pos = R2_pos(sess,:);
R2_vel = R2_vel(sess,:);
R2_ego = R2_ego(sess);

 save([savepath num2str(sss) '.mat'],...
     'num_units','R2_pos_paws','R2_vel_paws','R2_pos_nose','R2_vel_nose',...
     'R2_phase','R2_taps','R2_vel_paws_lag','R2_paws_lag',...
     'y_pred_pos','y_pred_vel','R2_pos','R2_vel', 'R2_ego');
 
return;
quit; 
%% version for 20+ units?
disp('edited so shouldnt be here');
if length(uUse)<20; return; end 

unitsPerGrp = 20;
u_grps = cell(length(unitsPerGrp),1);
for g = 1 : length(unitsPerGrp)
    if nchoosek(length(uUse),unitsPerGrp(g)) <= grpReps
        u_grps{g} = nchoosek(1:length(uUse),unitsPerGrp(g));
    else
        temp = zeros(grpReps,unitsPerGrp(g));
        comb_num = randperm(nchoosek(length(uUse),unitsPerGrp(g)),grpReps);
        for r = 1 : grpReps
           temp(r,:) = nchoosek_m(length(uUse),unitsPerGrp(g),comb_num(r)); 
        end
        u_grps{g} = temp;
    end    
end

disp('20 units now');
for uu = 1:length(u_grps)
    for uiter = 1:size(u_grps{uu},1)
        disp(uiter)
        trainVsTestInd = crossvalind('kfold', length(trSet), trainVsTest);
        [fitStats,y_pred] = myDecoder_Kfold_R2wts...
                                (X(:,u_grps{uu}(uiter,:)), y, trainVsTestInd, dec, R2_sets, R2_wts); 
                                %(u_dat(:,:,n,t), kin_sets{k}, trainVsTestInd, dec, R2_sets, R2_wts);
                                
        [fitStats_vel,y_pred_vel_temp] = myDecoder_Kfold_R2wts...
                                (X(:,u_grps{uu}(uiter,:)), y_vel, trainVsTestInd, dec, R2_sets, R2_wts); 
                            
        [fitStats_phase,y_pred_phase] = myDecoder_Kfold_R2wts...
                                (X(:,u_grps{uu}(uiter,:)), y_phase, trainVsTestInd, dec, {[1]}, {[1]}); 
                            
        [fitStats_taps,y_pred_taps] = myDecoder_Kfold_R2wts...
                                (X(:,u_grps{uu}(uiter,:)), y_taps, trainVsTestInd, dec, {[1,2,3]}, {[1,1,1]});           

                            
        % store data
       % R2_pos20{sess}(:,end+1) = [fitStats.R2_train, fitStats.R2_test];
       % R2_vel20{sess}(:,end+1) = [fitStats_vel.R2_train, fitStats_vel.R2_test];
        R2_pos_paws20{sess}(:,end+1) = [fitStats(1).R2_train, fitStats(1).R2_test];
        R2_vel_paws20{sess}(:,end+1) = [fitStats_vel(1).R2_train, fitStats_vel(1).R2_test];
        R2_pos_nose20{sess}(:,end+1) = [fitStats(2).R2_train, fitStats(2).R2_test];
        R2_vel_nose20{sess}(:,end+1) = [fitStats_vel(2).R2_train, fitStats_vel(2).R2_test];
        num_units20{sess}(:,end+1) = uUse(u_grps{uu}(uiter,:));
        R2_phase20{sess}(:,end+1) = [fitStats_phase(:).R2_train, fitStats_phase(:).R2_test];
        R2_taps20{sess}(:,end+1) = [fitStats_taps(:).R2_train, fitStats_taps(:).R2_test];
        
        
         % only take last iteration? no k-fold, just need for a plot
        y_pred_pos20{uiter,1} = y_pred;
        y_pred_pos20{uiter,2} = y;
        y_pred_vel20{uiter,1} = y_pred_vel_temp;
        y_pred_vel20{uiter,2} = y_vel;
        
        
        % store for lag
        for iii = 1:length(tuse)
            [fitStats_lag,y_pred_lag] = myDecoder_Kfold_R2wts...
                      (Xlag{iii}(:,u_grps{uu}(uiter,:)),ylag{iii},...
                      trainVsTestInd, dec, R2_sets, R2_wts); 
            [fitStats_lag_pos,y_pred_lag_pos] = myDecoder_Kfold_R2wts...
                      (Xlag{iii}(:,u_grps{uu}(uiter,:)),ylag_pos{iii},...
                      trainVsTestInd, dec, R2_sets, R2_wts); 
                  
                  
            R2_vel_paws_lag20{sess,iii}(:,end+1) = ...
                [fitStats_lag(1).R2_train, fitStats_lag(1).R2_test,... % paws
                 fitStats_lag(2).R2_train, fitStats_lag(2).R2_test]; % nose
            R2_paws_lag20{sess,iii}(:,end+1) = ...
                 [fitStats_lag_pos(1).R2_train, fitStats_lag_pos(1).R2_test,... % paws
                 fitStats_lag_pos(2).R2_train, fitStats_lag_pos(2).R2_test];
             
        end
        
    end
end
                    
%% save stuff
% R2_pos(sess,:) = [fitStats.R2_train, fitStats.R2_test];
% R2_vel(sess,:) = [fitStats_vel.R2_train, fitStats_vel.R2_test]; % test x train
% num_units{sess} = uUse;% track units using

% uncell things
R2_pos_paws20 = R2_pos_paws20{sess};
R2_vel_paws20 = R2_vel_paws20{sess};
R2_pos_nose20 = R2_pos_nose20{sess};
R2_vel_nose20 = R2_vel_nose20{sess};
num_units20 = num_units20{sess};
R2_phase20 = R2_phase20{sess};
R2_taps20 = R2_taps20{sess};
R2_vel_paws_lag20 = R2_vel_paws_lag20(sess,:);
R2_paws_lag20 = R2_paws_lag20(sess,:);


 
%save([savepath num2str(sess) '.mat'], 'R2_pos_paws20','R2_vel_paws20','R2_pos_nose20','R2_vel_nose20','num_units',...
%     'R2_pos_paws','R2_vel_paws','R2_pos_nose','R2_vel_nose','num_units20');
 
 
save([savepath num2str(sess) '.mat'],...
     'num_units','R2_pos_paws','R2_vel_paws','R2_pos_nose','R2_vel_nose',...
     'R2_phase','R2_taps','R2_vel_paws_lag','R2_paws_lag',...
     'y_pred_pos','y_pred_vel',...
     'num_units20','R2_pos_paws20','R2_vel_paws20','R2_pos_nose20','R2_vel_nose20',...
     'R2_phase20','R2_taps20','R2_vel_paws_lag20','R2_paws_lag20',...
     'y_pred_pos20','y_pred_vel20');
                    
%end % end loop over sess...
