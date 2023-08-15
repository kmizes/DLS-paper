%% encoding models plot?
%  - i never saved what i used before to make the bar plots i think...

% set up paths


% - is this data for paws or velocity?

%% pull data


nn_vals = [];
xgb_vals = [];
unitvals = [];

fpath = 'X:\Kevin\Analysis\Data\Rat33';
fpath = 'X:\Kevin\Analysis\Data\Rat45';
fpath = 'X:\Kevin\Analysis\Data\Rat47';
%fpath = 'X:\Kevin\Analysis\Data\Rat56';
fpath = 'X:\Kevin\Analysis\Data\Rat67';

%* for MC rats?


fname = dir(fullfile(fpath, '*context.mat'));
for f = 1:length(fname)
    load(fullfile(fpath, fname(f).name));
    u = strsplit(fname(f).name,'_');
    unitvals(f) = str2num(u{1});
    
    nn_vals(f,:) = [nn_cued, nn_ot, nn_test, nn_train];
    xgb_vals(f,:) = [xgboost_cued, xgboost_ot, xgboost_test, xgboost_train];
    % load for many iterations of shuffled dataset?
end

% load all the previous encoding model results too...see if stuff is
% consistent?
lstm_cv=[];
nn_cv = [];
xgb_cv = [];
xgbhist_cv = [];
unit_cv = [];
fpath = 'X:\Kevin\Analysis\Data\Rat33';
fpath = 'X:\Kevin\Analysis\Data\Rat45';
fpath = 'X:\Kevin\Analysis\Data\Rat47';

fname = dir(fullfile(fpath, '*psr2.mat'));
for f = 1:length(fname)
    load(fullfile(fpath, fname(f).name));
    u = strsplit(fname(f).name,'_');
    unit_cv(f) = str2num(u{1});
    nn_cv(f,:) = nn;
    xgb_cv(f,:) = xgboost;
    xgbhist_cv(f,:) = xgboost_hist;
    lstm_cv(f,:) = lstm;
end

%% pull data for 50 iterations
nn_vals = [];
xgb_vals = [];
unitvals = [];

%fpath = 'X:\Kevin\Analysis\Data\Rat33';
%fpath = 'X:\Kevin\Analysis\Data\Rat47';
%fpath = 'X:\Kevin\Analysis\Data\Rat45';
fpath = 'X:\Kevin\Analysis\Data\Rat35';

fpath = 'X:\Kevin\Analysis\Data\Rat35_NOTQUALITY';


%fpath = 'X:\Kevin\Analysis\Data\Rat67';

fname = dir(fullfile(fpath, '*context.mat'));
for f = 1:length(fname)
    load(fullfile(fpath, fname(f).name));
    u = strsplit(fname(f).name,'_');
    unitvals(f) = str2num(u{1});
    
    nn_vals(f,:) = [mean(nn_cued), mean(nn_ot), mean(nn_test), mean(nn_train)];
    xgb_vals(f,:) = [mean(xgboost_cued), mean(xgboost_ot), mean(xgboost_test),mean(xgboost_train)];
    % load for many iterations of shuffled dataset?
end

fname = dir(fullfile(fpath, '*context_vel_lag.mat')); % pull velocity encoding data
fname = dir(fullfile(fpath, '*context_vel.mat')); % pull velocity encoding data
%fname = dir(fullfile(fpath, '*context_vel_lag_new.mat')); % pull velocity encoding data
% - nope doing the time warp thing screws something up with the data...

nn_vals_vel = [];
xgb_vals_vel = [];
usave = [];

for f = 1:length(fname)
    load(fullfile(fpath, fname(f).name));
    u = strsplit(fname(f).name,'_'); u = str2num(u{1});
    
    %usave(f) = u;
    nn_vals_vel(f,:) = [mean(nn_cued), mean(nn_ot), mean(nn_test), mean(nn_train)];
    xgb_vals_vel(f,:) = [mean(xgboost_cued), mean(xgboost_ot), mean(xgboost_test),mean(xgboost_train)];
    % load for many iterations of shuffled dataset?
end

%% plot
% one rat, both xgb and nn
% plot all
edges = -1:.05:1;
figure; hold on; for j = 1:4; subplot(2,2,j);  histogram(nn_vals_vel(:,j),edges,'DisplayStyle','stairs'); end



%
edges = -.5:.05:1;
figure;histogram(xgb_vals(:,3),edges)
hold on;histogram(xgb_vals(:,4),edges)
hold on;histogram(xgb_vals(:,1),edges)

%* useful??
figure; hold on;
b = bar(1:3, mean(xgb_vals(:,1:3))); b.FaceAlpha=.4;
errorbar(1:3, mean(xgb_vals(:,1:3)),std(xgb_vals(:,1:3))/sqrt(size(xgb_vals,1)),'.')

% maybe look at difference (should be around 0)
 % or filter out those that dont fit as well?

id = find(xgb_vals(:,3)>.05);
figure; hold on;
plot(xgb_vals(id,3) - xgb_vals(id,1))

%** of the units that had psr2 explainatory power, which of those did it
%matter if were cued or ot?
% so redo, for a better shuffled version (cv partition? or use from
% earlier?) of what is an explainatory neuron
% - then use those in the analysis?

cutoff = -1;
id3 = find(xgb_vals(:,3)>cutoff & xgb_vals(:,4) > cutoff);

% bar plot
figure; hold on;
b = bar(1:3,mean(xgb_vals(id3,1:3))); b.FaceAlpha=.4;
errorbar(1:3,mean(xgb_vals(id3,1:3)),std(xgb_vals(id3,1:3))/sqrt(length(id3)),'k.')
xticks(1:3)
xticklabels({'train ot','train cued','train both'})
ylabel('pseudo r^2')
title('Relevant units (psr2 > 0)')

%% plot vels next to position
%** for some reason time lag screws up both?



%
figure; hold on;

% cutoff = -0;
% id3 = find(nn_vals(:,3)>cutoff & nn_vals(:,4) > cutoff);
% b = bar(1:3,mean(nn_vals(id3,1:3))); b.FaceAlpha=.4;
% errorbar(1:3,mean(nn_vals(id3,1:3)),std(nn_vals(id3,1:3))/sqrt(length(id3)),'k.')
% 
% cutoff = -0;
% id3 = find(nn_vals_vel(:,3)>cutoff & nn_vals_vel(:,4) > cutoff);
% b = bar(5:7,mean(nn_vals_vel(id3,1:3))); b.FaceAlpha=.4;
% errorbar(5:7,mean(nn_vals_vel(id3,1:3)),std(nn_vals_vel(id3,1:3))/sqrt(length(id3)),'k.')

% version with all?

cutoff = -0;
id3 = find( nn_vals(:,4) > cutoff);
b = bar(1:3,mean(nn_vals(id3,1:3))); b.FaceAlpha=.4;
errorbar(1:3,mean(nn_vals(id3,1:3)),std(nn_vals(id3,1:3))/sqrt(length(id3)),'k.')

cutoff = -0;
id3 = find( nn_vals_vel(:,4) > cutoff);
b = bar(5:7,mean(nn_vals_vel(id3,1:3))); b.FaceAlpha=.4;
errorbar(5:7,mean(nn_vals_vel(id3,1:3)),std(nn_vals_vel(id3,1:3))/sqrt(length(id3)),'k.')


%% DLS matlab script
fpath = 'X:\Kevin\Analysis\Data\Rat33';

files = dir(fullfile(fpath, '*matlab*.mat'));



R2_cued = []; R2_ot = []; R2_both = [];
for f = 1:length(files)
    load(fullfile(files(f).folder, files(f).name));
    R2_cued(f) = R2_cued_train_ot_test(2);
    R2_ot(f) = R2_ot_train_cued_test(2);
    R2_both(f) = R2_both_test(2);
end

figure; hold on;
edges = -.5:.05:1;
histogram(R2_cued,edges)
histogram(R2_ot,edges)
histogram(R2_both,edges)

%% data for sequences
% - multi rat?

%fpath = 'X:\Kevin\Analysis\Data\Rat33';
fpath = 'X:\Kevin\Analysis\Data\Rat47';
fpath = 'X:\Kevin\Analysis\Data\Rat45';
%fpath = 'X:\Kevin\Analysis\Data\Rat35';

fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat45'};
%fpaths_all = {'X:\Kevin\Analysis\Data\Rat72';'X:\Kevin\Analysis\Data\Rat45'};
fpaths_all = {'X:\Kevin\Analysis\Data\Rat33';'X:\Kevin\Analysis\Data\Rat47'};

fpaths_all = {'X:\Kevin\Analysis\Data\Rat33';'X:\Kevin\Analysis\Data\Rat47';...
    'X:\Kevin\Analysis\Data\Rat35';'X:\Kevin\Analysis\Data\Rat18'};

nn_vals_rat = {};
xgb_vals_rat = {};
for fn = 1:length(fpaths_all)
fpath = fpaths_all{fn};
fname = dir(fullfile(fpath, '*sequence_vel.mat')); % pull velocity encoding data
nn_vals = [];
xgb_vals = [];

for f = 1:length(fname)
    load(fullfile(fpath, fname(f).name));
    u = strsplit(fname(f).name,'_');
    
    nn_vals(f,:) = [mean(nn_cued), mean(nn_ot), mean(nn_test), mean(nn_train)];
    xgb_vals(f,:) = [mean(xgboost_cued), mean(xgboost_ot), mean(xgboost_test),mean(xgboost_train)];
    % load for many iterations of shuffled dataset?
end
nn_vals_rat{fn} = nn_vals;
xgb_vals_rat{fn} = xgb_vals;

end

% histogram for reference 

edges = -.5:.05:1;
figure;
for fn = 1:length(fpaths_all)
    subplot(1,2,fn); hold on;
    histogram(nn_vals_rat{fn}(:,3),edges)
    histogram(nn_vals_rat{fn}(:,4),edges)
    histogram(nn_vals_rat{fn}(:,1),edges)
end


% plot to save (bar plot)
figure; hold on;
cutoff = -1;
idx = find(xgb_vals_rat{1}(:,3) > cutoff & xgb_vals_rat{1}(:,4) > cutoff);
b1=bar([1,2], mean(xgb_vals_rat{1}(idx,[1,3])),'FaceColor',[1,.3,.3]);
errorbar(1:2, mean(xgb_vals_rat{1}(idx,[1,3])),std(xgb_vals_rat{1}(idx,[1,3]))/sqrt(length(idx)),'k.');

idx = find(xgb_vals_rat{2}(:,3) > cutoff & xgb_vals_rat{2}(:,4) > cutoff);
b2=bar([4,5], mean(xgb_vals_rat{2}(idx,[1,3])),'FaceColor',[.3,.3,1]);
errorbar(4:5, mean(xgb_vals_rat{2}(idx,[1,3])),std(xgb_vals_rat{2}(idx,[1,3]))/sqrt(length(idx)),'k.');

xticks([1,2,4,5]); xticklabels({'Subset','Full','Subset','Full'});
ylabel('pseudo r^2'); legend([b1,b2],{'DLS','MC'});
title('Neural encoding from velocity');

%% plot across DLS rats
subset_psr2 = []; full_psr2 = [];
for rat = 1:3
    cutoff = -1;
    idx = find(xgb_vals_rat{rat}(:,3) > cutoff & xgb_vals_rat{rat}(:,4) > cutoff);
    subset_psr2(rat) = mean(xgb_vals_rat{rat}(idx,[1]));
    full_psr2(rat) = mean(xgb_vals_rat{rat}(idx,[3]));
end
% * need a lag control...

%% to redo
% - redo the subset full kinematic encoding across sequences
% - im not sure what bence's kinda thinking here
%   - I feel like he's just making me redo ashesh's analyses and this paper
%   is gonna end up looking exactly the same...

%% plot subset, full, lag encoding of sequences
% * new!
fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat33';...
    'X:\Kevin\Analysis\Data\Rat35_NOTQUALITY';'X:\Kevin\Analysis\Data\Rat18';}; 



%% all DLS rats vels context subset
fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat33';...
    'X:\Kevin\Analysis\Data\Rat35';'X:\Kevin\Analysis\Data\Rat18-part1'};
fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat33';...
    'X:\Kevin\Analysis\Data\Rat35'}; % no rat 18 since somethings up...
% 
fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat33_NOTQUALITY';...
    'X:\Kevin\Analysis\Data\Rat35_NOTQUALITY'}; % no rat 18 since somethings up...

fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat33';...
     'X:\Kevin\Analysis\Data\Rat35'}; % no rat 18 since somethings up...
% fpaths_all = {'X:\Kevin\Analysis\Data\Rat47'};
% 
% % MC rats
% %   % - need to rerun rat 67 without lag
% fpaths_all = {'X:\Kevin\Analysis\Data\Rat45';...
%     'X:\Kevin\Analysis\Data\Rat25';...
%     'X:\Kevin\Analysis\Data\Rat63';...
%    % 'X:\Kevin\Analysis\Data\Rat67';... % these are from anterior MC
%     'X:\Kevin\Analysis\Data\Rat72'};

% - for context
%fpaths_all = {'X:\Kevin\Analysis\Data\Rat45'};

% files_to_pull = '*context_vel.mat'; % lag?
% - so have few units for quality, more for not quality (a lot more), but
% r2 not as good, though mean still ok? (prob still around ashesh+steffen)
% - redo rat 47 for not quality?
% - also! need to put a limit on num trials for encoding (if only matching
% with <10 trials, probably shouldnt use!)

% - for subseq
 files_to_pull = '*sequence_vel.mat';
    
nn_vals_vel = cell(1,length(fpaths_all));
xgb_vals_vel = cell(1,length(fpaths_all));
unit_use = cell(1,length(fpaths_all));
for rat = 1:length(fpaths_all)
    fpath = fpaths_all{rat};
   % if contains(fpath, 'Rat67'); files_to_pull = '*context_vel_lag.mat'; else; files_to_pull = '*context_vel.mat'; end
   % fname = dir(fullfile(fpath, '*context_vel.mat')); % pull velocity encoding data
    fname = dir(fullfile(fpath, files_to_pull)); % pull velocity encoding data
    if isempty(fname); disp(fpath); end;
    
    for f = 1:length(fname)
        load(fullfile(fpath, fname(f).name));
        u = strsplit(fname(f).name,'_');
        if contains(u{1}, 'Rat'); u{1} = u{1}(6:end); end
        unit_use{rat}(f) = str2num(u{1});
        
        nn_vals_vel{rat}(f,:) = [mean(nn_cued), mean(nn_ot), mean(nn_test), mean(nn_train)];
        xgb_vals_vel{rat}(f,:) = [mean(xgboost_cued), mean(xgboost_ot), mean(xgboost_test),mean(xgboost_train)];
        % load for many iterations of shuffled dataset?
    end

end

%% plots across rats?
% xgb_vals_vel_MC = xgb_vals_vel;
% nn_vals_vel_MC = nn_vals_vel;
% xgb_vals_vel_DLS = xgb_vals_vel;
% nn_vals_vel_DLS = nn_vals_vel;

% - plot to make:
%  - how much of neuronal activity can be explained by behavior
%  - how much of neuronal activity can generalize across contexts

% histogram of across?
valMC = cell2mat(cellfun(@(v) v(:,1),xgb_vals_vel_MC,'un',0)');
valDLS = cell2mat(cellfun(@(v) v(:,1),xgb_vals_vel_DLS,'un',0)');
% 
% valMC(valMC<0)=0;

figure; hold on;
edges = -.5:.05:1;
histogram(valMC,edges,'Normalization','probability'); 
histogram(valDLS,edges,'Normalization','probability'); 
legend('MC','DLS')


% plot togheter
figure; hold on;

r2_vals = [];
for f = 1:length(xgb_vals_vel_MC)
    r2_vals = [r2_vals; xgb_vals_vel_MC{f}(:,1:3)];
end

cued_vel = cellfun(@(v) v(:,1)',xgb_vals_vel_MC,'un',0); cued_vel = [cued_vel{:}];
ot_vel = cellfun(@(v) v(:,2)',xgb_vals_vel_MC,'un',0); ot_vel = [ot_vel{:}];
both_vel = cellfun(@(v) v(:,3)',xgb_vals_vel_MC,'un',0); both_vel = [both_vel{:}];
plotSpread([cued_vel;ot_vel;both_vel]','categoryIdx',repelem([0,1,2],length(cued_vel)),'categoryColors',{'k','k','k'})
bar(1,mean(r2_vals(:,1)),'FaceColor',[1,.3,.3]);
bar(2,mean(r2_vals(:,2)),'FaceColor',[.3,1,.3]);
bar(3,mean(r2_vals(:,3)),'FaceColor',[1,1,.3]);

r2_vals = [];
for f = 1:length(xgb_vals_vel_DLS)
    r2_vals = [r2_vals; xgb_vals_vel_DLS{f}(:,1:3)];
end

cued_vel = cellfun(@(v) v(:,1)',xgb_vals_vel_DLS,'un',0); cued_vel = [cued_vel{:}];
ot_vel = cellfun(@(v) v(:,2)',xgb_vals_vel_DLS,'un',0); ot_vel = [ot_vel{:}];
both_vel = cellfun(@(v) v(:,3)',xgb_vals_vel_DLS,'un',0); both_vel = [both_vel{:}];
plotSpread([cued_vel;ot_vel;both_vel]','categoryIdx',repelem([0,1,2],length(cued_vel)),'categoryColors',{'k','k','k'},'xvalues',[5,6,7])
bar(5,mean(r2_vals(:,1)),'FaceColor',[1,.3,.3]);
bar(6,mean(r2_vals(:,2)),'FaceColor',[.3,1,.3]);
bar(7,mean(r2_vals(:,3)),'FaceColor',[1,1,.3]);

%% find bad units
% rat 35
 % - only one in top 5 i buy is #1, others look ok?
%  - compare to matlab code version
 
rat = 3;

[~,id] = sort(xgb_vals_vel{rat}(:,3),'ascend');
worseunits = unit_use{rat}(id);

xgb_vals_vel{rat}(id(1:5),3)
worseunits(1:5)

%% plot?
figure; hold on;
for f = 1:length(fpaths_all)
    subplot(2,3,f); hold on;
    edges = -.5:.05:1;
%     histogram(xgb_vals_vel{f}(:,3),edges)
%     histogram(xgb_vals_vel{f}(:,4),edges)
%     histogram(xgb_vals_vel{f}(:,1),edges)
    histogram(nn_vals_vel{f}(:,3),edges)
    histogram(nn_vals_vel{f}(:,4),edges)
    histogram(nn_vals_vel{f}(:,1),edges)
end

% plot mean for each rat and distribution as a whole?
cutoff = -1;
cutoff = 0;
r2_vals = []; rat_means = []; 
figure; hold on; colorsuse = 'rbgmc';
for f = 1:length(fpaths_all)
%     id3 = find(nn_vals_vel{f}(:,3)>cutoff & nn_vals_vel{f}(:,4) > cutoff);
%     r2_vals = [r2_vals; nn_vals_vel{f}(id3,1:3)];
%     rat_means(f,:) = mean(nn_vals_vel{f}(id3,1:3));
%     errorbar(1:3,mean(nn_vals_vel{f}(id3,1:3)),std(nn_vals_vel{f}(id3,1:3))/sqrt(length(id3)),'.','Color',colorsuse(f));
    
    id3 = find(xgb_vals_vel{f}(:,3)>cutoff & xgb_vals_vel{f}(:,4) > cutoff);
    r2_vals = [r2_vals; xgb_vals_vel{f}(id3,1:3)];
    rat_means(f,:) = mean(xgb_vals_vel{f}(id3,1:3));
    errorbar(1:3,mean(xgb_vals_vel{f}(id3,1:3)),std(xgb_vals_vel{f}(id3,1:3))/sqrt(length(id3)),'.','Color',colorsuse(f));


end

% use distribtuion plot instead for all avaerages?
cued_vel = cellfun(@(v) v(:,1)',xgb_vals_vel,'un',0); cued_vel = [cued_vel{:}]; 
ot_vel = cellfun(@(v) v(:,2)',xgb_vals_vel,'un',0); ot_vel = [ot_vel{:}];
both_vel = cellfun(@(v) v(:,3)',xgb_vals_vel,'un',0); both_vel = [both_vel{:}];
% thresh
cued_vel(cued_vel<cutoff)=[]; ot_vel(ot_vel<cutoff) = []; both_vel(both_vel<cutoff)=[];
figure; hold on;
distributionPlot({cued_vel, ot_vel, both_vel},'showMM',4);

figure;distributionPlot(r2_vals,'histOpt',0,'showMM',4)

% swarm plot maybe better?
% - histgoram for now out of laziness
figure; hold on;
edges = -.5:.05:1;
histogram(cued_vel, edges,'FaceColor',[1,.3,.3],'EdgeColor','none');
histogram(ot_vel, edges, 'FaceColor',[.3,1,.3],'EdgeColor','none');
histogram(both_vel,edges,'FaceColor',[.8,.8,.3],'EdgeColor','none');
% - yeah this needs to be a swarm plot to tell them apart...
figure; hold on;
plotSpread([cued_vel;ot_vel;both_vel]','categoryIdx',repelem([0,1,2],[length(cued_vel), length(ot_vel), length(both_vel)]),...
    'categoryColors',{'r','g','b'})

% swarm plot
figure; hold on;
swarmchart(repelem([1,2,3],[length(cued_vel), length(ot_vel), length(both_vel)]),[cued_vel,ot_vel,both_vel],'r.');
[~,pcuedot] = ttest2(cued_vel, ot_vel);
[~,pcuedboth] = ttest2(cued_vel, both_vel);
[~,potboth] = ttest2(both_vel, ot_vel);

% plot averaged all together...

if contains(files_to_pull,'context')
figure; hold on;
bar(1,mean(r2_vals(:,1)),'FaceColor',[1,.3,.3]);
bar(2,mean(r2_vals(:,2)),'FaceColor',[.3,1,.3]);
bar(3,mean(r2_vals(:,3)),'FaceColor',[1,1,.3]);

errorbar(1:3, mean(r2_vals(:,[1:3])),std(r2_vals(:,[1:3]))/sqrt(length(r2_vals)),'k.','linewidth',2);
xticks(1:3)
xticklabels({'Train Cued','Train OT','Train Both'});
% sig test
[~,p12] = ttest(r2_vals(:,1), r2_vals(:,2));
[~,p13] = ttest(r2_vals(:,1), r2_vals(:,3));
[~,p23] = ttest(r2_vals(:,2), r2_vals(:,3));
sigstar({[1,2],[1,3],[2,3]},[p12,p13,p23]);

ylabel('Pseudo R2'); title('All DLS neurons for context (encoding velocity)');

% new colors
figure; hold on;
bar(1,mean(r2_vals(:,1)),'FaceColor',[.9,.4,.4]);
bar(2,mean(r2_vals(:,2)),'FaceColor',[.4,.9,.4]);
bar(3,mean(r2_vals(:,3)),'FaceColor',[.7,.7,.7]);
    
end


% **** if using this plot in the paper, need to do hypothesis testing *** %
if contains(files_to_pull,'sequence')
figure; hold on;
bar(1,mean(r2_vals(:,1)),'FaceColor',[1,.2,.2]);
bar(2,mean(r2_vals(:,3)),'FaceColor',[1,.5,.5])
errorbar(1:2, mean(r2_vals(:,[1,3])),std(r2_vals(:,[1,3]))/sqrt(length(r2_vals)),'k.');
xticks(1:2)
xticklabels({'Subset','Full'});
ylabel('Pseudo R2'); title('All DLS neurons (encoding velocity)');

% - hypothesis testing
% - yeah not significant
[h,p,ks2stat] = kstest2(r2_vals(:,1), r2_vals(:,3));

end



% plot context for individual rats
cued_vel = cellfun(@(v) v(:,1)',xgb_vals_vel,'un',0);
ot_vel = cellfun(@(v) v(:,2)',xgb_vals_vel,'un',0); 
both_vel = cellfun(@(v) v(:,3)',xgb_vals_vel,'un',0); 
% - toss those < 0?
thresh = 0;
cued_vel_rat = cellfun(@(v) mean(v(v>thresh)),cued_vel);
ot_vel_rat = cellfun(@(v) mean(v(v>thresh)),ot_vel);
both_vel_rat = cellfun(@(v) mean(v(v>thresh)),both_vel);
figure; hold on;
bar(1:3, [mean(cued_vel_rat), mean(ot_vel_rat), mean(both_vel_rat)]);
for rat = 1:length(cued_vel_rat)
    plot(1:3, [cued_vel_rat(rat), ot_vel_rat(rat), both_vel_rat(rat)],'k');
end


%% all DLS rats vels sequence subsets


%% save for across rats
nn_valsDLS = nn_vals;
xgb_valsDLS = xgb_vals;

nn_valsMC =nn_vals;
xgb_valsMC = xgb_vals;

% velocities?
nn_valsDLS = nn_vals_vel;
xgb_valsDLS = xgb_vals_vel;


nn_valsMC =nn_vals_vel;
xgb_valsMC = xgb_vals_vel;


%% for across rats

figure; hold on;


% cutoff = -0;
% id3 = find(nn_valsDLS(:,3)>cutoff & nn_valsDLS(:,4) > cutoff);
% b = bar(1:3,mean(nn_valsDLS(id3,1:3))); b.FaceAlpha=.4;
% errorbar(1:3,mean(nn_valsDLS(id3,1:3)),std(nn_valsDLS(id3,1:3))/sqrt(length(id3)),'k.')
% 
% cutoff = -0;
% id3 = find(nn_valsMC(:,3)>cutoff & nn_valsMC(:,4) > cutoff);
% b = bar(5:7,mean(nn_valsMC(id3,1:3))); b.FaceAlpha=.4;
% errorbar(5:7,mean(nn_valsMC(id3,1:3)),std(nn_valsMC(id3,1:3))/sqrt(length(id3)),'k.')


cutoff = -0;
id3 = find(xgb_valsDLS(:,3)>cutoff & xgb_valsDLS(:,4) > cutoff);
b = bar(1:3,mean(xgb_valsDLS(id3,1:3))); b.FaceAlpha=.4;
errorbar(1:3,mean(xgb_valsDLS(id3,1:3)),std(xgb_valsDLS(id3,1:3))/sqrt(length(id3)),'k.')

cutoff = -0;
id3 = find(xgb_valsMC(:,3)>cutoff & xgb_valsMC(:,4) > cutoff);
b = bar(5:7,mean(xgb_valsMC(id3,1:3))); b.FaceAlpha=.4;
errorbar(5:7,mean(xgb_valsMC(id3,1:3)),std(xgb_valsMC(id3,1:3))/sqrt(length(id3)),'k.')







xticks([1:3, 5:7])
%xticklabels({'train ot','train cued','train both','train ot','train cued','train both'})
xticklabels({'OT','Cued','Both','OT','Cued','Both'});
ylabel('pseudo r^2')
title('Relevant units (psr2 > 0)')

%% to do
% - version with lag analysis
% as a control 
% - basically replaced train with a shuffled version
% ** i think the issue here was the lag screwed up everything
% - so have to pull from lag and not from lag...think already did?

%% pull units to crop
load('D:\Kevin\Sequence_tap_analysis\basic_beh_ephys_analysis\MC_DLS_comparisons_conglomerate\FR_stats\FRstats_new_nosmooth.mat');
ratname = {}; unitNum = [];
for u = 1:length(unitPathSave_DLS)
    ratname{end+1} = unitPathSave_DLS{u}(strfind(unitPathSave_DLS{u},'Rat') + (0:4));
    utemp= regexp( unitPathSave_DLS{u},'\d*','Match'); unitNum(end+1) = str2num(utemp{end});
end
% for u = 1:length(unitPathSave_MC)
%     ratname{end+1} = unitPathSave_MC{u}(strfind(unitPathSave_MC{u},'Rat') + (0:4));
%     utemp = regexp( unitPathSave_MC{u},'\d*','Match'); unitNum(end+1) = str2num(utemp{end});
% end

% %% super quick script to rename the messed up files in rat 33
% for f = 1:length(fname)
%     if strcmp(fname(f).name(1:5), 'Rat33')
%         curfilename = fullfile(fpath, fname(f).name);
%         newfilename = fullfile(fpath, fname(f).name(6:end));
%         movefile(curfilename, newfilename);
%     end
% end

%% pull data with and without lag
% - pull lag separate from cued, ot both
fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat33_NOTQUALITY';...
    'X:\Kevin\Analysis\Data\Rat35_NOTQUALITY'}; % no rat 18 since somethings up...
fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat33';...
    'X:\Kevin\Analysis\Data\Rat35_NOTQUALITY'}; % no rat 18 since somethings up...
% reran for rat 18? cant recall what was up maybe itll work now
fpaths_all = {'X:\Kevin\Analysis\Data\Rat47';'X:\Kevin\Analysis\Data\Rat33';...
    'X:\Kevin\Analysis\Data\Rat35_NOTQUALITY';'X:\Kevin\Analysis\Data\Rat18';}; 
% list units??
% 
% % for MC now
% fpaths_all = {'X:\Kevin\Analysis\Data\Rat45';'X:\Kevin\Analysis\Data\Rat25';...
%     'X:\Kevin\Analysis\Data\Rat63';'X:\Kevin\Analysis\Data\Rat67';'X:\Kevin\Analysis\Data\Rat72'}; 


%files_to_pull = '*context_vel.mat'; % lag?
files_to_pull = '*sequence_vel.mat';

nn_vals_vel = cell(1,length(fpaths_all));
xgb_vals_vel = cell(1,length(fpaths_all));
unit_use = cell(1,length(fpaths_all));
for rat = 1:length(fpaths_all)
    fpath = fpaths_all{rat};
   % if contains(fpath, 'Rat67'); files_to_pull = '*context_vel_lag.mat'; else; files_to_pull = '*context_vel.mat'; end
   % fname = dir(fullfile(fpath, '*context_vel.mat')); % pull velocity encoding data
    fname = dir(fullfile(fpath, files_to_pull)); % pull velocity encoding data
    if isempty(fname); disp(fpath); end
    ratmatch = fpath(strfind(fpath, 'Rat') + (0:4));
    % better getting rat
    ratmatch  =strsplit(fpath, 'Rat'); ratmatch = ['Rat' ratmatch{2}(1:2)];
    
    count = 1;
    for f = 1:length(fname)
        % check if unit is on use path
        uuse = regexp( fname(f).name,'\d*','Match'); uuse = str2num(uuse{1});
        if ~ismember(uuse, unitNum(contains(ratname, ratmatch))); continue; end
        
        load(fullfile(fpath, fname(f).name));
        u = strsplit(fname(f).name,'_');
        if contains(u{1}, 'Rat'); u{1} = u{1}(6:end); end
        unit_use{rat}(count) = str2num(u{1});
        
        nn_vals_vel{rat}(count,:) = [mean(nn_cued), mean(nn_ot), mean(nn_test), mean(nn_train)];
        xgb_vals_vel{rat}(count,:) = [mean(xgboost_cued), mean(xgboost_ot), mean(xgboost_test),mean(xgboost_train)];
        % load for many iterations of shuffled dataset?
        count = count+1;
    end

end

%files_to_pull = '*context_vel_lag.mat'; % lag?
files_to_pull = '*sequence_vel_lag.mat';
nn_vals_vel_lag = cell(1,length(fpaths_all));
xgb_vals_vel_lag = cell(1,length(fpaths_all));
unit_use = cell(1,length(fpaths_all));
for rat = 1:length(fpaths_all)
    fpath = fpaths_all{rat};
   % if contains(fpath, 'Rat67'); files_to_pull = '*context_vel_lag.mat'; else; files_to_pull = '*context_vel.mat'; end
   % fname = dir(fullfile(fpath, '*context_vel.mat')); % pull velocity encoding data
    fname = dir(fullfile(fpath, files_to_pull)); % pull velocity encoding data
    if isempty(fname); disp(fpath); end
    ratmatch = fpath(end-4:end);
    ratmatch  =strsplit(fpath, 'Rat'); ratmatch = ['Rat' ratmatch{2}(1:2)];

    
    count=1;
    for f = 1:length(fname)
        % check if unit is on use path
        uuse = regexp( fname(f).name,'\d*','Match'); uuse = str2num(uuse{1});
        if ~ismember(uuse, unitNum(contains(ratname, ratmatch))); continue; end
        
        load(fullfile(fpath, fname(f).name));
        u = strsplit(fname(f).name,'_');
        if contains(u{1}, 'Rat'); u{1} = u{1}(6:end); end
        unit_use{rat}(count) = str2num(u{1});
        
        nn_vals_vel_lag{rat}(count,:) = [mean(nn_cued), mean(nn_ot), mean(nn_test), mean(nn_train)];
        xgb_vals_vel_lag{rat}(count,:) = [mean(xgboost_cued), mean(xgboost_ot), mean(xgboost_test),mean(xgboost_train)];
        % load for many iterations of shuffled dataset?
        count = count+1;
    end

end

%% make plot from above

cued_vel = cellfun(@(v) v(:,1)',xgb_vals_vel,'un',0);
ot_vel = cellfun(@(v) v(:,2)',xgb_vals_vel,'un',0); 
both_vel = cellfun(@(v) v(:,3)',xgb_vals_vel,'un',0); 
lag_vel = cellfun(@(v) v(:,3)', xgb_vals_vel_lag,'un',0);
% - toss those < 0?
thresh = -1.05;
cued_vel_rat = cellfun(@(v) mean(v(v>thresh)),cued_vel);
ot_vel_rat = cellfun(@(v) mean(v(v>thresh)),ot_vel);
both_vel_rat = cellfun(@(v) mean(v(v>thresh)),both_vel);
lag_vel_rat = cellfun(@(v) mean(v(v>thresh)),lag_vel);
figure; hold on;
bar(1:4, [mean(cued_vel_rat), mean(ot_vel_rat), mean(both_vel_rat), mean(lag_vel_rat)]);
for rat = 1:length(cued_vel_rat)
    plot(1:4, [cued_vel_rat(rat), ot_vel_rat(rat), both_vel_rat(rat), lag_vel_rat(rat)],'k');
end
[~,pcuedot] = ttest(cued_vel_rat, ot_vel_rat);
[~,pcuedboth] = ttest(cued_vel_rat, both_vel_rat);
[~,potboth] = ttest(ot_vel_rat, both_vel_rat);
% [pcuedot] = signrank(cued_vel_rat, ot_vel_rat);
% [pcuedboth] = signrank(cued_vel_rat, both_vel_rat);
% [potboth] = signrank(ot_vel_rat, both_vel_rat);
sigstar({[1,2],[1,3],[2,3]},[pcuedot, pcuedboth,potboth],0,1);

% add in a bar for lag
[~,p14] = ttest(cued_vel_rat, lag_vel_rat);
[~,p24] = ttest(ot_vel_rat, lag_vel_rat);
[~,p34] = ttest(both_vel_rat, lag_vel_rat);
% [p14] = signrank(cued_vel_rat, lag_vel_rat);
% [p24] = signrank(ot_vel_rat, lag_vel_rat);
% [p34] = signrank(both_vel_rat, lag_vel_rat);
sigstar({[1,4],[2,4],[3,4]},[p14,p24,p34],0,1);

xticks(1:4);
xticklabels({'Train cued','Train ot','train both','lag control'})
% for sequences
xticklabels({'Train half','Train half','Train all','lag control'});

%% save for across brain regions
xgb_vals_vel_DLS = xgb_vals_vel;
nn_vals_vel_DLS = nn_vals_vel;
xgb_vals_vel_lag_DLS = xgb_vals_vel_lag;

xgb_vals_vel_MC = xgb_vals_vel;
nn_vals_vel_MC = nn_vals_vel;
xgb_vals_vel_lag_MC = xgb_vals_vel_lag;
%%
% first remake effect for rat 47, and rat 45
% - of generalization
% - looks good
% when add all rats, DLS rat encodes much worse...

% do above, but compare for individual units?


cued_vel = cellfun(@(v) v(:,1)',xgb_vals_vel_DLS(1),'un',0); cued_vel = [cued_vel{:}];
ot_vel = cellfun(@(v) v(:,2)',xgb_vals_vel_DLS(1),'un',0); ot_vel = [ot_vel{:}];
both_vel = cellfun(@(v) v(:,3)',xgb_vals_vel_DLS(1),'un',0); both_vel = [both_vel{:}];
lag_vel = cellfun(@(v) v(:,3)', xgb_vals_vel_lag_DLS([1,2,4]),'un',0); lag_vel = [lag_vel{:}];
% toss thresh?
thresh = 0;%-5.1;
cued_vel(cued_vel<thresh) = 0; cued_vel_DLS = cued_vel;
ot_vel(ot_vel<thresh) = 0; ot_vel_DLS = ot_vel;
both_vel(both_vel<thresh) = 0; both_vel_DLS = both_vel;
lag_vel(lag_vel<thresh) = 0; lag_vel_DLS = lag_vel;

figure; hold on;
% histogram(cued_vel);
% histogram(ot_vel)
% histogram(both_vel);
pcuedot = ranksum(cued_vel, ot_vel);
pcuedboth = ranksum(cued_vel, both_vel);
potboth = ranksum(ot_vel, both_vel); % this turns up significant...
swarmchart(repelem(1:4, [length(cued_vel), length(ot_vel), length(both_vel), length(lag_vel)]),...
    [cued_vel, ot_vel, both_vel, lag_vel],'r.');
% boxplot([cued_vel, ot_vel, both_vel, lag_vel],...
%     repelem(1:4, [length(cued_vel), length(ot_vel), length(both_vel), length(lag_vel)]))
boxchart(repelem(1:4, [length(cued_vel), length(ot_vel), length(both_vel), length(lag_vel)]),...
    [cued_vel, ot_vel, both_vel, lag_vel],'MarkerStyle','none')
sigstar({[1,2],[1,3],[2,3]},[pcuedot, pcuedboth, potboth]);

% go through and filter units

cued_vel = cellfun(@(v) v(:,1)',xgb_vals_vel_MC(1),'un',0); cued_vel = [cued_vel{:}];
ot_vel = cellfun(@(v) v(:,2)',xgb_vals_vel_MC(1),'un',0); ot_vel = [ot_vel{:}];
both_vel = cellfun(@(v) v(:,3)',xgb_vals_vel_MC(1),'un',0); both_vel = [both_vel{:}];
lag_vel = cellfun(@(v) v(:,3)', xgb_vals_vel_lag_MC(1),'un',0); lag_vel = [lag_vel{:}];
% toss thresh?
thresh = 0;%-5.01;
cued_vel(cued_vel<thresh) = 0;
ot_vel(ot_vel<thresh) = 0;
both_vel(both_vel<thresh) = 0;
lag_vel(lag_vel<thresh)=0;

%figure; hold on;
% histogram(cued_vel);
% histogram(ot_vel)
% histogram(both_vel);
pcuedot = ranksum(cued_vel, ot_vel);
pcuedboth = ranksum(cued_vel, both_vel);
potboth = ranksum(ot_vel, both_vel);
% boxplot([cued_vel, ot_vel, both_vel],...
%     repelem(1:3, [length(cued_vel), length(ot_vel), length(both_vel)]))
boxchart(repelem(5:8, [length(cued_vel), length(ot_vel), length(both_vel), length(lag_vel)]),...
    [cued_vel, ot_vel, both_vel,lag_vel],'MarkerStyle','none')
swarmchart(repelem(5:8, [length(cued_vel), length(ot_vel), length(both_vel), length(lag_vel)]),...
    [cued_vel, ot_vel, both_vel, lag_vel],'b.');

sigstar({[5,6],[5,7],[6,7]},[pcuedot, pcuedboth, potboth]);

xlim([0,9])

% compaer across?
pcued = ranksum(cued_vel, cued_vel_DLS);
pot = ranksum(ot_vel, ot_vel_DLS);
pboth = ranksum(both_vel, both_vel_DLS);

sigstar({[1,5],[2,6],[3,7]},[pcued,pot,pboth]);

xticks(1:8); xticklabels({'Cued','OT','Both','Lag','Cued','OT','Both','Lag'});
ylabel('Pseudo r2'); title('Encoding');

%% new for full lag analysis 2022
% - goal is to look at encoding values for different neurons at different
% lags
% - different neurons should turn up different lags

% * rerunning for neurons that werent here long enough
% - make sure they are overwriting!!
% some: 133, 141, 191, 286, 593, 432

tlag = [-20,-16,-12,-8,-4,0,4,8,12,16,20];
tlag = tlag * .025 ; % 25 ms bins
num_repeats = 10;

%*** todo, rerun these with lag on the OT only sequences
%  - see if that's having the different on lead lag
%  - ok during DLS and MC comparisons, but not for FT and AO???
% - especially since DLS OT only is way different than DLS FT

% *** to do, update these to take only the OT sequence, or the one with
% everything
% - in order to compare to OT only animals, need to use only OT


fpaths_all = {'X:\Kevin\Analysis\Data\Rat47_alllag';'X:\Kevin\Analysis\Data\Rat33_alllag';...
    'X:\Kevin\Analysis\Data\Rat35_NOTQUALITY_alllag';'X:\Kevin\Analysis\Data\Rat18_alllag';}; 
 fpaths_all = {'X:\Kevin\Analysis\Data\Rat45_alllag'; 'X:\Kevin\Analysis\Data\Rat25_alllag';...
      'X:\Kevin\Analysis\Data\Rat63_alllag'}; 
 % fpaths_all = {'X:\Kevin\Analysis\Data\Rat81_quality_alllag'}; 
  fpaths_all = {'X:\Kevin\Analysis\Data\Rat56_alllag'}; 

%MC
fpaths_all = {'Z:\Kevin\Data\Rat45_alllag'; 'Z:\Kevin\Data\Rat25_alllag'; 'Z:\Kevin\Data\Rat63_alllag'};
%fpaths_all = {'Z:\Kevin\Data\Rat81_quality_alllag'};
%DLS
fpaths_all = {'Z:\Kevin\Data\Rat56_alllag'}; 
%fpaths_all = {'Z:\Kevin\Data\Rat47_alllag';'Z:\Kevin\Data\Rat33_alllag';'Z:\Kevin\Data\Rat35_NOTQUALITY_alllag';'Z:\Kevin\Data\Rat18_alllag';}; 
  
baduuse = cell(length(fpaths_all),1);

% for OT only, just use sequence_vel.mat
% for full task, need to specify
files_to_pull = '*sequence_vel*.mat'; % uncomment for OT only rats
%files_to_pull = '*sequence_vel*OTseqonly.mat'; % uncomment for FT rats
nn_vals_vel_lag = cell(1,length(fpaths_all));
xgboost_vals_vel_lag = cell(1,length(fpaths_all));
unit_use = cell(1,length(fpaths_all));

xgboost_vals_vel_all = cell(1,length(fpaths_all));

for rat = 1:length(fpaths_all)
    fpath = fpaths_all{rat};

    fname = dir(fullfile(fpath, files_to_pull)); % pull velocity encoding data
    if isempty(fname); disp(fpath); end
    ratmatch = fpath(end-4:end);
    ratmatch  =strsplit(fpath, 'Rat'); ratmatch = ['Rat' ratmatch{2}(1:2)];

    
    count=1;
    for f = 1:length(fname)
        % check if unit is on use path
        uuse = regexp( fname(f).name,'\d*','Match'); uuse = str2num(uuse{1});
       % if ~ismember(uuse, unitNum(contains(ratname, ratmatch))); continue; end
        unit_use{rat}(count) = uuse;
        
        load(fullfile(fpath, fname(f).name));
        
        % save even if didnt finish running because I need something
        xgboost_train = reshape(xgboost_train, [num_repeats, length(xgboost_train)/num_repeats]);
        xgboost_test = reshape(xgboost_test, [num_repeats, length(xgboost_test)/num_repeats]);
        xgboost_val =  mean([xgboost_train; xgboost_test]);
        xgboost_val = [xgboost_val, NaN(1,length(tlag) - length(xgboost_val))];
        xgboost_vals_vel_all{rat}(end+1,:) = xgboost_val;
        
        load(fullfile(fpath, fname(f).name));

        
        if length(xgboost_train)~=num_repeats*length(tlag)
            baduuse{rat}(end+1) = uuse;
            %disp(uuse);
            continue
        end
        
        % reshape data?
        % cued, ot, test, train
        % - this wont make sense since not testing on a lagged dataset...
        % - but want to keep train and test datasets, not split seqs
        
        nn_train = reshape(nn_train, [num_repeats, length(tlag)]);
        xgboost_train = reshape(xgboost_train, [num_repeats, length(tlag)]);
        xgboost_test = reshape(xgboost_test, [num_repeats, length(tlag)]);
        nn_test = reshape(nn_test, [num_repeats, length(tlag)]);
        
        % save
        xgboost_vals_vel_lag{rat}(count,:) = mean([xgboost_train; xgboost_test]);
        nn_vals_vel_lag{rat}(count,:) = mean([nn_train; nn_test]);

        count = count+1;
    end

    
end

%xgboost_vals_FT = xgboost_vals_vel_all;
%xgboost_vals_AO = xgboost_vals_vel_all;

%xgboost_vals_FT = xgboost_vals_vel_lag; unit_use_MC_FT = unit_use;
%xgboost_vals_AO = xgboost_vals_vel_lag;


%xgboost_vals_FT_DLS = xgboost_vals_vel_lag; unit_use_DLS_FT = unit_use;
xgboost_vals_AO_DLS = xgboost_vals_vel_lag;

%% plot a sample of neurons

% negative is future kinematics predicting past spikes
% positive is past kinematics predicting future spikes

% so negative is lead (spikes leading kinematics)
% and positive is lag (spikes following kinematics)

rat = 2;
rat = 1:3;
%rat = 1;
% random
sampuse = randsample(length(xgboost_vals_vel_lag{rat}),16);
% max encoders
[~,sampuse] = sort( mean( xgboost_vals_vel_lag{rat} ,2) ,'descend');
sampuse = sampuse(1:16);

figure;
for j = 1:16
    subplot(4,4,j); hold on;
    plot(tlag, xgboost_vals_vel_lag{rat}(sampuse(j),:));
    title(num2str(unit_use{rat}(sampuse(j))))
end

% make histogram of all lead lag
%* include a filter that only takes decoding value > 0?
peakid = []; peakval = [];
xgboost_vsls_vel_lag_all = vertcat(xgboost_vals_vel_lag{rat});
for u = 1:length(xgboost_vsls_vel_lag_all)
    [peakval(u),peakid(u)] = max(xgboost_vsls_vel_lag_all(u,:));
end
peakid(peakval<0) = [];
lagval = tlag(peakid);
dt = mode(diff(tlag));
edges = linspace(tlag(1)-dt/2, tlag(end)+dt/2, length(tlag)+1);
figure;histogram(lagval,edges)
xlabel('lead --- lag');

% all rats?

%% look across rats

%** one difference is that I'm looking only at a single sequence here. so
%the decoding could be completely different? 
% - maybe want to run AO encoding models on the FT rats (on their OT only
% sequence) as a better comparison?

%lagval_full = lagval;
%lagval_only = lagval;
%lagval_DLS = lagval;
%lagval_DLSonly = lagval;

figure; hold on;
histogram(lagval_DLS, edges,'Normalization','probability','DisplayStyle','stairs','linewidth',2)
histogram(lagval_DLSonly, edges,'Normalization','probability','DisplayStyle','stairs','linewidth',2)

histogram(lagval_full,edges,'Normalization','probability','DisplayStyle','stairs','linewidth',2)
histogram(lagval_only,edges,'Normalization','probability','DisplayStyle','stairs','linewidth',2)
legend('DLS - FT','DLS - AO','MC - FT','MC - AO')
xlabel('lead --- lag');

p1 = ranksum(lagval_DLS, lagval_full);
p2 = ranksum(lagval_only, lagval_full);

[~,p1] = kstest2(lagval_DLS, lagval_full);
[~,p2] = kstest2(lagval_only, lagval_full);


%% plot max FT and AO encoding
%*** this is from all neurons, including some who didnt run on all lags!!
% - but who did run on lag = 0

%xgboost_vals_AO; xgboost_vals_FT;

xgboost_vals_max_FT = cellfun(@(v) max(v')', xgboost_vals_FT,'un',0);
xgboost_vals_max_AO = cellfun(@(v) max(v')', xgboost_vals_AO,'un',0);


figure; hold on;
histogram(vertcat(xgboost_vals_max_FT{:}), 'Normalization','probability');
histogram(vertcat(xgboost_vals_max_AO{:}), 'Normalization','probability');

temp = cellfun(@mean, xgboost_vals_max_FT);
arrayfun(@(v) plot([v,v],ylim,'b--'),temp)
temp = cellfun(@mean, xgboost_vals_max_AO);
arrayfun(@(v) plot([v,v],ylim,'r--'),temp)

legend('FT','AO');
title('Encoding accuracy (kinematics -> spikes) of automatic sequence (xgboost)')
ylabel('Pseudo-R2'); xlabel('Fraction');

% plot lags for FT and AO
% - using only the fully ran neurons

peakidFT = {}; lagvalFT = {};
for r = 1:length(xgboost_vals_FT)
    peakval = [];
    for u = 1:length(xgboost_vals_FT{r})
        [peakval(u),peakidFT{r}(u)] = max(xgboost_vals_FT{r}(u,:));
    end
    %lagvalFT{r} = tlag(peakidFT{r});
    lagvalFT{r} = tlag(peakidFT{r}(peakval>=0)); % do I want this?
end

peakidAO = {}; lagvalAO = {};
for r = 1:length(xgboost_vals_AO)
    peakval = [];
    for u = 1:length(xgboost_vals_AO{r})
        [peakval(u),peakidAO{r}(u)] = max(xgboost_vals_AO{r}(u,:));
    end
    %lagvalAO{r} = tlag(peakidAO{r});
    lagvalAO{r} = tlag(peakidAO{r}(peakval>=0));
end

figure; hold on;
dt = mode(diff(tlag));
edges = linspace(tlag(1)-dt/2, tlag(end)+dt/2, length(tlag)+1);
histogram([lagvalFT{:}],edges, 'Normalization','probability');
histogram([lagvalAO{:}],edges, 'Normalization','probability');

temp = cellfun(@mean, lagvalFT);
arrayfun(@(v) plot([v,v],ylim,'b--'),temp)
temp = cellfun(@mean, lagvalAO);
arrayfun(@(v) plot([v,v],ylim,'r--'),temp)

legend('FT','AO');

title('Optimal lag for encoding (seconds)');
xlabel('lead --- lag'); ylabel('Fraction');

%% same as above but with DLS rats

%*** this is from all neurons, including some who didnt run on all lags!!
% - but who did run on lag = 0

%xgboost_vals_AO; xgboost_vals_FT;

% take max encoding
xgboost_vals_max_FT = cellfun(@(v) max(v')', xgboost_vals_FT,'un',0);
xgboost_vals_max_AO = cellfun(@(v) max(v')', xgboost_vals_AO,'un',0);
xgboost_vals_max_FT_DLS = cellfun(@(v) max(v')', xgboost_vals_FT_DLS,'un',0);
xgboost_vals_max_AO_DLS = cellfun(@(v) max(v')', xgboost_vals_AO_DLS,'un',0);

% take encoding only at 0
id_tuse = find(tlag==0);
xgboost_vals_max_FT = cellfun(@(v) v(:,id_tuse), xgboost_vals_FT,'un',0);
xgboost_vals_max_AO = cellfun(@(v) v(:,id_tuse), xgboost_vals_AO,'un',0);
xgboost_vals_max_FT_DLS = cellfun(@(v) v(:,id_tuse), xgboost_vals_FT_DLS,'un',0);
xgboost_vals_max_AO_DLS = cellfun(@(v) v(:,id_tuse), xgboost_vals_AO_DLS,'un',0);


figure; hold on;
edges = -.3:.1:1;
histogram(vertcat(xgboost_vals_max_FT{:}),edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',1.5);
histogram(vertcat(xgboost_vals_max_AO{:}),edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor','c','LineWidth',1.5);
histogram(vertcat(xgboost_vals_max_FT_DLS{:}),edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',1.5);
histogram(vertcat(xgboost_vals_max_AO_DLS{:}),edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor','m','LineWidth',1.5);


temp = cellfun(@mean, xgboost_vals_max_FT);arrayfun(@(v) plot([v,v],ylim,'b--'),temp)
temp = cellfun(@mean, xgboost_vals_max_AO);arrayfun(@(v) plot([v,v],ylim,'c--'),temp)
temp = cellfun(@mean, xgboost_vals_max_FT_DLS);arrayfun(@(v) plot([v,v],ylim,'r--'),temp)
temp = cellfun(@mean, xgboost_vals_max_AO_DLS);arrayfun(@(v) plot([v,v],ylim,'m--'),temp)

legend('FT - MC','AO - MC','FT - DLS','AO - DLS');

title('Encoding accuracy (kinematics -> spikes) of automatic sequence (xgboost)')
ylabel('Pseudo-R2'); xlabel('Fraction');

% plot lags for FT and AO
% - using only the fully ran neurons

peakidFT = {}; lagvalFT = {};
for r = 1:length(xgboost_vals_FT)
    peakval = [];
    for u = 1:length(xgboost_vals_FT{r})
        [peakval(u),peakidFT{r}(u)] = max(xgboost_vals_FT{r}(u,:));
    end
    %lagvalFT{r} = tlag(peakidFT{r});
    lagvalFT{r} = tlag(peakidFT{r}(peakval>=0)); % do I want this?
end
peakidAO = {}; lagvalAO = {};
for r = 1:length(xgboost_vals_AO)
    peakval = [];
    for u = 1:length(xgboost_vals_AO{r})
        [peakval(u),peakidAO{r}(u)] = max(xgboost_vals_AO{r}(u,:));
    end
    %lagvalAO{r} = tlag(peakidAO{r});
    lagvalAO{r} = tlag(peakidAO{r}(peakval>=0));
end

peakidFT = {}; lagvalFT_DLS = {};
for r = 1:length(xgboost_vals_FT_DLS)
    peakval = [];
    for u = 1:length(xgboost_vals_FT_DLS{r})
        [peakval(u),peakidFT{r}(u)] = max(xgboost_vals_FT_DLS{r}(u,:));
    end
    %lagvalFT{r} = tlag(peakidFT{r});
    lagvalFT_DLS{r} = tlag(peakidFT{r}(peakval>=0)); % do I want this?
end
peakidAO = {}; lagvalAO_DLS = {};
for r = 1:length(xgboost_vals_AO_DLS)
    peakval = [];
    for u = 1:length(xgboost_vals_AO_DLS{r})
        [peakval(u),peakidAO{r}(u)] = max(xgboost_vals_AO_DLS{r}(u,:));
    end
    %lagvalAO{r} = tlag(peakidAO{r});
    lagvalAO_DLS{r} = tlag(peakidAO{r}(peakval>=0));
end


figure; hold on;
dt = mode(diff(tlag));
edges = linspace(tlag(1)-dt/2, tlag(end)+dt/2, length(tlag)+1);
histogram([lagvalFT{:}],edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',1.5);
histogram([lagvalAO{:}],edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor','c','LineWidth',1.5);
histogram([lagvalFT_DLS{:}],edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor','r','LineWidth',1.5);
histogram([lagvalAO_DLS{:}],edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor','m','LineWidth',1.5);


temp = cellfun(@mean, lagvalFT);arrayfun(@(v) plot([v,v],ylim,'b--'),temp)
temp = cellfun(@mean, lagvalAO);arrayfun(@(v) plot([v,v],ylim,'c--'),temp)
temp = cellfun(@mean, lagvalFT_DLS);arrayfun(@(v) plot([v,v],ylim,'r--'),temp)
temp = cellfun(@mean, lagvalAO_DLS);arrayfun(@(v) plot([v,v],ylim,'m--'),temp)

legend('FT - MC','AO - MC','FT - DLS','AO - DLS');


title('Optimal lag for encoding (seconds)');
xlabel('lead --- lag'); ylabel('Fraction');


%% just make plots for DLS

%% 0) plot sample neurons AO 
rat = 1;
% random
sampuse = randsample(length(xgboost_vals_AO_DLS{rat}),16);
% max encoders
[~,sampuse] = sort( mean( xgboost_vals_AO_DLS{rat} ,2) ,'descend');
sampuse = sampuse(1:16);

figure;
for j = 1:16
    subplot(4,4,j); hold on;
    plot(tlag, xgboost_vals_AO_DLS{rat}(sampuse(j),:));
    %title(num2str(unit_use_DLS_FT{rat}(sampuse(j))))
end

%% 1) plot sample neurons


rat = 1;
% random
sampuse = randsample(length(xgboost_vals_FT_DLS{rat}),16);
% max encoders
[~,sampuse] = sort( mean( xgboost_vals_FT_DLS{rat} ,2) ,'descend');
sampuse = sampuse(1:16);

figure;
for j = 1:16
    subplot(4,4,j); hold on;
    plot(tlag, xgboost_vals_FT_DLS{rat}(sampuse(j),:));
    title(num2str(unit_use_DLS_FT{rat}(sampuse(j))))
end

% make histogram of all lead lag
%* include a filter that only takes decoding value > 0?
rat  =1:4; % all dls rats
peakid = []; peakval = [];
xgboost_vsls_vel_lag_all = vertcat(xgboost_vals_FT_DLS{rat});
for u = 1:length(xgboost_vsls_vel_lag_all)
    [peakval(u),peakid(u)] = max(xgboost_vsls_vel_lag_all(u,:));
end
peakid(peakval<0) = [];
lagval = tlag(peakid);
dt = mode(diff(tlag));
edges = linspace(tlag(1)-dt/2, tlag(end)+dt/2, length(tlag)+1);
figure;histogram(lagval,edges)
xlabel('lead --- lag');



%% 2)  plot max FT and AO encoding
%*** this is from all neurons, including some who didnt run on all lags!!
% - but who did run on lag = 0

%xgboost_vals_AO; xgboost_vals_FT;

xgboost_vals_max_FT = cellfun(@(v) max(v')', xgboost_vals_FT_DLS,'un',0);
xgboost_vals_max_AO = cellfun(@(v) max(v')', xgboost_vals_AO_DLS,'un',0);

edges = -0.3:.1:.9;
figure; hold on;
histogram(vertcat(xgboost_vals_max_FT{:}),edges, 'Normalization','probability');
histogram(vertcat(xgboost_vals_max_AO{:}),edges, 'Normalization','probability');

temp = cellfun(@mean, xgboost_vals_max_FT);
arrayfun(@(v) plot([v,v],ylim,'b--'),temp)
temp = cellfun(@mean, xgboost_vals_max_AO);
arrayfun(@(v) plot([v,v],ylim,'r--'),temp)

legend('FT','AO');
title('Encoding accuracy (kinematics -> spikes) of automatic sequence (xgboost)')
ylabel('Pseudo-R2'); xlabel('Fraction');

% plot lags for FT and AO
% - using only the fully ran neurons

peakidFT = {}; lagvalFT = {};
for r = 1:length(xgboost_vals_FT_DLS)
    peakval = [];
    for u = 1:length(xgboost_vals_FT_DLS{r})
        [peakval(u),peakidFT{r}(u)] = max(xgboost_vals_FT_DLS{r}(u,:));
    end
    %lagvalFT{r} = tlag(peakidFT{r});
    lagvalFT{r} = tlag(peakidFT{r}(peakval>=0)); % do I want this?
end

peakidAO = {}; lagvalAO = {};
for r = 1:length(xgboost_vals_AO_DLS)
    peakval = [];
    for u = 1:length(xgboost_vals_AO_DLS{r})
        [peakval(u),peakidAO{r}(u)] = max(xgboost_vals_AO_DLS{r}(u,:));
    end
    %lagvalAO{r} = tlag(peakidAO{r});
    lagvalAO{r} = tlag(peakidAO{r}(peakval>=0));
end

figure; hold on;
dt = mode(diff(tlag));
edges = linspace(tlag(1)-dt/2, tlag(end)+dt/2, length(tlag)+1);
histogram([lagvalFT{:}],edges, 'Normalization','probability');
histogram([lagvalAO{:}],edges, 'Normalization','probability');

temp = cellfun(@mean, lagvalFT);
arrayfun(@(v) plot([v,v],ylim,'b--'),temp)
temp = cellfun(@mean, lagvalAO);
arrayfun(@(v) plot([v,v],ylim,'r--'),temp)

legend('FT','AO');

title('Optimal lag for encoding (seconds)');
xlabel('lead --- lag'); ylabel('Fraction');
