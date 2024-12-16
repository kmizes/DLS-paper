
%% Setup data
uniRetrain = [1,0,0,1,0,0,0];
biRetrain = [0,0,0,0,0,0,0];
uniRetrainOT = [1,0,0,1,0,0,0];
biRetrainOT = [0,0,0,0,0,1,0];


ratPath = {'Y:\Users\kmizes\D_backup\Kevin\Sequence_tap\L2_output\Results-L2-Rat70\ratBEHstruct.mat';...
   'Y:\Users\kmizes\D_backup\Kevin\Sequence_tap\D7_output\Results-D7-Rat74\ratBEHstruct.mat';...
   'Y:\Users\kmizes\D_backup\Kevin\Sequence_tap\J1_output\Results-J1-Rat71\ratBEHstruct.mat';...
   'Y:\Users\kmizes\D_backup\Kevin\Sequence_tap\L3_output\Results-L3-Rat65\ratBEHstruct.mat';...
   'Y:\Users\kmizes\D_backup\Kevin\Sequence_tap\J7_output\Results-J7-Rat78\ratBEHstruct.mat';...
   'Y:\Users\kmizes\D_backup\Kevin\Sequence_tap\F5_output\Results-F5-Rat86\ratBEHstruct.mat';...
   'Y:\Users\kmizes\D_backup\Kevin\Sequence_tap\L3_output\Results-L3-Rat117\ratBEHstruct.mat'};

OTseqall = {'CLC','CLR','RCR','LRC','LRC','CRL','RCL'};



%% take 2

% days to actually sample from
tBef = 14; tAft = 14;
tBefOT = 14; tAftOT = 14;
num_actual = 7;

% tBef = 21; tAft = 21;
% tBefOT = 21; tAftOT = 21;
% num_actual = 14;

addOneWeek = 0;


num_avg = 3; % sessions to average cued and OT trials over to equal 1 day
num_avg_OT = 1; 

tBef = tBef * num_avg;
tAft = tAft * num_avg;

ratMockCued_rat = [];
ratUniCued_rat = [];
ratBiCued_rat = [];

ratMockWM_rat = [];
ratUniWM_rat = [];
ratBiWM_rat = [];

ratMockOT_rat = [];
ratUniOT_rat = [];
ratBiOT_rat = [];
ratInf_rat = [];
ratInfOT_rat = [];
ratInfWM_rat = [];

for rat = 1:length(ratPath)
    
try
    load(ratPath{rat})
catch
    load(['F' ratPath{rat}(2:end)]);
end
disp(ratPath{rat})


% J778
if strcmp(ratBEHstruct(1).name, 'J7-Rat78')
sMock = [539, 567];
sUni = [788, 815];
sBi = [943, 971];

% L2 rat 70
elseif strcmp(ratBEHstruct(1).name, 'L2-Rat70')
sMock = [890, 918];
sUni = [1062, 1090];
sBi = [1294, 1325];

% d7 rat 74
elseif strcmp(ratBEHstruct(1).name, 'D7-Rat74')
sMock = [621, 648];
sUni = [793, 813];
sBi = [982, 1000];
ratBEHstruct(end-50:end) = []; % box was broken

% e1 rat 75
elseif strcmp(ratBEHstruct(1).name, 'E1-Rat75')
sMock = [718, 745];
sUni = [802, 829];
sBi = [957, 984];

% j1 rat 71
elseif strcmp(ratBEHstruct(1).name, 'J1-Rat71')
sMock = [757, 785];
sUni = [930, 957];
sUni = [930-8,957];
sBi = [1133, 1164];

% j4 rat 76
elseif strcmp(ratBEHstruct(1).name, 'J4-Rat76')
sMock = [709, 736];
sUni = [793, 820];
sBi = [921, 952];

% L3 rat 65
elseif strcmp(ratBEHstruct(1).name, 'L3-Rat65')
sMock = [921, 949];
sUni = [1170, 1197];
sBi = [1363, 1381];

elseif strcmp(ratBEHstruct(1).name, 'L3-Rat34')
    disp('this is a MC rat for plotting purposes only')
    ratBEHstruct =  ratBEHstruct(890:end); % lesion
    sMock = [522, 565];
    sUni = [682,773 ];
    sBi = [921, 952];
    
elseif strcmp(ratBEHstruct(1).name, 'F5-Rat86')
sMock = [799, 827];
sUni = [1235, 1263];
sBi = [1360, 1386];
sBi = [1360,1416];
% new including the pre-wm task times!
sMock = [921, 949];
sUni = [1357,1385];
sBi = [1482,1538];


elseif strcmp(ratBEHstruct(1).name, 'L3-Rat117')
sMock = [335, 366];
sUni = [591, 618];
sBi = [703, 730];

elseif strcmp(ratBEHstruct(1).name, 'L5-Rat118')
sMock = [332, 363];
sUni = [588, 615];
sBi = [700, 726];

elseif strcmp(ratBEHstruct(1).name, 'F4-Rat100')
sMock = [571,603];
sUni = [815, 842];
sBi = [939, 967];

elseif strcmp(ratBEHstruct(1).name, 'E1-Rat110')
sMock = [518, 545];
sUni = [630, 657];
sBi = [754, 781];

else
    disp('somethings wrong')
end

% change to smooth days over!
if addOneWeek
    sUni(2) = sUni(2)+ num_actual * (num_avg + num_avg_OT);
    sBi(2) = sBi(2) + num_actual * (num_avg+num_avg_OT);
end

% pull data clean for cued
% - before
count = 0; s = sMock(1); accPreMock = []; accPreMockWM = [];
while count < tBef
    Hit = ratBEHstruct(s).Hit;
    blocknumRepair = ratBEHstruct(s).blocknumRepair;
    if isempty(blocknumRepair); blocknumRepair = ratBEHstruct(s).blocknum; end
    if length(Hit) < 10; s = s-1; continue; end
    if ratBEHstruct(s).protocol~=7; s=s-1; continue; end;
    % get cued accuracy
    idxCued = find(blocknumRepair < 3);
    idxWM = find(blocknumRepair>=3);
    if length(Hit) < 6; idxCued = ones(1,length(Hit)); end
    accPreMock(count+1) = sum(Hit(idxCued))/length(idxCued);
    accPreMockWM(count+1) = sum(Hit(idxWM))/length(idxWM);
    count = count+1; s = s-1;
end
count = 0; s = sMock(2); accPostMock = []; accPostMockWM = [];
while count < tAft
    Hit = ratBEHstruct(s).Hit;
    blocknumRepair = ratBEHstruct(s).blocknumRepair;
    if isempty(blocknumRepair); blocknumRepair = ratBEHstruct(s).blocknum; end
    if length(Hit) < 10; s = s+1; continue; end
    if ratBEHstruct(s).protocol~=7; s=s+1; continue; end;
    % get cued accuracy
    idxCued = find(blocknumRepair < 3); 
    idxWM = find(blocknumRepair>=3);
    if length(Hit) < 6; idxCued = ones(1,length(Hit)); end
    accPostMock(count+1) = sum(Hit(idxCued))/length(idxCued);
    accPostMockWM(count+1) = sum(Hit(idxWM))/length(idxWM);
    count = count+1; s = s+1;
end
% uni
count = 0; s = sUni(1); accPreUni = []; accPreUniWM = [];
while count < tBef
    Hit = ratBEHstruct(s).Hit;
    blocknumRepair = ratBEHstruct(s).blocknumRepair;
    if isempty(blocknumRepair); blocknumRepair = ratBEHstruct(s).blocknum; end
    if length(Hit) < 10; s = s-1; continue; end
    if ratBEHstruct(s).protocol~=7; s=s-1; continue; end;
    % get cued accuracy
    idxCued = find(blocknumRepair < 3);
    idxWM = find(blocknumRepair>=3);
    if length(Hit) < 6; idxCued = ones(1,length(Hit)); end
    accPreUni(count+1) = sum(Hit(idxCued))/length(idxCued);
    accPreUniWM(count+1) = sum(Hit(idxWM)) / length(idxWM);
    count = count+1; s = s-1;
end
count = 0; s = sUni(2); accPostUni = []; ssave = []; accPostUniWM = [];
while count < tAft
    Hit = ratBEHstruct(s).Hit;
    blocknumRepair = ratBEHstruct(s).blocknumRepair;
    if isempty(blocknumRepair); blocknumRepair = ratBEHstruct(s).blocknum; end
    if length(Hit) < 2; s = s+1; continue; end
    if ratBEHstruct(s).protocol~=7; s=s+1; continue; end;
    % get cued accuracy
    idxCued = find(blocknumRepair < 3); 
    idxWM = find(blocknumRepair>=3);
    if length(Hit) < 6; idxCued = ones(1,length(Hit)); end
    if length(Hit)<6; s=s+1; continue; end
   % if isempty(idxWM); disp('here'); end
    accPostUni(count+1) = sum(Hit(idxCued))/length(idxCued); ssave(count+1) = s;
    accPostUniWM(count+1) = sum(Hit(idxWM))/length(idxWM);
    count = count+1; s = s+1;
end
accPostUniWM(isnan(accPostUniWM))=0; % no trials, so 0
accPostUni(isnan(accPostUni))=0; % no trials, so 0
% bi
count = 0; s = sBi(1); accPreBi = [];  ssave = []; accPreBiWM = [];
while count < tBef
    Hit = ratBEHstruct(s).Hit;
    blocknumRepair = ratBEHstruct(s).blocknumRepair;
    if isempty(blocknumRepair); blocknumRepair = ratBEHstruct(s).blocknum; end
    if length(Hit) < 2; s = s-1; continue; end
    if ratBEHstruct(s).protocol~=7; s=s-1; continue; end;
    % get cued accuracy
    idxCued = find(blocknumRepair < 3); 
    idxWM = find(blocknumRepair>=3);
    if length(Hit) < 6; idxCued = ones(1,length(Hit)); end
    %if length(Hit)<6; s=s-1; continue; end
    accPreBi(count+1) = sum(Hit(idxCued))/length(idxCued); ssave(count+1) = s;
    accPreBiWM(count+1) = sum(Hit(idxWM)) / length(idxWM); 
    %if isnan(sum(Hit(idxWM)) / length(idxWM)); figure; waitforbuttonpress; end
    count = count+1; s = s-1; 
end
accPreBiWM(isnan(accPreBiWM))=0; % no trials, so 0

count = 0; s = sBi(2); accPostBi = []; accPostBiWM = [];
while count < tAft
    Hit = ratBEHstruct(s).Hit;
    blocknumRepair = ratBEHstruct(s).blocknumRepair;
    if isempty(blocknumRepair); blocknumRepair = ratBEHstruct(s).blocknum; end
    if length(Hit) < 2; s = s+1; continue; end
    if ratBEHstruct(s).protocol~=7; s=s+1; continue; end;
    % get cued accuracy
    idxCued = find(blocknumRepair < 3); 
    idxWM = find(blocknumRepair>=3);
    if length(Hit) < 6; idxCued = ones(1,length(Hit)); end
    accPostBi(count+1) = sum(Hit(idxCued))/length(idxCued);
    accPostBiWM(count+1) = sum(Hit(idxWM))/length(idxWM);
    count = count+1; s = s+1;
end
accPostBiWM(isnan(accPostBiWM))=0; % no trials, so 0
accPostBi(isnan(accPostBi))=0; % no trials, so 0

count = 0; s = length(ratBEHstruct); accInf = []; accInfWM = [];
% do 1 month insteod of length?
s = min(sBi(2)+3*28, length(ratBEHstruct)); % 3 sess * 7 days * 5 weeks
disp('doing 1 month for inf');
while count < tBef
    Hit = ratBEHstruct(s).Hit;
    blocknumRepair = ratBEHstruct(s).blocknumRepair;
    if isempty(blocknumRepair); blocknumRepair = ratBEHstruct(s).blocknum; end
    if length(Hit) < 2; s = s-1; continue; end
    if ratBEHstruct(s).protocol~=7; s=s-1; continue; end;
    % get cued accuracy
    idxCued = find(blocknumRepair < 3); 
    idxWM = find(blocknumRepair>=3);
    if length(Hit) < 6; idxCued = ones(1,length(Hit)); end
    accInf(count+1) = sum(Hit(idxCued))/length(idxCued);
    accInfWM(count+1) = sum(Hit(idxWM))/ length(idxWM);
    count = count+1; s = s-1;
end

% same as above but for OT
%* need to take unassist
count = 0; s = sMock(1); accPreMockOT = []; ssave = [];
while count < tBefOT
    Hit = ratBEHstruct(s).Hit;
    if length(Hit) < 10; s = s-1; continue; end
    if ratBEHstruct(s).protocol~=8 && ratBEHstruct(s).protocol~=6; s=s-1; continue; end;
    % get cued accuracy
    idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 1;
    if strcmp(ratBEHstruct(s).name, 'J7-Rat78') % fixing some bug
        idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 3; 
    end
    accPreMockOT(count+1) = sum(Hit(idx))/length(idx); ssave(count+1) = s;
    count = count+1; s = s-1;
end
count = 0; s = sMock(2); accPostMockOT = [];
while count < tAftOT
    Hit = ratBEHstruct(s).Hit;
    if length(Hit) < 10; s = s+1; continue; end
    if ratBEHstruct(s).protocol~=8 && ratBEHstruct(s).protocol~=6; s=s+1; continue; end;
    % get cued accuracy
    idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 1;
    if strcmp(ratBEHstruct(s).name, 'J7-Rat78')
        idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 3;
    end
    accPostMockOT(count+1) =  sum(Hit(idx))/length(idx);
    count = count+1; s = s+1;
end
% uni
count = 0; s = sUni(1); accPreUniOT = []; ssave = [];
while count < tBefOT
    Hit = ratBEHstruct(s).Hit;
    if length(Hit) < 10; s = s-1; continue; end
    if ratBEHstruct(s).protocol~=8 && ratBEHstruct(s).protocol~=6; s=s-1; continue; end;
    % get cued accuracy
    idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 1;
    if strcmp(ratBEHstruct(s).name, 'J7-Rat78')
        idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 3;
    end
    accPreUniOT(count+1) = sum(Hit(idx))/length(idx); ssave(count+1) = s;
    count = count+1; s = s-1;
end
count = 0; s = sUni(2); accPostUniOT = [];
while count < tAftOT
    Hit = ratBEHstruct(s).Hit;
    if length(Hit) < 2; s = s+1; continue; end
    if ratBEHstruct(s).protocol~=8 && ratBEHstruct(s).protocol~=6; s=s+1; continue; end;
    % get cued accuracy
    idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 1;
    if strcmp(ratBEHstruct(s).name, 'J7-Rat78')
        idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 3;
    end
    accPostUniOT(count+1) =  sum(Hit(idx))/length(idx);
    count = count+1; s = s+1;
end
accPostUniOT(isnan(accPostUniOT))=0; % no trials, so 0

% bi
count = 0; s = sBi(1); accPreBiOT = [];
while count < tBefOT
    Hit = ratBEHstruct(s).Hit;
    if length(Hit) < 2; s = s-1; continue; end
    if ratBEHstruct(s).protocol~=8 && ratBEHstruct(s).protocol~=6; s=s-1; continue; end;
    % get cued accuracy
    idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 1;
    if strcmp(ratBEHstruct(s).name, 'J7-Rat78')
        idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 3;
    end
    accPreBiOT(count+1) = sum(Hit(idx))/length(idx);
    count = count+1; s = s-1;
end
count = 0; s = sBi(2); accPostBiOT = [];
while count < tAftOT
    Hit = ratBEHstruct(s).Hit;
    if length(Hit) < 2; s = s+1; continue; end
    if ratBEHstruct(s).protocol~=8 && ratBEHstruct(s).protocol~=6; s=s+1; continue; end;
    % get cued accuracy
    idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 1;
    if strcmp(ratBEHstruct(s).name, 'J7-Rat78')
        idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 3;
    end
    accPostBiOT(count+1) =  sum(Hit(idx))/length(idx);
    count = count+1; s = s+1;
end
count = 0; s = length(ratBEHstruct); accInfOT = [];
s = min(sBi(2)+1*28, length(ratBEHstruct));
disp('doing month for inf for OT also');
while count < tAftOT
    Hit = ratBEHstruct(s).Hit;
    if length(Hit) < 2; s = s-1; continue; end
    if ratBEHstruct(s).protocol~=8 && ratBEHstruct(s).protocol~=6; s=s-1; continue; end;
    % get cued accuracy
    idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 1;
    if strcmp(ratBEHstruct(s).name, 'J7-Rat78')
        idx = cellfun(@length, ratBEHstruct(s).cuedNames) <= 3;
    end
    accInfOT(count+1) =  sum(Hit(idx))/length(idx);
    count = count+1; s = s-1;
end
accPostBiOT(isnan(accPostBiOT))=0; % no trials, so 0

% average from sessions to days
accPreMockOT = mean(reshape(accPreMockOT, [num_avg_OT, length(accPreMockOT)/num_avg_OT]),1);
accPostMockOT = mean(reshape(accPostMockOT, [num_avg_OT, length(accPostMockOT)/num_avg_OT]),1);
accPreUniOT = mean(reshape(accPreUniOT, [num_avg_OT, length(accPreUniOT)/num_avg_OT]),1);
accPostUniOT = mean(reshape(accPostUniOT, [num_avg_OT, length(accPostUniOT)/num_avg_OT]),1);
accPreBiOT = mean(reshape(accPreBiOT, [num_avg_OT, length(accPreBiOT)/num_avg_OT]),1);
accPostBiOT = mean(reshape(accPostBiOT, [num_avg_OT, length(accPostBiOT)/num_avg_OT]),1);
accInfOT = mean(reshape(accInfOT, [num_avg_OT, length(accInfOT)/num_avg_OT]),1);

accPreMock = mean(reshape(accPreMock, [num_avg, length(accPreMock)/num_avg]),1);
accPostMock = mean(reshape(accPostMock, [num_avg, length(accPostMock)/num_avg]),1);
accPreUni = mean(reshape(accPreUni, [num_avg, length(accPreUni)/num_avg]),1);
accPostUni = mean(reshape(accPostUni, [num_avg, length(accPostUni)/num_avg]),1);
accPreBi = mean(reshape(accPreBi, [num_avg, length(accPreBi)/num_avg]),1);
accPostBi = mean(reshape(accPostBi, [num_avg, length(accPostBi)/num_avg]),1);
accInf = mean(reshape(accInf, [num_avg, length(accInf)/num_avg]),1);

accPreMockWM = mean(reshape(accPreMockWM, [num_avg, length(accPreMockWM)/num_avg]),1);
accPostMockWM = mean(reshape(accPostMockWM, [num_avg, length(accPostMockWM)/num_avg]),1);
accPreUniWM = mean(reshape(accPreUniWM, [num_avg, length(accPreUniWM)/num_avg]),1);
accPostUniWM = mean(reshape(accPostUniWM, [num_avg, length(accPostUniWM)/num_avg]),1);
accPreBiWM = mean(reshape(accPreBiWM, [num_avg, length(accPreBiWM)/num_avg]),1);
accPostBiWM = mean(reshape(accPostBiWM, [num_avg, length(accPostBiWM)/num_avg]),1);
accInfWM = mean(reshape(accInfWM, [num_avg, length(accInfWM)/num_avg]),1);


% get accuracy at infinity?
dayoff = 0;
% clean up and smooth, as requested by bence
% - not a fan of this stuff....
 accPreMock(accPreMock==1) = []; accPreMock = smooth(accPreMock,3); accPreMock = accPreMock((1:num_actual)+dayoff)'; accPreMock = accPreMock(end:-1:1);
accPostMock(accPostMock==1) = []; accPostMock = smooth(accPostMock,3); accPostMock = accPostMock((1:num_actual)+dayoff)';
 accPreUni(accPreUni==1) = []; accPreUni = smooth(accPreUni,3); accPreUni = accPreUni((1:num_actual)+dayoff)'; accPreUni = accPreUni(end:-1:1);
accPostUni(accPostUni==1) = []; accPostUni = smooth(accPostUni,3); accPostUni = accPostUni((1:num_actual)+dayoff)';
 accPreBi(accPreBi==1) = []; accPreBi = smooth(accPreBi,3); accPreBi = accPreBi((1:num_actual)+dayoff)'; accPreBi = accPreBi(end:-1:1);
accPostBi(accPostBi==1) = []; accPostBi = smooth(accPostBi,3); accPostBi = accPostBi((1:num_actual)+dayoff)';
accInf(accInf==1) = []; accInf = smooth(accInf(1:num_actual),3);

 accPreMockOT(accPreMockOT==1 | accPreMockOT==0) = []; accPreMockOT = smooth(accPreMockOT,3); accPreMockOT = accPreMockOT((1:num_actual)+dayoff)'; accPreMockOT = accPreMockOT(end:-1:1);
accPostMockOT(accPostMockOT==1 | accPostMockOT==0) = []; accPostMockOT = smooth(accPostMockOT,3); accPostMockOT = accPostMockOT((1:num_actual)+dayoff)';
 accPreUniOT(accPreUniOT==1 | accPreUniOT==0) = []; accPreUniOT = smooth(accPreUniOT,3); accPreUniOT = accPreUniOT((1:num_actual)+dayoff)'; accPreUniOT = accPreUniOT(end:-1:1);
accPostUniOT(accPostUniOT==1) = []; accPostUniOT = smooth(accPostUniOT,3); accPostUniOT = accPostUniOT((1:num_actual)+dayoff)';
accPreBiOT(accPreBiOT==1) = [];accPreBiOT = smooth(accPreBiOT,3); accPreBiOT = accPreBiOT((1:num_actual)+dayoff)'; accPreBiOT = accPreBiOT(end:-1:1); 
accPostBiOT(accPostBiOT==1) = []; accPostBiOT = smooth(accPostBiOT,3); accPostBiOT = accPostBiOT((1:num_actual)+dayoff)';
accInfOT(accInfOT==1) = []; accInfOT = accInfOT(1:num_actual);

 accPreMockWM(accPreMockWM==1) = []; accPreMockWM = smooth(accPreMockWM,3); accPreMockWM = accPreMockWM((1:num_actual)+dayoff)'; accPreMockWM = accPreMockWM(end:-1:1);
accPostMockWM(accPostMockWM==3) = []; accPostMockWM = smooth(accPostMockWM,3); accPostMockWM = accPostMockWM((1:num_actual)+dayoff)';
 accPreUniWM(accPreUniWM==1 | isnan(accPreUniWM)) = []; accPreUniWM = smooth(accPreUniWM,3); accPreUniWM = accPreUniWM((1:num_actual)+dayoff)';accPreUniWM = accPreUniWM(end:-1:1);
accPostUniWM(accPostUniWM==1) = []; accPostUniWM = smooth(accPostUniWM,1); accPostUniWM = accPostUniWM((1:num_actual)+dayoff)';
 accPreBiWM(accPreBiWM==1) = []; accPreBiWM = smooth(accPreBiWM,3); accPreBiWM = accPreBiWM((1:num_actual)+dayoff)';accPreBiWM = accPreBiWM(end:-1:1);
accPostBiWM(accPostBiWM==1) = []; accPostBiWM = smooth(accPostBiWM,1); accPostBiWM = accPostBiWM((1:num_actual)+dayoff)';
accInfWM(accInfWM==1) = []; accInfWM = smooth(accInfWM(1:num_actual),3);

% normalize?
% accPreMock = accPreMock ./ mean([accPreMock, accPostMock]);
% accPostMock = accPreMock ./ mean([accPreMock, accPostMock]);

% make plot and save data?
ratMockCued_rat(:,rat) = [accPreMock, accPostMock];
ratUniCued_rat(:,rat) = [accPreUni, accPostUni];
ratBiCued_rat(:,rat) = [accPreBi, accPostBi];
ratInf_rat(:,rat) = mean(accInf);

ratMockOT_rat(:,rat) = [accPreMockOT, accPostMockOT];
ratUniOT_rat(:,rat) = [accPreUniOT, accPostUniOT];
ratBiOT_rat(:,rat) = [accPreBiOT, accPostBiOT];
ratInfOT_rat(:,rat) = mean(accInfOT);

ratMockWM_rat(:,rat) = [accPreMockWM, accPostMockWM];
ratUniWM_rat(:,rat) = [accPreUniWM, accPostUniWM];
ratBiWM_rat(:,rat) = [accPreBiWM, accPostBiWM];
ratInfWM_rat(:,rat) = mean(accInfWM);



end% rat

% normalize?


% figure
% figure; hold on;
% errorbar(-tBef:-1, mean(ratMockCued_rat(1:tBef, :),2), std(ratMockCued_rat(1:tBef, :),[],2),'r')
% errorbar(1:tAft, mean(ratMockCued_rat(end-tAft+1:end, :),2), std(ratMockCued_rat(end-tAft+1:end, :),[],2),'r')
% 

% prep for plotting in the other script
days_before = 7;
num_days_plot = 7;
acc_mock_aligned = ratMockCued_rat;
acc_uni_aligned = ratUniCued_rat;
acc_bi_aligned = ratBiCued_rat; 
acc_inf = ratInf_rat;

acc_mock_aligned = ratMockOT_rat;
acc_uni_aligned = ratUniOT_rat;
acc_bi_aligned = ratBiOT_rat; 
acc_inf = ratInfOT_rat;


acc_mock_aligned = ratMockWM_rat;
acc_uni_aligned = ratUniWM_rat;
acc_bi_aligned = ratBiWM_rat; 
acc_inf = ratInfWM_rat;

% set retrianing rats to 0
% ratUniCued_rat((days_before+1):end,:) = ratUniCued_rat((days_before+1):end,:).*~logical(uniRetrain);
% ratBiCued_rat((days_before+1):end,:) = ratBiCued_rat((days_before+1):end,:).*~logical(biRetrain);
% ratUniWM_rat((days_before+1):end,:) = ratUniWM_rat((days_before+1):end,:).*~logical(uniRetrain);
% ratBiWM_rat((days_before+1):end,:) = ratBiWM_rat((days_before+1):end,:).*~logical(biRetrain);
% ratUniOT_rat((days_before+1):end,:) = ratUniOT_rat((days_before+1):end,:).*~logical(uniRetrainOT);
% ratBiOT_rat((days_before+1):end,:) = ratBiOT_rat((days_before+1):end,:).*~logical(biRetrainOT);
% 
% 



%% plot script no uni
% - only pre, post, late (maybe not late)
% - normalized to pre lesion?
 % - i really think this wont look any better...
% - bar plot

% - this is just pre-lesion, post lesion, late

% redo this for the OT only rats?

days_before = 7;
num_days_plot = 7;


% first try plotting all 3...
figure('Position',[100,100,1000,400]); hold on;
offset_lesion = num_days_plot+days_before+2;
cbehuse = {[.2,.8,.2]; [.8,.2,.2];[.2,.2,.8]};
nrats = length(ratPath);

noretrain = 0; % remove rats that had to be retrained (set to 0)
donorm = 0;
dodrop = 0;
showindividuals = 1;

doMock = 0;

p_all = [];
valBarSave = [];

for context = [2,3,1]%[1,2,3];%[2,3,1]
    
    if context == 3
        acc_uni_aligned = ratUniWM_rat;
        acc_bi_aligned = ratBiWM_rat; 
        acc_inf = ratInfWM_rat;
        %acc_early = ratEarlyWM_rat;
    elseif context==2
        acc_uni_aligned = ratUniCued_rat;
        acc_bi_aligned = ratBiCued_rat; 
        acc_inf = ratInf_rat;
        %acc_early = ratEarly_rat;
    elseif context==1
        acc_uni_aligned = ratUniOT_rat;
        acc_bi_aligned = ratBiOT_rat; 
        acc_inf = ratInfOT_rat;
        %acc_early = ratEarlyOT_rat;
    end
    if doMock
        if context==3
            acc_uni_aligned = ratMockWM_rat;
            acc_bi_aligned = ratMockWM_rat;
        elseif context==2
            acc_uni_aligned = ratMockCued_rat; 
            acc_bi_aligned = ratMockCued_rat;
        elseif context==1
            acc_uni_aligned = ratMockOT_rat;
            acc_bi_aligned = ratMockOT_rat;
        end
    end
    
    usecolor = cbehuse{context}; 
    if noretrain
        acc_uni_aligned(acc_uni_aligned==0) = nan;
        acc_bi_aligned(acc_bi_aligned==0) = nan;
    end
    if donorm
        val = mean(acc_uni_aligned(1:days_before,:));
        acc_uni_aligned = acc_uni_aligned ./val;
        acc_bi_aligned = acc_bi_aligned./val;
        acc_inf = acc_inf ./ val;
    end
    if dodrop % look at change in accuracy
        val = mean(acc_uni_aligned(1:days_before,:));
        acc_uni_aligned = acc_uni_aligned - val;
        acc_bi_aligned = acc_bi_aligned - val;
        acc_inf = acc_inf - val;
    end
    
% temp add in early
% errorbar(-days_before+offset_lesion-2, nanmean(acc_early(:)), nanstd(acc_early)./sqrt(nrats),...
%     '-','Color',usecolor,'linewidth',2,'CapSize',0,'MarkerFaceColor',usecolor);
% 
% unilateral
errorbar([-days_before:-1]+offset_lesion,nanmean(acc_uni_aligned(1:days_before,:),2),...
    nanstd(acc_uni_aligned(1:days_before,:),[],2)./sqrt(nrats),...
    '-','Color',usecolor,'linewidth',2,'CapSize',0,'MarkerFaceColor',usecolor);%'Color',usecolor,'linewidth',2);
errorbar([1:num_days_plot]+offset_lesion,nanmean(acc_bi_aligned((days_before+1):end,:),2),...
    nanstd(acc_bi_aligned((days_before+1):end,:),[],2)./sqrt(nrats),...
    '-','Color',usecolor,'linewidth',2,'CapSize',0,'MarkerFaceColor',usecolor);%'Color',usecolor,'linewidth',2);

errorbar(offset_lesion+num_days_plot+3, nanmean(acc_inf), nanstd(acc_inf)./sqrt(nrats),...
    '-','Color',usecolor,'linewidth',2,'CapSize',0,'MarkerFaceColor',usecolor);%'Color',usecolor,'linewidth',2);

for days_off_lesion = 1:length(nanmean(acc_bi_aligned((days_before+1):end,:),2))
    acc_pre = mean(acc_uni_aligned(1:days_before,:),1);
    acc_post = acc_bi_aligned((days_before+days_off_lesion),:);
    %p = signrank(acc_pre, acc_post);
    [~,p] = ttest(acc_pre, acc_post);
    if p<0.05
        plot((days_before+days_off_lesion)+num_days_plot+2, 1+context/10, '*','Color',usecolor);
    end
    p_all(context,days_off_lesion) = p;
end

% add in individual rats
% xvals = [-days_before:-1,1:num_days_plot,[-days_before:-1]+offset_lesion,...
%     [1:num_days_plot]+offset_lesion,[-days_before:-1]+2*offset_lesion,...
%     [1:num_days_plot]+2*offset_lesion,2*offset_lesion+num_days_plot+3 ];
% symbs_use = 'o+*xsd^v><ph';
if showindividuals;
usecolor = usecolor + .4; usecolor(usecolor>1)=1;
for rat = 1:size(acc_uni_aligned,2)%length(ratPath)
    
    %plot(-days_before:-1, acc_mock_aligned(1:days_before,rat),[ '-'],'Color',usecolor,'linewidth',.5)
    %plot(1:num_days_plot, acc_mock_aligned((days_before+1):end,rat),[ '-'],'Color',usecolor,'linewidth',.5)
    plot([-days_before:-1]+offset_lesion,acc_uni_aligned(1:days_before,rat),[ '-'],'Color',usecolor,'linewidth',.5)
   % plot([1:num_days_plot]+offset_lesion, acc_uni_aligned((days_before+1):end,rat),[ '-'],'Color',usecolor,'linewidth',.5)
    %plot([-days_before:-1]+2*offset_lesion,acc_bi_aligned(1:days_before,rat),[ '-'],'Color',usecolor,'linewidth',.5)
    plot([1:num_days_plot]+offset_lesion, acc_bi_aligned((days_before+1):end,rat),[ '-'],'Color',usecolor,'linewidth',.5)
    plot(offset_lesion+num_days_plot+3 ,acc_inf(rat),[ '-*'],'Color',usecolor,'linewidth',.5)
    
end
end

valBarSave(1,:,context) = [nanmean(acc_uni_aligned(1:days_before,:))];
valBarSave(2,:,context) = [nanmean(acc_bi_aligned((days_before+1):end,:))];

end
xlim([-days_before+offset_lesion-1, offset_lesion+num_days_plot+3+1])
 ylim([-.05,1]); ylabel('Accuracy');title('Unilateral'); xlabel('Days');
xticks([-days_before:-1,1:num_days_plot,[-days_before:-1]+offset_lesion,...
    [1:num_days_plot]+offset_lesion,[-days_before:-1]+2*offset_lesion,...
    [1:num_days_plot]+2*offset_lesion,2*offset_lesion+num_days_plot+3 ])
xticklabels([num2str(repmat([-days_before:-1,1:num_days_plot],1,3)')])
val = xticklabels; val{end} = 'Late';
if donorm; ylim([-.05, 1.3]) ; end
if dodrop; ylim([-1,0.3]); ylabel('accuracy drop'); end;

% significance test on a day to day basis?
pWMOTpre = []; pWMOTpost = [];
for j = 1:days_before
    [~,pWMOTpre(j)] = ttest(ratUniWM_rat(j, :), ratUniOT_rat(j,:));
    [~,pWMOTpost(j)] = ttest(ratBiWM_rat(days_before+j, :), ratBiOT_rat(days_before+j,:));
end
% norm
val = mean(ratUniWM_rat(1:days_before,:));
ratUniWM_norm = ratUniWM_rat((days_before+1):end,:)./val;
ratBiWM_norm = ratBiWM_rat((days_before+1):end,:)./val;
val = mean(ratUniOT_rat(1:days_before,:));
ratUniOT_norm = ratUniOT_rat((days_before+1):end,:)./val;
ratBiOT_norm = ratBiOT_rat((days_before+1):end,:)./val;

pWMOTpre = []; pWMOTpost = [];
for j = 1:days_before
    [~,pWMOTpre(j)] = ttest(ratUniWM_norm(j, :), ratUniOT_norm(j,:));
    [~,pWMOTpost(j)] = ttest(ratBiWM_norm(j, :), ratBiOT_norm(j,:));
end

% put sig stars on everyboyd?

%% Bar plot


accCued = [mean(ratUniCued_rat(1:days_before,:)); mean(ratBiCued_rat((days_before+1):end,:)); ratInf_rat];
accWM = [mean(ratUniWM_rat(1:days_before,:)); mean(ratBiWM_rat((days_before+1):end,:)); ratInfWM_rat];
accOT = [mean(ratUniOT_rat(1:days_before,:)); mean(ratBiOT_rat((days_before+1):end,:)); ratInfOT_rat];

% do mock
accCued = [mean(ratUniCued_rat(1:days_before,:)); mean(ratBiCued_rat((days_before+1):end,:)); ratInf_rat];
accWM = [mean(ratUniWM_rat(1:days_before,:)); mean(ratBiWM_rat((days_before+1):end,:)); ratInfWM_rat];
accOT = [mean(ratUniOT_rat(1:days_before,:)); mean(ratBiOT_rat((days_before+1):end,:)); ratInfOT_rat];

% do crop differently as requested by Bence
accCued = [mean(ratUniCued_rat(1:days_before,:)); mean(ratBiCued_rat((3:7)+days_before,:)); ratInf_rat];
accWM = [mean(ratUniWM_rat(1:days_before,:)); mean(ratBiWM_rat((3:7)+days_before,:)); ratInfWM_rat];
accOT = [mean(ratUniOT_rat(1:days_before,:)); mean(ratBiOT_rat((3:7)+days_before,:)); ratInfOT_rat];



figure; hold on;

bar(1, mean(accCued(1,:)),'FaceColor',[1,.2,.2]);
bar(2, mean(accCued(2,:)),'FaceColor',[1,.4,.4]);
bar(3, mean(accCued(3,:)),'FaceColor',[1,.7,.7]);

bar(5, mean(accWM(1,:)),'FaceColor',[.2,.2,1]);
bar(6, mean(accWM(2,:)),'FaceColor',[.4,.4,1]);
bar(7, mean(accWM(3,:)),'FaceColor',[.7,.7,1]);

bar(9,  mean(accOT(1,:)),'FaceColor',[.2,1,.2]);
bar(10, mean(accOT(2,:)),'FaceColor',[.4,1,.4]);
bar(11, mean(accOT(3,:)),'FaceColor',[.7,1,.7]);

for rat = 1:length(ratPath)
    plot(1:3, accCued(:,rat),'Color',[0,0,0,.4]);
    plot(5:7, accWM(:,rat),'Color',[.0,0,0,.4]);
    plot(9:11, accOT(:,rat),'Color',[0,0,0,.4]);
end

% ttest
[~,ppreunicued] = ttest(accCued(1,:),accCued(2,:));
[~,pprebicued] =  ttest(accCued(1,:),accCued(3,:));
[~,ppreuniwm] = ttest(accWM(1,:),accWM(2,:));
[~,pprebiwm] =  ttest(accWM(1,:),accWM(3,:));
[~,ppreuniot] = ttest(accOT(1,:),accOT(2,:));
[~,pprebiot] =  ttest(accOT(1,:),accOT(3,:));
% wilcoxon
ppreunicued = signrank(accCued(1,:),accCued(2,:));
pprebicued =  signrank(accCued(1,:),accCued(3,:));
ppreuniwm = signrank(accWM(1,:),accWM(2,:));
pprebiwm =  signrank(accWM(1,:),accWM(3,:));
ppreuniot = signrank(accOT(1,:),accOT(2,:));
pprebiot =  signrank(accOT(1,:),accOT(3,:));

sigstar({[1,2],[1,3],[5,6],[5,7],[9,10],[9,11]},...
    [ppreunicued,pprebicued,ppreuniwm,pprebiwm,ppreuniot,pprebiot],0,1);

% compare across conditions?
[~,p1] = ttest(accCued(1,:), accWM(1,:));
[~,p2] = ttest(accCued(1,:), accOT(1,:));
[~,p3] = ttest(accWM(1,:), accOT(1,:));

[p1] = signrank(accCued(1,:), accWM(1,:));
[p2] = signrank(accCued(1,:), accOT(1,:));
[p3] = signrank(accWM(1,:), accOT(1,:));
sigstar({[1,5], [1,9],[5,9]},[p1,p2,p3],0,1);
[~,p1] = ttest(accCued(2,:), accWM(2,:));
[~,p2] = ttest(accCued(2,:), accOT(2,:));
[~,p3] = ttest(accWM(2,:), accOT(2,:));
% 
[p1] = signrank(accCued(2,:), accWM(2,:));
[p2] = signrank(accCued(2,:), accOT(2,:));
[p3] = signrank(accWM(2,:), accOT(2,:));
sigstar({[2,6], [2,10],[6,10]},[p1,p2,p3],0,1);
[~,ptest] = ttest(accWM(1,:), accOT(2,:));

xticks([2,6,10]); ylim([0,1.5])
xticklabels({'Cued','WM','OT'});

