%% Script for calibrating multiple cameras in Matlab
% Jesse Marshall 2019
% - Kevin Mizes 2022

%% enter parent paths here
% path name
box = 'L2';
ratname = 'L2-Rat70';
vidpath = 'Z:\Kevin\Video\';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};
savepath  = 'D:\Kevin\Sequence_tap_analysis\basic_beh_traj_analysis\';

% path name
box = 'E5';
ratname = 'E5-Rat123';
vidpath = 'Z:\Kevin\Video\';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','T','L'};
savepath  = 'D:\Kevin\Sequence_tap_analysis\basic_beh_traj_analysis\';

% J2-Rat17
box = 'J2'; ratname = 'J2-Rat17';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','T','L'};

% F4-Rat43
box = 'F4'; ratname = 'F4-Rat43';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','T','L'};

% E4-Rat39
box = 'E4'; ratname = 'E4-Rat39';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'L','R','T'};

% L3-Rat34
box = 'L3'; ratname = 'L3-Rat34';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% J1-Rat71
box = 'J1'; ratname = 'J1-Rat71';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% L3-Rat65
box = 'L3'; ratname = 'L3-Rat65';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% L5-Rat6
box = 'L5'; ratname = 'L5-Rat6';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% J5-Rat31
% no images for this box??
box = 'J5'; ratname = 'J5-Rat31';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','T','L'};

% E2-Rat21
box = 'E2'; ratname = 'E2-Rat21';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% E8-Rat25
box = 'E8'; ratname = 'E8-Rat25';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','T','L'};

% J7-Rat45
box = 'J7'; ratname = 'J7-Rat45';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% - additional rats

% D7 rat 74 early in learning
box = 'D7'; ratname = 'D7-Rat74';
cams = {'Master','Slave','Slave2'};
cam_angles = {'R','L','T'};

% J7 rat 78 early in learning
box = 'J7'; ratname = 'J7-Rat78';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% E1 75 mayb
box = 'E1'; ratname = 'E1-Rat75';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','T','L'};

% l2 rat 116
box = 'L2'; ratname = 'L2-Rat116';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% J3 rat 113
box = 'J3'; ratname = 'J3-Rat113';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% d8 rat 119
box = 'D8'; ratname = 'D8-Rat119';
cams = {'Master','Slave','Slave2'};
cam_angles = {'R','L','T'};

% J6 rat 105 - no master tracking?

% e4 rat 122
box = 'E4'; ratname = 'E4-Rat122';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'L','R','T'};

% f3 rat 115
box = 'F3'; ratname = 'F3-Rat115';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','T','L'};

% l2 rat 132
box = 'L2'; ratname = 'L2-Rat132';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% f5 rat 86
box = 'F5'; ratname = 'F5-Rat86';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% lets try an ephys camera which has the best tracking
box = 'EphysE7'; ratname = 'EphysE7-Rat47';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

% other ephys boxes out of sheer curiousity?
% - to see if can decode the forelimb above chance
% - honestly probably if it just gives up on the z-dimension
box = 'EphysE6'; ratname = 'EphysE6-Rat35';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

box = 'EphysE7'; ratname = 'EphysE7-Rat33';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

box = 'EphysE7'; ratname = 'EphysE7-Rat18';
cams = {'Master','Slave1','Slave2'};
cam_angles = {'R','L','T'};

%% generate intrinsics
% - only use if recorded intrinsics on!
% over each cam
num2use = 100;
for cam = cams
    calibpath = fullfile(vidpath, box, cam{1}, 'calib');
    % if calibration exists, skip for now
    if exist(fullfile(calibpath, 'intrinsics.mat'))
        disp('intrinsics already found!');
        continue; 
    end
    
    % generate set of calibration images from the video
    disp('creating imgs from video of calibration checkerboard');
    if ~exist(fullfile(calibpath, 'imgFiles'))
        mkdir(fullfile(calibpath, 'imgFiles'));
        vidnames = dir(fullfile(calibpath, '*.mp4'));
        imageFiles = {}; count = 1;
        for j = 1:length(vidnames)
            v = VideoReader([calibpath filesep vidnames(j).name]);
            while hasFrame(v)
                img = readFrame(v);
                imwrite(img,[calibpath filesep 'imgFiles' filesep sprintf('img%05d.png',count)]);
                count = count+1;
            end
        end
    end
    disp('calculating intrinsics')
    % subsamples frames
    imFiles1 = imageDatastore(fullfile(calibpath, 'imgFiles'));
    imUse1 = randsample(length(imFiles1.Files), num2use);
    % detect calibration pattern
    [imagePoints, boardSize, imagesUsed] = detectCheckerboardPoints(imFiles1.Files(imUse1));
    % generate world coordinates
    squareSize = 6.1;  % mm
    worldPoints = generateCheckerboardPoints(boardSize,squareSize);
    % calibrate the camera
    I = readimage(imFiles1,1); 
    imageSize = [size(I,1),size(I,2)];
    [params,~,estimationErrors] = estimateCameraParameters(imagePoints,worldPoints, ...
                                  'ImageSize',imageSize);
    % plot errors (extrinsics and reprojection errors)
    h100=figure; showExtrinsics(params, 'CameraCentric');
    saveas(gca,fullfile(calibpath, 'extrinsics_errors.png'),'png'); close(h100);
    h100=figure; showReprojectionErrors(params);
    saveas(gca,fullfile(calibpath, 'reprojection_errors.png'),'png'); close(h100);
    % save intrinsics
    save(fullfile(calibpath, 'intrinsics.mat'),'params','imUse1','estimationErrors');
    
end

%% ENTER SIZE OF LFRAME HERE
% intrinsics = load('Z:\Kevin\Video\F3\CalibImgs\1-2.mat');
% cam1_intrinsics = intrinsics.params.CameraParameters1.Intrinsics;
% cam2_intrinsics = intrinsics.params.CameraParameters2.Intrinsics;

% skip, use only the L@ intrinsics...
%cam1_intrinsics = load(fullfile(vidpath, box, cams{1}, 'calib','intrinsics.mat'));
cam1_intrinsics = load(fullfile(vidpath, 'L2', cams{1}, 'calib','intrinsics.mat'));
if strcmp(cam_angles{2},'L')
    %cam2_intrinsics = load(fullfile(vidpath, box, cams{2}, 'calib','intrinsics.mat'));
    cam2_intrinsics = load(fullfile(vidpath, 'L2', cams{2}, 'calib','intrinsics.mat'));
else
    %cam2_intrinsics = load(fullfile(vidpath, box, cams{3}, 'calib','intrinsics.mat'));
    cam2_intrinsics = load(fullfile(vidpath, 'L2', cams{3}, 'calib','intrinsics.mat'));
end

% for ephys
cam1_intrinsics = load(fullfile(vidpath, 'EphysE7', 'Right_view','intrinsics.mat'));
cam2_intrinsics = load(fullfile(vidpath, 'EphysE7', 'Left_view','intrinsics.mat'));

cam1_intrinsics = cam1_intrinsics.params.Intrinsics;
cam2_intrinsics = cam2_intrinsics.params.Intrinsics;

in2mm = @(x) x*25.4;
% origin is bottom left lever hole
% 2 is bottom right lever hole
% 3 is bottom left handle hole
% 4 is bottom right handle hole
LFrame_coordinates = in2mm([0 0 0; 0 0 3.5;  0 -.621 .75; 0 -.621 2.72]);
intrinsics = {cam1_intrinsics, cam2_intrinsics};

%% STEP 4 - label lframe
% cam1 = "Z:\Kevin\Video\L2\Master\imDirRat70";
% cam2 = "Z:\Kevin\Video\L2\Slave1\imDirRat70";
% cam1 = 'Z:\Kevin\Video\E5\Master\imDir-Rat123\';
% cam2 = 'Z:\Kevin\Video\E5\Slave2\imDir-Rat123\';

 cam1 = fullfile(vidpath, box, cams{contains(cam_angles, 'R')},['imDir',ratname(3:end)]);
 cam2 = fullfile(vidpath, box, cams{contains(cam_angles, 'L')},['imDir',ratname(3:end)]);

%cam1 = fullfile(vidpath, box, cams{contains(cam_angles, 'R')},['imDir',ratname(8:end)]);
%cam2 = fullfile(vidpath, box, cams{contains(cam_angles, 'L')},['imDir',ratname(8:end)]);


% programatically get path from the info

im1 = imread(fullfile(cam1, 'img00089.png'));
im2 = imread(fullfile(cam2, 'img00088.png'));


images = {im1, im2};
numcams = numel(images);
point_coordinates = cell(1,numcams);

for kk = [1:numcams]%1:numcams
    camparams = intrinsics{kk};
    camparams  = cameraParameters('IntrinsicMatrix', camparams.IntrinsicMatrix);
    figure(8077)
    LFrame_image{kk} = undistortImage(images{kk},camparams);
    image(LFrame_image{kk})
    [xi,yi] = getpts ;
    point_coordinates{kk} = [xi yi];
end
for kk = 1:numcams
    point_coordinates{kk} = reshape(point_coordinates{kk},[],2);
    point_coordinates{kk} = point_coordinates{kk}(1:size(LFrame_coordinates,1),:);
end

if ~exist(fullfile(savepath, ratname, 'calib'))
    mkdir(fullfile(savepath, ratname, 'calib'))
end

save(fullfile(savepath, ratname, 'calib','calib.mat'),...
    'images','LFrame_image','point_coordinates','intrinsics');%,'estimationErrors');
%save('calib','images','LFrame_image','point_coordinates','intrinsics','estimationErrors');


%% Step 5 compute the world orientation
worldOrientation = cell(1,numcams);
worldLocation = cell(1,numcams);
rotationMatrix = cell(1,numcams);
translationVector = cell(1,numcams);
% get the orienation and location of cameras
% if ~exist('calib', 'dir')
%     mkdir('calib')
% end
%saveFolder = 'calib';
saveFolder = fullfile(savepath, ratname, 'calib');
for kk = [1:numcams]%:numcams
    
    [worldOrientation{kk},worldLocation{kk}] = estimateWorldCameraPose(double(point_coordinates{kk}),...
        double(LFrame_coordinates),intrinsics{kk}, "MaxReprojectionError", 50);
%     [worldOrientation{kk},worldLocation{kk}] = estimateWorldCameraPose(double(point_coordinates{kk}),...
%         double(LFrame_coordinates),intrinsics{kk},'Confidence',99,'MaxReprojectionError',4,'MaxNumTrials',5000);
    [rotationMatrix{kk},translationVector{kk}] = cameraPoseToExtrinsics(worldOrientation{kk},worldLocation{kk});
    
    figure(222)
    plotCamera('Location',worldLocation{kk},'Orientation',worldOrientation{kk},'Size',50,'Label',num2str(kk));
    hold on
    if (kk == numcams)
        print('-dpng',strcat(saveFolder,filesep,'cameraarrangement.png'));
    end
    %[R,t] = extrinsics(double(point_coordinates{kk}),double(LFrame_coordinates(:,1:2)),camparams);
    
    % look at reprojection in image from estimated world pose
    figure(233+kk)
    image( LFrame_image{kk})
    hold on
    
    imagePoints = worldToImage(intrinsics{kk},rotationMatrix{kk},translationVector{kk},double(LFrame_coordinates));
    
    colorarray = {'r','g','b','y'};
    for llll = 1:size(imagePoints,1)
        plot(imagePoints(llll,1),imagePoints(llll,2),'o','MarkerSize',4,...
            'MarkerEdgeColor',colorarray{llll},'MarkerFaceColor',colorarray{llll})
    end
    hold off
    print('-dpng',strcat(saveFolder,filesep,'camerareproject',num2str(kk),'.png'));
    
end
fprintf('saving world coordinates \n')
save(fullfile(savepath, ratname, 'calib','extrinsics.mat'),'intrinsics',...
    'worldOrientation','worldLocation','rotationMatrix','translationVector');

%% validate reprojection error on a large set of forelimbs
% - this might need to be a new scrip
% - for each rat, flip through and load the kinematics from both sides (not
% just one side)
% - then load the extrniscs to calculate world coordinates
% - then compute the reprojection error

%% load some matched points to test

load('D:\Kevin\Sequence_tap\L2_output\Results-L2-Rat70\ratBehTrajSess\513.mat')

nose_master_ex1 = ratBehTrajStructSess.trajMaster(10).nose{1}(200:320,:);
nose_slave1_ex1 = ratBehTrajStructSess.trajSlave1(10).nose{1}(200:320,:);
matchedPoints1 = nose_master_ex1;
matchedPoints2 = nose_slave1_ex1;

% paw_master_ex1 = ratBehTrajStructSess.trajMaster(10).pawL{1}(200:320,:);
% paw_slave1_ex1 = ratBehTrajStructSess.trajSlave1(10).pawL{1}(200:320,:);
% matchedPoints1 = paw_master_ex1;
% matchedPoints2 = paw_slave1_ex1;

% traingulate
camMatrix1 = cameraMatrix(intrinsics{1},rotationMatrix{1},translationVector{1});
camMatrix2 = cameraMatrix(intrinsics{2},rotationMatrix{2},translationVector{2});

worldPoints = triangulate(matchedPoints1,matchedPoints2,camMatrix1,camMatrix2);

figure;plot3(worldPoints(:,1), worldPoints(:,2), worldPoints(:,3))


%% calculate reprojectino error error

matchedPoints1_reproject = worldToImage(intrinsics{1},rotationMatrix{1},translationVector{1},worldPoints);
matchedPoints2_reproject = worldToImage(intrinsics{2},rotationMatrix{2},translationVector{2},worldPoints);

% calculate diff?
figure; hold on;
plot(matchedPoints1,'r'); plot(matchedPoints2,'b');
plot(matchedPoints1_reproject,'m');
plot(matchedPoints2_reproject,'c');

% compute average pixel distance off
errDist1 = mean(sqrt(sum((matchedPoints1 - matchedPoints1_reproject).^2,2)));
errDist2 = mean(sqrt(sum((matchedPoints2 - matchedPoints2_reproject).^2,2)));


err1 = mean((matchedPoints1 - matchedPoints1_reproject).^2);
err2 = mean((matchedPoints2 - matchedPoints2_reproject).^2);
err = (sum((matchedPoints1 - matchedPoints1_reproject).^2) + sum((matchedPoints2 - matchedPoints2_reproject).^2))/4;

% ok, need to do this for every single rat, and put into their trajectory
% folder
% - then use that to conver all trajectories to 3d and compute stuff on
% that??

%** pixel error got worse when I had the new intrinsics...




%% STEP 6 - refine with checkerboard
fprintf('starting checkerboard analysis:\n MAKE SURE ALL POINTS ARE TRACKED AND RED POINT IS OPPOSITE THE BEND OF THE L AND BLUE POINT IS ON SHORT AXIS\n If this is a problem, know that you dont need to click on the L exactly, but can exaggerate the positions. In the code, consider making shortdim equal to the opposite of longdim')
%cam1 =num2str(basecam');
%cam2 = {'2','2','3','4','5','6'};
%params_comp = cell(1,numel(cam1));
%params_comp{kk} = params;

load(lframename_checkerboard)
allimagepoints = cell(1,numcams);
for kk = [1:numcams]%1:numcams
    camparams =intrinsics{kk};
    
    undistorted_checkerboard =  undistortImage(checkerboard_images{kk},camparams);
    
    fprintf('finding checkerboard points \n')
    [imagePoints_ch, boardSize_ch, imagesUsed] = ...
        detectCheckerboardPoints(undistorted_checkerboard,'MinCornerMetric',0.15);
    if size(imagePoints_ch,1)<40
        [imagePoints_ch, boardSize_ch, imagesUsed] = ...
            detectCheckerboardPoints(undistorted_checkerboard,'MinCornerMetric',0.15);
    end
    squareSize_ch = 16.57; %mm
    worldPoints = generateCheckerboardPoints(boardSize_ch,squareSize_ch);
    figure(8078)
    imshow(undistorted_checkerboard);
    hold on;
    [xi,yi] = getpts ;
    point_coordinates_ch{kk} = [xi yi];
    point_coordinates_ch{kk} = reshape(point_coordinates_ch{kk},[],2);
    point_coordinates_ch{kk} =point_coordinates_ch{kk}(1:4,:);
    origin = point_coordinates_ch{kk}(1,:);
    long_axis = -(point_coordinates_ch{kk}(3,:)-point_coordinates_ch{kk}(1,:));
    short_axis = -(point_coordinates_ch{kk}(4,:)-point_coordinates_ch{kk}(1,:));
    
    detected_long_axis = -(imagePoints_ch(1,:)-imagePoints_ch(min(boardSize_ch),:));
    detected_short_axis = -(imagePoints_ch(1,:)-imagePoints_ch(min(boardSize_ch)-1,:));
    
    long_axis_dot = dot(long_axis,detected_long_axis)./(norm(long_axis)*norm(detected_long_axis));
    short_axis_dot =dot(short_axis,detected_short_axis)./(norm(short_axis)*norm(detected_short_axis));
    disp(long_axis_dot)
    disp(short_axis_dot)
    
    [~,longdim] = max(abs(long_axis));
    [~,shortdim] = max(abs(short_axis));
    flipped_undistorted_checkerboard=undistorted_checkerboard;
    pointindex = 1:((boardSize_ch(1)-1)*((boardSize_ch(2)-1)));
    pointindex = reshape(pointindex,min(boardSize_ch)-1,[]);
    
    if long_axis_dot<0
        fprintf('flipped long \n')
        pointindex = flip(pointindex,longdim);
    end
    
    if short_axis_dot<0
        fprintf('flipped short \n')
        pointindex = flip(pointindex,shortdim);
    end
    
    pointindex_unflipped = reshape(pointindex,[],1);
    %   [imagePoints_ch, boardSize_ch, imagesUsed] = ...
    %         detectCheckerboardPoints(flipped_undistorted_checkerboard,'MinCornerMetric',0.15);
    imagePoints_ch = imagePoints_ch(pointindex_unflipped,:);
    
    
    figure(8078+kk)
    imshow(undistorted_checkerboard);
    hold on;
    plot(imagePoints_ch(1,1), imagePoints_ch(1,2),'ro');
    plot(imagePoints_ch(2,1), imagePoints_ch(2,2),'bo');
    plot(imagePoints_ch(3:end,1), imagePoints_ch(3:end,2),'go');
    hold off
    %plot(imagePoints_ch(:,1), imagePoints_ch(:,2),'go');
    allimagepoints{kk} = imagePoints_ch;
    checkerboard_images_undistorted{kk} = undistorted_checkerboard;
    boardSize_ch_full{kk} = boardSize_ch;
end
save(lframename_checkerboard,'checkerboard_images_undistorted','allimagepoints',...
    'intrinsics','estimationErrors','boardSize_ch_full');



%% STEP 7 compute the corrected world coordinates
for kk = 1:numcams
    camMatrix{kk} = cameraMatrix(intrinsics{kk},rotationMatrix{kk},translationVector{kk});
end
% camMatrix{1} = cameraMatrix(intrinsics{1},rotationMatrix_2{1},translationVector_2{1});

minpts = min(size(allimagepoints{1},1),size(allimagepoints{2},1));
point3d = triangulate(allimagepoints{2}(1:minpts,:), allimagepoints{3}(1:minpts,:),...
    camMatrix{2},camMatrix{3});

for    kk=1:numcams
    [X,Y] = meshgrid(1:(max(boardSize_ch_full{kk})-1),1:(min(boardSize_ch_full{kk})-1));
    y_step = squareSize;
    x_step = squareSize;
    location3d = [(reshape(Y,[],1)-1)*y_step (reshape(X,[],1)-1)*x_step  zeros(numel(X),1)];
    
    [worldOrientation_2{kk},worldLocation_2{kk}] = estimateWorldCameraPose(double(allimagepoints{kk}),...
        double(location3d),intrinsics{kk},'Confidence',95,'MaxReprojectionError',4,'MaxNumTrials',5000);
    
    [rotationMatrix_2{kk},translationVector_2{kk}] = cameraPoseToExtrinsics(worldOrientation_2{kk},worldLocation_2{kk});
    
    figure(333+kk)
    image( checkerboard_images_undistorted{kk})
    hold on
    imagePoints = worldToImage(intrinsics{kk},rotationMatrix{kk},...
        translationVector{kk},double(point3d));
    imagePoints2 = worldToImage(intrinsics{kk},rotationMatrix_2{kk},...
        translationVector_2{kk},double(location3d));
    plot(allimagepoints{kk}(:,1),allimagepoints{kk}(:,2),'or')
    plot(imagePoints2(:,1),imagePoints2(:,2),'og')
    plot(imagePoints(:,1),imagePoints(:,2),'ok')
    
    hold off
    legend('Ground truth','Checkerboard','LFrame')
    print('-dpng',strcat(saveFolder,filesep,'camerareproject_checkerboard',num2str(kk),'.png'));
    
    
    
    figure(227)
    plotCamera('Location',worldLocation_2{kk},'Orientation',worldOrientation_2{kk},'Size',50,'Label',num2str(kk));
    hold on
    view([-91 84])
    if (kk == numcams)
        print('-dpng',strcat(saveFolder,filesep,'cameraarrangement_2.png'));
    end
    %
    
end

fprintf('saving world coordinates \n')
save( worldcoordinates_framename_2,'intrinsics','worldOrientation_2','worldLocation_2','rotationMatrix_2','translationVector_2');

%% plot L frame on top of images
ratname = 'L2-Rat70';
numcams = 2;
savepath = 'D:\Kevin\Sequence_tap_analysis\basic_beh_traj_analysis\';
load(fullfile(savepath, ratname, 'calib','calib.mat'))
saveFolder = fullfile(savepath, ratname, 'calib');

for kk = [1:numcams]
    figure(233+kk)
    %image( LFrame_image{kk})
    imshow(LFrame_image{kk})
    hold on
    
    imagePoints = point_coordinates{kk};    
    colorarray = {'r','g','b','y','m','c'};
    for llll = 1:size(imagePoints,1)
        plot(imagePoints(llll,1),imagePoints(llll,2),'o','MarkerSize',8,...
            'MarkerEdgeColor',colorarray{llll},'MarkerFaceColor',colorarray{llll})
    end
    print('-dpng',strcat(saveFolder,filesep,'camerareproject_lframe',num2str(kk),'.png'));
end

