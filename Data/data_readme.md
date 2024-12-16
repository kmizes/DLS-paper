# Behavioral data files

## Behavioral data (Fig 2, Fig 3, Fig 6, Fig 7)

Files are matlab structures. Each element in the array is a separate session, with fieldnames consisting of high-level experiment data. Vector fields are split by trials, consisting of information on the lever names, lever times, and cue information.


Data format: .mat files contain the following variables.
-	Hit : 1 x trial array of whether trial was rewarded (1) or not (0)
-	WM : 1 x trial array of trial type (-1 OT, 0 CUE, 1 WM)
-	paws : tracked joint x trial cell array of the kinematics of every tracked joint. Each cell array contains a 1 x frame vector of the position of the tracked joint in camera coordinates (pixels). The camera, joint ID, and dimension (x or y) are all saved in variables below
-	vels : tracked joint x trial cell array of the velocity of every tracked joint.
-	joints : 1 x tracked joint cell array of which joint was tracked
-	cams : 1 x tracked joint cell array of which cam the joint was recorded from
-	dimension : 1 x tracked joint cell array of which dimension (x or y, relative to video) of the tracked joint
-	lever : 1 x trial cell array of which levers were pressed
-	protocol : 1 x trial array of the current protocol (7 = flexible, 8 = automatic)
-	sessID : 1 x trial array of which session this trial came from
-	spTimesAll : 1 x trial cell array of the spike counts in each bin (resampled to line up with frames, which is at 40hz). Samples are taken 0.75 seconds before the first, and after the last tap. 
-	tapTimes : 1 x trial cell array of the indices of the lever taps
-	tapTimesVid : joint x trial cell array of the frames in the video of the lever taps.
-	target : 1 x trial cell array of the target sequence
-	trialID : 1 x trial array of the index of the trial in the session
-	trialTimes : 1 x trial cell array of the time in the trial, in seconds. Each cell should contain an array of length bins in the trial.
-	vidFrames : 1 x trial cell array containing the video frames the taps occurred in. 
-	vidNames : 3 x trial cell array containing the video names the kinematic data is from (3 for each camera, Master, Slave1, Slave2).


## Kinematic data (Fig 6, Fig 7)

Kinematic data is preprocessed and aligned to behavioral data. Only trials in CUE and WM that match the automatic sequence are selected. For raw data across all recorded trials and different sequences, please contact corresponding author.


## Ephys Data (Fig 4, Fig 5)

Trial averaged ephys data is uploaded here:

Individual sessions (with spike times, trials, and kinematics) are uploaded here: 
https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/4M2NWF&version=DRAFT
