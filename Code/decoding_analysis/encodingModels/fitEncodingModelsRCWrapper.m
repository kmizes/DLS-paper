function fitEncodingModelsRCWrapper(ratName, sessNum, encPath)
%fitEncodingModelsRCWrapper load encParams file and run encoding models on
%requested session

addpath('./HelperFunctions/');

% create a local cluster object
pc = parcluster('local');
 
% explicitly set the JobStorageLocation to the temp directory that was created in your sbatch script
pc.JobStorageLocation = strcat('/scratch/dhawale/', getenv('SLURM_JOB_ID'));

% load encoding params
encParams = loadSingle([encPath ratName '/encParams']);

% load session
sess = loadSingle(['../Data/SingleSess/' ratName '/' ratName '_sess' num2str(sessNum)]);

% initialize parallel pool
parpool(pc, str2num(getenv('SLURM_CPUS_PER_TASK')));

encFits = fitEncodingModels(sess, encParams);

save([encPath ratName '/encFits_sess' num2str(sessNum)], 'encFits');

end

