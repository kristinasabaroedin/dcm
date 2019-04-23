% Runs a loop over subjects that need to be analysed
scriptdir = '/scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/first_level_fc';
addpath(scriptdir)

spmdir = '/projects/kg98/kristina/spm12';
addpath(spmdir)

fileID = fopen('/projects/kg98/kristina/GenofCog/scripts/sublists/PLEs_n353.txt');
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};
% compute numsubs
numSubs = length(ParticipantIDs);

TR = 0.754;
N = 616;

projectdir = '/projects/kg98/kristina/GenofCog/';

projectdir_scratch = '/scratch/kg98/kristina/Projects/GenofCog/';

%datadir = '/projects/kg98/kristina/GenofCog/datadir/derivatives/';


for i = 64:64
    subject = ParticipantIDs{i};
    derivativesdir = [projectdir_scratch, 'derivatives/',subject,'/DCM_project/']; % parent folder of subject's derivatives files
    firstleveldir = [derivativesdir, '/firstlevel_fc/'];

    preprodir = [projectdir,'datadir/derivatives/',subject,'/prepro.feat/preprocessed/']; % where preprocessed epi is
    tsdir = [derivativesdir,'/timeseries/']; % where ts file is
	epi = [subject,'_filtered_func_data_clean_mni_smooth4mm.nii'];
	data = [subject,'_filtered_func_data_clean_mni_smooth4mm_00000'];
	
	volumedir = [preprodir,'VolumeSplit/nogsr_smooth4mm/'];


    cd (derivativesdir)
    if exist(firstleveldir) == 0;
        sprintf ('%s: Initialising first-level folder\n', subject)
        mkdir (firstleveldir)
    else
        display ('First-level folder already exists')
    end


    first_level_subject_dcm(subject, 'freesurfer_leftthalamus_smooth4mm_nogsr', TR, N, firstleveldir, volumedir, tsdir, data)
    first_level_subject_dcm(subject, 'freesurfer_rightthalamus_smooth4mm_nogsr', TR, N, firstleveldir, volumedir, tsdir, data)
   


end    


