% Runs a loop over subjects that need to be analysed
scriptdir = '/scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/';
addpath(genpath(scriptdir))

spmdir = '/scratch/kg98/kristina/spm12_r7487';
addpath(spmdir)

fileID = fopen('/projects/kg98/kristina/GenofCog/scripts/sublists/PLEs_n353.txt');
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};
% compute numsubs
numSubs = length(ParticipantIDs);

TR = 0.754;
N = 616;
TE = 0.021;
units = 'scans';

datadir = '/projects/kg98/kristina/GenofCog/datadir/derivatives/';

Denoise = {'nogsr', 'gsr4'};
WhichDenoise= Denoise{1};

% Nodes = {'SubjectSpecific', 'NonSubjectSpecific'};
% WhichNodes = Nodes{2}


for i = 201:300
    subject = ParticipantIDs{i};
    derivativesdir = [datadir, subject, '/']; % parent folder of subject's derivatives files
    scratchsub = ['/scratch/kg98/kristina/Projects/GenofCog/derivatives/',subject,'/'];
    DCMfirstleveldir = [scratchsub,'DCM_project/firstlevel_dcm/',WhichDenoise,'_CST_fullyconn/'];

% 	switch WhichNodes
% 		case 'SubjectSpecific'
%     		nodedir = [derivativesdir, 'dcm/nodes/from_2ndlvl_peaks_',noise,'/'];
% 		case 'NonSubjectSpecific'
% 			nodedir = '/scratch/kg98/kristina/Projects/GenofCog/ROIs/dcm/sphere_3.5rad/';
% 	end

	switch WhichDenoise
		case 'nogsr'
            preprodir = [derivativesdir,'/prepro.feat/preprocessed/']; % where preprocessed epi is
            epi = [subject,'_filtered_func_data_clean_mni_smooth4mm.nii'];
            data = [subject,'_filtered_func_data_clean_mni_smooth4mm_00000'];
            VolDir=[preprodir,'VolumeSplit/',WhichDenoise,'_smooth4mm/'];
        case 'gsr4'
            preprodir = ['/projects/kg98/kristina/GenofCog/datadir/derivatives/',subject,'/gsr/']; % where preprocessed epi is
            epi = [subject,'_filtered_func_data_clean_gsr4_mni_smooth4mm.nii'];
            data = [subject,'_filtered_func_data_clean_gsr4_mni_smooth4mm_00000'];
            VolDir=[preprodir,'volumesplit/',WhichDenoise,'_smooth4mm/'];
    end
	
    cd (derivativesdir)
    if exist(DCMfirstleveldir) == 0;
        sprintf ('%s: Initialising first-level folder\n', subject)
        mkdir(DCMfirstleveldir)
    else
        fprintf(1,'Cleaning and re-initialising outdir\n')
        rmdir(DCMfirstleveldir,'s')
        mkdir(DCMfirstleveldir)
    end


    %spDCM_FirstLevel_PLEsFC(VolDir,data,units,N,TR,TE,DCMfirstleveldir)
    spDCM_FirstLevel_PsychosisCST(VolDir,data,units,N,TR,TE,DCMfirstleveldir)

end    
