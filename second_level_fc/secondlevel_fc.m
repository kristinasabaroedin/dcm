% This script runs a 1x2 multifactorial design for second level analysis (main effect of a midbrain seed for both hemispheres)


spmdir = '/scratch/kg98/kristina/spm12_r7487';
addpath(spmdir)


seed = 'freesurfer_rightthalamus';
datadir = '/scratch/kg98/kristina/Projects/GenofCog/derivatives/';
firstlevel_dir = ['/DCM_project/firstlevel_fc/',seed,'_smooth4mm_nogsr/'];
outdir = [datadir,'secondlevel/dcm_fc/',seed,'/'];

% make outdir
if exist(outdir) == 0
    fprintf(1,'Initialising outdir\n')
    mkdir(outdir)
elseif exist(outdir) == 7
    fprintf(1,'Output directory already exists\n')
end

% ------------------------------------------------------------------------------
% need to get a list of subs
% ------------------------------------------------------------------------------
sublist = '/projects/kg98/kristina/GenofCog/scripts/sublists/PLEs_n353.txt'
fileID = fopen(sublist);
ParticipantIDs = textscan(fileID,'%s');
subs = ParticipantIDs{1};
% number of subjects
numSubs = length(subs)


% ------------------------------------------------------------------------------
% Count number of subjects (subs variable in matlab)
% ------------------------------------------------------------------------------
% number of subjects
numSubs = length(subs);

spm

%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.factorial_design.dir = {outdir};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'seed';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = 1;

                                                                   
for i = 1:numSubs
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans(i,1) = {[datadir,subs{i},firstlevel_dir,'con_0001.nii,1']}';
end


matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/projects/kg98/kristina/templates/MNI152_T1_2mm_brain_mask.nii,1'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


% ------------------------------------------------------------------------------
% Estimate
% ------------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[outdir,'SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% ------------------------------------------------------------------------------
% Contrasts
% ------------------------------------------------------------------------------
matlabbatch{3}.spm.stats.con.spmmat = {[outdir,'SPM.mat']};
matlabbatch{3}.spm.stats.con.delete = 1;

matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = [seed,'+'];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = [seed,'-'];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';


% ------------------------------------------------------------------------------
% Run
% ------------------------------------------------------------------------------
spm_jobman('run',matlabbatch);
clear matlabbatch

cd(outdir)
