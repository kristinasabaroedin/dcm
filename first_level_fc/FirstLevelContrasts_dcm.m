function [] = FirstLevelContrasts_dcm(seed,matfile)
%
%
% ------
% INPUTS
% ------
%
% matfile 	- path and file name of first level estimation output .mat file (typically SPM.mat)
%
% -------
% OUTPUTS
% -------
%
% SPM.mat file containing contrats for ROIs for 1st level analysis. located in spmdir.
%
% =========================================================================

spm('defaults','fmri');
spm_jobman('initcfg')
spm_get_defaults('mask.thresh',-Inf)

% run contrasts

matlabbatch{1}.spm.stats.con.spmmat = {matfile};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = [seed,'+'];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = [seed,'-'];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-1 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;


spm_jobman('run',matlabbatch);

clear matlabbatch
