%-----------------------------------------------------------------------
% Job saved on 28-Mar-2019 15:04:51 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.util.voi.spmmat = {'/scratch/kg98/kristina/Projects/GenofCog/derivatives/sub-034/DCM_project/firstlevel_fc/freesurfer_leftamyg_smooth4mm_nogsr/SPM.mat'};
matlabbatch{1}.spm.util.voi.adjust = NaN;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = 'left_amygdala';
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-24 -6 -18];
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/templates/freesurfer/aparcANDaseg_2mm_leftamyg.nii,1'};
matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
