% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'/scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/roi_definition/roi_definition_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
