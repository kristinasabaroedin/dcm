function [] = spDCM_FirstLevel_PsychosisCST(datadir,data,units,N,TR,TE,outdir)
	% This script is designed to replicate the DCM for resting state tutorial in the
	% SPM12 manual.
    % data: processed epi data - epi in split vols
	% I left out the parts where noise time series are extracted and then used to 
	% adjust the SPM. Thus, the data input to this function has to have been cleaned
	% already.

	rng('default')

	estimateSPM = 1;
	createVOIs = 1;

%  	spmdir = '/scratch/kg98/kristina/spm12_r7487';
%  	addpath(spmdir)

	% ------------------------------------------------------------------------------
	% Initiate SPM
	% ------------------------------------------------------------------------------
	spm('defaults','fmri');
	spm_jobman('initcfg')
	spm_get_defaults('mask.thresh',-Inf)

	if estimateSPM == 1

		%-----------------------------------------------------------------------
		% Estimate SPM
		%-----------------------------------------------------------------------
		matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir};
		matlabbatch{1}.spm.stats.fmri_spec.timing.units = units;
		matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
		matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
		matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1; %double check what this value shld be
		%%

		for i = 1:N
		    if i < 10
		        ZeroPad = '0000';
		    elseif i >= 10 & i < 100
		        ZeroPad = '000';
		    elseif i >= 100
		        ZeroPad = '00';
		    end
		    matlabbatch{1}.spm.stats.fmri_spec.sess.scans{i,1} = [datadir,data(1:end-5),ZeroPad,num2str(i),'.nii'];
		end

		%%
		matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
		matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
		matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
		matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''};
		matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = -Inf;
		matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
		matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
		matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
		matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
		matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
		matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
		matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';%none
		matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
		matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
		matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

		spm_jobman('run',matlabbatch);
		clear matlabbatch

	end

	if createVOIs == 1
		%-----------------------------------------------------------------------
		% Generate VOIs
		%-----------------------------------------------------------------------

		%DLPFC	
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_DLPFC';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-45 23 37];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 6;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/Projects/GenofCog/masks/dcm/HarvardOxford_left_DLPFC.nii,1'};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch

		%VMPFC	
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_VMPFC';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-6 45 -9];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 6;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/Projects/GenofCog/masks/dcm/raal_left_vmpfc.nii,1'};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch


		%thalamus	
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_thalamus';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-6 -12 14];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/Projects/GenofCog/masks/dcm/aparcANDaseg_2mm_leftthalamus.nii,1'};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch

		%amygdala	
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_amygdala';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-24 -6 -18];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/Projects/GenofCog/masks/dcm/aparcANDaseg_2mm_leftamyg.nii,1'};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch


		%ant hippocampus	
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_anthipp';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-26 -16 -18];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/templates/freesurfer/aparcANDaseg_2mm_leftanthipp.nii,1'};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch


		%midbrain	
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_midbrain';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-12 -20 -10];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/Projects/GenofCog/masks/dcm/left_midbrain_ero.nii,1'};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch


		%dorsal caudate
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_dorsalcaudate';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-13 15 9];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
        matlabbatch{1}.spm.util.voi.expression = 'i1';
	

		spm_jobman('run',matlabbatch);
		clear matlabbatch

		%nucleus accumbens
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_nucleusacc';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-9 9 -8];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
        matlabbatch{1}.spm.util.voi.expression = 'i1';
	

		spm_jobman('run',matlabbatch);
		clear matlabbatch

	end


	cd(outdir)

	ROIs(:,1) = {'left_DLPFC','left_VMPFC','left_thalamus','left_amygdala','left_anthipp','left_midbrain','left_dorsalcaudate','left_nucleusacc'};
	numROIs = length(ROIs);

	% ------------------------------------------------------------------------------
	% Specify DCM template (calling this lower-case dcm)
	% There wasn't any obvious way to do this with the SPM command line tools that 
	% didn't require interaction with the SPM GUI (e.g., spm_dcm_specify) and the SPM12
	% manual only gives instructions for specifying via the GUI not via the batch editor.
	% 
	% So, I made a template DCM.mat file using the GUI and then wrote the code below
	% to replicate that template. This way the code can run without being dependent 
	% on a .mat file
	% 
	% Note that the notation of 'dcm' in lower case refers to a template upon which
	% different models are specified. the 'DCM' represents the structure I actually
	% invert.
	% ------------------------------------------------------------------------------
	fprintf(1, 'Specifying DCM \n\n');

	load([outdir,'SPM.mat'])

	clear dcm

	% Load regions of interest
	for i = 1:numROIs
		load(fullfile([outdir,'VOI_',ROIs{i},'_1.mat']),'xY');
		dcm.xY(i) = xY;
	end
	clear xY

	dcm.n = length(dcm.xY); % number of regions. I could have pulled this from numROIs
	dcm.v = length(dcm.xY(1).u); % number of time points. Could also pull this from function input N

	% Experimental inputs (which there are none)
	dcm.U.u = zeros(dcm.v,1);
	dcm.U.name = {'null'};

	% Time series
	dcm.Y.dt  = SPM.xY.RT;
	dcm.Y.X0  = dcm.xY(1).X0;
	for i = 1:dcm.n
	    dcm.Y.y(:,i) = dcm.xY(i).u;
	    dcm.Y.name{i} = dcm.xY(i).name;
	end

	dcm.Y.Q = spm_Ce(ones(1,dcm.n)*dcm.v);

	% dcm parameters and options
	dcm.delays = repmat(SPM.xY.RT/2,dcm.n,1);
	dcm.TE = TE;

	dcm.options.nonlinear = 0;
	dcm.options.two_state = 0;
	dcm.options.stochastic = 0;
	dcm.options.centre = 1;
	dcm.options.induced = 1;
	dcm.options.maxnodes = numROIs; 
    
	dcm.b = zeros(dcm.n,dcm.n);
	dcm.c = zeros(dcm.n,1);
	dcm.d = double.empty(dcm.n,dcm.n,0);

	% ------------------------------------------------------------------------------
	% Define model space
	% ------------------------------------------------------------------------------	
	%find indices of ROI names
    idx_DLPFC = find(strcmp(ROIs,'left_DLPFC'));
	idx_VMPFC = find(strcmp(ROIs,'left_VMPFC'));
	idx_Amygdala = find(strcmp(ROIs,'left_amygdala'));
	idx_AntHipp = find(strcmp(ROIs,'left_anthipp'));
	idx_Thalamus = find(strcmp(ROIs,'left_thalamus'));
	idx_Midbrain = find(strcmp(ROIs,'left_midbrain'));
	idx_DorsalStr = find(strcmp(ROIs,'left_dorsalcaudate'));
	idx_VentralStr = find(strcmp(ROIs,'left_nucleusacc'));
    
	% Build parent model
    parent_sparse = zeros(dcm.n); % initialise
	parent_sparse(eye(dcm.n) == 1) = 1; % self cons
    
    %From VMPFC
	parent_sparse(idx_AntHipp,idx_VMPFC) = 1; % VMPFC to Hipp
	parent_sparse(idx_Amygdala,idx_VMPFC) = 1; % VMPFC to Amyg
    parent_sparse(idx_Thalamus,idx_VMPFC) = 1; % VMPFC to Thal
    parent_sparse(idx_VentralStr,idx_VMPFC) = 1; % VMPFC to VSi
	parent_sparse(idx_Midbrain,idx_VMPFC) = 1; % VMPFC to midbrain
	parent_sparse(idx_DLPFC,idx_VMPFC) = 1; % VMPFC to DLPFC
    %From Hipp
	parent_sparse(idx_VentralStr,idx_AntHipp) = 1; % Hipp to VSi
	parent_sparse(idx_Amygdala,idx_AntHipp) = 1; % Hipp to Amyg
	parent_sparse(idx_Thalamus,idx_AntHipp) = 1; % Hipp to Thal
	parent_sparse(idx_VMPFC,idx_AntHipp) = 1; % Hipp to VMPFC
	%From Amyg
  	parent_sparse(idx_VentralStr,idx_Amygdala) = 1; % Amyg to VSi
	parent_sparse(idx_AntHipp,idx_Amygdala) = 1; % Amyg to Hipp
	parent_sparse(idx_Thalamus,idx_Amygdala) = 1; % Amyg to Thal
	parent_sparse(idx_VMPFC,idx_Amygdala) = 1; % Amyg to VMPFC
	parent_sparse(idx_Midbrain,idx_Amygdala) = 1; % Amyg to Midbrain
    %From Thal
    parent_sparse(idx_VMPFC,idx_Thalamus) = 1; % Thal to VMPFC  
    parent_sparse(idx_VentralStr,idx_Thalamus) = 1; % Thal to VSi  
    parent_sparse(idx_DorsalStr,idx_Thalamus) = 1; % Thal to dorsal caudate 
    parent_sparse(idx_DLPFC,idx_Thalamus) = 1; % Thal to DLPFC 
	parent_sparse(idx_AntHipp,idx_Thalamus) = 1;  % Thal to Hipp
	parent_sparse(idx_Amygdala,idx_Thalamus) = 1; % Thal to Amyg
    %From VSi
    parent_sparse(idx_Thalamus,idx_VentralStr) = 1; % VSi to thal
    parent_sparse(idx_Midbrain,idx_VentralStr) = 1; % VSi to midbrain
    %From DC
    parent_sparse(idx_Thalamus,idx_DorsalStr) = 1; % dorsal caudate to thal
    parent_sparse(idx_Midbrain,idx_DorsalStr) = 1; % dorsal caudate to midbrain
    %From Midbrain
    parent_sparse(idx_VentralStr,idx_Midbrain) = 1; % midbrain to VSi
    parent_sparse(idx_DorsalStr,idx_Midbrain) = 1; % midbrain to dorsal caudate
	parent_sparse(idx_DLPFC,idx_Midbrain) = 1; % midbrain to DLPFC
	parent_sparse(idx_VMPFC,idx_Midbrain) = 1; % midbrain to VMPFC
	parent_sparse(idx_Amygdala,idx_Midbrain) = 1; % midbrain to Amyg
	parent_sparse(idx_AntHipp,idx_Midbrain) = 1; % midbrain to Hipp
    %From DLPFC
    parent_sparse(idx_Thalamus,idx_DLPFC) = 1; % DLPFC to Thal
    parent_sparse(idx_DorsalStr,idx_DLPFC) = 1; % DLPFC to Dorsal Caudate
	parent_sparse(idx_VMPFC,idx_DLPFC) = 1; % DLPFC to VMPFC
	parent_sparse(idx_Midbrain,idx_DLPFC) = 1; % DLPFC to midbrain
    
	

	% ------------------------------------------------------------------------------
	% Specify DCM proper
	% ------------------------------------------------------------------------------
	
	DCM = dcm;
	DCM.a = parent_sparse;
	

	cd(outdir)
	% ------------------------------------------------------------------------------
	% Invert the parent model
	% ------------------------------------------------------------------------------
	fprintf(1, '\nInverting parent model... \n');
	DCM = spm_dcm_fit(DCM); 
    [DCM] = spm_dcm_fmri_check(DCM);
    DCM = DCM{1};
	% ------------------------------------------------------------------------------
	% Save first-level outputs
	% ------------------------------------------------------------------------------
	save('DCM.mat','DCM')
	

    fprintf(1, '\t\t Done\n');
end
