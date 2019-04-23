function [] = spDCM_FirstLevel_GoC(datadir,data,units,N,TR,TE,nodedir,outdir)
	% This script is designed to replicate the DCM for resting state tutorial in the
	% SPM12 manual.
    % data: processed epi data - epi in split vols?
	% I left out the parts where noise time series are extracted and then used to 
	% adjust the SPM. Thus, the data input to this function has to have been cleaned
	% already.

	% Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,
	rng('default')


	% ------------------------------------------------------------------------------
	% Initiate SPM
	% ------------------------------------------------------------------------------
	spm('defaults','fmri');
	spm_jobman('initcfg')
	spm_get_defaults('mask.thresh',-Inf)

	%-----------------------------------------------------------------------
	% Estimate GLM
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

	%-----------------------------------------------------------------------
	% Generate VOIs
	%-----------------------------------------------------------------------

 	ROIs = dir([nodedir,'*.nii']);
 	numROIs = size(ROIs,1);
	origorder = [1:numROIs];

	Hemi = {'left','right'};
	WhichHemi = Hemi{1};

	WhichNode = 'NonSubjectSpecific';
	switch WhichNode
		case 'NonSubjectSpecific'
		% find ROI indices
			if strcmp(WhichHemi, 'left')
			% change to if left hemi
				idx_DorsalCaudate = find(strcmp({ROIs.name},'DC_sphere_left.nii'));
				idx_VentStri = find(strcmp({ROIs.name},'VSi_sphere_left.nii'));
				idx_Thal = find(strcmp({ROIs.name},'thal_dandash_sphere_left.nii'));
				idx_Hipp = find(strcmp({ROIs.name},'anthipp_freesurfer_sphere_left.nii'));
			    idx_Amyg = find(strcmp({ROIs.name},'amyg_freesurfer_sphere_left.nii'));
			    idx_Midbrain = find(strcmp({ROIs.name},'midbrain_murty_sphere_left.nii'));
				idx_VMPFC = find(strcmp({ROIs.name},'VMPFC_zhang_sphere_left.nii'));
			    idx_DLPFC = find(strcmp({ROIs.name},'DLPFC_fornito_sphere_left.nii'));
			elseif strcmp(WhichHemi, 'right')
				idx_DorsalCaudate = find(strcmp({ROIs.name},'DC_sphere_right.nii'));
				idx_VentStri = find(strcmp({ROIs.name},'VSi_sphere_right.nii'));
				idx_Thal = find(strcmp({ROIs.name},'thal_right.nii'));
				idx_Hipp = find(strcmp({ROIs.name},'hipp_right.nii'));
			    idx_Amyg = find(strcmp({ROIs.name},'amyg_right.nii'));
			    idx_Midbrain = find(strcmp({ROIs.name},'midbrain_right.nii'));
				idx_VMPFC = find(strcmp({ROIs.name},'VMPFC_right.nii'));
			    idx_DLPFC = find(strcmp({ROIs.name},'DLPFC_right.nii'));
			end
		case 'SubjectSpecific'
		% find ROI indices
			if strcmp(WhicHemi, 'left')
			% change to if left hemi
				idx_DorsalCaudate = find(strcmp({ROIs.name},'DC_left.nii'));
				idx_VentStri = find(strcmp({ROIs.name},'VSi_left.nii'));
				idx_Thal = find(strcmp({ROIs.name},'thal_left.nii'));
				idx_Hipp = find(strcmp({ROIs.name},'hipp_left.nii'));
			    idx_Amyg = find(strcmp({ROIs.name},'amyg_left.nii'));
			    idx_Midbrain = find(strcmp({ROIs.name},'midbrain_left.nii'));
				idx_VMPFC = find(strcmp({ROIs.name},'VMPFC_left.nii'));
			    idx_DLPFC = find(strcmp({ROIs.name},'DLPFC_left.nii'));
			elseif strcmp(WhicHemi, 'left')
				idx_DorsalCaudate = find(strcmp({ROIs.name},'DC_sphere_right.nii'));
				idx_VentStri = find(strcmp({ROIs.name},'VSi_sphere_right.nii'));
				idx_Thal = find(strcmp({ROIs.name},'thal_right.nii'));
				idx_Hipp = find(strcmp({ROIs.name},'hipp_right.nii'));
			    idx_Amyg = find(strcmp({ROIs.name},'amyg_right.nii'));
			    idx_Midbrain = find(strcmp({ROIs.name},'midbrain_right.nii'));
				idx_VMPFC = find(strcmp({ROIs.name},'VMPFC_right.nii'));
			    idx_DLPFC = find(strcmp({ROIs.name},'DLPFC_right.nii'));
			end

	end
    
	perm = [idx_DorsalCaudate idx_VentStri idx_Thal idx_Hipp idx_Amyg idx_Midbrain idx_VMPFC idx_DLPFC];

	% permute
	ROIs = ROIs(perm);

	% print ROI names to command line for visual confirmation
	{ROIs.name}'

	for i = 1:numROIs

		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = ROIs(i).name(1:end-4);
		matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {[nodedir,ROIs(i).name,',1']};
		matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {[outdir,'mask.nii,1']};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch
	end

	cd(nodedir); cd('../')
	if exist('temp') == 7
		rmdir('temp','s')
	end
	cd(outdir)

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
		load(fullfile(outdir,['VOI_',ROIs(i).name(1:end-4),'_1.mat']),'xY');
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
	dcm.options.maxnodes = numROIs; % bastard!
    
	dcm.b = zeros(dcm.n,dcm.n);
	dcm.c = zeros(dcm.n,1);
	dcm.d = double.empty(dcm.n,dcm.n,0);

	% ------------------------------------------------------------------------------
	% Define model space
	% ------------------------------------------------------------------------------	

    if ismember(ROIs(1).name(end-4),'left')
        idx_DorsalCaudate = find(strcmp({ROIs.name},'DC_sphere_left.nii'));
        idx_VentStri = find(strcmp({ROIs.name},'VSi_sphere_left.nii'));
        idx_Thal = find(strcmp({ROIs.name},'thal_dandash_sphere_left.nii'));
        idx_Hipp = find(strcmp({ROIs.name},'anthipp_freesurfer_sphere_left.nii'));
        idx_Amyg = find(strcmp({ROIs.name},'amyg_freesurfer_sphere_left.nii'));
        idx_Midbrain = find(strcmp({ROIs.name},'midbrain_murty_sphere_left.nii'));
        idx_VMPFC = find(strcmp({ROIs.name},'VMPFC_zhang_sphere_left.nii'));
        idx_DLPFC = find(strcmp({ROIs.name},'DLPFC_fornito_sphere_left.nii'));
	elseif ismember(ROIs(1).name(end-4),'right')
		idx_DorsalCaudate = find(strcmp({ROIs.name},'DC_sphere_right.nii'));
		idx_VentStri = find(strcmp({ROIs.name},'VSi_sphere_right.nii'));
		idx_Thal = find(strcmp({ROIs.name},'thal_right.nii'));
		idx_Hipp = find(strcmp({ROIs.name},'hipp_right.nii'));
        idx_Amyg = find(strcmp({ROIs.name},'amyg_right.nii'));
        idx_Midbrain = find(strcmp({ROIs.name},'midbrain_right.nii'));
		idx_VMPFC = find(strcmp({ROIs.name},'VMPFC_right.nii'));
        idx_DLPFC = find(strcmp({ROIs.name},'DLPFC_right.nii'));
    end
    
	% Build parent model
    par_sparse = zeros(dcm.n); % initialise
	par_sparse(eye(dcm.n) == 1) = 1; % self cons
    
    %From VMPFC
	par_sparse(idx_Hipp,idx_VMPFC) = 1; % VMPFC to Hipp
	par_sparse(idx_Amyg,idx_VMPFC) = 1; % VMPFC to Amyg
    par_sparse(idx_Thal,idx_VMPFC) = 1; % VMPFC to Thal
    par_sparse(idx_VentStri,idx_VMPFC) = 1; % VMPFC to VSi
    %From Hipp
	par_sparse(idx_VentStri,idx_Hipp) = 1; % Hipp to VSi
	%From Amyg
  	par_sparse(idx_VentStri,idx_Amyg) = 1; % Amyg to VSi
    %From Thal
    par_sparse(idx_VMPFC,idx_Thal) = 1; % thal to VMPFC  
    par_sparse(idx_VentStri,idx_Thal) = 1; % thal to VSi  
    par_sparse(idx_DorsalCaudate,idx_Thal) = 1; % thal to dorsal caudate 
    par_sparse(idx_DLPFC,idx_Thal) = 1; % thal to DLPFC 
    %From VSi
    par_sparse(idx_Thal,idx_VentStri) = 1; % VSi to thal
    par_sparse(idx_Midbrain,idx_VentStri) = 1; % VSi to midbrain
    %From DC
    par_sparse(idx_Thal,idx_DorsalCaudate) = 1; % dorsal caudateto thal
    par_sparse(idx_Midbrain,idx_DorsalCaudate) = 1; % dorsal caudate to midbrain
    %From Midbrain
    par_sparse(idx_VentStri,idx_Midbrain) = 1; % midbrain to VSi
    par_sparse(idx_DorsalCaudate,idx_Midbrain) = 1; % midbrain to dorsal caudate
    %From DLPFC
    par_sparse(idx_Thal,idx_DLPFC) = 1; % DLPFC to Thal
    par_sparse(idx_DorsalCaudate,idx_DLPFC) = 1; % DLPFC to Dorsal Caudate
    
	% Build nested models
	WhichNest = 'Nested';
	switch WhichNest
		case 'Nested'
			par_sparse(eye(dcm.n) == 1) = 2; % self cons (set to 2 to make them immutable. they are set back to 1 below)
			% flatten
			x = dcm.n * (dcm.n - 1) + dcm.n;
			par_sparse = reshape(par_sparse,x,1);
			idx = find(par_sparse == 1); % index of connections to turn off
			% models
			models = [];
			% systematically switch off one connection at a time (not including self cons)
			for i = 1:length(idx)
				models(:,i) = par_sparse;
				models(idx(i),i) = 0;
			end
			models = [par_sparse,models];

			% reorganise into NxN matrices, remove 2 values from self cons and store in cell
			TheModels = cell(1,size(models,2));
			for i = 1:size(models,2)
				TheModels{i} = reshape(models(:,i),dcm.n,dcm.n) - eye(dcm.n);
			end
		case 'Manual'
			TheModels = cell(1,4);
			TheModels{1} = par_sparse;

			% model 2: actual model
			model2 = zeros(dcm.n); % initialise
            model2(eye(dcm.n) == 1) = 1; % self cons
            model2(idx_Hipp,idx_VMPFC) = 1; % VMPFC to Hipp
            model2(idx_Amyg,idx_VMPFC) = 1; % VMPFC to Amyg
            model2(idx_VentStri,idx_VMPFC) = 1; % VMPFC to VSi
            model2(idx_VentStri,idx_Hipp) = 1; % Hipp to VSi
            model2(idx_VentStri,idx_Amyg) = 1; % Amyg to VSi
            model2(idx_Midbrain,idx_VentStri) = 1; % VSi to midbrain
            model2(idx_VentStri,idx_Midbrain) = 1; % midbrain to VSi
            model2(idx_DorsalCaudate,idx_Midbrain) = 1; % midbrain to dorsal caudate
            model2(idx_Midbrain,idx_DorsalCaudate) = 1; % dorsal caudate to midbrain
            model2(idx_Thal,idx_DorsalCaudate) = 1; % dorsal caudate to thal
            model2(idx_Thal,idx_VentStri) = 1; % VSi to thal
            model2(idx_DLPFC,idx_Thal) = 1; % thal to DLPFC  
            model2(idx_VMPFC,idx_Thal) = 1; % thal to VMPFC  
            model2(idx_VentStri,idx_Thal) = 1; % thal to VSi  
            model2(idx_DorsalCaudate,idx_Thal) = 1; % thal to dorsal caudate 
            model2(idx_DorsalCaudate,idx_DLPFC) = 1; % DLPFC to Dorsal Caudate
            model2(idx_DLPFC,idx_VMPFC) = 1; % VMPFC to DLPFC
			
            TheModels{2} = model2;

			% model 3: no amygdala or hippocampus
			model3 = zeros(dcm.n); % initialise
            model3(eye(dcm.n) == 1) = 1; % self cons
            model3(idx_VentStri,idx_VMPFC) = 1; % VMPFC to VSi
            model3(idx_Midbrain,idx_VentStri) = 1; % VSi to midbrain
            model3(idx_VentStri,idx_Midbrain) = 1; % midbrain to VSi
            model3(idx_DorsalCaudate,idx_Midbrain) = 1; % midbrain to dorsal caudate
            model3(idx_Midbrain,idx_DorsalCaudate) = 1; % dorsal caudate to midbrain
            model3(idx_Thal,idx_DorsalCaudate) = 1; % dorsal caudate to thal
            model3(idx_Thal,idx_VentStri) = 1; % VSi to thal
            model3(idx_DLPFC,idx_Thal) = 1; % thal to DLPFC  
            model3(idx_VMPFC,idx_Thal) = 1; % thal to VMPFC  
            model3(idx_VentStri,idx_Thal) = 1; % thal to VSi  
            model3(idx_DorsalCaudate,idx_Thal) = 1; % thal to dorsal caudate 
            model3(idx_DorsalCaudate,idx_DLPFC) = 1; % DLPFC to Dorsal Caudate
            model3(idx_DLPFC,idx_VMPFC) = 1; % VMPFC to DLPFC
            
			TheModels{3} = model3;

			% model 4: no top-down conn from thalamus to striatum
			model4 = zeros(dcm.n); % initialise
            model4(eye(dcm.n) == 1) = 1; % self cons
            model4(idx_VentStri,idx_VMPFC) = 1; % VMPFC to VSi
            model4(idx_Midbrain,idx_VentStri) = 1; % VSi to midbrain
            model4(idx_VentStri,idx_Midbrain) = 1; % midbrain to VSi
            model4(idx_DorsalCaudate,idx_Midbrain) = 1; % midbrain to dorsal caudate
            model4(idx_Midbrain,idx_DorsalCaudate) = 1; % dorsal caudate to midbrain
            model4(idx_Thal,idx_DorsalCaudate) = 1; % dorsal caudate to thal
            model4(idx_Thal,idx_VentStri) = 1; % VSi to thal
            model4(idx_DLPFC,idx_Thal) = 1; % thal to DLPFC  
            model4(idx_VMPFC,idx_Thal) = 1; % thal to VMPFC  
            model4(idx_DorsalCaudate,idx_DLPFC) = 1; % DLPFC to Dorsal Caudate
            model4(idx_DLPFC,idx_VMPFC) = 1; % VMPFC to DLPFC
            
			TheModels{4} = model4;
	end

	% ------------------------------------------------------------------------------
	% Specify DCM proper
	% ------------------------------------------------------------------------------
	DCM = cell(1,length(TheModels));
	for i = 1:length(TheModels)
		DCM{i} = dcm;
		DCM{i}.a = TheModels{i};
	end

	cd(outdir)
	% ------------------------------------------------------------------------------
	% Invert the parent model
	% ------------------------------------------------------------------------------
	fprintf(1, '\nInverting parent model... \n');
	DCM{1} = spm_dcm_fit(DCM{1}); 
	%DCM{1} = spm_dcm_fmri_csd(DCM{1});
    DCM{1} = DCM{1}{1};

	% ------------------------------------------------------------------------------
	% Invert the nested model using Bayesian model reduction
	% ------------------------------------------------------------------------------
	%fprintf(1, '\nInverting nested models using BMR... \n');
	[RCM,BMC,BMA] = spm_dcm_bmr(DCM);

	% ------------------------------------------------------------------------------
	% Save first-level outputs
	% ------------------------------------------------------------------------------
	save('DCM.mat','DCM','RCM','BMC','BMA')
	%save('DCM.mat','DCM')

    fprintf(1, '\t\t Done\n');
end
