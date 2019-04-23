function [] = spDCM_FirstLevel_PLEsFC(datadir,data,units,N,TR,TE,outdir)
	% This script is designed to replicate the DCM for resting state tutorial in the
	% SPM12 manual.
    % data: processed epi data - epi in split vols
	% I left out the parts where noise time series are extracted and then used to 
	% adjust the SPM. Thus, the data input to this function has to have been cleaned
	% already.

	rng('default')

	estimateSPM = 1;
	createVOIs = 1;


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

        %left dorsal ACC
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_dacc';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-8 28 26];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
        matlabbatch{1}.spm.util.voi.expression = 'i1';
	

		spm_jobman('run',matlabbatch);
		clear matlabbatch


        %left thalamus	
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'left_thalamus';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [-16 -20 14];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/Projects/GenofCog/masks/dcm/aparcANDaseg_2mm_leftthalamus.nii,1'};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch
		
        
        %right thalamus	
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'right_thalamus';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [8 -18 0];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
		matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {'/scratch/kg98/kristina/Projects/GenofCog/masks/dcm/aparcANDaseg_2mm_rightthalamus.nii,1'};
		matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
		matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

		spm_jobman('run',matlabbatch);
		clear matlabbatch

        
		%left dorsal caudate
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
        
        
        %right dorsal caudate
		matlabbatch{1}.spm.util.voi.spmmat = {[outdir,'SPM.mat']};
		matlabbatch{1}.spm.util.voi.adjust = NaN;
		matlabbatch{1}.spm.util.voi.session = 1;
		matlabbatch{1}.spm.util.voi.name = 'right_dorsalcaudate';
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [13 15 9];
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 3.5;
		matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
        matlabbatch{1}.spm.util.voi.expression = 'i1';
	

		spm_jobman('run',matlabbatch);
		clear matlabbatch
        


	end


	cd(outdir)

	ROIs(:,1) = {'left_dacc','left_thalamus','right_thalamus','left_dorsalcaudate','right_dorsalcaudate'};
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

	dcm.n = length(dcm.xY); % number of region
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
    idx_leftdacc = find(strcmp(ROIs,'left_dacc'));
	idx_leftThalamus = find(strcmp(ROIs,'left_thalamus'));
    idx_rightThalamus = find(strcmp(ROIs,'right_thalamus'));
	idx_leftDorsalStr = find(strcmp(ROIs,'left_dorsalcaudate'));
    idx_rightDorsalStr = find(strcmp(ROIs,'right_dorsalcaudate'));
    
	% Build parent model
    parent_sparse = zeros(dcm.n); % initialise
	parent_sparse(eye(dcm.n) == 1) = 1; % self cons
    


    %From left dACC
	parent_sparse(idx_leftDorsalStr,idx_leftdacc) = 1; % left dacc to left dorsal str
	parent_sparse(idx_rightDorsalStr,idx_leftdacc) = 1; % left dacc to right dorsal str
    
    %From dorsal striatum
    parent_sparse(idx_leftThalamus,idx_leftDorsalStr) = 1; %left dorsal striatum to left thalamus
	parent_sparse(idx_rightThalamus,idx_rightDorsalStr) = 1; %right dorsal striatum to right thalamus
    parent_sparse(idx_leftdacc,idx_leftDorsalStr) = 1; %left dorsal striatum to left dacc
	parent_sparse(idx_leftdacc,idx_rightDorsalStr) = 1; %right dorsal striatum to left dacc
    
    %From thalamus
    parent_sparse(idx_leftdacc,idx_leftThalamus) = 1; %left thalamus to left dacc
	parent_sparse(idx_leftdacc,idx_rightThalamus) = 1; %right thalamus to left dacc
    parent_sparse(idx_leftDorsalStr,idx_leftThalamus) = 1; %left thalamus to left dorsal str
	parent_sparse(idx_rightDorsalStr,idx_rightThalamus) = 1; %right thalamus to right dorsal str
%     
    
%     % ------------------------------------------------------------------------------
% 	% Build nested models
% 	% ------------------------------------------------------------------------------
%     AllModels = cell(1,3);
%     AllModels{1} = parent_sparse;
%     
%     % Model 2 - left hemisphere only
%     model2 = zeros(dcm.n); % initialise
%     model2(eye(dcm.n) == 1) = 1; % self cons
%     model2(idx_leftDorsalStr,idx_leftdacc) = 1; % left dacc to left dorsal str
%     model2(idx_leftThalamus,idx_leftDorsalStr) = 1; %left dorsal str to left thalamus
%     model2(idx_leftDorsalStr,idx_leftThalamus) = 1; %left thalamus to left dorsal str
%     model2(idx_leftdacc,idx_leftThalamus) = 1; %left thalamus to left dACC
%     
%     AllModels{2} = model2;
%     
%     % Model 3 - no thalamus
%     model3 = zeros(dcm.n);
%     model3(eye(dcm.n) == 1) = 1;
%     model3(idx_leftdacc, idx_leftDorsalStr) = 1; % left dorsal striatum to left dacc
%     model3(idx_leftdacc, idx_rightDorsalStr) = 1; % right dorsal striatum to left dacc
%     model3(idx_leftDorsalStr, idx_leftdacc) = 1; % left dacc to left dorsal striatum
%     model3(idx_rightDorsalStr, idx_leftdacc) = 1; % left dacc to right dorsal striatum
%     
%     AllModels{3} = model3;
    
%     % ------------------------------------------------------------------------------
% 	% Specify DCM proper
% 	% ------------------------------------------------------------------------------
% 	DCM = cell(1,length(AllModels));
% 	for i = 1:length(AllModels)
% 		DCM{i} = dcm;
% 		DCM{i}.a = AllModels{i};
%     end
%     
% 	
% 	cd(outdir)
% 
%     % ------------------------------------------------------------------------------
% 	% Invert the parent model
% 	% ------------------------------------------------------------------------------
% 	fprintf(1, '\nInverting parent model... \n');
% 	DCM{1} = spm_dcm_fit(DCM{1}); 
%     DCM{1} = DCM{1}{1};
%     save('DCM_spmr7487.mat','DCM')
%     
% 	% ------------------------------------------------------------------------------
% 	% Invert the nested model using Bayesian model reduction
% 	% ------------------------------------------------------------------------------
% 	%fprintf(1, '\nInverting nested models using BMR... \n');
% 	[RCM,BMC,BMA] = spm_dcm_bmr(DCM);
% 
% 	% ------------------------------------------------------------------------------
% 	% Save first-level outputs
% 	% ------------------------------------------------------------------------------
% 	save('DCM_spmr7487.mat','DCM','RCM','BMC','BMA')
	

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
	save('DCM_model2.mat','DCM')

    fprintf(1, '\t\t Done\n');
end
