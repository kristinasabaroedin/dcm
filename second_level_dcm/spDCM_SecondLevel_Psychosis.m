% Script to run second level DCM
clc

createTable = 1;

% ------------------------------------------------------------------------------
% Parent dir
% ------------------------------------------------------------------------------
parentdir = '/projects/kg98/kristina/';
parentdir_scratch = '/scratch/kg98/kristina/';

% ------------------------------------------------------------------------------
% Add paths - edit this section
% ------------------------------------------------------------------------------
% where the scripts are
scriptsdir = [parentdir_scratch,'Projects/GenofCog/scripts/dcm_project/'];
addpath(genpath(scriptsdir))
clear scriptsdir

% where spm12 is
spmdir = '/projects/kg98/kristina/spm12';
addpath(spmdir)



% ------------------------------------------------------------------------------
% Set options
% ------------------------------------------------------------------------------
Projects = {'GenofCog', 'UCLA', 'STAGES'};
WhichProject = Projects{1}


Hemi = {'Left','Right'}
WhichHemi = Hemi{1};


Denoise = {'gsr4', 'nogsr', 'dbscan', '2p+AROMAaggr'};
WhichDenoise = Denoise{2}

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
    case 'GenofCog'
        projdir = [parentdir_scratch,'Projects/GenofCog/derivatives/'];
        datafile = [projdir,'secondlevel/DCM/GoC_metadata_n352.mat'];
        TR = 0.754;

    case 'UCLA'
        projdir = [parentdir_scratch,'projects/UCLA/'];
        datadir = [projdir,'datadir/derivatives/'];
        preprostr = '/func/prepro/';

        datafile = [projdir,'SecondLevel/SPM/Factorial/',WhichNoise,'TriStri/metadata.mat'];
        TR = 2;

    case 'STAGES'
        projdir = [parentdir_scratch,'projects/STAGES/'];
        datadir = [projdir,'datadir/derivatives/'];
        preprostr = '/func/prepro/';

        datafile = [projdir,'SecondLevel/SPM/Factorial/',WhichNoise,'TriStri/metadata.mat'];
        TR = 2;
end

% ------------------------------------------------------------------------------
% Setup output directory
% ------------------------------------------------------------------------------
% new
outdir = [parentdir_scratch,'/Projects/',WhichProject,'/derivatives/secondlevel/DCM/',WhichDenoise,'/18_conns/'];

if exist(outdir) == 0
    fprintf(1,'Initialising outdir\n')
    mkdir(outdir)
elseif exist(outdir) == 7
    fprintf(1,'Cleaning and re-initialising outdir\n')
    rmdir(outdir,'s')
    mkdir(outdir)
end
cd(outdir)

% ------------------------------------------------------------------------------
% Non-imaging metadata
% ------------------------------------------------------------------------------
load(datafile)
numSubs = size(metadata.SubjectID,1);
numGroups = numel(unique(metadata.Diagnosis));

% ------------------------------------------------------------------------------
% Create demographic table
% ------------------------------------------------------------------------------
if createTable == 1

	switch WhichProject
    	case 'GenofCog'
            Age = zeros(numGroups,1); AgeSquared = zeros(numGroups,1); Sex = zeros(numGroups,1); 
            IQ = zeros(numGroups,1); MeanFD = zeros(numGroups,1); 
            PC1 = zeros(numGroups,1); PC2 = zeros(numGroups,1); 
            for i = 1:numGroups
                logi = metadata.Diagnosis == i;
                Age(i) = mean(metadata(logi,:).Age); Age_SD(i) = std(metadata(logi,:).Age); AgeSquared(i) = mean(metadata(logi,:).AgeSquared);
                Sex(i) = sum(metadata(logi,:).Sex == 1); Sex_percent(i) = Sex(i) / sum(logi);
                IQ(i) = mean(metadata(logi,:).IQ); IQ_SD(i) = std(metadata(logi,:).IQ);
                MeanFD(i) = mean(metadata(logi,:).MeanFD); MeanFD_SD(i) = std(metadata(logi,:).MeanFD);
            end

    	demo_table = table(Age,Age_SD,AgeSquared,Sex,Sex_percent,IQ,IQ_SD,MeanFD,MeanFD_SD)
    	clear Age* Sex* IQ* PC* MeanFD* 

	end
end
% ------------------------------------------------------------------------------
% Create a GDCMs
% i.e., a group of reduced DCMs generate at first level using BMR
% ------------------------------------------------------------------------------
% load subject 1 to find number of models
DCMfirstleveldir = [parentdir_scratch,'Projects/GenofCog/derivatives/',metadata.SubjectID{1},'/DCM_project/firstlevel_dcm/',WhichDenoise,'_18conns/'];
load([DCMfirstleveldir,'DCM.mat'],'RCM')
numModels = length(RCM); 
clear RCM 

% Load DCMs for each subject into a cell array
GRCM = cell(numSubs,numModels);
GBMC.F = zeros(numSubs,numModels);
GBMC.P = zeros(numSubs,numModels);
%modelChosen=1;
fprintf(1, 'Loading DCMs...\n');
for i = 1:numSubs
    % new
    load([parentdir_scratch,'Projects/GenofCog/derivatives/',metadata.SubjectID{i},'/DCM_project/firstlevel_dcm/',WhichDenoise,'_18conns/DCM.mat'],'RCM','BMC')

    GRCM(i,:) = RCM;
    % GRCM(i,:) = RCM(1,1:3);
    GBMC.F(i,:) = BMC.F;
    GBMC.P(i,:) = BMC.P;
    clear RCM BMC
    
end


fprintf(1, 'Done\n');

% ------------------------------------------------------------------------------
% 1) primary results
% BMA using PEB routine, but modeling group effects at the second level
% ------------------------------------------------------------------------------
    vars = {'PC2','PC1'};
    nuisance_vars = {'Age','AgeSquared','Sex','IQ','MeanFD'};
    vars = [vars nuisance_vars];

    numVars = length(vars);
    
    Variable = {'detrend','zscore','none'};
    
    WhichVariable = Variable{2}
    
    switch WhichVariable
        case Variable{1}
            covs = detrend(metadata{:,vars});
        case Variable{2}
            covs = zscore(metadata{:,vars});
        case Variable{3}
            covs = metadata{:,vars};
    end


    % generate diagnosis dummy variables for f-test
    addContrasts = 0;
    if addContrasts == 1
        % exhaustive t-tests (case-control & case-case)
        cons = [kron([1 1 -1 -1 0 0],ones(sum(metadata.Diagnosis == 1),1));...
                kron([-1 0 1 0 1 -1],ones(sum(metadata.Diagnosis == 2),1));...
                kron([0 -1 0 1 -1 1],ones(sum(metadata.Diagnosis == 3),1))];

        str = [];
        for j = 1:size(cons,2)
            str{j} = ['C ',num2str(j)];
        end

        covs = [cons covs];
        vars = [str vars];
        numVars = length(vars);
    end

    vars

    % Specify PEB model settings (see batch editor for help on each setting)
    M = struct();
    M.alpha = 1;
    M.beta  = 16;
    M.hE    = 0;
    M.hC    = 1/16;
    M.Q     = 'single';

    % Specify design matrix for N subjects. It should start with a constant column
    %M.X = [ones(numSubs,1) covs_demean];
    M.X = [ones(numSubs,1) covs];
    M.Xnames = ['Group Mean', vars];

    rng('default') % for reproducibility
    clear PEB BMA BMR
    close all
    field = {'A'};
    PEB = spm_dcm_peb(GRCM,M,field);
    %[BMC,PEB] = spm_dcm_bmc_peb(GRCM,M,field); % This one isn't usually run
%     [BMA,BMR] = spm_dcm_peb_bmc(PEB(1), GRCM(1,:));
% 
%     % ------------------------------------------------------------------------------
%     % Review results
%     % ------------------------------------------------------------------------------
%     pause(2)
%     GRCMs = GRCM(1,:);
%     spm_dcm_peb_review(BMA,GRCMs)
    % spm_dcm_peb_review(PEB,GRCMs)

    % ------------------------------------------------------------------------------
    % Save DCM outputs
    % ------------------------------------------------------------------------------
%     save('output_bayesian.mat','metadata','GRCMs','PEB','BMA','BMR','-v7.3')

% % ------------------------------------------------------------------------------
% % 2) secondary analysis excluding HCs and incorporating aggregate severity
% % Correlations to severity for clinical patients
% % ------------------------------------------------------------------------------
%     vars = {'Disinhibition_Pheno_exc','Compulsivity_Pheno_exc','Impulsivity_Pheno_exc'};
%     nuisance_vars = {'Age','Gender','IQ','fdJenk_m','Medicated'};
%     vars = [vars nuisance_vars];
% 
%     numVars = length(vars);
%     covs = zscore(metadata{metadata.Diagnosis == 2 | metadata.Diagnosis == 3,vars});
% 
%     % get severity ratings, normalise, and combine
%     sev_ocd = metadata(metadata.Diagnosis == 2,:).OCI_R_Total; sev_ocd = zscore(sev_ocd);
%     sev_gd = metadata(metadata.Diagnosis == 3,:).PGSI; sev_gd = zscore(sev_gd);
%     sev = [sev_ocd; sev_gd]; sev = zscore(sev);
% 
%     vars = ['Severity',vars];
%     covs = [sev,covs];
% 
%     % vars
%     [R,p] = corr(covs,sev)
% 
%     % Specify PEB model settings (see batch editor for help on each setting)
%     M = struct();
%     M.alpha = 1;
%     M.beta  = 16;
%     M.hE    = 0;
%     M.hC    = 1/16;
%     M.Q     = 'single';
% 
%     % Specify design matrix for N subjects. It should start with a constant column
%     M.X = [ones(size(covs,1),1) covs];
%     M.Xnames = ['Group Mean', vars];
% 
%     % Choose field
%     rng('default') % for reproducibility
%     clear PEB BMA BMR
%     close all
%     GRCM_clinicals = GRCM(metadata.Diagnosis == 2 | metadata.Diagnosis == 3,:);
%     field = {'A'};
%     PEB = spm_dcm_peb(GRCM_clinicals,M,field);
%     [BMA,BMR] = spm_dcm_peb_bmc(PEB(1), GRCM_clinicals(1,:));
% 
%     % ------------------------------------------------------------------------------
%     % Review results
%     % ------------------------------------------------------------------------------
%     pause(2)
%     GRCMs = GRCM_clinicals(1,:);
%     spm_dcm_peb_review(BMA,GRCMs)
%     % spm_dcm_peb_review(PEB,GRCMs)
% 
% % ------------------------------------------------------------------------------
% % Connection mask and names for below func. conn. analyses
% % ------------------------------------------------------------------------------
%     idx = [12,53];
%     % number of connections less the self connections
%     numNodes = size(GRCM{1,1}.Ep.A,1);
% 
%     % create a node-node connectivity string matrix for labels
%     nodeNames = {GRCM{1,1}.xY.name};
%     nodeNames = strrep(nodeNames,[WhichHemi,'_TriStri_',num2str(WhichStri),'-'],'');
%     connectionNames = cell(numNodes,numNodes);
%     for i = 1:numNodes
%         for j = 1:numNodes
%             connectionNames{i,j} = [nodeNames{i},'<-',nodeNames{j}];
%         end
%     end
% 
%     connectionNames = connectionNames(idx);
%     numConnections = length(connectionNames);
% 
% % ------------------------------------------------------------------------------
% % Run functional connectivity analysis
% % ------------------------------------------------------------------------------
%     % get FC
%     FC = zeros(numNodes,numNodes,numSubs);
%     FCv = zeros(numSubs,length(idx));
%     for i = 1:numSubs
%         fc = corr(GRCM{i,1}.Y.y);
%         fc = atanh(fc);
%         fc(eye(numNodes)==1) = 1;
%         FC(:,:,i) = fc;
%         FCv(i,:) = fc(idx);
%     end
%     clear fc
% 
%     % ------------------------------------------------------------------------------
%     % Fit linear model
%     % ------------------------------------------------------------------------------
%     vars = {'Disinhibition_Pheno_exc','Compulsivity_Pheno_exc','Impulsivity_Pheno_exc'};
%     nuisance_vars = {'Age','Gender','IQ','fdJenk_m','Medicated'};
%     vars = [vars nuisance_vars];
% 
%     numVars = length(vars);
%     covs = zscore(metadata{:,vars});
% 
%     betas = zeros(numVars,numConnections);
%     ps = zeros(numVars,numConnections);
%     for i = 1:numConnections
%         t2 = table();
%         for j = 1:numVars
%             t2 = [t2,table(covs(:,j),'VariableNames',{vars{j}});];
%         end
%         t = [t2,table(FCv(:,i),'VariableNames',{'Connection'})];
% 
%         mdl = fitglm(t);
% 
%         b = mdl.Coefficients.Estimate;
%         b(1) = [];
%         betas(:,i) = b;
%         
%         p = mdl.Coefficients.pValue;
%         p(1) = [];
%         ps(:,i) = p;
%     end
% 
%     % ------------------------------------------------------------------------------
%     % Plot
%     % ------------------------------------------------------------------------------
%     % plotMat = R;
%     plotMat = betas;
% 
%     figure('color','w'); box('on'); hold on;
%     imagesc(plotMat)
%     title('Correlations between connections and Jeg''s phenotypes');
%     xlabel('Connections')
%     ylabel('Phenotypes')
%     colormap([flipud(BF_getcmap('blues',6,0));1,1,1;BF_getcmap('reds',6,0)])
%     caxis([-0.2 0.2])
%     axis square
%     axis tight
%     colorbar
% 
%     ax = gca;
%     % Set axis stuff
%     ax.XTick = 1:length(connectionNames);
%     % connectionNames_temp = strrep(connectionNames,'_','\_');
%     connectionNames_temp = strrep(connectionNames,'<-','<->');
%     ax.XTickLabel = strrep(connectionNames_temp,'_','\_');
%     ax.XTickLabelRotation = 45;
%     ax.YTick = 1:length(vars);
%     ax.YTickLabel = strrep(vars,'_','\_');
% 
%     % plot values
%     for i = 1:size(plotMat,1)
%         for j = 1:size(plotMat,2)
%             text(j,i,num2str(plotMat(i,j),'%0.2f'),'HorizontalAlignment','center',...
%                 'Color','k','FontSize',15,'FontWeight','normal');
%         end
%     end
