clear all;
close all;

parentdir_scratch = '/scratch/kg98/kristina/';
%------------------------------------------------------------------------------

scriptdir = [parentdir_scratch,'Projects/GenofCog/scripts/dcm_project/'];
addpath(genpath(scriptdir))

spmdir = [parentdir_scratch,'spm12_r7487'];
addpath(spmdir)

% ------------------------------------------------------------------------------
% Set options
% ------------------------------------------------------------------------------
Projects = {'GenofCog', 'UCLA', 'STAGES'};
WhichProject = Projects{1}

Hemi = {'Left','Right'}
WhichHemi = Hemi{1};

Denoise = {'gsr4', 'nogsr', 'dbscan', '2p+AROMAaggr'};
WhichDenoise = Denoise{1}

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
switch WhichProject
    case 'GenofCog'
        projdir = [parentdir_scratch,'Projects/GenofCog/derivatives/'];
        datafile = [projdir,'secondlevel/DCM/GoC_metadata_n353.mat'];
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
outdir = [parentdir_scratch,'/Projects/',WhichProject,'/derivatives/secondlevel/DCM/',WhichDenoise,'/cst_fullyconn/'];

if exist(outdir) == 0
    fprintf(1,'Initialising outdir\n')
    mkdir(outdir)
elseif exist(outdir) == 7
    fprintf(1,'Cleaning and re-initialising outdir\n')
%   rmdir(outdir,'s')
%   mkdir(outdir)
end


% ------------------------------------------------------------------------------
% Non-imaging metadata
% ------------------------------------------------------------------------------
load(datafile)
numSubs = size(metadata.SubjectID,1);
numGroups = numel(unique(metadata.Diagnosis));

% ------------------------------------------------------------------------------
% Create demographic table
% ------------------------------------------------------------------------------
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

%------------------------------------------------------------------------------
% First provide estimated DCM as GCM cell array. Individual DCMs can be
% estimated by using spm_dcm_fit.m

for i = 1:120
    GCM{i,1} = [projdir,metadata.SubjectID{i},'/DCM_project/firstlevel_dcm/',WhichDenoise,'_CST_fullyconn/DCM.mat'];
end

% ------------------------------------------------------------------------------
% Define covariates
% ------------------------------------------------------------------------------
vars = {'PC1','PC2'};
nuisance_vars = {'Age','Sex','IQ','MeanFD'};
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
    
M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';



% Specify design matrix for N subjects. It should start with a constant column
M.X = [ones(numSubs,1) covs];
M.Xnames = ['Group Mean', vars];

% Choose field
field = {'A'};

% Estimate model
PEB = spm_dcm_peb(GCM,M,field);


BMA = spm_dcm_peb_bmc(PEB);

cd(outdir)
save('PEB_fullyconn_CST.mat','BMA','PEB','GCM');

spm_dcm_peb_review(BMA,GCM)

% % Perform leave-one-out cross validation (GCM,M,field are as before)
% spm_dcm_loo(GCM,M,field);
