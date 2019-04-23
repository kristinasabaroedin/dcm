clear all
clc

fileID = fopen('/projects/kg98/kristina/GenofCog/scripts/sublists/PLEs_n352.txt');
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};
% compute numsubs
numSubs = length(ParticipantIDs);

WhichDCM = 1;
Denoise = {'nogsr', 'gsr4'};
WhichDenoise= Denoise{1};

spmdir = '/projects/kg98/kristina/spm12';
addpath(spmdir)

for i=1:numSubs
    subject = ParticipantIDs{i};
    scratchsub = ['/scratch/kg98/kristina/Projects/GenofCog/derivatives/',subject,'/'];
    DCMfirstleveldir = [scratchsub,'DCM_project/firstlevel_dcm/',WhichDenoise,'_18conns/'];
    
    cd(DCMfirstleveldir)
    load('DCM.mat')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Number of estimable parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    qE    = spm_vec(DCM{WhichDCM}.Ep);
    pE    = spm_vec(DCM{WhichDCM}.M.pE);
    qC    = DCM{WhichDCM}.Cp;
    pC    = DCM{WhichDCM}.M.pC;
    k     = rank(full(pC));
    pC    = pinv(pC);
    
    D  = trace(pC*qC) + (pE - qE)'*pC*(pE - qE) - spm_logdet(qC*pC) - k;
    N_est_par(i,1)  = full(D/log(DCM{WhichDCM}.v));
    
    clear qE pE qC pC k D;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %maximum extrinsic connection strength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Max_conn(i,1)=max(max(abs(DCM{WhichDCM}.Ep.A - diag(diag(DCM{WhichDCM}.Ep.A)))));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Proportion explained variance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PSS   = sum(sum(sum(abs(DCM{WhichDCM}.Hc).^2)));
    RSS   = sum(sum(sum(abs(DCM{WhichDCM}.Rc).^2)));
    Expl_var(i,1) = 100*PSS/(PSS + RSS);
    
    clear RSS PSS;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Free energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Free_energy(i,1)=DCM{WhichDCM}.F;

    
    %%%%%%%%
    %Motion
    %%%%%%%%

%     rp_name='rp_afunctional_disc_4D.txt';
%     rp = load(deblank(rp_name));
%     radius=50;
%     
%     %Change rotation in radials to mm (arch=radial*radius)
%     rp(:,4:6)=rp(:,4:6)*radius;
%     
%     %Calculate change in position between current and previous volume
%     rp_diff=diff(rp);
%     
%     %add zero as first element (see Power et al., 2014)
%     rp_zero=[zeros(1,size(rp_diff,2)); rp_diff];
%     
%     %compute framewise displacement (FD)
%     FD(:,session)=sum(abs(rp_zero),2);
% 
%     clear rp rp_diff rp_zero
end