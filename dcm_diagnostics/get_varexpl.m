fileID = fopen('/projects/kg98/kristina/GenofCog/scripts/sublists/PLEs_n353.txt');
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};
% compute numsubs
numSubs = length(ParticipantIDs);

for i = 1:numSubs
    load(['/scratch/kg98/kristina/Projects/GenofCog/derivatives/',ParticipantIDs{i},'/DCM_project/firstlevel_dcm/nogsr_CST_fullyconn/DCM.mat'])
    var_expl(i,1) = DCM.diagnostics(1)
    clear DCM
    load(['/scratch/kg98/kristina/Projects/GenofCog/derivatives/',ParticipantIDs{i},'/DCM_project/firstlevel_dcm/gsr4_CST_fullyconn/DCM.mat'])
    var_expl(i,2) = DCM.diagnostics(1)
    clear DCM
end 