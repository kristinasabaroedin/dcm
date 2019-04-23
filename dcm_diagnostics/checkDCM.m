fileID = fopen('/projects/kg98/kristina/GenofCog/scripts/sublists/PLEs_n353.txt');
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};
% compute numsubs
numSubs = length(ParticipantIDs);

connections = {'11','21','31','41','51','12','22','32','42','52','13','23','33','43','53','14','24','34','44','54',...
'15','25','35','45','55'};    
for i = 1:numSubs
    subdir1 = ['/scratch/kg98/kristina/Projects/GenofCog/derivatives/',ParticipantIDs{i},'/DCM_project/firstlevel_dcm/nogsr_DCleftdACC_02/'];
    load([subdir1,'DCM.mat'])
    A = DCM.Ep.A;
    A2_mat_reshape(i,:) = reshape(A,1,25);
    clear DCM A
    subdir2 = ['/scratch/kg98/kristina/Projects/GenofCog/derivatives/',ParticipantIDs{i},'/DCM_project/firstlevel_dcm/nogsr_DCleftdACC/'];
    load([subdir2,'DCM.mat'])
    A = DCM.Ep.A;
    A1_mat_reshape(i,:) = reshape(A,1,25);
    clear DCM A
end

for k = 1:25;
    correlations(k) = corr(A1_mat_reshape(:,k),A2_mat_reshape(:,k));
end