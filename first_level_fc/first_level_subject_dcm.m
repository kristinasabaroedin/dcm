function [] = first_level_subject_dcm(subject, seed, TR, N, firstleveldir, volumedir, tsdir, data)
% Runs first level analysis that loops over a subject list
%     addpath('/scratch/kg98/kristina/spm12_r7487')
%     addpath('/scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/first_level_fc')
    
	cfg.mat = 0; % if a cfg.mat containing time-series exists

	display(subject)
	projdir=['/projects/kg98/GenofCog/datadir/derivatives/',subject,'/']

    outdir = [firstleveldir,seed,'/'];

    
    cd (firstleveldir)


 	if exist(outdir) == 0;
 		mkdir (outdir)
    elseif exist(outdir) == 7
 		fprintf(1,'Cleaning and re-initialising outdir\n')
 		rmdir(outdir,'s')
 		mkdir(outdir)
    end
 
	cd(tsdir)
	ts = dlmread([seed,'.txt']);
	R = ts; % get ts, save it to a matrix called R   
	save([seed,'_ts.mat'],'R')
	cd (outdir)
	movefile ([tsdir,seed,'_ts.mat'], outdir)

	clear R
	clear ts


	% Estimate first-level analysis for hemisphere
	first_level_3D(outdir,'scans',TR,[volumedir,data],N,[outdir,seed,'_ts.mat'],'')
	% Run contrasts
	FirstLevelContrasts_dcm(seed,[outdir,'/SPM.mat'])


	display ('First-level analysis done')

end
