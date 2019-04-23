#!/bin/bash

module load  matlab/r2018a

cd /scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/first_level_dcm

units=scans
subject=sub-494

WhichDenoise=nogsr
derivativesdir="/projects/kg98/kristina/GenofCog/datadir/derivatives/"$subject
scratchsub="/scratch/kg98/kristina/Projects/GenofCog/derivatives/"$subject
DCMfirstleveldir=$scratchsub"/DCM_project/firstlevel_dcm/"$WhichDenoise"_CST_fullyconn/"

preprodir=$derivativesdir"/"$subject"/prepro.feat/preprocessed"
epi=$subject"_filtered_func_data_clean_mni_smooth4mm.nii"
data=$subject"_filtered_func_data_clean_mni_smooth4mm_00000"
VolDir=$preprodir"/VolumeSplit/"$WhichDenoise"_smooth4mm/"



if [ ! -d $DCMfirstleveldir ]; then mkdir -p $DCMfirstleveldir; echo "$subject - making directory"; fi

matlab -nodisplay -r "spDCM_FirstLevel_PsychosisCST('$VolDir', '$data', '$units', 616, 0.754, 0.021, '$DCMfirstleveldir'); exit"

echo -e "\t\t\t ----- DONE ----- "
