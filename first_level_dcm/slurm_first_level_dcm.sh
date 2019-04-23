#!/bin/env bash

#SBATCH --job-name=firstlevelDCM
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=11:00:00
#SBATCH --mail-user=kristina.sabaroedin@monash.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem=56384
#SBATCH --qos=normal
#SBATCH -A kg98
#SBATCH --array=1-342%19


declare -i ID=$SLURM_ARRAY_TASK_ID
echo -e $SLURM_ARRAY_TASK_ID

module load  matlab/r2018a


SUBJECT_LIST="/scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/do_fc_firstlevel.txt"


cd /scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/first_level_dcm


subject=$(sed -n "${ID}p" ${SUBJECT_LIST})
echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${ID} ${subject} ----- "
echo -e "\t\t\t --------------------------- \n"


units='scans'


WhichDenoise=nogsr
derivativesdir='/projects/kg98/kristina/GenofCog/datadir/derivatives/'${subject}
scratchsub='/scratch/kg98/kristina/Projects/GenofCog/derivatives/'${subject}
DCMfirstleveldir=$scratchsub'/DCM_project/firstlevel_dcm/'$WhichDenoise'_CST_fullyconn/'

preprodir=$derivativesdir'/prepro.feat/preprocessed/'
epi=${subject}'_filtered_func_data_clean_mni_smooth4mm.nii'
data=${subject}'_filtered_func_data_clean_mni_smooth4mm_00000'
VolDir=$preprodir'VolumeSplit/'$WhichDenoise'_smooth4mm/'



if [ ! -d $DCMfirstleveldir ]; then mkdir -p $DCMfirstleveldir; echo "$subject - making directory"; fi

matlab -nodisplay -r "spDCM_FirstLevel_PsychosisCST('$VolDir', '$data', '$units', 616, 0.754, 0.021, '$DCMfirstleveldir'); exit"

echo -e "\t\t\t ----- DONE ----- "
