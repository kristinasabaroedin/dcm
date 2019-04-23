#!/bin/env bash

#SBATCH --job-name=firstlevel
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
#SBATCH --array=1-352%44


declare -i ID=$SLURM_ARRAY_TASK_ID
echo -e $SLURM_ARRAY_TASK_ID

module load  matlab/r2018a


SUBJECT_LIST="/scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/do_fc_firstlevel.txt"


cd /scratch/kg98/kristina/Projects/GenofCog/scripts/dcm_project/first_level_fc


subject=$(sed -n "${ID}p" ${SUBJECT_LIST})
echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${ID} ${subject} ----- "
echo -e "\t\t\t --------------------------- \n"

projectdir=/projects/kg98/kristina/GenofCog
projectdir_scratch=/scratch/kg98/kristina/Projects/GenofCog
derivativesdir=$projectdir_scratch'/derivatives/'${subject}'/DCM_project' 
firstleveldir=$derivativesdir'/firstlevel_fc/'

preprodir=$projectdir'/datadir/derivatives/'$subject'/prepro.feat/preprocessed/'
tsdir=$derivativesdir'/timeseries/'
epi=${subject}'_filtered_func_data_clean_mni_smooth4mm.nii'
data=${subject}'_filtered_func_data_clean_mni_smooth4mm_00000'
volumedir=$preprodir'VolumeSplit/nogsr_smooth4mm/'



#if [ ! -d $firstleveldir ]; then mkdir -p $firstleveldir; echo "$subject - making directory"; fi

matlab -nodisplay -r "first_level_subject_dcm('$subject', 'freesurfer_leftthalamus_smooth4mm_nogsr', 0.754, 616, '$firstleveldir', '$volumedir', '$tsdir', '$data'); first_level_subject_dcm('$subject', 'freesurfer_rightthalamus_smooth4mm_nogsr', 0.754, 616, '$firstleveldir', '$volumedir', '$tsdir', '$data'); exit"


echo -e "\t\t\t ----- DONE ----- "
