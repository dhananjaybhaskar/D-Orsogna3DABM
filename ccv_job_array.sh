#!/bin/bash
#SBATCH -J MATLAB
#SBATCH -t 12:00:00
#SBATCH --array=1-25

# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID
#SBATCH -e arrayjob-%a.err
#SBATCH -o arrayjob-%a.out

echo "Starting job $SLURM_ARRAY_TASK_ID on $HOSTNAME"
matlab -nodisplay -r "sim_3D_model($SLURM_ARRAY_TASK_ID,3); exit"
