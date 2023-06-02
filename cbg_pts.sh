#!/bin/bash -l
#SBATCH --job-name=pts
# speficity number of nodes 
#SBATCH -N 20
#SBATCH --array=1-3

# specify number of tasks/cores per node required
#SBATCH --ntasks-per-node 1

# specify the walltime e.g 20 mins
#SBATCH -t 2:00:00

# set to email at start,end and failed jobs
#SBATCH --mail-type=ALL
#SBATCH --mail-user=abel.rassat@ucdconnect.ie

module load openmpi/3.1.4
source ~/cbg-environment/bin/activate
export PYTHONUNBUFFERED=yes
export PYNN_OUTPUT_DIRNAME=~22213094/scratch/PTS_test1/PTS-$(date +"%Y%m%d%H%M%S")-${SLURM_ARRAY_TASK_ID}
configfile=~22213094/CBG_Model_Fleming_PTS/configs/pts_test_${SLURM_ARRAY_TASK_ID}.yml

cd ~22213094/CBG_Model_Fleming_PTS/Cortex_BasalGanglia_DBS_model

# command to use
mpirun -n ${SLURM_NPROCS} ./run_model.py -o ${PYNN_OUTPUT_DIRNAME} ${configfile}
cp ~22213094/CBG_Model_Fleming_PTS/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${PYNN_OUTPUT_DIRNAME}