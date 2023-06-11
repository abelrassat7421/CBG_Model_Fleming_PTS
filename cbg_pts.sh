#!/bin/bash -l
#SBATCH --job-name=pts
#SBATCH -N 10
#SBATCH --array=401-600
#SBATCH --ntasks-per-node=1
#SBATCH -t 3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=abel.rassat@ucdconnect.ie

module load openmpi/3.1.4
source ~/cbg-environment/bin/activate
export PYTHONUNBUFFERED=yes
export PYNN_OUTPUT_DIRNAME=~22213094/scratch/PTS_1_3mA/PTS-$((SLURM_ARRAY_TASK_ID+600))
configfile=~22213094/CBG_Model_Fleming_PTS/configs/pts_test_${SLURM_ARRAY_TASK_ID}.yml

cd ~22213094/CBG_Model_Fleming_PTS/Cortex_BasalGanglia_DBS_model

mpirun -n ${SLURM_NPROCS} ./run_model.py -o ${PYNN_OUTPUT_DIRNAME} ${configfile}
cp ~22213094/CBG_Model_Fleming_PTS/slurm-${SLURM_ARRAY_JOB_ID}_$((SLURM_ARRAY_TASK_ID+600)).out ${PYNN_OUTPUT_DIRNAME}

