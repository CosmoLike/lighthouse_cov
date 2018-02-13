#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -q longq
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=48:00:00
#PBS -J 1-10000
#PBS -N DEScov_1-1035k
#PBS -e /aurora_nobackup/sunglass/teifler/output/
#PBS -o /aurora_nobackup/sunglass/teifler/output/

cd $PBS_O_WORKDIR
/home/teifler/CosmoLike/lighthouse_cov/./multi_covariances_real_mpp $PBS_ARRAY_INDEX >& /aurora_nobackup/sunglass/teifler/job_output.log


