#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -q mediumq
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=12:00:00
#PBS -J 1-5886
#PBS -N HSC
#PBS -e /nobackup0/sunglass/hironaom/output/ 
#PBS -o /nobackup0/sunglass/hironaom/output/ 

cd $PBS_O_WORKDIR
/home/hironaom/mg/cosmolike/top-level/miyatake/./compute_covariances_fourier $PBS_ARRAY_INDEX >& /nobackup0/sunglass/hironaom/job_output.log 
