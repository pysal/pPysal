#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=2:ppn=8
#PBS -N Cooperative Max-P
#PBS -o out
#PBS -e err
#PBS -A Jay
#PBS -l walltime=00:20:00
use openmpi-1.8.1
cd $PBS_O_WORKDIR
/packages/openmpi-1.8.1/bin/mpiexec -n 16 python mpicluster.py
