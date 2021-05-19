#!/bin/bash

#PBS -l walltime=00:01:00,nodes=2:ppn=4
#PBS -N lab_job
#PBS -q batch

cd $PBS_O_WORKDIR
mpirun --hostfile $PBS_NODEFILE -np 8 ./lab 1000 1000 0.1 0.1
