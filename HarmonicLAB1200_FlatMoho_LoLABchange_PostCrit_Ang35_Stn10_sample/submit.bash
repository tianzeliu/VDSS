#!/bin/bash
#SBATCH -J run.output
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=16
#SBATCH --threads-per-core=1
#SBATCH --mem=4GB
#SBATCH --time=08:00:00
#SBATCH --output run.output
#SBATCH --partition=serc
#SBATCH --mail-type=END

module load openmpi

dir_bin=$HOME/specfem2d/bin


mkdir -p OUTPUT_FILES
mkdir -p DATA

rm -rf OUTPUT_FILES/*
rm -rf DATA/*

cd DATA/
rm -f Par_file SOURCE
ln -s ../Par_file Par_file
ln -s ../SOURCE SOURCE
cd ../



rm -f xmeshfem2D xspecfem2D
ln -s $dir_bin/xmeshfem2D
ln -s $dir_bin/xspecfem2D

srun --mpi=pmi2 -n $SLURM_NTASKS ./xmeshfem2D
srun --mpi=pmi2 -n $SLURM_NTASKS ./xspecfem2D

cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES
