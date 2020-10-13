#!/bin/bash
#SBATCH -A oxwasp                       
#SBATCH --partition=grey-fast
#SBATCH --nodelist=grey01.cpu.stats.ox.ac.uk
#SBATCH -J inla01                         
#SBATCH --time=3-00:00:00                   
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40000
#SBATCH --mail-user=juha.kreula@stats.ox.ac.uk     
#SBATCH --mail-type=ALL                   
#SBATCH --chdir="/data/localhost/not-backed-up/jukreula/inla/"
#SBATCH --output="/data/localhost/not-backed-up/jukreula/inla/Output/slurm-%u-%A-std-output"
#SBATCH --error="/data/localhost/not-backed-up/jukreula/inla/Output/slurm-%u-%A-err-output"

Rscript inla_main.R
