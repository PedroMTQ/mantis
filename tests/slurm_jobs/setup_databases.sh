#!/bin/bash -l
#SBATCH -J setup_mantis
#SBATCH --time=1-00:00:00
#SBATCH --mem=32GB
#SBATCH -p batch
#SBATCH --qos=qos-batch
conda activate mantis_env
python mantis setup_databases
