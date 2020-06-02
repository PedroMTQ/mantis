#!/bin/bash -l
#SBATCH -J 25c250w
#SBATCH --time=1-00:00:00
#SBATCH --mem=100GB
#SBATCH -p batch
#SBATCH --qos=qos-batch
conda activate mantis_env
python mantis run_mantis -t mantis/tests/test_sample.faa
python mantis run_mantis -t mantis/tests/test_file.tsv