#!/bin/bash -l
#SBATCH -J 25c250w
#SBATCH --time=0-00:10:00
#SBATCH --mem=100GB
#SBATCH -p batch
#SBATCH --qos=qos-batch

#change the time accordingly. a genome should be analyzed in around 10-30 minutes. A metagenome will take longer, can be around 15-20 hours but of course it depends a lot on the size of the metagenome.
#If you want faster results consider removing EGGNOG
conda activate mantis_env
python mantis run_mantis -t mantis/tests/test_sample.faa
python mantis run_mantis -t mantis/tests/test_file.tsv