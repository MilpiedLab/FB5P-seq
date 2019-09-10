#!/bin/bash

module load userspace/all
module load python3/3.6.3
source /home/cdong/snakemake/bin/activate

snakemake --jobs 480 --use-singularity --singularity-args "-B /scratch/cdong:/scratch/cdong" --snakefile snakefile/snakefile_all_BCR.yml --configfile config/config_BCR.yml --cluster-config config/cluster.json --cluster 'sbatch -A {cluster.project} --job-name {cluster.job-name} --ntasks {cluster.ntasks} --cpus-per-task {threads} --mem-per-cpu {cluster.mem-per-cpu} --partition {cluster.partition} --time {cluster.time} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --error {cluster.error} --output {cluster.output}'

