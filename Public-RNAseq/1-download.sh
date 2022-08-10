#!/bin/bash
# First init the folder structure with 1-prep-schlange.R

conda activate ~/mysnake

# init donwload procedure of all RNA-seq libraries
for i in 1-RNA-Schlange/* ; do
    echo $i
    cd $i 
    git pull
    rm samples.csv
    snakemake --unlock
    rm samples.csv
    sbatch --ntasks=1 --job-name=download --partition=bullseye \
        --cpus-per-task=8 --mem=8GB --time=20:00:00            \
        --output='log-download.out.txt'                        \
        --error='log-download.err.txt'                         \
        run_singularity.sh
    cd ../..
done 

echo 'jobs submitted.'
