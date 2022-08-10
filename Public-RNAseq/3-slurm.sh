#!/bin/bash
# Start the RNA-Browser mapping etc steps

conda activate ~/mysnake

for i in 3-RNA-Browser/* ; do
    echo $i
    cd $i 
    #git pull
    #snakemake --unlock
    bash run_slurm_singularity.sh
    cd ../..
    sleep 5
done 

echo 'jobs submitted.'
