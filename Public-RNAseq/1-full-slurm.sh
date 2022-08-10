#!/bin/bash
# After download run the full pipeline per genome

conda activate ~/mysnake

# init donwload procedure of all RNA-seq libraries
for i in 1-RNA-Schlange/* ; do
    echo $i
    cd $i 
    #git pull
    #snakemake --unlock
    bash run_slurm_singularity.sh
    cd ../..
    sleep 5
done 

echo 'jobs submitted.'
