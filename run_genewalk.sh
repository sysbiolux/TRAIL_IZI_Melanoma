#!/bin/bash -l
#SBATCH -J 200602
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=aurelien.ginolhac@uni.lu
#SBATCH --qos=normal
#SBATCH -c 6
#SBATCH --mem-per-cpu=4096
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-05:03:00
#SBATCH -p batch


module load bio/SAMtools
conda activate genewalk

genewalk --project IZI --base_folder /scratch/users/aginolhac/200602/  --genes same_pdirection.txt --id_type hgnc_symbol --nproc 6
