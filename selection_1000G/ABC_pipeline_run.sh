#!/bin/bash
# Run nextflow on sbatch
#SBATCH -J run_abc.nf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -p hbfraser,hns,normal
#SBATCH -t 47:59:00
#SBATCH -o sbatch.out
#SBATCH -e sbatch.err

module load devel java/11.0.11 conda/4.8.2
export NXF_OPTS='-Xms1G -Xmx4G'
export NXF_CONDA_CACHEDIR="/home/users/kpettie/conda_nf/"
alias sbatch="./sbatch_repeat.sh"

nextflow -c ABC_pipeline.config run ABC_pipeline.nf -w /scratch/users/kpettie/nextflow/abc_1000G/work -resume

echo "Completed with code $err!"
exit $err
