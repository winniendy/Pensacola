#!/usr/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=CandidaPB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --output=sanibelpb.%j.out
#SBATCH --error=sanibelpb.err
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=FAIL,END

module load nextflow
APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR


nextflow run candidapb.nf -params-file params.yaml

mv ./*.out ./output
mv ./*err ./output

dt=$(date "+%Y%m%d%H%M%S")
mv ./output ./output-$dt
rm -r ./work
