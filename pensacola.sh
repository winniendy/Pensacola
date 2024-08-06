#!/usr/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=pensacola
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300gb
#SBATCH --time=48:00:00
#SBATCH --output=pensacola.%j.out
#SBATCH --error=pensacola.err
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=FAIL,END

module load nextflow
module load longqc
APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR

nextflow run pensacola.nf -params-file params.yaml
mv ./*.out ./output
mv ./*err ./output

#gfa to fa
mkdir -p ./output/assemble
cp ./output/*/assemble/*.bp.p_ctg.gfa ./output/assemble
gfas=`ls ./output/assemble/*.gfa`
for eachfile in $gfas
do
  #echo $eachfile
  gawk '/^S/{print ">"$2"\n"$3}' $eachfile|fold > ${eachfile}.fa
done

# drug resistant detect
python drugres.py


dt=$(date "+%Y%m%d%H%M%S")
mv ./output ./output-$dt
#mv ./work ./work-$dt
rm -r ./work
rm -r ./cache