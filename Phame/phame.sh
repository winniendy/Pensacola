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

#rm -rf test/workdirs
export PATH=src:$PATH

############## #3 Test with just contigs using ANI based reference ##############
#echo "Test with just contigs ";
#mkdir -p ./test/workdirs/mytest
#cp ./test/data/*.contig test/workdirs/mytest/

src/phame test/ctl_files/candida.ctl
	