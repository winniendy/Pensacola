# Pensacola
A Nextflow pipeline to analyze *Candida auris* long-sequencing data from PacBio. The QC, species identification and species abundance of a sample, contig assembly, SNP calling and annotation, drug resistance detection, *etc* are included in the pipeline.  

## Prerequisites
Nextflow is needed. The detail of installation can be found in https://github.com/nextflow-io/nextflow. For HiPerGator users, its installation is not needed. 

Singularity/APPTAINER is needed. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/. For HiPerGator users, its installation is not needed.

SLURM is needed. For HiPerGator users, its installation is not needed.

Python3 is needed. The package "pandas" should be installed by ``` pip3 install pandas ``` if not included in your python3.

LongQC is needed. Please install it to your local computer from its github repository (https://github.com/yfukasawa/LongQC). For HiPerGator users, its installation is not needed.

PacBio SMRTLINK stand-alone tools are needed. About how to install them, please see the file "How_to_install_smrtlink_tools.txt" in the pipeline.

The Kraken2/Bracken Refseq index--PlusPF is needed. Please download PlusPF index (over 77 GB) from the link (https://benlangmead.github.io/aws-indexes/k2) to the "PlusPF" folder in your local computer. And then extract the tar.gz archive. For HiPerGator users, its downloading is not needed. It has been downloaded and configed in the pipeline.

## Recommended conda environment installation
   ```bash
   conda create -n PENSACOLA -c conda-forge python=3.10 pandas
   ```
   ```bash
   conda activate PENSACOLA
   ```

## Reference
1. In this pipeline, the latest genome annotation of *Candida.auris* strain B8441 is used as the reference. They are available from http://www.candidagenome.org/download/sequence/C_auris_B8441/current/.
2. If you want to use another *Candida.auris* strain as the reference in the pipeline, please refer to the guide "How to build custom snpeff database using snpeff docker image in the pipeline.txt" in the pipeline.

## How to run

1. Rename your data files and make them looks like "bc2024bc2024.bam.pbi" and "bc2024bc2024.bam". You can use to the script "rename.sh" in the pipeline to rename your data files.
2. put the renamed data files (*.bam and *.bam.pbi) into the directory /pbbams.
3. open file "params.yaml", set the full paths of the parameters.   
   **input** : the full path to pbbams dir of the pipeline in your computer. It looks like "/\<the full path to the pipeline dir in your computer\>/pbbams".    
   **output** : the full path to output dir of the pipeline in your computer. It looks like "/\<the full path to the pipeline dir in your computer\>/output".            
   **reference** : the full path to reference dir of the pipeline in your computer. It looks like "/\<the full path to the pipeline dir in your computer\>/reference".    
   **snpeffconfig** : the full path to configs dir of the pipeline in your computer. It looks like "/\<the full path to the pipeline dir in your computer\>/configs".      
          
   **db** : the full path to kraken/bracken-database (PlusPF) in your computer. It looks like "/\<the full path to the parent dir of PlusPF foler in your computer\>/PlusPF".    
   **qc** : the full path to LongQC dir in your computer. It looks like "/\<the full path to the parent dir of LongQC foler in your computer\>/LongQC".     
           
   **Note: For HiperGator users, the parameters "db" and "qc" do not need to be changed. Just keep default settings.**     

4. get to the top directory of the pipeline, run 
```bash
sbatch ./pensacola.sh
```
#### Note1: If you want to get email notification when the pipeline running ends, please input your email address in the line "#SBATCH --mail-user=<EMAIL>" in the batch file that you will run (namely, pensacola.sh). 
