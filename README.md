# Pensacola
A Nextflow pipeline to analyze Candida auris data from PacBio sequencing. 

## Prerequisites
Nextflow is needed. The detail of installation can be found in https://github.com/nextflow-io/nextflow. For HiPerGator users, its installation is not needed. 

Singularity/APPTAINER is needed. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/. For HiPerGator users, its installation is not needed.

SLURM is needed. For HiPerGator users, its installation is not needed.

Python3 is needed. The package "pandas" should be installed by ``` pip3 install pandas ``` if not included in your python3.

PacBio SMRTLINK stand-alone tools are needed. About how to install them, please see the file "How_to_install_smrtlink_tools.txt" in the pipeline.

## Recommended conda environment installation
   ```bash
   conda create -n PENSACOLA -c conda-forge python=3.10 pandas
   ```
   ```bash
   conda activate PENSACOLA
   ```
## How to run

1. Rename your data files and make them looks like "bc2024bc2024.bam.pbi" and "bc2024bc2024.bam". You can use to the script "rename.sh" in the pipeline to rename your data files.
2. put the renamed data files (*.bam and *.bam.pbi) into the directory /pbbams.
3. open file "params.yaml", set the three parameters absolute paths. They should be ".../.../pbbams",  ".../.../output", and ".../.../reference" . 
4. get to the top directory of the pipeline, run 
```bash
sbatch ./candidapb.sh
```
#### Note: If you want to get email notification when the pipeline running ends, please input your email address in the line "#SBATCH --mail-user=<EMAIL>" in the batch file that you will run (namely, candidapb.sh). 
