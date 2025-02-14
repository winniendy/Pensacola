#!/bin/bash

# Directory containing BAM files
BAM_DIR="./pbbams"
OUTPUT_DIR="./output"

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Loop through each BAM file in the directory
for BAM_FILE in ${BAM_DIR}/*.bam; do
    # Extract base name of the file (without extension)
    BASE_NAME=$(basename ${BAM_FILE} .bam)
    
    # Convert BAM to FASTQ and save output in the specified output directory
    bam2fastq -o ${OUTPUT_DIR}/${BASE_NAME} ${BAM_FILE}
done
