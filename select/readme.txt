How to select special gene annotations from snpeff vcf file:
1. Generate a gene name list in CSV format.
2. Copy the two files "ann_select.py" and "genelist.csv" to ".../output_*/variants/"
3. Get ".../output_*/variants/" directory and run the command:
         $python3 ann_select.py
4. Selected gene annotations can be found in the same directory.