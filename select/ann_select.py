#!/usr/bin/env python

import glob

file1 = open('./genelist.csv')
#file2 = open('./sequence.fasta')
filelist=glob.glob("*_snpeff.ann.vcf")

content1 = file1.readlines()
#content2 = file2.readlines()

#### read 1st line from the file
#print(content2[0])

for afile in filelist:
   file2=open(afile)
   content2 = file2.readlines()
   for aline in content2:
      if aline.startswith("##"):
         with open("selected_"+afile, "a") as f:
              f.write(aline)
      else:
         genelist=content1[0].strip().split(",")
         for agene in genelist:
              if agene in aline:
                 with open("selected_"+afile, "a") as f:
                      f.write(aline)


