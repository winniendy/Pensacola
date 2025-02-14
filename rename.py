#!/usr/bin/env python
import os

file1 = open('./barcoded_samples.csv')
list_files = os.listdir('./')
content1 = file1.readlines()
#content2 = file2.readlines()
#print(list_files)
#print(content1)

for i in range(len(list_files)):
   for pair in content1:
      items=pair.strip().split(",")
      #print(items[0])
      if items[0] in list_files[i]:
         if "pbi" in list_files[i]:
             os.rename(list_files[i],items[1]+".bam.pbi")
         else:
             os.rename(list_files[i],items[1]+".bam")
    


