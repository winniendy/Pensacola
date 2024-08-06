#!/usr/bin/env python

import glob

with open("./output/Drug_resistance_report.txt", "w") as f:
        f.write("########## Detected drug resistance ##########"+"\n\n")
        #f.write("Detected drug resistances:  "+"%\n\n")


file1 = open('./drugres/drugres_snpeff_list.csv')
filelist=glob.glob("./output/variants/*_snpeff.ann.vcf")

content1 = file1.readlines()
druglistlen = len(content1)
#content2 = file2.readlines()

#### read 1st line from the file
#print(content2[0])
filenum=0
for afile in filelist:
   filenum += 1
   #print(afile)
   paths = afile.strip().split("/")
   with open("./output/Drug_resistance_report.txt", "a") as f:
      f.write(str(filenum)+". "+paths[-1]+": \n\n")
      f.write(content1[0])
   file2=open(afile)
   content2 = file2.readlines()
   for aline in content2:
      if not aline.startswith("#"):
         tabitems=aline.strip().split("\t")
         subtabs=tabitems[7].strip().split(";")
         for cell in subtabs:
            if "ANN=" in cell:
                #anns = subtabs[10].strip().split(",")
                anns = cell.strip().split(",")
                for oneann in anns:
                    for indx in range(1,druglistlen):
                        cols=content1[indx].strip().split(",")
                        muts=[]
                        if ";" in cols[3]:
                            cells=cols[3].strip().split(";")
                            for x in cells:
                                muts.append(x)        
                        else:
                            muts.append(cols[3])
                        if cols[0] in oneann:
                            #print("ID found: "+cols[0])
                            for amut in muts:
                                if amut in oneann:
                                   # print("Mutation found: "+amut)
                                   # print("Drug resistant found: "+cols[1])
                                   # print(content1[indx])
                                   with open("./output/Drug_resistance_report.txt", "a") as f:
                                        f.write(content1[indx]+"\n\n")
                        else:
                            continue


