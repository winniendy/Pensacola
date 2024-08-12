      refdir = ../data/candida_ref  # directory where reference files are located
      workdir = ../workdirs/mytest # directory where contigs/reads files are located and output is stored

    reference = 1  # 0:pick a random reference; 1:use given reference; 2: use ANI based reference
      reffile = GCA_003013715_clade2.fasta  # reference filename 

      project = candida  # main alignment file name

      cdsSNPS = 1  # 0:no cds SNPS; 1:cds SNPs

      buildSNPdb = 0 # 0: only align to reference 1: build SNP database of all complete genome
      #aligner = bwa

    FirstTime = 1  # 1:yes; 2:update existing SNP alignment

         data = 1  # *See below 0:only complete(F); 1:only contig(C); 2:only reads(R); 
                   # 3:combination F+C; 4:combination F+R; 5:combination C+R; 
                   # 6:combination F+C+R; 7:realignment  *See below 
        #reads = 1  # 1: single reads; 2: paired reads; 3: both types present;
        #SNPSfilter = 0.0 # Only use snps that are above given threshold

         tree = 2  # 0:no tree; 1:use FastTree; 2:use RAxML; 3:use IQ tree; 4:use all
    bootstrap = 1  # 0:no; 1:yes;  # Run bootstrapping  *See below
            N = 100  # Number of bootstraps to run *See below    
  
    PosSelect = 0  # 0:No; 1:use PAML; 2:use HyPhy; 3:use both

         code = 2  # 0:Bacteria; 1:Virus; 2: Eukarya

        clean = 1  # 0:no clean; 1:clean

      threads = 20  # Number of threads to use

       cutoff = 0.1  # Linear alignment (LA) coverage against reference - ignores SNPs from organism that have lower cutoff.

* When using data option 1,2,5 need a complete reference to align/map to. 
* Use data option 7 when need to extract SNPs using a sublist of already aligned genomes. 

