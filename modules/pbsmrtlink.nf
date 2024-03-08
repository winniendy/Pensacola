process pbsmrtlink {
   input:
      val x
   output:
      //stdout
      //path 'xfile.txt', emit: aLook
      val "${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """  
   
   mkdir -p ${params.output}/variants
   #mkdir -p ${params.output}/${x}
   #cp ${params.input}/${x}.bam ${params.output}/${x}
   #cp ${params.input}/${x}.bam.pbi ${params.output}/${x}

   
   #align reads to reference
   pbmm2 align --sort --preset HiFi --log-level INFO ${params.reference}/*.fasta ${params.output}/${x}/${x}.bam ${params.output}/${x}/${x}.mapped.bam

   """
}
