process statistics {
   input:
      val x
   output:
      //stdout
      //path 'xfile.txt', emit: aLook
      val "${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """  
   #calculate statistics
   samtools coverage ${params.output}/${x}/${x}.mapped.bam -o ${params.output}/${x}/${x}.coverage.txt
   samtools mpileup --min-BQ 1 -f ${params.reference}/*.fasta -s -o ${params.output}/${x}/${x}.mapped.bam.mpileup ${params.output}/${x}/${x}.mapped.bam
   samtools depth -q 0 -Q 0 -o ${params.output}/${x}/${x}.mapped.bam.depth ${params.output}/${x}/${x}.mapped.bam

   
   """
}
