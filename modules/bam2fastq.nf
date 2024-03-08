process bam2fastq {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   mkdir -p ${params.output}/${x}
   #cp ${params.input}/${x}_*.fastq.gz ${params.output}/${x}
   cp ${params.input}/${x}.bam ${params.output}/${x}
   cp ${params.input}/${x}.bam.pbi ${params.output}/${x}
   
   #QC
   python ${params.qc}/longQC.py sampleqc -x pb-sequel -o ${params.output}/${x}/qc ${params.output}/${x}/${x}.bam
   #fastqc ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz

   #bam to fastq
   bam2fastq -o ${params.output}/${x}/${x} ${params.output}/${x}/${x}.bam

   #mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_original_fastqc.html
   #mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_original_fastqc.zip
   #mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_original_fastqc.html
   #mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_original_fastqc.zip
   
   """
}
