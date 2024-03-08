process snpeff {
   input:
      val x
   output:
      //stdout
      //path 'xfile.txt', emit: aLook
      val "${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """  
   #call variants
   #singularity exec docker://staphb/bcftools:1.16 bcftools mpileup --op en-prob25 --indel-size 450 --gap-frac 0.01 --ext-prob 1 --min-ireads 3 --max-depth 1000 --max-idepth 5000 --seed 1984 -h 500 -B -a FORMAT/AD -f ${params.reference}/*.fasta ${params.output}/${x}/${x}.mapped.bam | singularity exec docker://staphb/bcftools:1.16 bcftools call -mv -Ov | singularity exec docker://staphb/bcftools:1.16 bcftools norm -f ${params.reference}/*.fasta - | singularity exec docker://staphb/bcftools:1.16 bcftools filter -e 'QUAL < 20' - > ${params.output}/${x}/${x}.variants_bcftools.vcf
   #bcftools mpileup --op en-prob25 --indel-size 450 --gap-frac 0.01 --ext-prob 1 --min-ireads 3 --max-depth 1000 --max-idepth 5000 --seed 1984 -h 500 -B -a FORMAT/AD -f ${params.reference}/*.fasta ${params.output}/${x}/${x}.mapped.bam | bcftools call -mv -Ov | bcftools norm -f ${params.reference}/*.fasta - | bcftools filter -e 'QUAL < 20' - > ${params.output}/${x}/${x}.variants_bcftools.vcf
   
   #annotate snp
   #singularity exec -B ./:/snpEff/data docker://staphb/snpeff:5.2a snpEff -c ${params.input}/snpEff.config -v C_auris_B8441 ${params.output}/${x}/${x}.variants_bcftools.vcf > ${params.output}/${x}/${x}.variants_bcftools.ann.vcf
   snpEff -c ${params.snpeffconfig}/snpEff.config -v C_auris_B8441 -stats ${params.output}/${x}/${x}_snpeff.html  ${params.output}/${x}/${x}.variants_bcftools.vcf > ${params.output}/${x}/${x}.variants_snpeff.ann.vcf

   #copy vcf files to variants folder
   cp ${params.output}/${x}/${x}.variants_bcftools.vcf ${params.output}/variants
   cp ${params.output}/${x}/${x}.variants_snpeff.ann.vcf ${params.output}/variants
   cp ${params.output}/${x}/${x}_snpeff.html ${params.output}/variants
   cp ${params.output}/${x}/${x}_snpeff.genes.txt ${params.output}/variants

   """
}
