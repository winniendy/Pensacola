process bracken {
    input:
        val x
    output:
        //stdout
        val x
        
    """   
    mkdir -p ${params.output}/species_abundance
    
    #kraken2 --db ${params.db} --use-names --report ${params.output}/${x}/kraken_out/${x}.report --output ${params.output}/${x}/kraken_out/${x}_kraken.out ${params.output}/${x}/${x}.fastq.gz
    bracken -d ${params.db} -i ${params.output}/${x}/kraken_out/${x}.report -o ${params.output}/${x}/kraken_out/${x}_bracken.species_abundance
    cp ${params.output}/${x}/kraken_out/${x}_bracken.species_abundance ${params.output}/species_abundance/
    """
}