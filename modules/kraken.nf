process kraken {
    input:
        val x
    output:
        //stdout
        val x
        
    """   
    mkdir -p ${params.output}/${x}/kraken_out/
    
    kraken2 --db ${params.db} --use-names --report ${params.output}/${x}/kraken_out/${x}.report --output ${params.output}/${x}/kraken_out/${x}_kraken.out ${params.output}/${x}/${x}.fastq.gz

    """
}