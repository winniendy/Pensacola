process assemble {
    input:
        val x
    output:
        //stdout
        val x
        
    """    
    mkdir -p ${params.output}/${x}/assemble/
    hifiasm -o ${params.output}/${x}/assemble/${x} -t 15 ${params.output}/${x}/${x}.fastq.gz
    """
}