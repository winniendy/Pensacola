for f in ./*.bam*; do
    mv "$f" "${f//*.hifi_reads.bc/bc}"
done

for f in ./*.bam*; do
    mv "$f" "${f//--/}"
done


