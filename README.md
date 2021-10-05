** Pull-a-long**
Python3, pysam, samtools and minimap2 are required.
Step 1 Aligne reads to reference genome
minimap2 -ax splice -t 7 -B 3 -O 3,20 --junc-bed  /path to transcriptome bed file /path to reference genome  /path to Nanopore data    > output.sam
samtools sort output.sam > output.bam
samtools index output.bam
