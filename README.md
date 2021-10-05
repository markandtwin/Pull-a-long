# Pull-a-long

Python3, pysam, samtools and minimap2 are required.

## Step 1 Aligne reads to reference genome

minimap2 -ax splice -t 7 -B 3 -O 3,20 --junc-bed  /path to transcriptome bed file /path to reference genome  /path to Nanopore data    > output.sam

samtools sort output.sam > output.bam

samtools index output.bam

## Step 2 Filter out reads covering required regions

python  exon_covering.py -i output.bam -b /pathto for file including upstream exon coordinates in bed4 format -o up.bam -t 7
	  
python  exon_covering.py -i up.bam -b /pathto for file including downstream exon coordinates in bed4 format -o down.bam -t 7 

python  exon_covering.py -i down.bam -b /pathto for file including universal 3'UTR coordinates in bed4 format -o UTR.bam -t 7
