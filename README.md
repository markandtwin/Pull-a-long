# Pull-a-long

Python3, pysam, samtools and minimap2 are required.

## Step 1 Aligne reads to reference genome

minimap2 -ax splice -t 7 -B 3 -O 3,20 --junc-bed  /path to transcriptome bed file /path to reference genome  /path to Nanopore data    > output.sam

samtools sort output.sam > output.bam

samtools index output.bam

## Step 2 Filter out reads covering required regions

python  exon_covering.py -i output.bam -b /path to for file including upstream exon coordinates in bed4 format -o up.bam -t 7
	  
python  exon_covering.py -i up.bam -b /path to for file including downstream exon coordinates in bed4 format -o down.bam -t 7 

python  exon_covering.py -i down.bam -b /path to for file including universal 3'UTR coordinates in bed4 format -o UTR.bam -t 7

## Step 3 Parse reads covering extended 3'UTR regions to long 3'UTR isoforms

python  exon_covering.py -i UTR.bam -b /path to for file including extended 3'UTR coordinates in bed4 format -o long.bam -t 7 

## Step 3 Parse reads that do not cover extended 3'UTR regions and have correct 3' end to short 3'UTR isoforms  

samtools view -h -b UTR.bam -L /path to for file including extended 3'UTR coordinates in bed4 format  -U short_1.bam > other.bam
	  
samtools index short_1.bam

python  polyA_filtering.py -i  short_1.bam  -o short.bam -t 7

## Step 4 Calculate PSI of cassette exon for each isoform

python calculate_PSI.py -i long.bam -b /path to for file including cassette exon coordinates in bed4 format -t 7  -o long_PSI.csv

python calculate_PSI.py -i short.bam -b /path to for file including cassette exon coordinates in bed4 format -t 7  -o short_PSI.csv
