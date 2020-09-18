# TEfinder
A bioinformatics tool for detecting novel transposable element insertions

Authors: Vista Sohrab & Dilay Hazal Ayhan 

TEfinder uses discordant reads to detect novel transposable element insertion events in paired-end sample sequencing data. 

Software dependencies:
* Bedtools 2.28.0
* Samtools 1.3 or later
* Picard 2.0.1 or later

Required inputs:
* Sample Alignment File (.BAM)
* Reference Genome FASTA Index (.FAI)
* Reference TE Annotation (.GFF)
* TEs of interest (.TXT)

Usage:

TEfinder -bam sample.bam -fai reference.fai -gff TEs.gff -te List_of_TEs.txt

Output files:
* TE_insertions.bed contains novel TE insertion events in sample 
* TE_insertions_inRepeatRegions.bed contains TE insertions in known repeat regions based on TE annotation of the reference genome
* DiscordantReads.bam contains all discordant reads that have been identified based on the TEs of interest that have been submitted to TEfinder
