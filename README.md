# TEfinder
A bioinformatics tool for detecting novel transposable element insertions

Authors: Vista Sohrab & Dilay Hazal Ayhan 

TEfinder uses discordant reads to detect novel transposable element insertion events in paired-end sample sequencing data. 

Software dependencies:
* Bedtools 2.28.0 or later
* Samtools 1.3 or later
* Picard 2.0.1 or later

Required inputs:
* Sample alignment file(.bam | .sam)
* Reference genome FASTA (.fa)
* Reference TE annotation (.gff | .gtf)
* TEs of interest in a single column text file (.txt)

Optional arguments with default values in brackets: 
* -fis:        short-read sequencing fragment insert size [400]
* -picard:     path to Picard Tools .jar file [picard.jar]
* -md:         maximum distance between reads for bedtools merge [150]
* -k:          maximum TE target site duplication (TSD) length [20]
* -maxHeapMem: java maximum heap memory allocation for picard in Mb [2000]
* -workingdir: working directory name [TEfinder_Date]
* -out:        output format as GTF [BED]
* -outname:    output name prefix added to file names [null]
* -threads:    number of threads for samtools multi-threading [1]
* -h:          prints help option
 
Usage:

TEfinder -alignment sample.bam -fa reference.fa -gtf TEs.gtf -te List_of_TEs.txt

Output files:
* TE_insertions.bed contains identified TE insertion events in sample (in the final column, FILTER attribute with "PASS" refers to high confidence insertion events while instances labeled as "in_repeat", "weak_evidence", "strand bias" or a combination of these three labels indicate less confident insertion events)
* TE_insertions.gtf is provided with the same information as the BED file if using -out GTF
* DiscordantReads.bam contains all discordant reads that have been identified based on the TEs of interest that have been submitted to TEfinder
