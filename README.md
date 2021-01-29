# TEfinder
A bioinformatics tool for detecting novel transposable element insertions

Authors: Vista Sohrab & Dilay Hazal Ayhan 

TEfinder uses discordant reads to detect novel transposable element insertion events in paired-end sample sequencing data. 

**Software dependencies:**
* Bedtools 2.28.0 or later
* Samtools 1.3 or later
* Picard 2.0.1 or later

**Required inputs:**
* Sample alignment file(.bam | .sam) 
* Reference genome FASTA (.fa)
* Reference TE annotation (.gff | .gtf)
* TEs of interest in a single column text file (.txt)

**Optional arguments with default values in brackets:**
* -fis:        short-read sequencing fragment insert size [400]
* -picard:     path to Picard Tools .jar file [picard.jar]
* -md:         maximum distance between reads for bedtools merge [150]
* -k:          maximum TE target site duplication (TSD) length [20]
* -maxHeapMem: java maximum heap memory allocation for picard in MB [2000]
* -workingdir: working directory name [TEfinder_Date]
* -out:        output format as GTF [BED]
* -outname:    output name prefix added to file names [null]
* -threads:    number of threads for samtools multi-threading [1]
* -intermed:   keep intermediate files created by pipeline [no]
* -h:          prints help option
 
**Usage:**

Default run:
```
TEfinder -alignment sample.bam -fa reference.fa -gtf TEs.gtf -te List_of_TEs.txt
```
Note: In the test_dataset folder, the alignment file used for testing is gzipped so please do the following prior to running TEfinder: 

      ```
      gunzip sample.bam.gz
      ```

      The expected test run result is available as sample_TEinsertions.bed for comparison.

Example command to change fragment insert size to 500, set maximum TSD length to 30, create GTF output, and include all intermediate files:
```
TEfinder -alignment sample.bam -fa reference.fa -gtf TEs.gtf -te List_of_TEs.txt \
         -fis 500 -k 30 -out GTF -intermed yes
```
**Output files:**
* TE_insertions.bed contains identified TE insertion events in sample (in the final column, FILTER attribute with "PASS" refers to high confidence insertion events while instances labeled as "in_repeat", "weak_evidence", "strand bias" or a combination of these three labels indicate less confident insertion events)
* TE_insertions.gtf is provided with the same information as the BED file if using -out GTF
* DiscordantReads.bam contains all discordant reads that have been identified based on the TEs of interest that have been submitted to TEfinder

**Notes:**
* The sample alignment file must include read group information in the header. An example command to include read group if using BWA for alignment would be:
  
  ```
  bwa index reference.fasta
  bwa mem -R ‘@RG\tID:sample\tPL:illumina\tLB:LIB\tSM:sample’ \
          reference.fasta sample_R1.fq sample_R2.fq > sample.bam
  ```

  If alignment has already been completed without specifying read group information, Picard's AddOrReplaceReadGroups function can generate read group information provided
  the binary alignment file (BAM):
  ```
  java -jar picard.jar AddOrReplaceReadGroups \
      I=input.bam \
      O=sample.bam \
      RGID=sample \
      RGLB=LIB \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=sample
   ```
* Please ensure that the TE names of interest are separated by newline character "\n" and not carriage return "\r". 
  If carriage return is present, please replace all instances of "\r\n" with "\n".
  *This issue is generally encountered by Windows users.*

* To remove simple repeats from the RepeatMasker TE annotation output file, an example command would be the following: 
  ```
  grep -v -iE '(Motif\:[ATGC]+\-rich)|(Motif\:\([ATGC]+\)n)' RepeatMasker_output.gff > TEs.gtf
  ```
* To create a list of TE names from the TE GTF file: 

  ```
   awk -F '\t' '{print $9}' TEs.gtf | awk -F '"' '{print $2}' | sort | uniq > List_of_TEs.txt
  ```
  *This command may not yield the desired format of TE names as it depends on the structure of the attribute column of the GTF file and would need to be modified for specific cases.* 

* Modifying the maximum TSD length (-k) could be useful if there is an unexpected number of insertion events identified with the default parameter.
  The optimal maximum TSD length can vary across datasets.

* Modifying the fragment insert size (-fis) based on the sequencing library preparation can be useful.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4446971.svg)](https://doi.org/10.5281/zenodo.4446971)
