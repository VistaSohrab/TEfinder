#!/bin/bash
##
##
## Authors: Vista Sohrab & Dilay Hazal Ayhan
## Date: September 16, 2020
## Description: TEfinder uses discordant reads to detect novel transposable element insertion events in paired-end sample sequencing data. 
##              Software dependencies include bedtools 2.28.0, samtools 1.3 or later, picard 2.0.1 or later
##              Required inputs include sample alignment file (.BAM), reference genome FASTA index (.FAI), reference TE annotation (.GFF),and TEs of interest (.TXT)
##
## University of Massachusetts Amherst
##
##
##
##

set -e

margs=4

# Functions
function example {
     echo -e "example: TEfinder -bam sample.bam -fai reference.fai -gff TEs.gff -te List_of_TEs.txt"
}

function help {
     echo -e "REQUIRED:"
     echo -e "  -bam,        --alignmentFile       STR  sample reads aligned to reference genome (BAM file)"
     echo -e "  -fai,        --FastaIndexFile      STR  reference genome FASTA index (FAI file)"
     echo -e "  -gff,        --TransposonsInGenome STR  reference genome TE annotation (GFF file)"
     echo -e "  -te,         --TransposonsToSearch STR  TE names\n (single column text file)"
     echo -e "OPTIONAL:"
     echo -e "  -fis,        --FragmentInsertSize  INT  short-read sequencing fragment insert size [400]"
     echo -e "  -picard,     --pathToPicardjar     STR  path to picardtools .jar file [picard.jar]"
     echo -e "  -md,         --MaxDistanceForMerge INT  maximum distance between reads for bedtools merge [150]"
     echo -e "  -k,          --MaxTSDLength	       INT  maximum TE target site duplication (TSD) length [20]"
     echo -e "  -maxHeapMem, --MaxHeapMemory       INT  java maximum heap memory allocation for picard [null]" 
     echo -e "  -workingdir, --WorkingDirectory    STR  working directory name [TEfinder_<Date>]"
     echo -e "  -h,          --help                     prints help\n"
     example
}

# check if mandatory args are empty
function margs_check {
     if [ $# -lt $margs ]; then
          echo -e "One or more required parameters are missing."
          example
          exit 1 # error
     fi
}

# main workflow
 #### : comment out
function pipeline() {
     mkdir ${workingdir}/${line}
     currdir=${workingdir}/${line}
     echo -e $(date) " Transposon analysis for "${line}" has started\n"
     
     grep -i '"Motif:'${line}'"' $gff > ${currdir}/${line}_TE.gff
     echo -e $(date) " Individual TE GFF has been created for "${line}"\n" ####
     
     bedtools intersect -abam ${workingdir}/Alignments.bam -b ${currdir}/${line}_TE.gff -wa > ${currdir}/${line}_MappedReadsToTE.bam
     echo -e $(date) " Mapped reads to TE via bedtools intersect has been completed for "${line}"\n" ####
     samtools view ${currdir}/${line}_MappedReadsToTE.bam | \
          awk -v Ins=`expr $fis \* 10` '{if (($7 != "=") || ($9 > Ins) || ($9 < -Ins)) print $1}' > ${currdir}/${line}_ReadID.txt
     echo -e $(date) " Identifying discordant read IDs has been completed for "${line}"\n" ####
     
	 # if discordant readID file exists, then continue with remainder of TE analysis
	 if  [[ -s  ${currdir}/${line}_ReadID.txt ]]
	 then
          java $maxHeapMem -jar $picard FilterSamReads I=${workingdir}/Alignments.bam O=${currdir}/${line}_DiscordantPairs.bam \
               READ_LIST_FILE=${currdir}/${line}_ReadID.txt FILTER=includeReadList WRITE_READS_FILES=false
          echo -e $(date) " Filtering original alignment based on discordant reads IDs is complete for "${line}"\n"  ####
          
          bedtools merge -d $md -S + -c 1 -o count -i ${currdir}/${line}_DiscordantPairs.bam | \
               awk '{if ($4 > 3) print $0}' > ${currdir}/${line}_plusCluster.bed 
          echo -e $(date) " Primary reads from the + strand have been merged if read count greater than 3 for "${line}"\n" ####
          
          bedtools merge -d $md -S - -c 1 -o count -i ${currdir}/${line}_DiscordantPairs.bam | \
               awk '{if ($4 > 3) print $0}' > ${currdir}/${line}_minusCluster.bed
          echo -e $(date) " Primary reads from the - strand have been merged if read count greater than 3 for "${line}"\n" ####
          
          # filtering edges piped into bedtools merge (keeping read counts greater than 3 in the line above)
          ## find the closest minus strand to the plus strand in the cluster
          ## filter by the distance between the plus and minus clusters - only retain pairs if reads are 0-100 bases away
          ## if plus strand start is less than minus strand start and plus strand end is less than minus strand end then in proper orientation 
          bedtools closest -d -g $fai -t first -a ${currdir}/${line}_plusCluster.bed -b ${currdir}/${line}_minusCluster.bed | \
               awk -v TSD=$k '{if ($9 <= TSD && $9 >= 0) print $0}' | \
               awk '{if ($2 < $6 && $3 < $7) print $0}' > ${currdir}/${line}_plusminus.bed
          echo -e $(date) " Filtration of clusters in proper orientation using bedtools closest has been completed for "${line}"\n" ####
          
          # if plus strand end is greater than minus strand start, then report the pair 
          awk '{if ($3 > $6) print $1"\t"$6"\t"$3"\t"$0}' ${currdir}/${line}_plusminus.bed > ${currdir}/${line}_plusminus_1.bed 
          echo -e $(date) " Overlapping reads TE insertions reported for "${line}"\n" ####
          
          #if plus  strand end is less than or equal to minus strand start and the region in between is less than a user-defined value k, report the pair
          awk -v TSD=$k '{if ($3 <= $6 && $6 - $3 < TSD) print $1"\t"$3 - 1"\t"$6 + 1"\t"$0}' ${currdir}/${line}_plusminus.bed  > \
               ${currdir}/${line}_plusminus_2.bed
          echo -e $(date) " Non-overlapping reads TE insertions reported for "${line}"\n" ####
          
          #combine reported TE insertions
          cat ${currdir}/${line}_plusminus_1.bed ${currdir}/${line}_plusminus_2.bed | \
               awk -v TEname=$line '{$0=TEname"\t"$0}1' | sort -k 1 | sort -k 2 > \
               ${currdir}/${line}_insertionRegion.txt
          echo -e $(date) " TE insertions of "${line}" not present in repeat regions of reference genome have been reported.\n" ####
          
          cat ${currdir}/${line}_insertionRegion.txt >> ${workingdir}/insertions.txt
          echo -e $(date) " TE insertions for "${line}" have been reported.\n" ####
          
          #remove intermediate files
          rm ${currdir}/${line}_MappedReadsToTE.bam ${currdir}/${line}_ReadID.txt 
          echo -e $(date) " Transposon named "${line}" is processed.\n"
	 else
	      echo -e $(date) " Transposon named "${line}" is processed. No discordant reads found.\n"
		  rm -r ${currdir}
	 fi
}

# functions end

# get arguments
fai=
bam=
gff=
te=
maxHeapMem=
fis=400
picard="picard.jar"
md=150
k=20
d=$(date +%Y%m%d%H%M%S)
workingdir=TEfinder_${d}
 
while [ "$1" != "" ];
do
     case $1 in
     -fai | --FastaIndexFile )
          shift
          fai=$1 ;;
     -bam | --alignmentFile )
	      shift
          bam=$1 ;;
     -gff | --TransposonsInGenome )
	      shift
          gff=$1 ;;
     -te | --TransposonsToSearch )
	      shift
          te=$1  ;;
     -fis | --FragmentInsertSize )
	      shift
          fis=$1 ;;
     -picard | --pathToPicardjar )
	      shift
          picard=$1 ;;
     -md | --MaxDistanceForMerge )
	      shift
          md=$1  ;;
     -k | --MaxTSDLength )
	      shift
          k=$1  ;;
     -maxHeapMem | --MaxHeapMemory )
	      shift
          maxHeapMem="-Xmx"$1"m" ;;
     -workingdir | --WorkingDirectory )
	      shift
          workingdir=$1 ;;
     -h   | --help )
	      help
          exit;;
     *) # error
          echo "TEfinder: illegal option $1"
          example
          exit 1 ;;
     esac
     shift
done
margs_check $fai $bam $gff $te

# main

mkdir ${workingdir}

# remove empty lines from user provided TE list if present
sed '/^$/d' ${te} > ${workingdir}"/userTE_noEmptyLines.txt"

# create output files
touch ${workingdir}/TEinsertions.bed
printf "%s\t"  "track name=TEfinder" "type=bedDetail" "description=FR:forward read, RR:reverse read, IS:insertion region start, IE:insertion region end" > ${workingdir}/TEinsertions.bed
printf "%s\n" >> ${workingdir}/TEinsertions.bed

touch ${workingdir}/TEinsertions_inRepeatRegions.bed
printf "%s\t"  "track name=TEfinder" "type=bedDetail" "description=FR:forward read, RR:reverse read, IS:insertion region start, IE:insertion region end" > ${workingdir}/TEinsertions_inRepeatRegions.bed
printf "%s\n" >> ${workingdir}/TEinsertions_inRepeatRegions.bed

samtools view -F 2304 -o ${workingdir}/Alignments.bam ${bam} 
echo -e $(date) " Alignments are filtered - secondary and supplementary alignments have been removed. \n"

# run pipeline for each TE
while IFS="" read -r line || [ -n "$line" ]
do
     pipeline &
done < ${workingdir}/userTE_noEmptyLines.txt
wait
echo -e $(date) " All transposons are processed. Finalizing...\n"

# combine bam files
samtools merge -r ${workingdir}/DiscordantReads.bam ${workingdir}/*/*_DiscordantPairs.bam

# convert files
awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$8+$12"\t.\t"$6"\t"$11"\tFR="$8";RR="$12";IS="$6";IE="$11}' ${workingdir}/insertions.txt >> ${workingdir}/TEinsertions_putative.bed
bedtools intersect -v -a ${workingdir}/TEinsertions_putative.bed -b $gff > ${workingdir}/TEinsertions_temp.bed
bedtools intersect -a ${workingdir}/TEinsertions_putative.bed -b $gff -u > ${workingdir}/TEinsertions_inRepeatRegions_temp.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9}' ${workingdir}/TEinsertions_temp.bed >> ${workingdir}/TEinsertions.bed
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$9}' ${workingdir}/TEinsertions_inRepeatRegions_temp.bed >> ${workingdir}/TEinsertions_inRepeatRegions.bed
rm ${workingdir}/TEinsertions_putative.bed ${workingdir}/TEinsertions_inRepeatRegions_temp.bed ${workingdir}/TEinsertions_temp.bed

echo -e $(date) " TE insertion BED files have been created. TEfinder completed successfully."
