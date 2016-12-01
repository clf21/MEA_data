#!/bin/bash

if [ -z "$1" ]; then
		echo -e "\n  usage: $0 fastq_directory genome_build \n"
	exit
fi

WD=$1
Build=$2
echo -e "Genome build set to: $Build"
echo -e "Processing .fastq files found in: $WD \n"


cd $WD || { echo -e "ERROR: Could not find $WD, please check if directory is correct! \n" ; exit 1; }

# check if sample had already been processed
if [ -d "./tophat_out" ]; then
	echo -e "ERROR: Looks like this sample has already been aligned, exiting... \n" ; exit 1
fi

mkdir raw || { echo -e "ERROR: Could not write in $WD, please check permissions on this directory! \n" ; exit 1; }
FASTQS=$(ls -m ./*.fastq.gz | tr -d '\n' | tr -d ' ') # variable to store all fastq files in comma-separated list
echo -e "These are the .fastq files to be aligned: $FASTQS \n"

echo -e "Producing QC report on raw fastq files... \n"
for f in ./*.fastq.gz ; do fastqc --threads 6 --noextract $f ; done
echo -e "Inspect raw read QC report in raw directory after run \n"

echo -e "Beginning alignment to reference genome... \n"
# this alignment strategy runs Tophat2 which calls bowtie2 engine to align single-end reads to reference transcriptome first, then aligns to rest of genome with reads that fail to align to reference transcriptome
# this alignment command has to be modified if we have paired-end data
tophat2 -p 6 --transcriptome-index ~/programs/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/transcriptome_data/refSeq ~/programs/Rattus_norvegicus/UCSC/rn6/Sequence/Bowtie2Index/genome $FASTQS
mv ./*fastq* ./raw/

echo -e "Indexing alignment... \n"
samtools index ./tophat_out/accepted_hits.bam

echo -e "Generating scaled bedGraph coverage... \n"
# First get total number of succesfully mapped reads from tophat2 log file
NUMRDS=$(sed '3q;d' ./tophat_out/align_summary.txt | sed 's/^.*:  //; s/(.*$//')
NRM=$(echo "scale=6; $NUMRDS/1000000" | bc) # divide by 1 million
SF=$(echo "scale=6; 1/$NRM" | bc) # take inverse as scaling factor (this value gets multiplied by raw count to get coverage per million mapped reads
echo -e "Total number of mapped reads is $NUMRDS"
echo -e "Scaling factor for genome coverage is therefore $SF \n"

# The -split command does not count coverage of exon junction spanning reads to avoid false intron coverage counts
bedtools genomecov -ibam ./tophat_out/accepted_hits.bam -bg -scale "$SF" -split > ./tophat_out/accepted_hits.bedGraph
echo -e "Use igv.sh to view alignment and genome coverage with .bam and .bedGraph files.\n"

echo -e "Annotating reads per gene for UCSC knownGenes and ENSEMBL genes with featureCounts... \n"
featureCounts -T 6 -t exon -g gene_id -a ~/programs/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf -o geneCounts_UCSC.txt ./tophat_out/accepted_hits.bam
featureCounts -T 6 -t exon -g gene_id -a ~/programs/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/Rattus_norvegicus.Rnor_6.0.85.gtf -o geneCounts_ENSEMBL.txt ./tophat_out/accepted_hits.bam
mkdir gene_counts ; mv ./geneCounts* ./gene_counts/
 

echo -e "\nInitial RNA-seq processing complete ; Check above for any error messages.\n"
