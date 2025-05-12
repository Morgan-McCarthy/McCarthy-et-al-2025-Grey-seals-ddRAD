module load bwa 
module load samtools
samtools faidx Halichoerus_grypus_HiC_RM.fasta
bwa index Halichoerus_grypus_HiC_RM.fasta

for file in *.fq

do

shortID=`echo ${file} | cut -f1 -d "."`

nuref=Halichoerus_grypus_HiC_RM.fasta

#align to full genome
bwa mem -t 28 ${nuref} ${file} > ${shortID}_nuDNA_aligned.sam

#convert to bam file
samtools view -@ 28 -T ${nuref} -b ${shortID}_nuDNA_aligned.sam | samtools sort -@ 28 -o ${shortID}_nuDNA_aligned_sorted.bam

done
