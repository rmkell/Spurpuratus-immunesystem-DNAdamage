#!/bin/bash
#Map trimmed RNA-seq reads to the S. purpuratus 5.0 genome using HISAT2. Paired reads (R1-qc.fq, R2-qc.fq) in addition to unpaired R1 reads (R1-unpaired.fq) were used as input.
#R. Kell

module load histat2/v2.2.1
module load samtools/v1.18

fileList="SP1_C0_S29 SP1_C1_S35 SP1_C24_S53 SP1_C3_S41 SP1_C6_S47 SP1_Z0_S30 SP1_Z1_S36 SP1_Z24_S54 SP1_Z6_S48 SP3_C0_S31 SP3_C1_S37 SP3_C24_S55 SP3_C3_S43 SP3_C6_S49 SP3_Z0_S32 SP3_Z1_S38 SP3_Z24_S56 SP3_Z3_S44 SP3_Z6_S50 SP4_C0_S33 SP4_C1_S39 SP4_C24_S57 SP4_C3_S45 SP4_C6_S51 SP4_Z0_S34 SP4_Z1_S40 SP4_Z24_S58 SP4_Z3_S46 SP4_Z6_S52"

for file in ${fileList}
do

exec &> hisat2Output/${file}_hisat2.log
hisat2 -p 8 -x /data/prj/urchin/Spurp10XSeq/BulkRNASeq/hisat2Index/Spur5.0_hisat2index \
	-1 ${file}_R1-qc.fq \
	-2 ${file}_R2-qc.fq \
	-U ${file}_R1-unpaired.fq | \
	samtools view -@ 8 -S -h -u - | \
	samtools sort -@ 8 -T ${file} - > hisat2Output/${file}.bam

done
