#!/bin/bash
#R. Kell

module load Trimmomatic/v0.39

fileList="SP1_C0_S29 SP1_C1_S35 SP1_C24_S53 SP1_C3_S41 SP1_C6_S47 SP1_Z0_S30 SP1_Z1_S36 SP1_Z24_S54 SP1_Z6_S48 SP3_C0_S31 SP3_C1_S37 SP3_C24_S55 SP3_C3_S43 SP3_C6_S49 SP3_Z0_S32 SP3_Z1_S38 SP3_Z24_S56 SP3_Z3_S44 SP3_Z6_S50 SP4_C0_S33 SP4_C1_S39 SP4_C24_S57 SP4_C3_S45 SP4_C6_S51 SP4_Z0_S34 SP4_Z1_S40 SP4_Z24_S58 SP4_Z3_S46 SP4_Z6_S52”


for file in ${fileList}
do

java -jar /data/resources/app_modules/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 -phred33 \
	${file}_R1_001.fastq.gz \
	${file}_R2_001.fastq.gz \
  ${file}_R1-qc.fq ${file}_R1-unpaired.fq \
  ${file}_R2-qc.fq ${file}_R2-unpaired.fq \
  ILLUMINACLIP:/data/resources/app_modules/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 \
  2>&1 | tee -a trimmomatic.log

done
