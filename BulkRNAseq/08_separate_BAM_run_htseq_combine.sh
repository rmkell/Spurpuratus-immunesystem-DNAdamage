module load htseq/v2.0.4

#Run htseq-count on paired-end reads and on single-end reads separately, then combine outputs into one counts file

fileList="SP1_C0_S29 SP1_C1_S35 SP1_C24_S53 SP1_C3_S41 SP1_C6_S47 SP1_Z0_S30 SP1_Z1_S36 SP1_Z24_S54 SP1_Z6_S48 SP3_C0_S31 SP3_C1_S37 SP3_C24_S55 SP3_C3_S43 SP3_C6_S49 SP3_Z0_S32 SP3_Z1_S38 SP3_Z24_S56 SP3_Z3_S44 SP3_Z6_S50 SP4_C0_S33 SP4_C1_S39 SP4_C24_S57 SP4_C3_S45 SP4_C6_S51 SP4_Z0_S34 SP4_Z1_S40 SP4_Z24_S58 SP4_Z3_S46 SP4_Z6_S52"

INDIR="/data/prj/urchin/Spurp10XSeq/BulkRNASeq/11_19_24_BulkRNASequencingRAWFastqs/trimmomatic_trimmed/hisat2Output"

OUTDIR="/data/prj/urchin/Spurp10XSeq/BulkRNASeq/11_19_24_BulkRNASequencingRAWFastqs/trimmomatic_trimmed/hisat2Output/htseq_counts/"

GTF="GCF_000002235.5_Spur_5.0_genomic.gtf.gz"

for file in ${fileList}
do
    # Separate paired-end and single-end reads
    samtools view -b -f 1 ${INDIR}/${file}.bam > ${INDIR}/${file}.paired.bam
    samtools view -b -F 1 ${INDIR}/${file}.bam > ${INDIR}/${file}.single.bam

    # Run htseq-count on paired-end reads
    htseq-count -f bam -r pos -s no -t gene -i gene_id \
        ${INDIR}/${file}.paired.bam \
        ${GTF} > ${OUTDIR}/${file}.paired.counts

    # Run htseq-count on single-end reads
    htseq-count -f bam -r pos -s no -t gene -i gene_id \
        ${INDIR}/${file}.single.bam \
        ${GTF} > ${OUTDIR}/${file}.single.counts

    # Combine the counts
    paste ${OUTDIR}/${file}.paired.counts ${OUTDIR}/${file}.single.counts | \
    awk '{print $1, $2 + $4}' > ${OUTDIR}/${file}.combined.counts
done
