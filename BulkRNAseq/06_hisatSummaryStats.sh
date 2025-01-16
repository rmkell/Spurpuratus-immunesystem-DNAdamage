module load samtools/v1.18

fileList="*.bam"

for file in ${fileList}
do

	prefix=$(echo ${file} | cut -d "." -f 1)
	samtools stats ${file} > ${prefix}.stats

done
