#!/bin/bash
#Make the S. purpuratus 5.0 genome into a hisat2 index
#R. Kell

module load histat2/v2.2.1

hisat2-build GCF_000002235.5_Spur_5.0_genomic.fna Spur5.0_hisat2index
