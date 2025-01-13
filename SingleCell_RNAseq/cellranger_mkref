#!/bin/bash
# 10x cellranger mkgtf and mkref
# Author: Riss Kell
# 10X tutorial: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr

#Filter GTF file using mkgtf
#cellranger mkgtf \
    #Strongylocentrotus_purpuratus.Spur_5.0.58.gtf \
    #Strongylocentrotus_purpuratus.Spur_5.0.58.filtered.gtf \
    #--attribute=gene_biotype:protein_coding
#This will output the file Strongylocentrotus_purpuratus.Spur_5.0.58.filtered.gtf, which will be used in the next step.

#Make the S.purpuratus v5.0 genome into a 10x index using mkref
#cellranger mkref \
#    --genome=Spurpv5_refgenome \
#    --fasta=Strongylocentrotus_purpuratus.Spur_5.0.dna.toplevel.fa \
#    --genes=Strongylocentrotus_purpuratus.Spur_5.0.58.filtered.gtf \



#make v3.0 L variegatus genome into a 10x index - redo after adding gene_name attribute and filtering GTF file
/data/app/cellranger-6.1.2/cellranger mkref \
--genome=Lv_genome_filteredGenes_addGeneName_10xindx \
--fasta=GCF_018143015.1_Lvar_3.0_genomic.fna \
--gene=GCF_018143015.1_Lvar_3.0_genomic_addGeneName.filter.gtf \
2>&1 | tee -a cellranger_mkref_filterGenes_addGeneName__log.txt
