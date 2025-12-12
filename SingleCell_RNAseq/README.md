# *S. purpuratus* immune cell (coelomocyte) single-cell analysis: comparison of n=3 UVB treated (1000 mJ/cm2) and n=3 control biological replicates

# 1. Convert *S. purpuratus* genome into a cellranger reference
  * Complete 10X tutorial here: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr
  * Program: cellranger v7.1.0
  * Genome: S. purpuratus v5.0 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002235.5/)
  * Script: cellranger_mkref.sh

# 2. Count reads using cellranger 'count' in 10x Cloud Analysis platform
  * Complete 10X tutorial here: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct
  * This was run on each sample individually
  * Program: cellranger v7.1.0
  * --expect cells = 10,000
  * 'Count' was performed using 10X Cloud Analysis.
  * Output Feature / cell matrix HDF5 (raw) files were downloaded from the 10X Cloud Analysis project and used as CellBender input

# 3. Remove background RNA using CellBender
  * Program: cellbender v.0.3.0 (https://github.com/broadinstitute/CellBender)
  * Script: run_cellbender.sh

# 4. Clustering and Differential Expression analysis in Python
  * Programs: SCVI tools v1.1.2 (https://scvi-tools.org/), Scanpy v1.10.1 (https://scanpy.readthedocs.io/en/stable/)
  * Script: scRNAseq-DataAnalysisPipeline.ipynb
