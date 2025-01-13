# *S. purpuratus* immune cell (coelomocyte) single-cell analysis

# 1. Convert *S. purpuratus* genome into a cellranger reference
  * Complete 10X tutorial here: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr
  * Program: cellranger-6.1.2
  * Genome: S. purpuratus v5.0 (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002235.5/)
  * Script: cellranger_mkref.sh

# 2. Count reads using cellranger count
  * Complete 10X tutorial here: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct
  * This was run on each sample individually
  * Program: cellranger-6.1.2
  * --expect cells = 10,000
  * Script: cellranger_count.sh

# 3. Remove background RNA using CellBender
  * Program: cellbender
  * Script:

# 4. Clustering and Differential Expression analysis
  * Programs: SCVI tools, Scanpy
  * Script:
