#Using cellbenderv0.3.0
#Running on Northeastern University HPC Discovery Cluster
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:p100:1
#SBATCH --job-name=LOOP_allfiles_pgu_p100               # Job name
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --mem=4GB
#SBATCH --output=output_%j.txt       # Standard output file
#SBATCH --error=error_%j.txt        # Standard error file

#Loop through .h5files and apply cellbender
conda init
conda activate cellbender

for file in /home/r.kellogg/RawH5Files/*.h5; do cellbender remove-background --input "$file" --output "$(basename "$file" .h5)denoised.h5" --total-droplets-included 50000
 --expected-cells 10000 --fpr 0.01 --epochs 150 --cuda; done
