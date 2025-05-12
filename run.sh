#!/bin/bash

# Load required module
module load matlab

# Define arrays
gammas=(0.10 0.20)
ms=(0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9)
rs=(0.15)
thetas=(5)
alphas=(0.5 0.6 0.7 0.8 0.9)

# Iterate over all combinations
for gamma in "${gammas[@]}"; do
  for m in "${ms[@]}"; do
    for r in "${rs[@]}"; do
      for theta in "${thetas[@]}"; do
        for alpha in "${alphas[@]}"; do

          jobname="g${gamma}_m${m}_r${r}_t${theta}_a${alpha}"
          args="$gamma $m $r $theta $alpha"

          sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=$jobname
#SBATCH --output=./job-outs/$jobname.out
#SBATCH --error=./job-outs/$jobname.err
#SBATCH --account=pi-lhansen
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-04:00:00

module load matlab
export MATLAB_ARGS="$args"
matlab -nodisplay -r "main_wrapper; exit;"
EOF

        done
      done
    done
  done
done
