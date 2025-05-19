#!/bin/bash


# Define arrays
gammas=(0.10 0.20)
# ms=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85)
# ms=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
ms=(0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85)
# rs=(0.15)
rs=(0.10)
# rs=(0.15 0.30)
phi_hs=(1.05)
thetas=(5)
alphas=(0.8 0.9)
betas=(2 8)
epsilons=(0 1)

mkdir -p job-outs
# Iterate over all combinations
for gamma in "${gammas[@]}"; do
  for m in "${ms[@]}"; do
    for r in "${rs[@]}"; do
      for theta in "${thetas[@]}"; do
        for alpha in "${alphas[@]}"; do
            for phi_h in "${phi_hs[@]}"; do
                for beta in "${betas[@]}"; do
                    for epsilon in "${epsilons[@]}"; do

                        jobname="g${gamma}_m${m}_r${r}_t${theta}_a${alpha}_ph_${phi_h}_b${beta}_e${epsilon}"
                        args="$gamma $m $r $theta $alpha $phi_h $beta $epsilon"

          sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=$jobname
#SBATCH --output=./job-outs/$jobname.out
#SBATCH --error=./job-outs/$jobname.err
#SBATCH --account=pi-lhansen
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-04:00:00
#SBATCH --exclude=mcn08,mcn10,mcn12,mcn15,mcn52,mcn53,mcn57,mcn67,mcn68,mcn69

module load matlab/R2024b
export MATLAB_ARGS="$args"
matlab -nodisplay -r "main_wrapper; exit;"
EOF

                    done
                done
            done
        done
      done
    done
  done
done
