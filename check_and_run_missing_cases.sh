#!/bin/bash

# --------------------------
# Define full parameter grid
# --------------------------
gammas_all=(0.10 0.20)
ms_all=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85)
rs_all=(0.10 0.15)
phi_hs_all=(1.05)
thetas_all=(5)
alphas_all=(0.8 0.9)
betas_all=(2 8)
epsilons_all=(0 1)

sigma=3
data_dir="./data_beta"
mkdir -p job-outs

# --------------------------
# Initialize missing arrays
# --------------------------
gammas=()
ms=()
rs=()
phi_hs=()
thetas=()
alphas=()
betas=()
epsilons=()

# --------------------------
# Loop through all combinations and check file existence
# --------------------------
for gamma in "${gammas_all[@]}"; do
  for m in "${ms_all[@]}"; do
    for r in "${rs_all[@]}"; do
      for theta in "${thetas_all[@]}"; do
        for alpha in "${alphas_all[@]}"; do
          for phi_h in "${phi_hs_all[@]}"; do
            for beta in "${betas_all[@]}"; do
              for epsilon in "${epsilons_all[@]}"; do

                file_name=$(printf "epsilon=%.2f, gamma=%.2f, sigma=%d, phi_h=%.2f, m_bar=%.2f, r=%.2f, alpha=%.2f, theta=%.2f, beta=%.2f.mat" \
                            $epsilon $gamma $sigma $phi_h $m $r $alpha $theta $beta)

                full_path="${data_dir}/${file_name}"

                if [[ ! -f "$full_path" ]]; then
                  # Append missing case to rerun arrays
                  gammas+=($gamma)
                  ms+=($m)
                  rs+=($r)
                  thetas+=($theta)
                  alphas+=($alpha)
                  phi_hs+=($phi_h)
                  betas+=($beta)
                  epsilons+=($epsilon)
                fi

              done
            done
          done
        done
      done
    done
  done
done

# --------------------------
# Submit missing cases
# --------------------------
echo "Submitting ${#gammas[@]} missing jobs..."

for ((i=0; i<${#gammas[@]}; i++)); do
  gamma=${gammas[$i]}
  m=${ms[$i]}
  r=${rs[$i]}
  theta=${thetas[$i]}
  alpha=${alphas[$i]}
  phi_h=${phi_hs[$i]}
  beta=${betas[$i]}
  epsilon=${epsilons[$i]}

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
