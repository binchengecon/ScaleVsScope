#!/bin/bash
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

data_dir="./data_beta"
sigma=3

# Initialize output arrays
gammas=()
ms=()
rs=()
phi_hs=()
thetas=()
alphas=()
betas=()
epsilons=()

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
                  echo "Missing: $file_name"

                  # Append to arrays
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

echo ""
echo "Total missing cases: ${#gammas[@]}"
echo "Arrays of missing parameters to rerun:"
echo "gammas=(${gammas[@]})"
echo "ms=(${ms[@]})"
echo "rs=(${rs[@]})"
echo "thetas=(${thetas[@]})"
echo "alphas=(${alphas[@]})"
echo "phi_hs=(${phi_hs[@]})"
echo "betas=(${betas[@]})"
echo "epsilons=(${epsilons[@]})"
