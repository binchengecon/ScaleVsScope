
# EGM for Neoclassical Growth Model (Consumption-Policy-Based)

## 📘 Model Setup

We solve the deterministic neoclassical growth model:

$$
V(k) = \max_{k'} \left\{ u(k^\alpha - k') + \beta V(k') \right\}
$$

subject to:

$$
c = k^\alpha - k'
$$

with CRRA utility:

$$
u(c) = \frac{c^{1-\gamma}}{1-\gamma}
$$

where $$ \gamma > 0 $$ and $$ 0 < \alpha < 1 $$.

---

## 🛠 Method: Directly Update Consumption Policy

Instead of operating on $$V$$ or $$V'$$, we work **directly with the consumption policy** using the envelope and Euler equations.

### Key Euler-Envelope Link:

$$
u_c(c(k)) = \beta u_c(c(k')) \alpha (k')^{\alpha-1}
$$

where $$ c(k') $$ is future consumption.

---

## ⚙️ Steps

1. **Fix a grid** over next-period capital $$k'$$.
2. **Guess initial consumption policy** $$c(k)$$.
3. For each $$k'$$:
   - Interpolate $$c'(k')$$ from previous iteration.
   - Apply the Euler condition to find today's marginal utility.
   - Invert marginal utility to recover today's $$c(k)$$.
4. **Budget constraint**:
   - Recover current $$k$$ from:
     $$
     k^\alpha = c(k) + k'
     \quad \Rightarrow \quad
     k = (c + k')^{1/\alpha}
     $$
5. **Interpolation**:
   - Interpolate policies back onto fixed $$k$$-grid.
6. **Update consumption policy** and iterate until convergence.

---

## 📜 MATLAB Code Sketch

```matlab
% Parameters
alpha = 0.3; beta = 0.95; gamma = 2;

% Grids
nk = 500;
k_min = 0.01; k_max = 5;
k_prime_grid = linspace(k_min, k_max, nk)';
k_grid = linspace(k_min, k_max, nk)';

% Initial guess
c_policy_old = k_grid.^alpha;

tol = 1e-6;
max_iter = 1000;

for iter = 1:max_iter
    c_prime = interp1(k_grid, c_policy_old, k_prime_grid, 'linear', 'extrap');
    uc_prime = c_prime.^(-gamma);
    uc_today = beta * uc_prime .* alpha .* k_prime_grid.^(alpha-1);
    c_today = uc_today.^(-1/gamma);
    k_now = (c_today + k_prime_grid).^(1/alpha);

    valid = (k_now >= k_min) & (k_now <= k_max);
    k_now = k_now(valid);
    k_prime_valid = k_prime_grid(valid);
    c_valid = c_today(valid);

    k_policy = interp1(k_now, k_prime_valid, k_grid, 'linear', 'extrap');
    c_policy = interp1(k_now, c_valid, k_grid, 'linear', 'extrap');

    if max(abs(c_policy - c_policy_old)) < tol
        break;
    end
    c_policy_old = c_policy;
end
```

---

## ✅ Summary

- No need to compute $$V$$ or $$V'$$ explicitly.
- Directly work with $$c(k)$$ and $$k'(k)$$.
- Euler equation plus envelope condition is sufficient.
- Very fast, simple, and robust.

---

*End of EGM (Consumption-Policy-Based) Notes.*
