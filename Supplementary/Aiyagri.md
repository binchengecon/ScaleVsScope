# Endogenous Grid Method (EGM) for the Aiyagari Model - Documentation

## Overview

This document describes how to solve the standard **Aiyagari (1994)** incomplete markets model using the  **Endogenous Grid Method (EGM)** .

---

## 1. The Aiyagari Model

### Environment

* **Time** : Discrete, infinite horizon t=0,1,2,…t = 0,1,2,\ldots
* **Agents** : Continuum of ex-ante identical households facing idiosyncratic income risk.
* **Markets** : Incomplete asset markets (only risk-free saving).

### Technologies and Preferences

* **Preferences** :
  ∑t=0∞βtu(ct)withu(c)=c1−γ1−γ\sum_{t=0}^{\infty} \beta^t u(c_t) \quad \text{with} \quad u(c) = \frac{c^{1-\gamma}}{1-\gamma}
* **Budget constraint** :
  c+a′≤Ra+yc + a' \leq Ra + y
* **Borrowing constraint** :
  a′≥a‾(usually a‾=0 or some negative number)a' \geq \underline{a} \quad (\text{usually } \underline{a}=0 \text{ or some negative number})
* **Income process** :
  * Idiosyncratic labor income yty_t follows a Markov process π(y′∣y)\pi(y'|y).

---

## 2. Bellman Equation

The household's problem is:

V(a,y)=max⁡a′{u(Ra+y−a′)+βE[V(a′,y′)]}V(a,y) = \max_{a'} \left\{ u(Ra + y - a') + \beta \mathbb{E}[V(a',y')] \right\}

subject to:

* Budget constraint
* Borrowing constraint a′≥a‾a' \geq \underline{a}

---

## 3. Endogenous Grid Method (EGM) Algorithm

The EGM exploits the Euler equation to avoid nested maximization.

### Steps

1. **Discretize** the grids:
   * Grid for a′a' (next period assets)
   * Grid for income yy
2. **Guess** an initial consumption policy c^0(a′,y)\hat{c}_0(a',y).
3. **Construct the RHS** of the Euler equation:
   B(a′,y)=βR∑y′π(y′∣y)uc(c^0(a′,y′))B(a',y) = \beta R \sum_{y'} \pi(y'|y) u_c(\hat{c}_0(a',y'))
4. **Invert the Euler equation** to solve for consumption today:
   uc(c(a′,y))=B(a′,y)⇒c(a′,y)=(B(a′,y))−1/γu_c(c(a',y)) = B(a',y) \quad \Rightarrow \quad c(a',y) = (B(a',y))^{-1/\gamma}
5. **Recover the current asset aa** from the budget constraint:
   a=1R(c(a′,y)+a′−y)a = \frac{1}{R}(c(a',y) + a' - y)
6. **Interpolation** :
   * Since aa is endogenous, interpolate c(a,y)c(a,y) onto a fixed asset grid.
   * For values where aa would imply binding borrowing constraint, use constraint conditions.
7. **Update policy functions** .
8. **Check convergence** .

---

## 4. MATLAB Code Sketch

```matlab
% Parameters
beta = 0.96; gamma = 2; R = 1.04;
ygrid = [1, 2];
P = [0.9, 0.1; 0.1, 0.9];
na = 100; a_min = 0; a_max = 20;
a_prime_grid = linspace(a_min, a_max, na)';

% Initial guess for c(a',y)
c_hat = zeros(na, 2);
for iy = 1:2
    c_hat(:,iy) = R * a_prime_grid + ygrid(iy);
end

tol = 1e-6; max_iter = 1000;

for iter = 1:max_iter
    c_new = zeros(na,2); a_star = zeros(na,2);
    for iy = 1:2
        % Step 3: Compute RHS of Euler equation
        B = zeros(na,1);
        for jy = 1:2
            B = B + beta * R * P(iy,jy) * (c_hat(:,jy)).^(-gamma);
        end
      
        % Step 4: Invert Euler equation
        c_tilde = B.^(-1/gamma);

        % Step 5: Recover a
        a_now = (c_tilde + a_prime_grid - ygrid(iy)) / R;

        % Store
        a_star(:,iy) = a_now;
        c_new(:,iy) = c_tilde;
    end

    % Step 6: Interpolate onto standard asset grid
    a_grid = linspace(a_min, a_max, na)';
    c_hat_interp = zeros(na,2);
    for iy = 1:2
        valid = (a_star(:,iy) >= a_min);
        c_hat_interp(:,iy) = interp1(a_star(valid,iy), c_new(valid,iy), a_grid, 'linear', 'extrap');
    end

    % Step 7: Convergence check
    if max(abs(c_hat_interp(:) - c_hat(:))) < tol
        break;
    end
    c_hat = c_hat_interp;
end
```

---

## 5. Summary Table


| Concept           | Description                                             |
| :------------------ | :-------------------------------------------------------- |
| State Variables   | Current assetsaa, current incomeyy                      |
| Control Variable  | Next-period assetsa′a'                                 |
| Budget Constraint | c+a′=Ra+yc + a' = Ra + y                               |
| Euler Equation    | uc(c)=βRE[uc(c′)]u_c(c) = \beta R \mathbb{E}[u_c(c')] |
| Market            | Incomplete (only saving, no insurance)                  |

---

## 6. Economic Intuition

* Households **smooth consumption** over time despite income shocks.
* **Precautionary savings** motive arises because income risk cannot be fully insured.
* Borrowing constraint can bind for low-income, low-asset households.
* Wealth distribution emerges endogenously.

---

*End of Notes.*
