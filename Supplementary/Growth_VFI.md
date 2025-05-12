
# Basic Neoclassical Growth Model (Discrete Time) - Documentation

## Overview

This document describes the setup, economic intuition, and numerical solution of a basic discrete-time neoclassical growth model using **Value Function Iteration (VFI)**.

---

## 1. Problem Setup

### Environment

- **Time**: Discrete, infinite horizon $t=0,1,2,\ldots$
- **Agent**: Representative household or planner
- **State variable**: Capital stock $k$
- **Control variable**: Next period capital $k'$ (or consumption $c$)

### Technologies and Preferences

- **Production function**:

  $$
  y_t = k_t^{\alpha}
  $$

  where $0 < \alpha < 1$ is the output elasticity of capital.

- **Resource constraint**:

  $$
  c_t + k_{t+1} = y_t = k_t^{\alpha}
  $$

- **Utility function**:

  $$
  u(c) = \frac{c^{1-\gamma}}{1-\gamma}
  $$

  where $\gamma$ is the coefficient of relative risk aversion.

- **Objective**:

  $$
  \max_{\{k_{t+1}\}_{t=0}^{\infty}} \sum_{t=0}^{\infty} \beta^t u(c_t)
  $$

  subject to the resource constraint.

---

## 2. Bellman Equation

The recursive formulation is:

$$
V(k) = \max_{k'} \left\{ u(k^\alpha - k') + \beta V(k') \right\}
$$

subject to the feasibility constraint on $k'$.

---

### Economic Tradeoff

- **Consume today**: enjoy immediate utility.
- **Save (invest) today**: enjoy higher future consumption.

---

## 3. Numerical Solution via Value Function Iteration (VFI)

### Steps

1. Discretize the capital stock $k$ over a finite grid.
2. Guess an initial value function $V_0(k)$.
3. At each grid point $k$:
   - Search over all possible $k'$ choices.
   - Calculate consumption $c = k^\alpha - k'$.
   - Compute utility $u(c)$.
   - Compute total return $u(c) + \beta V(k')$.
   - Pick $k'$ that maximizes the total return.
4. Update the value function $V_1(k)$.
5. Repeat until convergence.

---

### MATLAB Code Sketch

```matlab
% Parameters
alpha = 0.3; beta = 0.95; gamma = 2;

% Grids
nk = 1000;
k_min = 0.01; k_max = 5;
k_grid = linspace(k_min, k_max, nk)';

% Initialization
V = zeros(nk,1);
policy_k = zeros(nk,1);
policy_c = zeros(nk,1);
tol = 1e-6;
max_iter = 1000;

for iter = 1:max_iter
    V_new = zeros(nk,1);
    for i = 1:nk
        k = k_grid(i);
        y = k^alpha;
        c = y - k_grid;
        c(c <= 0) = NaN;
        u = (c.^(1-gamma)) / (1-gamma);
        RHS = u + beta * V;
        [V_new(i), idx] = max(RHS);
        policy_k(i) = k_grid(idx);
        policy_c(i) = c(idx);
    end
    if max(abs(V_new - V)) < tol
        break;
    end
    V = V_new;
end

% Plot
plot(k_grid, policy_k);
xlabel('k today'); ylabel('k tomorrow');
title('Policy Function for Capital');
```

---

## 4. Summary Table

| Concept            | Description                        |
|:-------------------|:-----------------------------------|
| State Variable     | Current capital stock $k$          |
| Control Variable   | Next-period capital $k'$           |
| Production         | $y = k^\alpha$                     |
| Budget Constraint  | $c + k' = y$                       |
| Objective          | Maximize discounted lifetime utility |

---

## 5. Extensions

- Add **stochastic productivity** $z_t$: 
  $$
  y_t = z_t k_t^{\alpha}
  $$
- Solve with **Endogenous Grid Method (EGM)** for faster convergence.
- Analyze **transition dynamics**.

---

## 6. Economic Intuition

- Higher current capital $k$ leads to higher production.
- The agent faces a trade-off between **current consumption** and **future consumption**.
- The solution balances the marginal utility today versus the discounted value of future utility.

---

*End of Notes.*
