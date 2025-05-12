
# EGM for Neoclassical Growth Model (Documentation)

## 📘 Model Setup

We solve the deterministic neoclassical growth model:

$$
V(k) = \max_{k'} \left\{ u(k^\alpha - k') + \beta V(k') \right\}
$$

subject to the resource constraint:

$$
c = k^\alpha - k'
$$

and utility:

$$
u(c) = \frac{c^{1-\gamma}}{1-\gamma}
$$

---

## ⚙️ Endogenous Grid Method (EGM) Steps

### 1. Fix a Grid on \(k'\)

- We define a grid over next period's capital: `k_prime_grid`.

### 2. Initialize \(V'(k')\)

- Start with a guess for the marginal value function \(V'(k')\).
- This drives the inversion of the Euler equation.

### 3. Invert the Euler Equation

$$
u'(c) = \beta V'(k') \quad \Rightarrow \quad c = (\beta V'(k'))^{-1/\gamma}
$$

### 4. Recover Implied Current Capital \(k\)

$$
k^\alpha = c + k' \quad \Rightarrow \quad k = (c + k')^{1/\alpha}
$$

### 5. Interpolate Back to Fixed Grid

- Interpolate \(k'(k)\) and \(c(k)\) onto a fixed grid over current capital \(k\).

### 6. Envelope Condition Update

We use the envelope condition:

$$
V'(k) = u'(c(k)) \cdot \alpha k^{\alpha - 1}
$$

to update the marginal value function.

### 7. Iterate to Convergence

Repeat until \(V'(k)\) converges.

---

## 🧮 MATLAB Code Highlights

- Uses vectorized operations to compute \(c\) and \(k\).
- Efficient due to lack of inner optimization loop.
- Interpolation used to link endogenous and exogenous grids.

---

## ✅ Final Output

- The policy function \(k'(k)\) is computed and plotted.
- Fast convergence with minimal tuning.

---

## 🔍 Why EGM is Efficient

- Avoids nested maximization.
- Only requires Euler inversion and interpolation.
- Very suitable for deterministic and stochastic dynamic models.

---

## 📌 Key Equations

Euler inversion:

$$
c = (\beta V'(k'))^{-1/\gamma}
$$

Envelope condition:

$$
V'(k) = u'(c(k)) \cdot \alpha k^{\alpha - 1}
$$

Budget constraint:

$$
k^\alpha = c + k'
$$

---

*End of EGM notes.*
