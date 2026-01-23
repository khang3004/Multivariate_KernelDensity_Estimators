# Mathematical Proof: Convergence of Leverage Scores ($W_{ii} \to 1$) in the Limit $h \to 0$


---

## 1. Problem Definition

The LOOCV shortcut formula is given by:
$$
CV(h) = \frac{1}{n} \sum_{i=1}^n \left( \frac{Y_i - \hat{m}(X_i)}{1 - W_{ii}} \right)^2
$$
**Objective:** Prove that $\lim_{h \to 0} W_{ii} = 1$, causing numerical singularity ($1-W_{ii} \to 0$).

**Notation:**
* $\hat{m}(X_i) = \mathbf{w}(X_i)^T \mathbf{Y}$ (Linear Smoother).
* $W_{ii}$: The $i$-th element of vector $\mathbf{w}(X_i)$.
* $Z$: Design matrix centered at $X_i$.

---

## 2. Asymptotic Analysis of Kernel Weights

Let $K(\cdot)$ be a kernel function with $K(0) > 0$ and $K(\infty) = 0$. As bandwidth $h \to 0$:

1.  **Target point ($j=i$):**
    $$W_{ii} = K\left(\frac{X_i - X_i}{h}\right) = K(0) \neq 0$$
2.  **Neighbors ($j \neq i$):**
    $$\lim_{h \to 0} W_{jj} = K\left(\frac{X_j - X_i}{h}\right) \to 0$$

$\implies \mathbf{W}$ converges to a sparse matrix with a single non-zero entry at $(i,i)$.

---

## 3. Method 1: Matrix Algebra (Moore-Penrose Pseudoinverse)

The estimator is $\hat{\beta} = (Z^T \mathbf{W} Z)^{-1} Z^T \mathbf{W} \mathbf{Y}$.
At the limit $h \to 0$, $A = Z^T \mathbf{W} Z$ becomes singular. We apply the **Moore-Penrose Pseudoinverse ($A^+$)**.

### 3.1. Limit of the Information Matrix
$$
A = Z^T \mathbf{W} Z \approx z_i^T W_{ii} z_i
$$
Since $z_i = [1, 0, \dots, 0]^T$ (intercept only):
$$
A \approx \text{diag}(K(0), 0, \dots, 0)
$$

### 3.2. Computing the Pseudoinverse ($A^+$)
For a diagonal matrix, $A^+$ inverts non-zero elements and preserves zeros:
$$
A^+ = \text{diag}\left(\frac{1}{K(0)}, 0, \dots, 0\right)
$$

### 3.3. Deriving the Estimator
$$
\hat{\beta} = A^+ (Z^T \mathbf{W} \mathbf{Y})
$$
Substituting terms:
$$
\hat{\beta} \approx \begin{bmatrix} 1/K(0) & 0 & \dots \\ 0 & 0 & \dots \end{bmatrix} \begin{bmatrix} K(0)Y_i \\ 0 \\ \vdots \end{bmatrix} = \begin{bmatrix} Y_i \\ 0 \\ \vdots \end{bmatrix}
$$

**Result:** $\hat{m}(X_i) = \hat{\beta}_0 = 1 \cdot Y_i$.
Comparing with $\hat{m}(X_i) = \sum W_{ij}Y_j$, we deduce:
$$
\boxed{\lim_{h \to 0} W_{ii} = 1}
$$

---

## 4. Method 2: Optimization (Normal Equations)

We analyze the **Weighted Least Squares** objective function directly.

### 4.1. The Objective Function
$$
J(\beta) = \sum_{j=1}^n W_{jj} \left( Y_j - \sum_{r=0}^p \beta_r (X_j - X_i)^r \right)^2
$$

### 4.2. Limit Behavior ($h \to 0$)
Since $W_{jj} \to 0$ for all $j \neq i$, the summation collapses to the $i$-th term:
$$
J(\beta) \approx W_{ii} (Y_i - \beta_0)^2 = K(0)(Y_i - \beta_0)^2
$$
*(Higher-order terms vanish since $X_i - X_i = 0$)*.

### 4.3. Minimization
$$
\frac{\partial J}{\partial \beta_0} = -2K(0)(Y_i - \beta_0) = 0
$$
$$
\implies \hat{\beta}_0 = Y_i
$$

### 4.4. Conclusion
The model acts as an **Exact Interpolator** ($\hat{Y}_i = Y_i$).
This implies the effective weight vector $\mathbf{w}(X_i)$ converges to the basis vector $e_i$:
$$
\mathbf{w}(X_i) \to [0, \dots, 1, \dots, 0]^T
$$
$$
\boxed{W_{ii} \to 1}
$$

---

## 5. Summary

Both derivations confirm that as $h \to 0$, Local Polynomial Regression degenerates to Nearest Neighbor Interpolation.
$$
\text{Singularity:} \quad \lim_{h \to 0} \frac{1}{(1 - W_{ii})^2} = \infty
$$

## 6. Numerical Stabilization Proposals

To prevent the singularity ($1 - W_{ii} \to 0$) in practical implementation, we propose three robust strategies.

### 6.1. Bandwidth Lower Bound Constraint
**Concept:** Enforce a hard constraint on $h$ based on the data geometry (nearest-neighbor distance).
$$
h_{min} = \frac{1}{2} \min_{i \neq j} |X_i - X_j|
$$
**Algorithm:**
$$
\text{Search Space: } \mathcal{H} = \{ h \in \mathbb{R}^+ \mid h > h_{min} \}
$$

### 6.2. Ridge Regularization (L2)
**Concept:** Add a damping factor $\lambda$ to the inversion to prevent rank deficiency.
$$
\hat{\beta}_{Ridge} = (Z^T \mathbf{W} Z + \lambda I)^{-1} Z^T \mathbf{W} \mathbf{Y}
$$
**Effect:**
$$
W_{ii}^{Ridge} \approx \frac{K(0)}{K(0) + \lambda} < 1 \quad (\text{Strictly})
$$
This ensures the denominator $1 - W_{ii}$ is always non-zero.

### 6.3. Numerical Clipping (Epsilon Trick)
**Concept:** Introduce a small perturbation $\epsilon$ (e.g., $10^{-6}$) to the denominator to bound the penalty.
$$
CV(h) \approx \frac{1}{n} \sum_{i=1}^n \left( \frac{e_i}{\max(1 - W_{ii}, \epsilon)} \right)^2
$$
**Behavior:**
* If $W_{ii} \to 1$, the term becomes $\frac{e_i}{\epsilon} \to \text{Large Penalty}$.
* The optimizer will naturally reject $h \to 0$ without crashing.