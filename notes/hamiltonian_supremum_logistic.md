# Computing the Hamiltonian Supremum for Logistic Intensity

## Problem Statement

In the HJB PDE (Eq. 1, [1]), the quoting Hamiltonian is

$$
H^{n,i,j}(z,p) = \sup_\delta \; f^{n,i,j}(z,\delta)\,(\delta - p),
$$

where the logistic intensity function is

$$
f(z,\delta) = \frac{1}{1 + e^{\alpha + \beta\delta}},
$$

and the amplitude $\lambda^{n,i,j}(z)$ factors out (it does not depend on $\delta$). The argument $p$ in the PDE context is

$$
p = \frac{\theta(t,Y) - \theta(t, Y + ze^i - ze^j)}{z},
$$

computed from finite differences of the value function on the grid.

The paper [1] states (referencing [2]) that the optimal quote $\bar\delta(z,p) = (f)^{-1}(-\partial_p H(z,p))$ "can easily be computed numerically in the logistic case." This characterisation follows from the envelope theorem ($\partial_p H = -f(\delta^*)$), but is circular as a computational recipe since evaluating $\partial_p H$ already requires knowing $\delta^*$.

## First-Order Condition

The objective $g(\delta) = f(\delta)(\delta - p)$ is strictly concave for the logistic function (it satisfies Guéant's sufficient condition $\sup_\delta \frac{\Lambda\Lambda''}{\Lambda'^2} < 2$, see Lemma 3.1 in [2]). Setting $g'(\delta) = 0$:

$$
f'(\delta)(\delta - p) + f(\delta) = 0.
$$

Using $f'(\delta) = -\beta f(\delta)(1 - f(\delta))$ and dividing by $f(\delta) > 0$:

$$
\delta - p = \frac{1}{\beta\,(1 - f(\delta))}. \tag{FOC}
$$

This is a single transcendental equation in one unknown. It has a unique solution $\delta^* > p$.

## Semi-Closed Form via Lambert W

### Derivation

Set $s = \frac{1}{1 - f(\delta^*)} > 1$. Then $f(\delta^*) = \frac{s-1}{s}$ and inverting the logistic gives $\delta^* = \frac{1}{\beta}[\ln\frac{1}{s-1} - \alpha]$. Substituting into the FOC:

$$
\frac{1}{\beta}\Bigl[\ln\frac{1}{s-1} - \alpha\Bigr] - p = \frac{s}{\beta}.
$$

Rearranging:

$$
\ln(s-1) = -(s + \alpha + \beta p),
$$

$$
(s-1)\,e^{s-1} = e^{-(1 + \alpha + \beta p)}.
$$

Setting $w = s - 1 > 0$, this is exactly the defining equation of the Lambert W function:

$$
\boxed{w = W_0\!\bigl(e^{-(1 + \alpha + \beta p)}\bigr)}
$$

where $W_0$ is the principal (real, non-negative) branch.

### Resulting Formulas

| Quantity | Expression |
|----------|------------|
| Optimal quote $\delta^*(p)$ | $p + \dfrac{W_0(e^{-(1+\alpha+\beta p)}) + 1}{\beta}$ |
| Fill probability $f(\delta^*)$ | $\dfrac{W_0(e^{-(1+\alpha+\beta p)})}{W_0(e^{-(1+\alpha+\beta p)}) + 1}$ |
| Hamiltonian $H(z,p)$ | $\lambda(z) \cdot \dfrac{W_0(e^{-(1+\alpha+\beta p)})}{\beta}$ |

**Verification.** $f(\delta^*)(\delta^* - p) = \frac{w}{w+1} \cdot \frac{w+1}{\beta} = \frac{w}{\beta}$. Consistent.

### Derivatives (useful for the quadratic approximation $\hat{H}$)

By the envelope theorem:

$$
\partial_p H(z,p) = -\lambda(z)\,f(\delta^*) = -\lambda(z)\,\frac{w}{w+1}.
$$

Differentiating once more and using $W_0'(x) = \frac{W_0(x)}{x(1+W_0(x))}$:

$$
\partial_{pp} H(z,p) = \lambda(z)\,\frac{\beta\,w}{(w+1)^3},
$$

where $w = W_0(e^{-(1+\alpha+\beta p)})$ throughout. These give the quadratic approximation coefficients (Eq. 4 in [1]):

$$
\alpha_0 = H(z,0), \quad \alpha_1 = \partial_p H(z,0), \quad \alpha_2 = \partial_{pp} H(z,0).
$$

## Numerical Methods

### Method 1: SciPy Lambert W

The most direct approach. `scipy.special.lambertw` computes $W_0$ to machine precision and accepts array inputs for vectorised evaluation.

```python
import numpy as np
from scipy.special import lambertw

def H_logistic_lambertw(p, alpha, beta):
    """Compute (H/lambda, delta_star, f_star) via Lambert W."""
    x = np.exp(-(1.0 + alpha + beta * p))
    w = np.real(lambertw(x))          # principal branch, real for x >= 0
    delta_star = p + (w + 1.0) / beta
    f_star = w / (w + 1.0)
    H = w / beta                      # multiply by lambda(z) externally
    return H, delta_star, f_star
```

- Fully vectorised over arrays of $p$ values (one call per grid slice).
- Adds `scipy` as a dependency.

### Method 2: Newton Iteration on $w e^w = x$

Avoids the scipy dependency and gives explicit control over precision. The Lambert W equation $w e^w = x$ is solved by Newton's method:

$$
w_{k+1} = w_k - \frac{w_k e^{w_k} - x}{e^{w_k}(w_k + 1)}.
$$

Starting from $w_0 = \max(0, \ln x)$ (or simply $w_0 = x$ for small $x$, since $W_0(x) \approx x$ when $x \to 0$), this converges in 3--5 iterations to machine precision due to the quadratic convergence of Newton's method.

```python
import numpy as np

def _lambert_w0_newton(x, tol=1e-12, max_iter=8):
    """Compute W_0(x) for x >= 0 via Newton iteration."""
    # Starting guess
    w = np.where(x < 1.0, x, np.log(1.0 + x))
    for _ in range(max_iter):
        ew = np.exp(w)
        wew = w * ew
        dw = (wew - x) / (ew * (w + 1.0))
        w = w - dw
        if np.all(np.abs(dw) < tol):
            break
    return w

def H_logistic_newton(p, alpha, beta):
    """Compute (H/lambda, delta_star, f_star) via Newton-based Lambert W."""
    x = np.exp(-(1.0 + alpha + beta * p))
    w = _lambert_w0_newton(x)
    delta_star = p + (w + 1.0) / beta
    f_star = w / (w + 1.0)
    H = w / beta
    return H, delta_star, f_star
```

- No external dependency beyond NumPy.
- Vectorised: `x`, `w`, `dw` can all be arrays.
- 3--5 iterations suffice; 8 is a conservative upper bound.
- Numerically stable: the argument $x = e^{-(1+\alpha+\beta p)}$ is always positive, and $W_0$ is smooth on $(0, \infty)$.

### Comparison

| Aspect | SciPy `lambertw` | Newton iteration |
|--------|-------------------|------------------|
| Precision | Machine precision | Machine precision (with 5--8 iters) |
| Dependencies | `scipy` | `numpy` only |
| Vectorisation | Native | Native |
| Overhead per call | One function call | 3--5 multiply/exp per iteration |
| Transparency | Black box | Explicit, easy to verify |

Both are far more efficient than the bisection approach (which requires ~50 iterations per scalar evaluation) and both vectorise naturally over the grid.

## Context in the PDE Solver

In the full HJB PDE (Eq. 1 of [1]), at each time step and each inventory grid point $Y$, for every client tier $n$, currency pair $(i,j)$, and trade size $z$:

1. Compute $p = \frac{\theta(t,Y) - \theta(t, Y + ze^i - ze^j)}{z}$ from the current value function grid.
2. Evaluate $H^{n,i,j}(z, p)$ using the Lambert W formula.
3. The integral $\int_{\mathbb{R}^*_+} z\, H^{n,i,j}(z, p)\, \lambda^{n,i,j}(z)\, dz$ becomes a weighted sum over discrete trade sizes.

Since this must be done at every grid point and every time step, vectorisation over $p$ (i.e., over the inventory grid) is the key performance lever. Both methods above support this naturally.

## References

- [1] A. Barzykin, P. Bergault, O. Guéant, "Dealing with multi-currency inventory risk in FX cash markets," arXiv:2207.04100v4, 2023.
- [2] O. Guéant, "Optimal market making," arXiv:1605.01862v5, 2017.
