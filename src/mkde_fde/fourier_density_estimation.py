from base_model import BaseDensityEstimator
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.linalg import fractional_matrix_power
from abc import ABC, abstractmethod
from typing import Optional, Tuple

class FourierDensityEstimator(BaseDensityEstimator):
    """
    Fourier Density Estimator based on the Fourier Integral Theorem.

    Unlike WhitedMKDE, this estimator does NOT explicitly model covariance.
    It relies on high-frequency sampling (Sinc kernel) to capture structure.

    Ref: Ho & Walker (2021). Multivariate Smoothing via the Fourier Integral Theorem.

    Analogy: Using a Retina Display (High R) so individual pixels are small enough
    to draw any diagonal line without aliasing, removing the need to rotate the image.
    """

    def __init__(self, radius: float = 1.0) -> None:
        """
        Args:
            radius (float): The 'Resolution' parameter R.
                            Analogous to PPI (Pixels Per Inch).
        """
        self.radius = radius
        self.X_train: Optional[np.ndarray] = None
        self.n_samples: int = 0
        self.d_features: int = 0

    def fit(self, X: np.ndarray) -> 'FourierDensityEstimator':
        """
        Store data for lazy evaluation. No matrix inversion happens here.
        """
        self.X_train = X
        self.n_samples, self.d_features = X.shape
        return self

    def score_samples(self, X_test: np.ndarray) -> np.ndarray:
        """
        Compute density using the product of Sinc kernels.

        Formula: f(x) = 1/(n * pi^d) * sum( prod( sin(R(x-Xi))/(x-Xi) ) )
        """
        if self.X_train is None:
            raise RuntimeError("Estimator not fitted.")

        n_queries = X_test.shape[0]
        densities = np.zeros(n_queries)

        # Constant factor from the paper: 1 / (n * pi^d)
        const_factor = 1.0 / (self.n_samples * (np.pi ** self.d_features))

        # Loop through query points (can be vectorized but kept explicit for clarity)
        for i in range(n_queries):
            # 1. Compute difference vector (x - X_i)
            # Shape: (n_samples, d_features)
            diff = self.X_train - X_test[i]

            # 2. Scale by Radius R
            u = self.radius * diff

            # 3. Apply Sinc Kernel: sin(u)/u (Dimension-wise)
            # Handle u=0 singularity: limit is 1.0
            # Note: The paper's kernel is K_R(u) = sin(Ru)/u, so we divide by diff, not u.
            # Actually, Eq 48 in paper: sin(R(y-x)) / (y-x).

            with np.errstate(divide='ignore', invalid='ignore'):
                kernel_vals = np.sin(u) / diff

            # Fix singularity where diff is 0. Limit sin(Rx)/x as x->0 is R.
            kernel_vals = np.where(diff == 0, self.radius, kernel_vals)

            # 4. Product over dimensions (The "Grid" logic)
            # This preserves dependence structure without a matrix!
            prod_kernel = np.prod(kernel_vals, axis=1)

            # 5. Sum over samples
            densities[i] = const_factor * np.sum(prod_kernel)

        return densities


