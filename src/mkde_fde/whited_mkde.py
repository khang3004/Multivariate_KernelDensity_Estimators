from base_model import BaseDensityEstimator
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.linalg import fractional_matrix_power
from abc import ABC, abstractmethod
from typing import Optional, Tuple

class WhitedMKDE(BaseDensityEstimator):
    """
    Multivariate Kernel Density Estimator with Data Whitening.

    This estimator explicitly models the covariance structure of the data.
    It 'whitens' (decorrelates) the data before applying a standard Kernel,
    effectively rotating and scaling the data space to fit the kernel.

    Analogy: Rotating the image so it aligns with the pixels of a low-res screen.
    """

    def __init__(self, bandwidth: str = 'scott') -> None:
        """
        Initialize the estimator.

        Args:
            bandwidth (str): Method to select bandwidth for the internal KDE.
        """
        self.bandwidth_method = bandwidth
        self.whitening_matrix: Optional[np.ndarray] = None
        self.kde_engine: Optional[gaussian_kde] = None
        self.data_mean: Optional[np.ndarray] = None

    def fit(self, X: np.ndarray) -> 'WhitedMKDE':
        """
        Fit the model by computing the covariance and whitening the data.

        Args:
            X (np.ndarray): Training data (n_samples, n_features).
        """
        # 1. Compute Mean and Covariance Matrix (The H Matrix equivalent)
        self.data_mean = np.mean(X, axis=0)
        cov_matrix = np.cov(X, rowvar=False)

        # 2. Compute Whitening Matrix: Sigma^(-1/2)
        # This step is computationally expensive O(d^3) due to matrix power/inverse
        self.whitening_matrix = fractional_matrix_power(cov_matrix, -0.5)

        # 3. Whiten the data (Transform to spherical/uncorrelated space)
        # Formula: Z = (X - mu) @ W
        X_centered = X - self.data_mean
        X_whited = X_centered @ self.whitening_matrix

        # 4. Fit standard KDE on uncorrelated data
        self.kde_engine = gaussian_kde(X_whited.T, bw_method=self.bandwidth_method)

        return self

    def score_samples(self, X_test: np.ndarray) -> np.ndarray:
        """
        Estimate density.

        Must adjust for the volume change caused by the whitening transformation.
        Density_X(x) = Density_Z(z) * |det(W)|
        """
        if self.whitening_matrix is None or self.kde_engine is None:
            raise RuntimeError("Estimator not fitted.")

        # 1. Transform query points to the same whitened space
        X_test_centered = X_test - self.data_mean
        X_test_whited = X_test_centered @ self.whitening_matrix

        # 2. Get density in whitened space
        # gaussian_kde expects (n_features, n_samples)
        density_z = self.kde_engine(X_test_whited.T)

        # 3. Jacobian adjustment: density increases if we compressed space
        # Determinant of the whitening matrix
        det_W = np.linalg.det(self.whitening_matrix)

        return density_z * det_W