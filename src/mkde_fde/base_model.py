import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.linalg import fractional_matrix_power
from abc import ABC, abstractmethod
from typing import Optional, Tuple

# --- Configuration for Professional Plotting ---
plt.style.use('seaborn-v0_8-whitegrid')
np.random.seed(42)


class BaseDensityEstimator(ABC):
    """
    Abstract Base Class for Density Estimators to ensure a standardized API.
    """

    @abstractmethod
    def fit(self, X: np.ndarray) -> 'BaseDensityEstimator':
        """Load and prepare data."""
        pass

    @abstractmethod
    def score_samples(self, X_test: np.ndarray) -> np.ndarray:
        """Estimate density for query points."""
        pass


def generate_correlated_data(n_samples: int = 1000) -> np.ndarray:
    """
    Generates a 2D dataset with strong correlation (diagonal shape).

    This simulates the 'diagonal line' problem that is hard for standard KDE.

    Returns:
        np.ndarray: Dataset of shape (n, 2).
    """
    mean = [0, 0]
    # Strong covariance: 0.9 correlation between X and Y
    cov = [[1, 0.9],
           [0.9, 1]]
    return np.random.multivariate_normal(mean, cov, n_samples)