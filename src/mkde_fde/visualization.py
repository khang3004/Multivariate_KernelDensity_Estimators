from whited_mkde import *
from fourier_density_estimation import *
from base_model import *
def visualize_comparison():
    """
    Visual benchmark between Whited MKDE and FDE.
    """
    # 1. Generate Correlated Data
    X = generate_correlated_data(n_samples=500)

    # 2. Define Models
    # Whited MKDE (The "Rotation" approach)
    mkde = WhitedMKDE()
    mkde.fit(X)

    # FDE (The "Retina" approach)
    # Heuristic for R: log(n) is a safe start, usually around 3-10 for normalized data.
    # Since our data variance is ~1, R=3 is decent.
    fde = FourierDensityEstimator(radius=1.0)
    fde.fit(X)

    # 3. Create Grid for Plotting
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, 100),
                         np.linspace(y_min, y_max, 100))
    grid_points = np.c_[xx.ravel(), yy.ravel()]

    # 4. Predict
    z_mkde = mkde.score_samples(grid_points).reshape(xx.shape)
    z_fde = fde.score_samples(grid_points).reshape(xx.shape)

    # 5. Plotting
    fig, ax = plt.subplots(1, 2, figsize=(16, 7))

    # Plot Whited MKDE
    ax[0].scatter(X[:, 0], X[:, 1], s=5, alpha=0.5, color='gray')
    ax[0].contour(xx, yy, z_mkde, levels=10, cmap='viridis')
    ax[0].set_title("Whited MKDE\n(Explicit Matrix Rotation)", fontsize=14)
    ax[0].set_xlabel("Feature 1")
    ax[0].set_ylabel("Feature 2")

    # Plot FDE
    ax[1].scatter(X[:, 0], X[:, 1], s=5, alpha=0.5, color='gray')
    # Use 'levels' to show density. Note: FDE can be negative in tails!
    # We plot positive contours primarily.
    ax[1].contour(xx, yy, z_fde, levels=10, cmap='plasma')
    ax[1].set_title(f"Fourier Density Estimator (R={fde.radius})\n(High Resolution Sampling)", fontsize=14)
    ax[1].set_xlabel("Feature 1")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    visualize_comparison()