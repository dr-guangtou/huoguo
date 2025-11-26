import numpy as np
import sys
import os

# Add parent directory to path to import confeti
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import confeti
import visualize

def run_benchmark():
    print("Running Benchmark: Synthetic Gaussian Profiles")
    
    # --- Example 0: Pure Gaussian Profile (Noiseless) ---
    print("\n--- Example 0: Pure Gaussian ---")
    xx, yy = np.meshgrid(np.linspace(-8, 8, 100, endpoint=True), np.linspace(-8, 8, 100, endpoint=True))
    h1 = 10. * np.exp(-0.5 * (xx**2 + yy**2))
    extent = [-8, 8, -8, 8]
    
    fitter0 = confeti.Confeti(h1, extent=extent)
    fitter0.fit(mode='simple', optimizer_method='L-BFGS-B')
    model0 = fitter0.model_image()
    
    residual0 = h1 - model0
    metric0 = np.sum(np.abs(residual0)) / np.sum(h1)
    print(f"Metric: {metric0:.4f}")
    
    # Save CSV
    fitter0.get_isophotes().write('benchmark_synthetic_gaussian_pure.csv', format='csv', overwrite=True)
    
    # QA Dashboard
    results0 = fitter0.results
    profile_data0 = {'r': results0['r'], 'val': results0['val']}
    
    visualize.plot_comprehensive_qa(
        original=h1,
        model=model0,
        residual=residual0,
        contours=fitter0.contours_out,
        profile_data=profile_data0,
        extent=extent,
        title_prefix="Pure Gaussian",
        output_filename="benchmark_synthetic_gaussian_pure.png"
    )

    # --- Example 1: Simple Gaussian Profile with Noise ---
    print("\n--- Example 1: Gaussian with Noise ---")
    np.random.seed(42)
    h_noise = h1 + np.random.normal(0, 1, h1.shape)
    
    fitter1 = confeti.Confeti(h_noise, extent=extent)
    fitter1.fit(mode='simple', optimizer_method='L-BFGS-B', smoothing_sigma=1.0)
    model1 = fitter1.model_image()
    
    residual1 = h_noise - model1
    metric1 = np.sum(np.abs(residual1)) / np.sum(h_noise)
    print(f"Metric: {metric1:.4f}")
    
    # Save CSV
    fitter1.get_isophotes().write('benchmark_synthetic_gaussian_noisy.csv', format='csv', overwrite=True)
    
    print("Fitted Isophote Levels:")
    print(fitter1.results['val'])
    
    # QA Dashboard
    results1 = fitter1.results
    profile_data1 = {'r': results1['r'], 'val': results1['val']}
    
    visualize.plot_comprehensive_qa(
        original=h_noise,
        model=model1,
        residual=residual1,
        contours=fitter1.contours_out,
        profile_data=profile_data1,
        extent=extent,
        title_prefix="Noisy Gaussian",
        output_filename="benchmark_synthetic_gaussian_noisy.png"
    )

if __name__ == "__main__":
    run_benchmark()
