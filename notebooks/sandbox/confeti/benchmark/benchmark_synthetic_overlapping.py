import numpy as np
import sys
import os
from scipy.stats import norm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import confeti
import visualize

def run_benchmark():
    print("Running Benchmark: Overlapping Profiles")
    
    # Create fake data
    xx, yy = np.meshgrid(np.linspace(-8, 8, 100, endpoint=True), np.linspace(-8, 8, 100, endpoint=True))
    h1 = 10. * np.exp(-0.5 * (xx**2 + yy**2))
    h2 = 5. * np.exp(-0.5 * ((xx - 3)**2 + (yy - 3)**2) / 0.5) # Smaller, offset companion
    
    np.random.seed(42)
    noise = np.random.normal(0, 1, h1.shape)
    H = h1 + h2 + noise
    extent = [-8, 8, -8, 8]
    
    # Mask the companion
    mask = np.ones_like(H)
    mask[(xx - 3)**2 + (yy - 3)**2 < 2] = 0
    
    # Fit Main Object
    print("Fitting Main Object...")
    fitter = confeti.Confeti(H, mask=mask, extent=extent)
    fitter.fit(mode='simple', optimizer_method='L-BFGS-B')
    model = fitter.model_image()
    
    # Fit Companion on Residual
    print("Fitting Companion...")
    residual = H - model
    # Mask the main object center for companion fit? Or just fit the residual peak?
    # Let's just fit the residual with a box around the companion
    
    # Simple approach: Fit the residual directly
    fitter2 = confeti.Confeti(residual, extent=extent)
    # Give it a hint for center
    fitter2.fit(mode='simple', center=[3, 3], optimizer_method='L-BFGS-B') 
    model2 = fitter2.model_image()
    
    final_model = model + model2
    final_residual = H - final_model
    
    metric = np.sum(np.abs(final_residual)) / np.sum(H)
    print(f"Metric: {metric:.4f}")
    
    # Save CSVs
    fitter.get_isophotes().write('benchmark_synthetic_overlapping_main.csv', format='csv', overwrite=True)
    fitter2.get_isophotes().write('benchmark_synthetic_overlapping_companion.csv', format='csv', overwrite=True)
    
    # QA Dashboard
    # We use the main object's contours and profile for visualization
    results = fitter.results
    profile_data = {'r': results['r'], 'val': results['val']}
    
    # Combine contours from both fits for visualization?
    # Or just show main object. Let's show main object contours.
    
    visualize.plot_comprehensive_qa(
        original=H,
        model=final_model,
        residual=final_residual,
        contours=fitter.contours_out,
        profile_data=profile_data,
        extent=extent,
        title_prefix="Overlapping",
        output_filename="benchmark_synthetic_overlapping.png"
    )

if __name__ == "__main__":
    run_benchmark()
