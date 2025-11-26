import numpy as np
import sys
import os
from astropy.modeling.models import Sersic2D

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import confeti
import confeti_viz

def run_benchmark():
    print("Running Benchmark: Sersic Model (n=3.5, ell=0.35)")
    
    # Create Sersic model
    x, y = np.meshgrid(np.arange(200), np.arange(200))
    mod = Sersic2D(amplitude=10, r_eff=20, n=3.5, x_0=100, y_0=100, ellip=0.35, theta=np.radians(30))
    H_sersic = mod(x, y)
    extent = [0, 200, 0, 200]
    
    # Fit
    fitter = confeti.Confeti(H_sersic, extent=extent)
    fitter.fit(mode='c0', optimizer_method='L-BFGS-B')
    model = fitter.model_image()
    
    # Metric
    residual = H_sersic - model
    metric = np.sum(np.abs(residual)) / np.sum(H_sersic)
    print(f"Metric: {metric:.4f}")
    
    # QA Dashboard
    results = fitter.results
    profile_data = {'r': results['r'], 'val': results['val']}
    import visualize
    visualize.plot_comprehensive_qa(
        original=H_sersic,
        model=model,
        residual=residual,
        contours=fitter.contours_out,
        profile_data=profile_data,
        extent=extent,
        title_prefix="Sersic",
        output_filename="benchmark_synthetic_sersic_qa.png"
    )
    
    # Export Table
    table = fitter.get_isophotes()
    table.write('benchmark_synthetic_sersic.csv', format='csv', overwrite=True)
    print("Isophote table saved.")

if __name__ == "__main__":
    run_benchmark()
