import numpy as np
import sys
import os
from astropy.io import fits

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import confeti
import confeti_viz

def run_benchmark():
    print("Running Benchmark: Real Galaxy M51")
    
    # Load Data
    fits_path = os.path.join(os.path.dirname(__file__), 'M51.fits')
    try:
        hdul = fits.open(fits_path)
        data = hdul[0].data
        hdul.close()
    except Exception as e:
        print(f"Error loading FITS file: {e}")
        return

    if data is None:
        print("No data found.")
        return

    # Preprocess
    data = np.nan_to_num(data)
    h, w = data.shape
    corners = np.concatenate([
        data[0:50, 0:50].ravel(),
        data[0:50, w-50:w].ravel(),
        data[h-50:h, 0:50].ravel(),
        data[h-50:h, w-50:w].ravel()
    ])
    bg = np.median(corners)
    data_sub = data - bg
    extent = [0, w, 0, h]

    # Fit
    fitter = confeti.Confeti(data_sub, extent=extent)
    fitter.fit(mode='c0', optimizer_method='L-BFGS-B')
    model = fitter.model_image()

    # Metric
    residual = data_sub - model
    metric = np.sum(np.abs(residual)) / np.sum(data_sub)
    print(f"Metric: {metric:.4f}")

    # Save CSV
    fitter.get_isophotes().write('benchmark_real_galaxy_m51.csv', format='csv', overwrite=True)

    # Plot
    # QA Dashboard
    results = fitter.results
    profile_data = {'r': results['r'], 'val': results['val']}
    
    visualize.plot_comprehensive_qa(
        original=data_sub,
        model=model,
        residual=residual,
        contours=fitter.contours_out,
        profile_data=profile_data,
        extent=extent,
        title_prefix="M51",
        output_filename="benchmark_real_galaxy_m51.png"
    )

if __name__ == "__main__":
    run_benchmark()
