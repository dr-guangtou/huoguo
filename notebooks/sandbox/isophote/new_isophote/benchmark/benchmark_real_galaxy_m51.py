"""
Benchmark: Real Galaxy M51 (FITS Image)

Tests the optimized isophote fitting on a real galaxy image (M51)
to assess performance on realistic data with complex structure and noise.
"""

import sys
import time
import numpy as np
sys.path.append('..')
from astropy.io import fits
from photutils.isophote import Ellipse, EllipseGeometry
from isophote_optimize import fit_image
from isophote_plotting import plot_isophote_comparison


def run_photutils(image, x0, y0, sma0=10.0, step=0.1):
    """Run photutils isophote fitting."""
    geometry = EllipseGeometry(x0=x0, y0=y0, sma=sma0, eps=0.1, pa=0.0)
    ellipse = Ellipse(image, geometry=geometry)
    
    t0 = time.time()
    isolist = ellipse.fit_image(sma0=sma0, minsma=0.0, maxsma=image.shape[0]/2, step=step)
    t1 = time.time()
    
    return isolist, t1 - t0


def run_optimized(image, x0, y0, sma0=10.0, step=0.1):
    """Run optimized isophote fitting."""
    config = {
        'x0': x0, 'y0': y0, 'sma0': sma0,
        'eps': 0.1, 'pa': 0.0,
        'minsma': 0.0, 'maxsma': image.shape[0]/2,
        'astep': step,
        'maxit': 50, 'conver': 0.05,
        'fix_center': False, 'fix_eps': False, 'fix_pa': False,
        'compute_errors': True
    }
    
    t0 = time.time()
    results = fit_image(image, None, config)
    t1 = time.time()
    
    # Filter valid results
    valid = [r for r in results if r['stop_code'] in [0, 1, 2]]
    return valid, t1 - t0


def main():
    print("=" * 70)
    print("Benchmark: Real Galaxy M51")
    print("=" * 70)
    
    # Load M51 FITS file
    try:
        print("\nLoading M51.fits...")
        with fits.open('../../M51.fits') as hdul:
            image = hdul[0].data
            
        print(f"  Image shape: {image.shape}")
        
        # Use image center as initial guess
        cy, cx = image.shape
        cx /= 2
        cy /= 2
        print(f"  Using center guess: ({cx:.1f}, {cy:.1f})")
        
        # Run photutils
        print("\nRunning Photutils...")
        p_results, p_time = run_photutils(image, cx, cy)
        print(f"  Time: {p_time:.4f}s")
        print(f"  Isophotes fitted: {len(p_results)}")
        
        # Run optimized
        print("\nRunning Optimized...")
        o_results, o_time = run_optimized(image, cx, cy)
        print(f"  Time: {o_time:.4f}s")
        print(f"  Isophotes fitted: {len(o_results)}")
        
        # Calculate speedup
        speedup = p_time / o_time
        print(f"\n  Speedup: {speedup:.2f}x")
        
        # Generate plot
        output_filename = "benchmark_real_galaxy_m51.png"
        print(f"\nGenerating comparison plot: {output_filename}")
        plot_isophote_comparison(
            "Real Galaxy: M51",
            p_results, o_results, output_filename, image
        )
        
        print("\n" + "=" * 70)
        print("Benchmark completed successfully!")
        print("=" * 70)
        
    except FileNotFoundError:
        print("\nERROR: M51.fits not found!")
        print("Please ensure M51.fits is in the parent directory.")
        print("Benchmark skipped.")
    except Exception as e:
        print(f"\nERROR: {e}")
        print("Benchmark failed.")


if __name__ == "__main__":
    main()
