"""
Benchmark: Real Galaxy M51 (Fixed Geometry)

Tests the optimized isophote fitting on M51 with fixed geometry parameters:
- Fixed Center
- Fixed Ellipticity (eps=0.1)
- Fixed Position Angle (pa=0.0)

This verifies that the 'fix_*' configuration parameters work correctly.
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
    """Run photutils isophote fitting with fixed geometry."""
    # Note: photutils allows fixing parameters via geometry object
    geometry = EllipseGeometry(x0=x0, y0=y0, sma=sma0, eps=0.1, pa=0.0)
    geometry.fix_center = True
    geometry.fix_eps = True
    geometry.fix_pa = True
    
    ellipse = Ellipse(image, geometry=geometry)
    
    t0 = time.time()
    isolist = ellipse.fit_image(sma0=sma0, minsma=0.0, maxsma=image.shape[0]/2, step=step)
    t1 = time.time()
    
    return isolist, t1 - t0


def run_optimized(image, x0, y0, sma0=10.0, step=0.1):
    """Run optimized isophote fitting with fixed geometry."""
    config = {
        'x0': x0, 'y0': y0, 'sma0': sma0,
        'eps': 0.1, 'pa': 0.0,
        'minsma': 0.0, 'maxsma': image.shape[0]/2,
        'astep': step,
        'maxit': 50, 'conver': 0.05,
        # FIXED GEOMETRY PARAMETERS
        'fix_center': True, 
        'fix_eps': True, 
        'fix_pa': True,
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
    print("Benchmark: Real Galaxy M51 (Fixed Geometry)")
    print("=" * 70)
    print("Constraints:")
    print("  - Center: Fixed")
    print("  - Ellipticity: Fixed at 0.1")
    print("  - Position Angle: Fixed at 0.0")
    
    # Load M51 FITS file
    try:
        print("\nLoading M51.fits...")
        with fits.open('../../M51.fits') as hdul:
            image = hdul[0].data
            
        print(f"  Image shape: {image.shape}")
        
        # Use image center
        cy, cx = image.shape
        cx /= 2
        cy /= 2
        print(f"  Center: ({cx:.1f}, {cy:.1f})")
        
        # Run photutils
        print("\nRunning Photutils (Fixed)...")
        p_results, p_time = run_photutils(image, cx, cy)
        print(f"  Time: {p_time:.4f}s")
        print(f"  Isophotes fitted: {len(p_results)}")
        
        # Run optimized
        print("\nRunning Optimized (Fixed)...")
        o_results, o_time = run_optimized(image, cx, cy)
        print(f"  Time: {o_time:.4f}s")
        print(f"  Isophotes fitted: {len(o_results)}")
        
        # Calculate speedup
        speedup = p_time / o_time
        print(f"\n  Speedup: {speedup:.2f}x")
        
        # Verify constraints (exclude SMA=0 which is central pixel)
        print("\nVerifying Constraints (SMA > 0):")
        
        # Filter out SMA=0
        o_results_valid = [r for r in o_results if r['sma'] > 0.0]
        
        # Check Center
        o_x0 = np.array([r['x0'] for r in o_results_valid])
        o_y0 = np.array([r['y0'] for r in o_results_valid])
        x0_dev = np.max(np.abs(o_x0 - cx))
        y0_dev = np.max(np.abs(o_y0 - cy))
        print(f"  Max Center Deviation: {max(x0_dev, y0_dev):.6f} pixels")
        if max(x0_dev, y0_dev) < 1e-5:
            print("  [PASS] Center is fixed.")
        else:
            print("  [FAIL] Center moved!")
            
        # Check Ellipticity
        o_eps = np.array([r['eps'] for r in o_results_valid])
        eps_dev = np.max(np.abs(o_eps - 0.1))
        print(f"  Max Eps Deviation:    {eps_dev:.6f}")
        if eps_dev < 1e-5:
            print("  [PASS] Ellipticity is fixed.")
        else:
            print("  [FAIL] Ellipticity changed!")
            
        # Check PA
        o_pa = np.array([r['pa'] for r in o_results_valid])
        pa_dev = np.max(np.abs(o_pa - 0.0))
        print(f"  Max PA Deviation:     {pa_dev:.6f} rad")
        if pa_dev < 1e-5:
            print("  [PASS] PA is fixed.")
        else:
            print("  [FAIL] PA changed!")
        
        # Generate plot
        output_filename = "benchmark_real_galaxy_m51_fixed.png"
        print(f"\nGenerating comparison plot: {output_filename}")
        plot_isophote_comparison(
            "M51 (Fixed Geometry: eps=0.1, pa=0.0)",
            p_results, o_results, output_filename, image
        )
        
        print("\n" + "=" * 70)
        print("Benchmark completed successfully!")
        print("=" * 70)
        
    except FileNotFoundError:
        print("\nERROR: M51.fits not found!")
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
