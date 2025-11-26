"""
Benchmark: Synthetic Sersic n=4 Profile with Noise

Tests the optimized isophote fitting on a synthetic Sersic n=4 profile
(de Vaucouleurs profile) with added noise. This is representative of
elliptical galaxy profiles.
"""

import sys
import time
import numpy as np
from scipy.special import gammaincinv
sys.path.append('..')
from photutils.isophote import Ellipse, EllipseGeometry
from isophote_optimize import fit_image, get_elliptical_coordinates
from isophote_plotting import plot_isophote_comparison, calculate_accuracy_metrics


def create_sersic_image(shape=(500, 500), x0=250, y0=250, eps=0.3, pa=1.0, 
                        n=4.0, re=50, Ie=1000, noise_sigma=1.0):
    """Create a synthetic Sersic profile image."""
    y, x = np.mgrid[:shape[0], :shape[1]]
    sma, _ = get_elliptical_coordinates(x, y, x0, y0, pa, eps)
    
    # Sersic profile: I(r) = Ie * exp(-bn * ((r/re)^(1/n) - 1))
    # bn is chosen such that re encloses half the total flux
    bn = gammaincinv(2.0 * n, 0.5)
    
    image = Ie * np.exp(-bn * ((sma / re)**(1.0/n) - 1.0))
    
    if noise_sigma > 0:
        image += np.random.normal(0, noise_sigma, shape)
    return image


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
    print("Benchmark: Synthetic Sersic n=4 Profile with Noise")
    print("=" * 70)
    
    # Create synthetic image
    print("\nGenerating synthetic Sersic n=4 image...")
    image = create_sersic_image(eps=0.3, pa=1.0, n=4.0, noise_sigma=1.0)
    true_x0, true_y0 = 250, 250
    true_eps, true_pa = 0.3, 1.0
    
    # Run photutils
    print("\nRunning Photutils...")
    p_results, p_time = run_photutils(image, true_x0, true_y0)
    print(f"  Time: {p_time:.4f}s")
    print(f"  Isophotes fitted: {len(p_results)}")
    
    # Run optimized
    print("\nRunning Optimized...")
    o_results, o_time = run_optimized(image, true_x0, true_y0)
    print(f"  Time: {o_time:.4f}s")
    print(f"  Isophotes fitted: {len(o_results)}")
    
    # Calculate speedup
    speedup = p_time / o_time
    print(f"\n  Speedup: {speedup:.2f}x")
    
    # Calculate accuracy metrics
    eps_diff, pa_diff, center_dev = calculate_accuracy_metrics(
        o_results, true_eps, true_pa, true_x0, true_y0, min_sma=2.0
    )
    print(f"\nAccuracy Metrics (SMA > 2 pixels):")
    print(f"  Eps fractional difference: {eps_diff:.4f} ({eps_diff*100:.2f}%)")
    print(f"  PA fractional difference:  {pa_diff:.4f} ({pa_diff*100:.2f}%)")
    print(f"  Center deviation:          {center_dev:.2f} pixels")
    
    # Generate plot
    output_filename = "benchmark_synthetic_sersic_n4_noisy.png"
    print(f"\nGenerating comparison plot: {output_filename}")
    plot_isophote_comparison(
        "Synthetic Sersic n=4 Profile (Noise sigma=1.0, eps=0.3)",
        p_results, o_results, output_filename, image
    )
    
    print("\n" + "=" * 70)
    print("Benchmark completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
