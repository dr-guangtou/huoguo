"""
Benchmark: Synthetic Gaussian Profile (Noise-Free)

Tests the optimized isophote fitting on a noise-free synthetic Gaussian
intensity profile. This provides a best-case accuracy assessment.
"""

import sys
import time
import numpy as np
sys.path.append('..')
from photutils.isophote import Ellipse, EllipseGeometry
from isophote_optimize import fit_image, get_elliptical_coordinates
from isophote_plotting import plot_isophote_comparison, calculate_accuracy_metrics


def create_synthetic_image(shape=(500, 500), x0=250, y0=250, eps=0.4, pa=1.0, sigma=50):
    """Create a noise-free synthetic Gaussian profile image."""
    y, x = np.mgrid[:shape[0], :shape[1]]
    sma, _ = get_elliptical_coordinates(x, y, x0, y0, pa, eps)
    image = 1000 * np.exp(-0.5 * (sma / sigma)**2)
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
    print("Benchmark: Synthetic Gaussian Profile (Noise-Free)")
    print("=" * 70)
    
    # Create synthetic image
    print("\nGenerating noise-free synthetic image...")
    image = create_synthetic_image()
    true_x0, true_y0 = 250, 250
    true_eps, true_pa = 0.4, 1.0
    
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
    output_filename = "benchmark_synthetic_gaussian_clean.png"
    print(f"\nGenerating comparison plot: {output_filename}")
    plot_isophote_comparison(
        "Synthetic Gaussian Profile (Noise-Free)",
        p_results, o_results, output_filename, image
    )
    
    print("\n" + "=" * 70)
    print("Benchmark completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
