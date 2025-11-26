"""
Profile: Optimized Isophote Fitting (Sersic Benchmark)

Runs cProfile on the optimized isophote fitting using the Sersic benchmark case.
"""

import sys
import cProfile
import pstats
import io
import numpy as np
from scipy.special import gammaincinv
sys.path.append('..')
from isophote_optimize import fit_image, get_elliptical_coordinates

def create_sersic_image(shape=(500, 500), x0=250, y0=250, eps=0.3, pa=1.0, 
                        n=4.0, re=50, Ie=1000, noise_sigma=1.0):
    """Create a synthetic Sersic profile image."""
    y, x = np.mgrid[:shape[0], :shape[1]]
    sma, _ = get_elliptical_coordinates(x, y, x0, y0, pa, eps)
    
    bn = gammaincinv(2.0 * n, 0.5)
    image = Ie * np.exp(-bn * ((sma / re)**(1.0/n) - 1.0))
    
    if noise_sigma > 0:
        image += np.random.normal(0, noise_sigma, shape)
    return image

def run_optimized(image, x0, y0, sma0=10.0, step=0.1):
    """Run optimized isophote fitting."""
    config = {
        'x0': x0, 'y0': y0, 'sma0': sma0,
        'eps': 0.1, 'pa': 0.0,
        'minsma': 0.0, 'maxsma': image.shape[0]/2,
        'astep': step,
        'maxit': 50, 'conver': 0.05,
        'fix_center': False, 'fix_eps': False, 'fix_pa': False,
        'compute_errors': True,
        'compute_deviations': True
    }
    
    fit_image(image, None, config)

def main():
    print("Generating synthetic Sersic n=4 image...")
    image = create_sersic_image(eps=0.3, pa=1.0, n=4.0, noise_sigma=1.0)
    true_x0, true_y0 = 250, 250
    
    print("Profiling Optimized Fitting...")
    pr = cProfile.Profile()
    pr.enable()
    
    run_optimized(image, true_x0, true_y0)
    
    pr.disable()
    
    s = io.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(30)
    print(s.getvalue())
    
    # Also print by total time (internal time)
    s = io.StringIO()
    sortby = 'tottime'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(30)
    print("\nSorted by Internal Time (tottime):")
    print(s.getvalue())

if __name__ == "__main__":
    main()
