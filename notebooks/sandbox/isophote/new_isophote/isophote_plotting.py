"""
Plotting utilities for isophote analysis visualization.

This module provides reusable plotting functions for comparing isophote
fitting results between photutils and the optimized implementation.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as MPLEllipse


def normalize_angle(angle):
    """Normalize angle to [0, pi)"""
    return np.mod(angle, np.pi)


def plot_isophote_comparison(title, photutils_res, optimized_res, filename, image=None):
    """
    Create a comprehensive comparison plot of isophote fitting results.
    
    Parameters
    ----------
    title : str
        Plot title
    photutils_res : list
        Photutils IsophoteList
    optimized_res : list
        Optimized fitting results (list of dicts)
   filename : str
        Output filename for the plot
    image : 2D array, optional
        Original image to display with overplotted isophotes
    """
    # Extract data from photutils results
    p_sma = np.array([iso.sma for iso in photutils_res])
    p_intens = np.array([iso.intens for iso in photutils_res])
    p_eps = np.array([iso.eps for iso in photutils_res])
    p_pa = np.array([normalize_angle(iso.pa) for iso in photutils_res])
    p_x0 = np.array([iso.x0 for iso in photutils_res])
    p_y0 = np.array([iso.y0 for iso in photutils_res])
    
    # Extract uncertainties from photutils
    p_intens_err = np.array([iso.int_err for iso in photutils_res])
    p_eps_err = np.array([iso.ellip_err for iso in photutils_res])
    p_pa_err = np.array([iso.pa_err for iso in photutils_res])
    p_x0_err = np.array([iso.x0_err for iso in photutils_res])
    p_y0_err = np.array([iso.y0_err for iso in photutils_res])
    
    # Extract data from optimized results
    o_sma = np.array([r['sma'] for r in optimized_res])
    o_intens = np.array([r['intens'] for r in optimized_res])
    o_eps = np.array([r['eps'] for r in optimized_res])
    o_pa = np.array([normalize_angle(r['pa']) for r in optimized_res])
    o_x0 = np.array([r['x0'] for r in optimized_res])
    o_y0 = np.array([r['y0'] for r in optimized_res])
    
    # Extract uncertainties if available
    o_intens_err = np.array([r.get('intens_err', 0.0) for r in optimized_res])
    o_eps_err = np.array([r.get('eps_err', 0.0) for r in optimized_res])
    o_pa_err = np.array([r.get('pa_err', 0.0) for r in optimized_res])
    o_x0_err = np.array([r.get('x0_err', 0.0) for r in optimized_res])
    o_y0_err = np.array([r.get('y0_err', 0.0) for r in optimized_res])
    
    # Create figure layout
    if image is not None:
        fig = plt.figure(figsize=(25, 5))
        gs = fig.add_gridspec(1, 5)
        ax_img = fig.add_subplot(gs[0])
        axes = [fig.add_subplot(gs[i]) for i in range(1, 5)]
        
        # Display image with arcsinh scaling
        vmin, vmax = np.percentile(image[~np.isnan(image)], [1, 99])
        ax_img.imshow(np.arcsinh((image - vmin) / (vmax - vmin)), origin='lower', cmap='gray')
        ax_img.set_title('Image + Isophotes')
        
        # Overplot isophotes (sample every 5th to avoid clutter)
        for r in optimized_res[::5]:
            if r['sma'] > 0.5:
                b_over_a = 1 - r['eps']
                ellipse = MPLEllipse(
                    (r['x0'], r['y0']), 
                    2 * r['sma'], 2 * r['sma'] * b_over_a,
                    angle=np.degrees(r['pa']),
                    fill=False, edgecolor='red', linewidth=0.5, alpha=0.7
                )
                ax_img.add_patch(ellipse)
        ax_img.set_xlim(0, image.shape[1])
        ax_img.set_ylim(0, image.shape[0])
    else:
        fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    
    fig.suptitle(title)
    
    # Plot 1: Intensity profile
    axes[0].plot(p_sma**0.25, p_intens, 'k-', label='Photutils', linewidth=1.5)
    axes[0].plot(o_sma**0.25, o_intens, 'r--', label='Optimized', linewidth=1.5)
    axes[0].set_xlabel('SMA$^{0.25}$')
    axes[0].set_ylabel('Intensity')
    axes[0].set_yscale('log')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Ellipticity profile
    axes[1].plot(p_sma**0.25, p_eps, 'k-', label='Photutils', linewidth=1.5)
    # Add photutils error bars
    if np.any(p_eps_err > 0):
        axes[1].errorbar(p_sma**0.25, p_eps, yerr=p_eps_err, fmt='none', 
                        ecolor='black', elinewidth=0.5, alpha=0.3, capsize=2)
                        
    axes[1].plot(o_sma**0.25, o_eps, 'r--', label='Optimized', linewidth=1.5)
    # Add error bars if available and non-zero
    if np.any(o_eps_err > 0):
        axes[1].errorbar(o_sma**0.25, o_eps, yerr=o_eps_err, fmt='none', 
                        ecolor='red', elinewidth=0.5, alpha=0.5, capsize=2)
    axes[1].set_xlabel('SMA$^{0.25}$')
    axes[1].set_ylabel('Ellipticity')
    axes[1].set_ylim(0, 1)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Plot 3: Position Angle profile
    axes[2].plot(p_sma**0.25, np.degrees(p_pa), 'k-', label='Photutils', linewidth=1.5)
    # Add photutils error bars
    if np.any(p_pa_err > 0):
        axes[2].errorbar(p_sma**0.25, np.degrees(p_pa), yerr=np.degrees(p_pa_err), 
                        fmt='none', ecolor='black', elinewidth=0.5, alpha=0.3, capsize=2)
                        
    axes[2].plot(o_sma**0.25, np.degrees(o_pa), 'r--', label='Optimized', linewidth=1.5)
    # Add error bars if available and non-zero
    if np.any(o_pa_err > 0):
        axes[2].errorbar(o_sma**0.25, np.degrees(o_pa), yerr=np.degrees(o_pa_err), 
                        fmt='none', ecolor='red', elinewidth=0.5, alpha=0.5, capsize=2)
    axes[2].set_xlabel('SMA$^{0.25}$')
    axes[2].set_ylabel('PA (deg)')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    # Plot 4: Center coordinates
    axes[3].plot(p_sma**0.25, p_x0, 'k-', label='Photutils X0', linewidth=1.5)
    axes[3].plot(p_sma**0.25, p_y0, 'k:', label='Photutils Y0', linewidth=1.5)
    axes[3].plot(o_sma**0.25, o_x0, 'r--', label='Optimized X0', linewidth=1.5)
    axes[3].plot(o_sma**0.25, o_y0, 'r:', label='Optimized Y0', linewidth=1.5)
    axes[3].set_xlabel('SMA$^{0.25}$')
    axes[3].set_ylabel('Center (pixels)')
    axes[3].legend(fontsize=8)
    axes[3].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved plot to {filename}")


def calculate_accuracy_metrics(results, true_eps, true_pa, true_x0=None, true_y0=None, min_sma=2.0):
    """
    Calculate accuracy metrics for isophote fitting.
    
    Parameters
    ----------
    results : list
        Fitting results (list of dicts)
    true_eps, true_pa : float
        True ellipticity and position angle
    true_x0, true_y0 : float or None
        True center (for center deviation metric)
    min_sma : float
        Minimum SMA to include in metrics
        
    Returns
    -------
    eps_diff, pa_diff, center_dev : float
        Mean fractional differences and mean center deviation (pixels)
    """
    sma = np.array([r['sma'] for r in results])
    eps = np.array([r['eps'] for r in results])
    pa = np.array([normalize_angle(r['pa']) for r in results])
    
    mask = sma > min_sma
    if not np.any(mask):
        return np.nan, np.nan, np.nan
        
    # Use mean instead of median
    mean_eps = np.mean(eps[mask])
    mean_pa = np.mean(pa[mask])
    
    eps_diff = abs(mean_eps - true_eps) / true_eps
    pa_diff = abs(mean_pa - true_pa) / true_pa
    
    # Center deviation
    if true_x0 is not None and true_y0 is not None:
        x0 = np.array([r['x0'] for r in results])
        y0 = np.array([r['y0'] for r in results])
        center_dev = np.mean(np.sqrt((x0[mask] - true_x0)**2 + (y0[mask] - true_y0)**2))
    else:
        center_dev = np.nan
    
    return eps_diff, pa_diff, center_dev
