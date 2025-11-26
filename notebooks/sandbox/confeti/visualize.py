import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from matplotlib.lines import Line2D
from typing import List, Optional

def plot_fitting_results(original: np.ndarray, model: np.ndarray, residual: np.ndarray, 
                         extent: Optional[List[float]] = None, 
                         title_prefix: str = "", 
                         output_filename: Optional[str] = None, 
                         metric: Optional[float] = None,
                         mask: Optional[np.ndarray] = None):
    """
    Plot Original, Model, and Residual images side-by-side.
    """
    plt.figure(figsize=(15, 5))
    
    # Use arcsinh stretching for visualization if dynamic range is high
    # Or just simple_norm
    norm = simple_norm(original, 'asinh', percent=99)
    
    # Handle mask for display
    display_orig = original.copy()
    display_resid = residual.copy()
    if mask is not None:
        display_orig[~mask] = np.nan
        display_resid[~mask] = np.nan

    # 1. Original
    plt.subplot(131)
    plt.imshow(display_orig, extent=extent, origin='lower', norm=norm, cmap='gray')
    plt.title(f'{title_prefix} Original')
    plt.colorbar()

    # 2. Model
    plt.subplot(132)
    plt.imshow(model, extent=extent, origin='lower', norm=norm, cmap='gray')
    plt.title(f'{title_prefix} Model')
    plt.colorbar()

    # 3. Residual
    plt.subplot(133)
    plt.imshow(display_resid, extent=extent, origin='lower', norm=norm, cmap='gray')
    metric_str = f' (Metric: {metric:.4f})' if metric is not None else ''
    plt.title(f'{title_prefix} Residual{metric_str}')
    plt.colorbar()

    plt.tight_layout()
    if output_filename:
        plt.savefig(output_filename, dpi=150)
        print(f"Saved figure to {output_filename}")
    else:
        plt.show()
    plt.close()

def plot_contours(image: np.ndarray, all_contours: List[List[np.ndarray]], selected_contours: List[np.ndarray],
                  extent: Optional[List[float]] = None, 
                  mask: Optional[np.ndarray] = None,
                  title: str = "Contour Selection",
                  output_filename: Optional[str] = None,
                  profile_data: Optional[dict] = None):
    """
    Plot the image with all candidate contours and selected contours overlayed.
    Optionally plot a 1D surface brightness profile.
    
    Args:
        image: 2D image array.
        all_contours: List of lists of candidate contours.
        selected_contours: List of selected contours.
        extent: Physical extent of the image.
        mask: Optional mask.
        title: Plot title.
        output_filename: Output filename.
        profile_data: Dictionary with 'r' and 'val' arrays for 1D profile.
    """
    if profile_data is not None:
        fig = plt.figure(figsize=(12, 6))
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 1])
        ax_img = fig.add_subplot(gs[0])
        ax_prof = fig.add_subplot(gs[1])
    else:
        fig = plt.figure(figsize=(10, 10))
        ax_img = fig.add_subplot(111)
    
    # --- Image Plot ---
    norm = simple_norm(image, 'asinh', percent=99)
    
    im = ax_img.imshow(image, extent=extent, origin='lower', norm=norm, cmap='gray_r')
    plt.colorbar(im, ax=ax_img, label='Intensity')
    
    # Overlay mask if present
    if mask is not None:
        ax_img.contour(mask, levels=[0.5], extent=extent, colors='yellow', linestyles='--')

    # Plot ALL contours in faint blue
    for level_contours in all_contours:
        for c in level_contours:
            ax_img.plot(c[:, 0], c[:, 1], color='cyan', alpha=0.1, lw=0.5)
            
    # Plot SELECTED contours in red
    for c in selected_contours:
        ax_img.plot(c[:, 0], c[:, 1], color='red', alpha=0.8, lw=1.0)
        
    legend_elements = [
        Line2D([0], [0], color='cyan', alpha=0.5, label='All Candidates'),
        Line2D([0], [0], color='red', label='Selected Isophotes')
    ]
    if mask is not None:
        legend_elements.append(Line2D([0], [0], color='yellow', linestyle='--', label='Mask Boundary'))
        
    ax_img.legend(handles=legend_elements, loc='upper right')
    ax_img.set_title(title)
    
    # --- 1D Profile Plot ---
    if profile_data is not None:
        r = profile_data['r']
        val = profile_data['val']
        
        ax_prof.plot(r, val, 'o-', color='red', markersize=3, label='Extracted Profile')
        ax_prof.set_xlabel('Generalized Radius (pixels)')
        ax_prof.set_ylabel('Intensity')
        ax_prof.set_yscale('log')
        ax_prof.set_title('1D Surface Brightness Profile')
        ax_prof.grid(True, which="both", ls="-", alpha=0.2)
        ax_prof.legend()

    plt.tight_layout()
    
    if output_filename:
        plt.savefig(output_filename, dpi=150)
        print(f"Saved contour QA to {output_filename}")
    else:
        plt.show()
    plt.close()

def plot_comprehensive_qa(original: np.ndarray, model: np.ndarray, residual: np.ndarray, 
                          contours: List[np.ndarray], profile_data: dict,
                          extent: Optional[List[float]] = None, 
                          mask: Optional[np.ndarray] = None,
                          title_prefix: str = "", 
                          output_filename: Optional[str] = None):
    """
    Generate a comprehensive 2x2 QA dashboard.
    
    Layout:
    - Top-Left: Original Image (arcsinh scaled)
    - Top-Right: Model Image (arcsinh scaled) + Selected Contours
    - Bottom-Left: Residual Image
    - Bottom-Right: 1D Surface Brightness Profile (r^0.25)
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Common font sizes
    TITLE_SIZE = 12
    LABEL_SIZE = 10
    
    # --- 1. Original Image (TL) ---
    ax_orig = axes[0, 0]
    # Use arcsinh stretching
    norm_orig = simple_norm(original, 'asinh', percent=99)
    im_orig = ax_orig.imshow(original, extent=extent, origin='lower', norm=norm_orig, cmap='gray')
    ax_orig.set_title(f'{title_prefix} Original', fontsize=TITLE_SIZE)
    plt.colorbar(im_orig, ax=ax_orig, fraction=0.046, pad=0.04)
    
    # --- 2. Model Image + Contours (TR) ---
    ax_model = axes[0, 1]
    norm_model = simple_norm(model, 'asinh', percent=99)
    im_model = ax_model.imshow(model, extent=extent, origin='lower', norm=norm_model, cmap='gray')
    ax_model.set_title(f'{title_prefix} Model + Contours', fontsize=TITLE_SIZE)
    plt.colorbar(im_model, ax=ax_model, fraction=0.046, pad=0.04)
    
    # Overplot contours
    for c in contours:
        ax_model.plot(c[:, 0], c[:, 1], color='red', alpha=0.5, lw=0.8)
    
    # Legend for contours
    legend_elements = [Line2D([0], [0], color='red', lw=1, label='Isophotes')]
    ax_model.legend(handles=legend_elements, loc='upper right', fontsize='small')

    # --- 3. Residual Image (BL) ---
    ax_resid = axes[1, 0]
    # Use arcsinh stretching and percentile clipping (99%) to handle outliers
    norm_resid = simple_norm(residual, 'asinh', percent=99)
    im_resid = ax_resid.imshow(residual, extent=extent, origin='lower', norm=norm_resid, cmap='gray')
    ax_resid.set_title(f'{title_prefix} Residual', fontsize=TITLE_SIZE)
    plt.colorbar(im_resid, ax=ax_resid, fraction=0.046, pad=0.04)
    
    # --- 4. 1D Profile (BR) ---
    ax_prof = axes[1, 1]
    r = profile_data['r']
    val = profile_data['val']
    
    # X-axis: r^0.25
    x_axis = r**0.25
    
    ax_prof.plot(x_axis, val, 'o-', color='blue', markersize=3, lw=1, label='Measured Profile')
    ax_prof.set_xlabel(r'Generalized Radius $R^{1/4}$ (pix$^{1/4}$)', fontsize=LABEL_SIZE)
    ax_prof.set_ylabel('Intensity', fontsize=LABEL_SIZE)
    ax_prof.set_yscale('log')
    
    # Limit X-axis range if extent is provided
    if extent is not None:
        # Calculate half-size of image
        half_size = max(abs(extent[1] - extent[0]), abs(extent[3] - extent[2])) / 2.0
        max_r = half_size * 1.3
        ax_prof.set_xlim(left=0, right=max_r**0.25)
        
    ax_prof.set_title('Surface Brightness Profile', fontsize=TITLE_SIZE)
    ax_prof.grid(True, which="both", ls="-", alpha=0.2)
    ax_prof.legend(fontsize='small')

    plt.tight_layout()
    
    if output_filename:
        plt.savefig(output_filename, dpi=150)
        print(f"Saved QA dashboard to {output_filename}")
    else:
        plt.show()
    plt.close()
