# Confeti: Contour Fitting for EllipTical Isophotes

**Confeti** is a Python tool for empirical isophote fitting of galaxy images. Unlike traditional methods that fit intensity profiles along fixed geometric paths, Confeti extracts actual iso-intensity contours from the image and fits generalized ellipse models to them. This allows it to handle complex, varying galaxy morphologies (boxy, disky, twisting isophotes) robustly.

## Algorithm Overview

The Confeti algorithm proceeds in four main steps:

### 1. Isophote Level Selection
The algorithm first determines which intensity levels to fit.
-   It analyzes the distribution of pixel values in the masked image.
-   It selects a set of levels (default 100) that logarithmically cover the dynamic range of the galaxy, from the bright core to the faint outskirts.
-   This ensures that the fitting captures the full structural detail of the object.

### 2. Contour Extraction
For each selected intensity level:
-   **Extraction**: Uses the marching squares algorithm (`skimage.measure.find_contours`) to find all contours at that level.
-   **Filtering**: 
    -   Removes contours that are too short (noise).
    -   Removes contours that touch the image boundaries.
    -   **Mask Handling**: Explicitly filters out contours that trace the boundaries of masked regions, ensuring that missing data does not corrupt the fit.
-   **Selection**: Selects the "best" contour, typically the one that encloses the center of the previous fit or is the longest.

### 3. Ellipse Fitting (Optimization)
For each selected contour, Confeti fits a generalized ellipse model.
-   **Model**: The generalized ellipse is defined by center $(x_0, y_0)$, axis ratio $q$, position angle $\theta$, and boxiness parameter $c_0$.
    $$ (|x'|)^{c_0+2} + \left(\frac{|y'|}{q}\right)^{c_0+2} = r^{c_0+2} $$
-   **Objective**: The optimizer minimizes the variance of the generalized radius $r$ for all points on the contour. If the contour is a perfect ellipse, $r$ is constant.
-   **Optimization**: Uses the Nelder-Mead simplex algorithm to find the parameters that best describe the contour shape.

### 4. Model Reconstruction
Once all levels are fitted, Confeti builds a continuous 2D model of the galaxy.
-   **Interpolation**: Creates continuous functions (interpolators) for each parameter ($x_0, y_0, q, \theta, c_0, Intensity$) as a function of radius $r$.
-   **Iterative Radius Mapping**: To generate the model image, we need to find the "elliptical radius" of every pixel. Since the geometry changes with radius, this is solved iteratively:
    1.  Guess $r$ (circular distance).
    2.  Look up geometry parameters for this $r$.
    3.  Calculate new generalized radius $r'$.
    4.  Repeat until convergence.
-   **Intensity Lookup**: The final pixel value is determined by looking up the intensity at the converged radius.

## Usage

```python
import confeti
from astropy.io import fits

# Load data
data = fits.getdata('galaxy.fits')

# Initialize Fitter
fitter = confeti.Confeti(data)

# Run Fit
# mode='c0' enables boxiness parameter fitting
fitter.fit(mode='c0')

# Generate Model
model = fitter.model_image()

# Get Results Table
table = fitter.get_isophotes()
table.write('isophotes.csv')
```

## Features

-   **Generalized Ellipses**: Supports boxy ($c_0 > 0$) and disky ($c_0 < 0$) shapes.
-   **Fourier Modes**: Can optionally fit higher-order Fourier perturbations ($a_3, b_3, a_4, b_4$).
-   **Robust Masking**: Handles masked regions gracefully without artifacts.
-   **QA Visualization**: Includes tools to visualize candidate contours and fitting quality.
