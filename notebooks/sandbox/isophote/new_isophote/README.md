# Isophote Optimization Walkthrough

I have successfully optimized the isophote analysis code by refactoring it into a streamlined, function-based structure and vectorizing the core integration steps.

## Changes Implemented

### 1. Vectorized Core (`isophote_optimize.py`)
- **`get_elliptical_coordinates`**: Helper to convert image coordinates to elliptical coordinates (used in testing).
- **`extract_isophote_data`**: Extracts pixels along the elliptical path using `scipy.ndimage.map_coordinates`. This scales linearly with the ellipse circumference $O(N)$, providing significant speedups over area-based or pixel-looping methods.
- **`fit_isophote`**: Implements the iterative harmonic fitting and geometry update logic in a single function.
- **`fit_image`**: Driver function that manages the outward and inward fitting loops.

### 2. CLI Interface (`run_isophote.py`)
- A command-line script to run the analysis.
- Accepts FITS images and YAML configuration files.
- Outputs results to a CSV or FITS table.

### 3. Configuration (`config.yaml`)
- A sample YAML file demonstrating how to configure the analysis (geometry guesses, fitting parameters).

## Verification & Benchmarking

I verified the correctness and performance using `benchmark_isophote.py`.

- **Speedup**:
    - **Synthetic Image (500x500)**: ~10x faster than `photutils` (path-based sampling).
    - **M51 Example**: ~2x faster than `photutils`.
- **Accuracy**:
    - **Noise-free data**: Near-perfect (< 0.001 fractional error)
    - **Noisy data (σ=1)**: Excellent (< 0.001 fractional error for SMA > 2 pixels)
    - Both implementations show instability at SMA < 2 pixels with noise where signal-to-noise is very low

## Known Limitations

- **Small SMA with Noise**: At SMA < 2 pixels with noisy data, the fitting becomes unstable (ellipticity corrections can drive values to clamp limits). This behavior matches `photutils` and appears inherent to the method rather than an implementation issue.
- **Central Pixel**: The SMA=0 fit is implemented via bilinear interpolation but provides limited geometric information.

## Usage

To run the analysis on an image:

```bash
python run_isophote.py image.fits --config config.yaml --output results.csv
```

To run the verification test:

```bash
python test_optimization.py
```
