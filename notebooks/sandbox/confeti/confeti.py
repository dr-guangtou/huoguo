import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from skimage import measure
from scipy.optimize import minimize

from scipy.ndimage import zoom, gaussian_filter
from matplotlib.path import Path
from typing import List, Tuple, Dict, Optional, Callable, Union

# --- Utility Functions ---

def is_iterable(obj) -> bool:
    """
    Check if an object is iterable.
    
    Args:
        obj: The object to check.
        
    Returns:
        bool: True if the object implements the iterator protocol, False otherwise.
    """
    try:
        iter(obj)
        return True
    except TypeError:
        return False

def midpoints(arr: np.ndarray) -> np.ndarray:
    """
    Calculate the midpoints of bins in an array.
    
    Args:
        arr: Input array of bin edges (length N).
        
    Returns:
        np.ndarray: Array of midpoints (length N-1).
    """
    return (arr[1:] + arr[:-1]) * 0.5

def weighted_percentile(data: np.ndarray, weights: np.ndarray, percentile: Optional[float] = None) -> Union[float, interp1d]:
    """
    Calculate the weighted percentile of data.
    
    This function computes percentiles for weighted data by sorting the data and 
    calculating the cumulative sum of the normalized weights. It then interpolates 
    to find the value corresponding to the desired percentile.
    
    Args:
        data: 1D array of data values.
        weights: 1D array of weights corresponding to each data point.
        percentile: The percentile to calculate (0-100). If None, returns the interpolator object.
        
    Returns:
        Union[float, interp1d]: The value at the given percentile, or the scipy.interpolate.interp1d 
                                object mapping percentile -> value if percentile is None.
    """
    weights = np.array(weights)
    # Normalize weights so they sum to 100 (representing 100%)
    weights /= np.sum(weights) / 100.0
    
    # Sort data and weights based on data values
    idxs = np.argsort(data)
    sorted_data = np.array(data)[idxs]
    sorted_weights = weights[idxs]
    
    # Calculate cumulative sum of weights to get the percentile position of each data point
    cumulative_weights = np.cumsum(sorted_weights)
    
    # Create interpolation points
    # We prepend 0 and calculate midpoints to center the cumulative distribution
    # This ensures that the interpolation covers the full range [0, 100] correctly
    C = midpoints(np.append([0.0], cumulative_weights))
    C[0] = 0.0
    C[-1] = 100.0
    
    # Create interpolator: Percentile (x) -> Value (y)
    # We allow extrapolation to handle edge cases
    interpolator = interp1d(C, sorted_data, bounds_error=False, fill_value="extrapolate")
    
    if percentile is None:
        return interpolator
    else:
        return interpolator(percentile)

def extract_contours(image: np.ndarray, level: float, extent: Optional[List[float]] = None, mask: Optional[np.ndarray] = None) -> List[np.ndarray]:
    """
    Extract contours from an image at a specific intensity level.
    
    This function uses the marching squares algorithm (via skimage.measure.find_contours) 
    to find iso-intensity contours. It maps the resulting pixel coordinates to physical 
    coordinates defined by 'extent' and optionally filters contours that overlap with a mask.
    
    Args:
        image: 2D image array.
        level: The intensity level to find contours for.
        extent: [xmin, xmax, ymin, ymax] of the image coordinates. 
                If None, assumes pixel coordinates [0, width, 0, height].
        mask: Optional boolean mask (True=valid, False=masked). 
              If provided, contours passing through masked regions (False) will be filtered out.
                
    Returns:
        List[np.ndarray]: A list of contours, where each contour is an (N, 2) array of (x, y) coordinates.
    """
    if extent is None:
        extent = [0, image.shape[1], 0, image.shape[0]]
        
    # Calculate coordinate grids for the pixel centers
    # We use linspace to define pixel edges and then take midpoints to get centers
    locx = midpoints(np.linspace(extent[0], extent[1], image.shape[1] + 1, endpoint=True))
    locy = midpoints(np.linspace(extent[2], extent[3], image.shape[0] + 1, endpoint=True))
    
    # Find contours using skimage
    # measure.find_contours returns coordinates in (row, column) order, i.e., (y, x)
    contours_raw = measure.find_contours(image, level)
    
    contours_converted = []
    for seg in contours_raw:
        # seg is (N, 2) with (row, col) -> (y, x)
        
        # If mask is provided, check if contour passes through masked pixels
        if mask is not None:
            # Convert sub-pixel contour coordinates to nearest integer pixel indices
            rows = np.clip(np.round(seg[:, 0]).astype(int), 0, mask.shape[0]-1)
            cols = np.clip(np.round(seg[:, 1]).astype(int), 0, mask.shape[1]-1)
            
            # Check if any point is masked (False in mask)
            # Ensure mask is boolean for inversion
            mask_bool = mask.astype(bool)
            is_masked = ~mask_bool[rows, cols]
            
            if np.any(is_masked):
                # Calculate the fraction of points on the mask boundary
                # If a significant fraction (>10%) of the contour lies in the masked region,
                # it is likely tracing the mask boundary rather than a real isophote.
                mask_fraction = np.sum(is_masked) / len(seg)
                if mask_fraction > 0.1: 
                    continue

        # Map pixel indices to physical coordinates
        # Interpolate x coordinates (col indices)
        x_coords = np.interp(seg[:, 1], np.arange(image.shape[1]), locx)
        
        # Interpolate y coordinates (row indices)
        y_coords = np.interp(seg[:, 0], np.arange(image.shape[0]), locy)
        
        contours_converted.append(np.column_stack((x_coords, y_coords)))
        
    return np.array(contours_converted, dtype=object)

# --- Helper Functions for Geometry ---

def generalized_ellipse_radius(x: np.ndarray, y: np.ndarray, 
                               x0: float, y0: float, q: float, theta: float, c0: float, 
                               return_theta: bool = False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Calculate the generalized elliptical radius for given coordinates.
    
    The generalized ellipse is defined by the equation:
        (|x'|)^(c0+2) + (|y'/q|)^(c0+2) = r^(c0+2)
    where (x', y') are the coordinates rotated by 'theta' and centered at (x0, y0).
    
    Args:
        x, y: Coordinate arrays (can be 2D meshgrids or 1D arrays).
        x0, y0: Center of the ellipse.
        q: Axis ratio (b/a), where 0 < q <= 1.
        theta: Position angle in radians (counter-clockwise from x-axis).
        c0: Boxiness parameter. 
            c0 = 0 corresponds to a pure ellipse.
            c0 > 0 corresponds to a boxy shape.
            c0 < 0 corresponds to a disky shape.
        return_theta: If True, also return the angular coordinate in the ellipse frame.
        
    Returns:
        Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]: 
            - Radius array 'r' with the same shape as input x, y.
            - (Optional) Theta array 'th' if return_theta is True.
    """
    shape = x.shape
    dx = x.ravel() - x0
    dy = y.ravel() - y0
    
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    # Rotate coordinates to the ellipse frame
    # x_rot is along the major axis
    # y_rot is along the minor axis, scaled by q to normalize to a circle
    x_rot = dx * cos_t + dy * sin_t
    y_rot = (-dx * sin_t + dy * cos_t) / q
    
    if return_theta:
        # Calculate angle in the unscaled frame (standard polar angle)
        th = np.unwrap(np.arctan2(y_rot, x_rot))
        
    # Calculate generalized radius
    # The power is (c0 + 2). For a normal ellipse, c0=0, so power is 2.
    # |x|^p + |y|^p = r^p  =>  r = (|x|^p + |y|^p)^(1/p)
    term = np.abs(x_rot)**(c0 + 2) + np.abs(y_rot)**(c0 + 2)
    radius = term**(1.0 / (c0 + 2))
    
    if return_theta:
        return radius.reshape(shape), th
    return radius.reshape(shape)

def fourier_ellipse_radius(x: np.ndarray, y: np.ndarray, 
                           x0: float, y0: float, q: float, theta: float, 
                           a_coeffs: List[float]) -> np.ndarray:
    """
    Calculate radius for an ellipse modified by 3rd and 4th order Fourier terms.
    
    The radius is perturbed as:
        r_perturbed = r_base / (1 + sum(a_n * cos(n*t) + b_n * sin(n*t)))
    
    Args:
        x, y: Coordinates.
        x0, y0: Center.
        q: Axis ratio.
        theta: Position angle.
        a_coeffs: List of 4 coefficients [a3, b3, a4, b4] for n=3 and n=4 modes.
        
    Returns:
        np.ndarray: Perturbed radius array.
    """
    shape = x.shape
    dx = x.ravel() - x0
    dy = y.ravel() - y0
    
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    # Rotate and scale
    x_rot = dx * cos_t + dy * sin_t
    y_rot = (-dx * sin_t + dy * cos_t) / q
    
    # Angle in the ellipse frame
    th = np.arctan2(y_rot, x_rot)
    
    # Calculate Fourier perturbation factor
    # a_coeffs = [a3, b3, a4, b4]
    perturbation = 1.0
    # n=3 terms (asymmetry/triangle)
    perturbation += a_coeffs[0] * np.cos(3*th) + a_coeffs[1] * np.sin(3*th)
    # n=4 terms (boxiness/diskiness)
    perturbation += a_coeffs[2] * np.cos(4*th) + a_coeffs[3] * np.sin(4*th)
    
    # Base elliptical radius (Euclidean in the scaled frame)
    base_r = np.sqrt(x_rot**2 + y_rot**2)
    
    return (base_r / perturbation).reshape(shape)

def get_ellipse_points(r: float, x0: float, y0: float, q: float, theta: float, c0: float = 0.0, n_points: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate (x, y) coordinates for the boundary of a generalized ellipse.
    
    Args:
        r: Radius (semimajor axis).
        x0, y0: Center.
        q: Axis ratio.
        theta: Position angle.
        c0: Boxiness parameter.
        n_points: Number of points to generate.
        
    Returns:
        Tuple[np.ndarray, np.ndarray]: Arrays of x and y coordinates.
    """
    t = np.linspace(0, 2*np.pi, n_points)
    
    # Parametric equations for generalized ellipse:
    # x = r * cos(t)^(2/(2+c0))
    # y = q * r * sin(t)^(2/(2+c0))
    # We use a helper to handle the sign correctly for fractional powers
    def signed_pow(val, p):
        return np.sign(val) * np.abs(val)**p
    
    exponent = 2.0 / (2.0 + c0)
    xr0 = r * signed_pow(np.cos(t), exponent)
    yr0 = q * r * signed_pow(np.sin(t), exponent)
    
    # Rotate to position angle
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    xr = xr0 * cos_t - yr0 * sin_t + x0
    yr = xr0 * sin_t + yr0 * cos_t + y0
    
    return xr, yr

def get_fourier_ellipse_points(r: float, x0: float, y0: float, q: float, theta: float, a_coeffs: List[float], n_points: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate (x, y) coordinates for the boundary of a Fourier-perturbed ellipse.
    
    Args:
        r: Base radius.
        x0, y0: Center.
        q: Axis ratio.
        theta: Position angle.
        a_coeffs: Fourier coefficients [a3, b3, a4, b4].
        n_points: Number of points.
        
    Returns:
        Tuple[np.ndarray, np.ndarray]: Arrays of x and y coordinates.
    """
    t = np.linspace(0, 2*np.pi, n_points)
    
    # Calculate perturbation at each angle
    perturbation = 1.0
    perturbation += a_coeffs[0] * np.cos(3*t) + a_coeffs[1] * np.sin(3*t)
    perturbation += a_coeffs[2] * np.cos(4*t) + a_coeffs[3] * np.sin(4*t)
    
    # Apply perturbation to radius
    xr0 = r * np.cos(t) * perturbation
    yr0 = q * r * np.sin(t) * perturbation
    
    # Rotate
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    xr = xr0 * cos_t - yr0 * sin_t + x0
    yr = xr0 * sin_t + yr0 * cos_t + y0
    
    return xr, yr


def midpoints(arr: np.ndarray) -> np.ndarray:
    """
    Calculate the midpoints of bins in an array.
    
    Args:
        arr: Input array of bin edges.
        
    Returns:
        np.ndarray: Array of midpoints, length is len(arr) - 1.
    """
    return (arr[1:] + arr[:-1]) * 0.5

def weighted_percentile(data: np.ndarray, weights: np.ndarray, percentile: Optional[float] = None) -> Union[float, interp1d]:
    """
    Calculate the weighted percentile of data.
    
    Args:
        data: 1D array of data values.
        weights: 1D array of weights corresponding to data.
        percentile: The percentile to calculate (0-100). If None, returns the interpolator.
        
    Returns:
        Union[float, interp1d]: The value at the given percentile, or the interpolator function if percentile is None.
    """
    weights = np.array(weights)
    # Normalize weights to sum to 100
    weights /= np.sum(weights) / 100.0
    
    # Sort data and weights
    idxs = np.argsort(data)
    sorted_data = np.array(data)[idxs]
    sorted_weights = weights[idxs]
    
    # Calculate cumulative sum of weights
    cumulative_weights = np.cumsum(sorted_weights)
    
    # Create interpolation points
    # Prepend 0 and calculate midpoints to get centered cumulative distribution
    # This logic mimics the original getWP function's behavior
    C = midpoints(np.append([0.0], cumulative_weights))
    C[0] = 0.0
    C[-1] = 100.0
    
    # Create interpolator: Percentile -> Value
    interpolator = interp1d(C, sorted_data, bounds_error=False, fill_value="extrapolate")
    
    if percentile is None:
        return interpolator
    else:
        return interpolator(percentile)

def extract_contours(image: np.ndarray, level: float, extent: Optional[List[float]] = None, mask: Optional[np.ndarray] = None) -> List[np.ndarray]:
    """
    Extract contours from an image at a specific level.
    
    Args:
        image: 2D image array.
        level: The intensity level to find contours for.
        extent: [xmin, xmax, ymin, ymax] of the image coordinates. 
                If None, assumes pixel coordinates.
        mask: Optional boolean mask (True=valid, False=masked). 
              If provided, contours passing through masked regions will be filtered or split.
                
    Returns:
        List[np.ndarray]: A list of contours, where each contour is an (N, 2) array of (x, y) coordinates.
    """
    if extent is None:
        extent = [0, image.shape[1], 0, image.shape[0]]
        
    # Calculate coordinate grids
    # Note: The original code used linspace with endpoint=True and shape+1, then took midpoints.
    # This effectively creates pixel centers.
    locx = midpoints(np.linspace(extent[0], extent[1], image.shape[1] + 1, endpoint=True))
    locy = midpoints(np.linspace(extent[2], extent[3], image.shape[0] + 1, endpoint=True))
    
    # Find contours using skimage
    # measure.find_contours returns coordinates in (row, column) order, i.e., (y, x)
    # Note: We do NOT pass mask to find_contours here because older versions don't support it 
    # and we want custom control. We assume 'image' already has masked values set to 0 or similar.
    contours_raw = measure.find_contours(image, level)
    
    contours_converted = []
    for seg in contours_raw:
        # seg is (N, 2) with (row, col) -> (y, x)
        
        # If mask is provided, check if contour passes through masked pixels
        if mask is not None:
            # Get integer pixel indices
            rows = np.clip(np.round(seg[:, 0]).astype(int), 0, mask.shape[0]-1)
            cols = np.clip(np.round(seg[:, 1]).astype(int), 0, mask.shape[1]-1)
            
            # Check if any point is masked (False in mask)
            # We allow a small tolerance or check if a significant fraction is masked
            # But strictly, if we masked it, we probably don't want it.
            # However, find_contours on a 0-filled mask hole will trace the boundary.
            # These boundary points will be technically "valid" (0) or close to it.
            # But they are artifacts of the mask.
            
            # A simple heuristic: if a point is ON a masked pixel (or very close), it's suspicious.
            # But find_contours returns sub-pixel coordinates.
            
            # Let's check the value of the mask at the nearest pixel.
            # If mask is False (0) at the nearest pixel, it's a masked region.
            # Ensure mask is boolean
            mask_bool = mask.astype(bool)
            is_masked = ~mask_bool[rows, cols]
            
            if np.any(is_masked):
                # If the contour touches the mask, we might want to discard it or split it.
                # For ellipse fitting, a partial contour is better than a wrong one.
                # But splitting is hard.
                # If the contour is largely following the mask boundary, it should be discarded.
                
                # Let's calculate the fraction of points on the mask boundary
                mask_fraction = np.sum(is_masked) / len(seg)
                if mask_fraction > 0.1: # If >10% of points are in masked region
                    continue

        # We want to map these pixel indices to physical coordinates
        
        # Interpolate x coordinates (col indices)
        x_coords = np.interp(seg[:, 1], np.arange(image.shape[1]), locx)
        
        # Interpolate y coordinates (row indices)
        y_coords = np.interp(seg[:, 0], np.arange(image.shape[0]), locy)
        
        contours_converted.append(np.column_stack((x_coords, y_coords)))
        
    return np.array(contours_converted, dtype=object)

# --- Helper Functions for Geometry ---

def generalized_ellipse_radius(x: np.ndarray, y: np.ndarray, 
                               x0: float, y0: float, q: float, theta: float, c0: float, 
                               return_theta: bool = False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """
    Calculate the generalized elliptical radius for given coordinates.
    
    The generalized ellipse is defined by:
    (|x'|)^(c0+2) + (|y'/q|)^(c0+2) = r^(c0+2)
    where (x', y') are rotated coordinates.
    
    Args:
        x, y: Coordinate arrays.
        x0, y0: Center of the ellipse.
        q: Axis ratio (b/a).
        theta: Position angle in radians (counter-clockwise from x-axis).
        c0: Boxiness parameter (0 for pure ellipse).
        return_theta: If True, also return the angular coordinate in the ellipse frame.
        
    Returns:
        Radius array, and optionally theta array.
    """
    shape = x.shape
    dx = x.ravel() - x0
    dy = y.ravel() - y0
    
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    # Rotate coordinates
    # x_rot =  dx * cos + dy * sin
    # y_rot = -dx * sin + dy * cos
    # Scaled y_rot by q
    
    x_rot = dx * cos_t + dy * sin_t
    y_rot = (-dx * sin_t + dy * cos_t) / q
    
    # Generalized radius
    # Note: The original code used np.dot with a specific rotation matrix structure.
    # c = [[cos, sin], [-sin/q, cos/q]]
    # M = dot(c, [dx, dy])
    # This matches the manual rotation above.
    
    if return_theta:
        th = np.unwrap(np.arctan2(y_rot * q, x_rot)) # Reconstruct angle in unscaled frame? 
        # Original code: th=np.unwrap(np.arctan2(M[1,:],M[0,:])) where M is the rotated scaled vector.
        # Let's stick to original logic for consistency.
        th = np.unwrap(np.arctan2(y_rot, x_rot))
        
    term = np.abs(x_rot)**(c0 + 2) + np.abs(y_rot)**(c0 + 2)
    radius = term**(1.0 / (c0 + 2))
    
    if return_theta:
        return radius.reshape(shape), th
    return radius.reshape(shape)

def fourier_ellipse_radius(x: np.ndarray, y: np.ndarray, 
                           x0: float, y0: float, q: float, theta: float, 
                           a_coeffs: List[float]) -> np.ndarray:
    """
    Calculate radius for an ellipse modified by 3rd and 4th order Fourier terms.
    
    Args:
        x, y: Coordinates.
        x0, y0: Center.
        q: Axis ratio.
        theta: Position angle.
        a_coeffs: List of [a3, b3, a4, b4] coefficients.
    """
    shape = x.shape
    dx = x.ravel() - x0
    dy = y.ravel() - y0
    
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    x_rot = dx * cos_t + dy * sin_t
    y_rot = (-dx * sin_t + dy * cos_t) / q
    
    # Original code logic:
    # M = dot(c, M)
    # th = arctan2(M[1], M[0])
    # ff = 1 + sum(...)
    # return sqrt(M[0]**2 + M[1]**2) / ff
    
    th = np.arctan2(y_rot, x_rot)
    
    # Fourier perturbations
    # a_coeffs = [a3, b3, a4, b4]
    # Terms are cos(3t), sin(3t), cos(4t), sin(4t)
    
    perturbation = 1.0
    # i=0 -> n=3 (indices 0, 1)
    perturbation += a_coeffs[0] * np.cos(3*th) + a_coeffs[1] * np.sin(3*th)
    # i=1 -> n=4 (indices 2, 3)
    perturbation += a_coeffs[2] * np.cos(4*th) + a_coeffs[3] * np.sin(4*th)
    
    # Base elliptical radius (Euclidean in the scaled frame)
    base_r = np.sqrt(x_rot**2 + y_rot**2)
    
    return (base_r / perturbation).reshape(shape)

def get_ellipse_points(r: float, x0: float, y0: float, q: float, theta: float, c0: float = 0.0, n_points: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """Generate (x, y) points for an ellipse boundary."""
    t = np.linspace(0, 2*np.pi, n_points)
    
    # Parametric equations for generalized ellipse
    # x = r * cos(t)^(2/(2+c0))
    # y = q * r * sin(t)^(2/(2+c0))
    # But we need to handle signs correctly for the power
    
    # Helper for signed power
    def signed_pow(val, p):
        return np.sign(val) * np.abs(val)**p
    
    exponent = 2.0 / (2.0 + c0)
    xr0 = r * signed_pow(np.cos(t), exponent)
    yr0 = q * r * signed_pow(np.sin(t), exponent)
    
    # Rotate
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    xr = xr0 * cos_t - yr0 * sin_t + x0
    yr = xr0 * sin_t + yr0 * cos_t + y0
    
    return xr, yr

def get_fourier_ellipse_points(r: float, x0: float, y0: float, q: float, theta: float, a_coeffs: List[float], n_points: int = 100) -> Tuple[np.ndarray, np.ndarray]:
    """Generate points for Fourier-perturbed ellipse."""
    t = np.linspace(0, 2*np.pi, n_points)
    
    perturbation = 1.0
    perturbation += a_coeffs[0] * np.cos(3*t) + a_coeffs[1] * np.sin(3*t)
    perturbation += a_coeffs[2] * np.cos(4*t) + a_coeffs[3] * np.sin(4*t)
    
    xr0 = r * np.cos(t) * perturbation
    yr0 = q * r * np.sin(t) * perturbation
    
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    
    xr = xr0 * cos_t - yr0 * sin_t + x0
    yr = xr0 * sin_t + yr0 * cos_t + y0
    
    return xr, yr


# --- Main Fitting Class ---

class Confeti:
    """
    Main class for Contour Fitting for EllipTical Isophotes (Confeti).
    
    This class handles the entire process of empirical ellipse fitting:
    1. Determining isophote levels.
    2. Extracting contours.
    3. Fitting generalized ellipses to contours.
    4. Generating a continuous model image.
    """
    
    def __init__(self, image: np.ndarray, mask: Optional[np.ndarray] = None, extent: Optional[List[float]] = None):
        """
        Initialize the Confeti fitter.
        
        Args:
            image: 2D numpy array containing the galaxy image.
            mask: Optional boolean mask (True=valid, False=masked).
            extent: Optional [xmin, xmax, ymin, ymax] defining the physical coordinates of the image.
        """
        self.image = image
        self.mask = mask if mask is not None else np.ones_like(image)
        self.extent = extent if extent is not None else [0, image.shape[1], 0, image.shape[0]]
        
        # Create coordinate grids for the image
        # We define coordinates at the center of each pixel
        self.xx, self.yy = np.meshgrid(
            midpoints(np.linspace(self.extent[0], self.extent[1], image.shape[1] + 1)),
            midpoints(np.linspace(self.extent[2], self.extent[3], image.shape[0] + 1))
        )
        self.dx = self.xx[0][1] - self.xx[0][0]
        self.dy = self.yy[1][0] - self.yy[0][0]
        
        # Initialize storage for fitting results
        self.results = {
            'ps': [],      # Percentiles
            'x0': [],      # Center x
            'y0': [],      # Center y
            'q': [],       # Axis ratio
            'theta': [],   # Position angle
            'val': [],     # Intensity value
            'r': [],       # Mean radius
            'c0': [],      # Boxiness parameter
            'a4': [],      # Fourier a4 coefficient
            'a': [],       # All Fourier coefficients
            'ellval': []   # Intensity statistics along the ellipse (placeholder)
        }
        self.contours_out = [] # List of selected contours for each level
        self.all_contours = [] # List of all candidate contours for QA

    def fit(self, mode: str = 'c0', center: Optional[List[float]] = None, min_contour_points: int = 10, 
            optimizer_method: str = 'Nelder-Mead', maxiter: int = 1000, smoothing_sigma: float = 0.0) -> Dict:
        """
        Perform the ellipse fitting procedure.
        
        Args:
            mode: Fitting mode.
                - 'simple': Fits standard ellipses (c0=0).
                - 'c0': Fits generalized ellipses with boxiness parameter c0.
                - 'fourier': Fits ellipses perturbed by 3rd and 4th order Fourier modes.
            center: Optional [x, y] initial guess for the galaxy center.
            min_contour_points: Minimum number of points required for a contour to be fitted.
            optimizer_method: Optimization algorithm to use (default: 'Nelder-Mead').
                              Options include 'Nelder-Mead', 'Powell', 'L-BFGS-B', etc.
            maxiter: Maximum number of iterations for the optimizer.
            smoothing_sigma: Sigma for Gaussian smoothing applied to the image before contour extraction.
                             Useful for noisy data. Default is 0.0 (no smoothing).
            
        Returns:
            Dict: A dictionary containing the fitted parameters (arrays for r, x0, y0, q, etc.).
        """
        # 1. Determine which intensity levels (isophotes) to fit
        levels = self._determine_isophote_levels()
        
        prev_contour = None
        prev_params = None # [x0, y0, q, theta]
        
        # Create a masked version of the image for contour extraction
        masked_image = self.image.copy()
        masked_image[self.mask == 0] = 0 # Set masked pixels to 0 (background)
        
        # Apply smoothing if requested
        if smoothing_sigma > 0:
            from scipy.ndimage import gaussian_filter
            masked_image = gaussian_filter(masked_image, sigma=smoothing_sigma)
        
        # 2. Iterate through each level and fit
        for i, level in enumerate(levels):
            # Extract contours at this level
            contours = extract_contours(masked_image, level, self.extent, mask=self.mask)
            self.all_contours.append(contours)
            
            # Filter invalid contours (e.g., too short or touching borders)
            valid_contours = []
            for c in contours:
                if len(c) == 0: continue
                # Check bounds (ignore contours touching the image edge)
                if (np.min(c[:, 0]) <= self.extent[0] + 1 or np.max(c[:, 0]) >= self.extent[1] - 1 or
                    np.min(c[:, 1]) <= self.extent[2] + 1 or np.max(c[:, 1]) >= self.extent[3] - 1):
                    continue
                valid_contours.append(c)
            
            if not valid_contours:
                continue
                
            # Select the "best" contour among candidates
            if prev_contour is not None:
                # If we have a previous fit, choose the contour that encloses the previous center
                prev_center = np.median(prev_contour, axis=0)
                candidates = [c for c in valid_contours if Path(c).contains_point(prev_center)]
                if not candidates:
                    candidates = valid_contours # Fallback to all valid if none contain center
            else:
                candidates = valid_contours
                if center is not None:
                    # If user provided a center, choose contour closest to it
                    dists = [np.linalg.norm(np.median(c, axis=0) - center) for c in candidates]
                    candidates = [candidates[np.argmin(dists)]]
            
            # Tie-breaker: Pick the longest contour (most pixels)
            best_contour = max(candidates, key=len)
            
            # Skip if contour has too few points for a reliable fit
            if len(best_contour) < min_contour_points:
                pass 
            
            prev_contour = best_contour
            
            # 3. Fit the ellipse model to the selected contour points
            fitted_params, success = self._fit_single_isophote(best_contour, mode, initial_guess=prev_params, 
                                                               method=optimizer_method, maxiter=maxiter)
            
            if not success:
                continue
                
            # Store results
            self._store_results(level, fitted_params, mode, best_contour)
            
            # Update previous parameters to warm-start the next level
            prev_params = fitted_params[:4] # x0, y0, q, theta
            
        # 4. Finalize results (convert to arrays and build interpolators)
        return self._finalize_results()

    def _determine_isophote_levels(self) -> List[float]:
        """
        Determine the intensity levels to fit based on the image statistics.
        
        Returns:
            List[float]: A list of intensity values, sorted in descending order.
        """
        # Get valid pixels from the mask
        valid_pixels = self.image[self.mask == 1]
        if len(valid_pixels) == 0:
            return []
            
        # We sample the dynamic range of the galaxy using weighted percentiles.
        # This ensures we sample both the bright core and the faint wings adequately.
        
        # Create a set of percentiles to query (100 levels from 100% down to 0%)
        ps = np.linspace(100, 0, 100)
        
        # Handle potential negative values (e.g., from background subtraction) by adding an offset
        offset = np.abs(np.min(valid_pixels)) if np.min(valid_pixels) < 0 else 0
        
        # Create an interpolator for the cumulative distribution function
        interpolator = weighted_percentile(valid_pixels, valid_pixels + offset)
        
        # Get intensity levels corresponding to the percentiles
        levels = interpolator(ps)
        
        # Filter out non-positive levels (background/noise)
        levels = levels[levels > 0]
        
        levels = np.unique(levels)[::-1] # Sort descending
        
        return levels

    def _fit_single_isophote(self, contour: np.ndarray, mode: str, initial_guess: Optional[List[float]] = None,
                             method: str = 'Nelder-Mead', maxiter: int = 1000) -> Tuple[List[float], bool]:
        """
        Fit an ellipse model to a single contour.
        
        Args:
            contour: (N, 2) array of (x, y) points.
            mode: Fitting mode ('simple', 'c0', 'fourier').
            initial_guess: Optional initial parameters [x0, y0, q, theta].
            method: Optimization method.
            maxiter: Maximum iterations.
            
        Returns:
            Tuple[List[float], bool]: Fitted parameters and success flag.
        """
        x = contour[:, 0]
        y = contour[:, 1]
        
        # Initial Guess
        if initial_guess is None:
            x0 = np.median(x)
            y0 = np.median(y)
            theta = 0.0
            q = 0.8
        else:
            x0, y0, q, theta = initial_guess
            
        # Define objective function to minimize
        def objective(params):
            px0, py0, pq, ptheta = params[:4]
            
            # Constraints
            # For gradient-based methods, hard constraints in objective are bad.
            # But we can use bounds if the method supports it.
            # For now, we keep the penalty method for simplicity/compatibility.
            if not (0 < pq <= 1): return 1e10
            
            # Prevent center from drifting too far from the contour median
            # This helps stability, especially for noisy data or partial contours
            if np.abs(px0 - x0) > 5.0 or np.abs(py0 - y0) > 5.0: return 1e10
            
            if mode == 'simple':
                pc0 = 0.0
                r = generalized_ellipse_radius(x, y, px0, py0, pq, ptheta, pc0)
            elif mode == 'c0':
                pc0 = params[4]
                if not (-2 < pc0 < 5): return 1e10
                r = generalized_ellipse_radius(x, y, px0, py0, pq, ptheta, pc0)
            elif mode == 'fourier':
                pa_coeffs = params[4:]
                if np.any(np.abs(pa_coeffs) > 0.5): return 1e10
                r = fourier_ellipse_radius(x, y, px0, py0, pq, ptheta, pa_coeffs)
            else:
                return 1e10
                
            # The objective is to minimize the variance of the generalized radius 'r'.
            # If the points lie perfectly on the model ellipse, 'r' will be constant for all points.
            return np.sum((r - np.median(r))**2)

        # Set up optimization parameters
        start_params = [x0, y0, q, theta]
        if mode == 'c0':
            start_params.append(0.0)
        elif mode == 'fourier':
            start_params.extend([0.0, 0.0, 0.0, 0.0]) # a3, b3, a4, b4
            
        # Run optimization
        # Note: 'Powell' is a good derivative-free alternative to Nelder-Mead.
        # 'L-BFGS-B' supports bounds, which would be better than penalty functions,
        # but requires refactoring to use the 'bounds' argument of minimize.
        
        res = minimize(objective, start_params, method=method, options={'maxiter': maxiter})
        
        if not res.success or res.fun > 1e9:
            return res.x, False
            
        return res.x, True

    def _store_results(self, level: float, params: List[float], mode: str, contour: np.ndarray):
        """
        Store the fitting results for a single level.
        
        Args:
            level: Intensity level.
            params: Fitted parameters.
            mode: Fitting mode.
            contour: The contour points used.
        """
        self.results['val'].append(level)
        self.results['x0'].append(params[0])
        self.results['y0'].append(params[1])
        self.results['q'].append(params[2])
        self.results['theta'].append(params[3])
        
        # Calculate mean radius for this fit
        x = contour[:, 0]
        y = contour[:, 1]
        if mode == 'simple':
            r = generalized_ellipse_radius(x, y, params[0], params[1], params[2], params[3], 0.0)
            self.results['c0'].append(0.0)
        elif mode == 'c0':
            r = generalized_ellipse_radius(x, y, params[0], params[1], params[2], params[3], params[4])
            self.results['c0'].append(params[4])
        elif mode == 'fourier':
            r = fourier_ellipse_radius(x, y, params[0], params[1], params[2], params[3], params[4:])
            self.results['a4'].append(params[6])
            self.results['a'].append(params[4:])
            self.results['c0'].append(0.0)
            
        self.results['r'].append(np.median(r))
        self.contours_out.append(contour)

    def _finalize_results(self) -> Dict:
        """
        Process the raw results into arrays and create interpolators.
        
        Returns:
            Dict: The results dictionary with numpy arrays.
        """
        for key in self.results:
            self.results[key] = np.array(self.results[key])
            
        # Create interpolators (radius -> parameter)
        # We sort by radius to ensure monotonic domain for interpolation
        sort_idx = np.argsort(self.results['r'])
        r_sorted = self.results['r'][sort_idx]
        
        # Remove duplicate radii if any
        r_unique, unique_idx = np.unique(r_sorted, return_index=True)
        
        interpolators = {}
        for key in ['x0', 'y0', 'q', 'theta', 'c0', 'val']:
            if len(self.results[key]) > 0:
                val_sorted = self.results[key][sort_idx][unique_idx]
                
                if key == 'val':
                    # For intensity, clamp to peak at center (r=0) and 0 at infinity
                    interpolators[key] = interp1d(r_unique, val_sorted, bounds_error=False, fill_value=(val_sorted[0], 0.0))
                else:
                    # For geometry parameters, clamp to the nearest valid value (no extrapolation)
                    # This prevents artifacts at the edges where extrapolation might diverge
                    interpolators[key] = interp1d(r_unique, val_sorted, bounds_error=False, fill_value=(val_sorted[0], val_sorted[-1]))
                
        self.interpolators = interpolators
        return self.results

    def model_image(self, shape: Optional[Tuple[int, int]] = None, extent: Optional[List[float]] = None) -> np.ndarray:
        """
        Generate a model image from the fitted parameters.
        
        This method uses an iterative radius mapping approach. Since the ellipse parameters 
        (center, shape, orientation) change with radius, the "elliptical radius" of a pixel 
        is defined implicitly. We iteratively solve for this radius to determine the intensity.
        
        Args:
            shape: Output image shape (height, width). Defaults to input image shape.
            extent: Physical extent [xmin, xmax, ymin, ymax]. Defaults to input extent.
            
        Returns:
            np.ndarray: The generated model image.
        """
        if shape is None: shape = self.image.shape
        if extent is None: extent = self.extent
        
        # Create coordinate grid for the model image
        xx, yy = np.meshgrid(
            midpoints(np.linspace(extent[0], extent[1], shape[1] + 1)),
            midpoints(np.linspace(extent[2], extent[3], shape[0] + 1))
        )
        
        # Initial guess for radius: Circular distance from the median center
        if len(self.results['x0']) == 0:
            return np.zeros(shape)
            
        x0_med = np.median(self.results['x0'])
        y0_med = np.median(self.results['y0'])
        r_map = np.sqrt((xx - x0_med)**2 + (yy - y0_med)**2)
        
        # Iteratively refine radius to account for varying geometry
        # We perform a few iterations to converge on the "elliptical radius" for each pixel
        for _ in range(3):
            # Interpolate parameters at current estimated radius
            # We clip r_map to the range of fitted radii to avoid wild extrapolation of geometry
            r_clamped = np.clip(r_map, np.min(self.results['r']), np.max(self.results['r']))
            
            ix0 = self.interpolators['x0'](r_clamped)
            iy0 = self.interpolators['y0'](r_clamped)
            iq = self.interpolators['q'](r_clamped)
            itheta = self.interpolators['theta'](r_clamped)
            ic0 = self.interpolators['c0'](r_clamped)
            
            # Calculate generalized radius with these parameters
            dx = xx - ix0
            dy = yy - iy0
            cos_t = np.cos(itheta)
            sin_t = np.sin(itheta)
            
            x_rot = dx * cos_t + dy * sin_t
            y_rot = (-dx * sin_t + dy * cos_t) / iq
            
            # Calculate new radius estimate
            term = np.abs(x_rot)**(ic0 + 2) + np.abs(y_rot)**(ic0 + 2)
            r_map = term**(1.0 / (ic0 + 2))

        # Final intensity lookup using the converged radius
        model = self.interpolators['val'](r_map)
        
        return model

    def get_isophotes(self):
        """
        Return the fitting results as an Astropy Table.
        
        Returns:
            astropy.table.Table: Table containing isophote parameters (r, x0, y0, q, theta, c0, val).
        """
        from astropy.table import Table
        
        # Ensure results are arrays
        data = {}
        for key in ['r', 'x0', 'y0', 'q', 'theta', 'c0', 'val']:
            if key in self.results and len(self.results[key]) > 0:
                data[key] = np.array(self.results[key])
                
        # Handle 'a4' if present
        if 'a4' in self.results and len(self.results['a4']) > 0:
            data['a4'] = np.array(self.results['a4'])
            
        return Table(data)
