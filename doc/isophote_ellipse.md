# Analysis of the `isophote.ellipse` IRAF code

----

## Basic Information

* The native IRAF language is **SPP (Subset Pre-Processor)**, a portable
pre-processor language which **resembles C**.  Tasks written in SPP will run on
all supported IRAF platforms without change.
* In addition, the **CL** provides a scripting language which can be used to
write high level tasks.
* [The **`CVOS`**](http://www.stsci.edu/institute/software_hardware/stsdas/cvos_stsci) is the C Interface to the Image Reduction and Analysis Facility (IRAF) Virtual Operating System libraries (VOS). The CVOS was developed by the Science Software Branch at the Space Telescope Science Institute.
    - [CVOS User's Guide and Reference Manual](http://www.stsci.edu/institute/software_hardware/stsdas/cvos/CVOS_UsersGuide.pdf)
    - [A C Interface to IRAF's VOS Library](http://www.stsci.edu/institute/software_hardware/stsdas/cvos/CVOSRefMan.pdf)
* The algorithm is closely based on the description given by [Jedrzejewski (Mon.Not.R.Astr.Soc., 226, 747, 1987)](https://ui.adsabs.harvard.edu/abs/1987MNRAS.226..747J/abstract).

## Long Term Goals:

* Improve the readability of the code.
* Use multiprocessing to speed things up.
    - At either the C or Python level.
* Use multiprocessing to fit multiple objects at the same time
* Use GPU to speed things up.
* Integrated profile extrapolation function.

## Contents of All Files

* `ellipse` takes FITS image as input, and return a STSDAS table.

### Top level

* `x_isophote.x`: defines tasks: `t_ellipse`, `t_isoexam`, `t_map`, `t_model`
    - These are the functions available in `x_isophote.e`
* `isophote.cl`: defines the `isophote` procedure and the list of tasks within.
    - Including the functions in `x_isophote.e`
    - And `isoplot.cl`, `isopall.cl`, `isomap.cl`, `isoimap.cl`, `bmodel.cl`
    - Also defint the list of parameters used: `geompar`, `controlpar`, `samplepar`, `magpar`.
* `isophote.hlp`: brief introduction of the `isophote` package.
    - `isophote.hd`: index of all documents.
* `mkpkg`: script to build `isophote` package.
    - `link	x_isophote.o isophote.a -lxtools -lstxtools -ltbtables -liminterp -o  xx_isophote.e`
    - `move	xx_isophote.e stsdasbin$x_isophote.e`

### Sub-tasks

* `isoexam.cl`: script to interactively plot ellipse on the image.
    - Based on `t_isoexam.x`.
    - Has its own parameter file `isoexam.par`.
        * Basic grammer: `isoexam(table, image, color, **interactive_mode)`
    - **Note**: Interactive fitting is an interesting feature, but it is not absolutely necessary. Can ignore at first.
* `isomap.cl`: script to draw ellipses taken from a `isophote` table over image contour.
    - Draw contours based on the `ellipse` result.
    - `map.par` is the parameter file.
* `isoimap.cl`: script to plot ellipses superposed over the gray-scale image display
* `isoplot.cl`: script to make plot between radius and one `ellipse` parameter.
* `isopall.cl`: script to make a summary plot of the `ellipese` results: radial plot of all important parameters.
    - **Note**: Above functions can be handled solely on the Python side as visualization functions.
* `bmodel.cl`: script to build a model image based on the isophotal analysis result from `ellipse`
    - Basic grammer: `bmodel(table, output, parent_image, minsma=1, maxsma=1, background=0., interp='spline', add_harmonics=False, verbose=False)`.
    - Support `nearest`, `linear`, `poly3`, and `spline` interpolation algorithms.
    - Use a lot of functions from the `Tables` package: `tfile`, `tcopy`, `tsort` (sort the table by `SMA`), `tinfo`, `tabpar` (get column), `tcalc` (normalize and correct PA values), `trebin` (interpolation).
    - [`trebin`](http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?trebin) is a key function here
    - **Note**: This is actually an important function. It involves interpolate the 1-D profile to a 2-D image.

### Parameter files

* **Note**: In Python, a `yaml` file can contains all these parameters. They can be converted into a multi-layer configuration dictionary in Python.
    - See `default.yaml`

* `controlpar.par`: sets the algorithm control parameters for the `ellipse` task
* `geompar.par`: sets the geometry of the isophote.
* `magpar.par`: sets the magnitude scale parameters.
* `samplepar.par`: sets the image sampling control parameters.

### `src` directory

#### Top-level scripts

* `mkpkg`: script that builds the `isophote` tasks.
* `t_ellipse.x`: the core of the `ellipse` procedure.
  - Depends on `imset.h`, `imhdr.h`, `imio.h`, `ctype.h`, `tbset.h`, `error.h`, `ellipse.h`, `eldisplay.h`, `elcolnames.h`
* `t_model.x`: the core of the `bmodel` procdure.
  - It takes a isophote table with equaly spaced SMA column with step < 1. pixel as input. The `SMA` column must be sorted too.
  - It works together with the `bmodel.cl` script that takes care of the sorting and interpolation.
* `t_isoexam.x`: the core of `isoexam` procedure.
* `t_map.x`: the core of `isomap` procedure.

#### Header files

* `ellipse.h`: defines the structure that holds the isophote parameters, all the control and sampling parameters.
    - Maximum harmonic level is 4.
    - Also defines stopping criteria. **Note**: put in `const.h`
* `elcolnames.h`: defines the column names, units, and formats of an output isophote table.
  - There are at least 40 columns (depending on whether to fit harmonics and how many levels are included). Maximum string length for column name is 10.
  - **Note**: In Python, this can be taken care of a pre-defined `astropy.table` object or a structured `numpy` array. But need to think about how to deal with multiple harmonic terms.
* `eldisplay.h`: defines parameters and files related to display image and interactive fitting.

#### Individual source files

##### Core

* `elfind.x`: finds the object center.
    - `el_find (img, is, al, sec, recenter, list, thresh)`
    - **Note**: this can be replaced by other algorithm

* `elfit.x`: fits one elliptical isophote on the image.
    - **Note**: core function including details of the fitting algorithm.
    - `el_fit (im, sec, is, al, file)`

* `elget.x`: samples the image array over an elliptical path and returns the azimutal angles and the sampled intensities at the corresponding angles with the mean value subtracted.
    - `el_get (im, sec, x0, y0, a, eps, theta, bufx, bufy, nbuffer, npoint, ndata, mean, sigma, astep, linear, integrmode, usclip, lsclip, nclip, aarea)`
    - `el_sarea (a, eps, phi, r)` computes the area of an elliptical sector.

* `elgetsec.x`: gets an image sub-raster.
    - `el_getsec (im, sec, i1, i2, j1, j2)`.

* `elharmonics.x`: decomposes data vectors in harmonic components.
    - `el_harmonics (bufx, bufy, ndata, coef, cerror, norder, rms, err)` 

* `elintegr.x`: computes integrated flux inside current ellipse.
    - `el_integr (im, sec, is)`

* `elcopy.x`: routines that handle isophote structure.
    - `el_copy (is, iso)` copies data from one isophote structure to another.
    - `el_alloc (nharm)` allocates structure for holding isophote parameters.
    - `el_free (is)` frees structure that holds isophote parameters.

##### Math

* `elbilinear.x`: bi-linear interpolation
    - `el_bilnear (image, pixel_pointer, sma, ip, jp, fx, fy)`
        - `ip, jp`: pixel address; `fx, fy`: fractional pixel position on pixel array.
        - Depends on `el_subpix()`
    
* `elmatrix.x`: algorithms dealing with matrix.
    - `el_matinv (array, norder, det)` inverts a symmetric matrix and calculates the determinant.
    - `el_matmul (array, column, norder)` multiplies a square matrix by a columnar matrix.
    - **Note**: maybe there is better way to do this.

* `elclip.x`: sigma-clipping function
    - `el_clip (data, usclip, lsclip, nclip, npoint, ndata, mean, sigma)`

* `elpolar.x`: converts image coordinate (x, y) into polar coordinate.
    - `el_polar (x, y, x0, y0, theta, radius, angle)`

* `elsubpix.x`: integrates a bi-linear interpolation on a 2x2 grid using sub-pixel resolution.
    - `el_subpix (z1, z2, z3, z4, fx, fy)`; `z1-z4` are the pixel contents; `fx, fy` are the fractional positions.

* `elqsortr.x`: re-orders the elements of the real array.
    - `el_qsortr (x, nelem, el_comparer)`

* `reorder.x`: re-orders table rows according to an index array. `reorder (tp, nindex, index)`
    - The algorithm used is taken from Knuth's Sorting and Searching p.595.
    - The index array comes from `tsort1` or `tsortm`.

##### IO and other helpers

* `elopens.x`: opens the main output table.
    - `el_opens (is, im, outname, tp, column, extracols, list)`

* `elopm.x`: opens the associated `.pl` pixel mask.
    - ` el_opm (im, sec, mode)`

* `elpflag.x`: flag or unflag a pixel.
    - `el_pflag (im, sec, i, j, interact)` flags a single pixel.
    - `el_punf (img, sec, i, j)` un-flags a single pixel.
    - `el_pmem (im, sec, i, j, value)` sets a pixel in a memory array.

* `elgetpix.x`: gets a single pixel from image
    - `el_getpix (im, sec, ip, jp)`

* `eldump.x`: writes isophote data to the main output table.
    - `el_dump (is, al, im, outname, tp, column, extracols, rownumber)`
* `eldsample.x`: dumps the azimutal sample of intensity as a binary table.
    - `el_dsample (file, is, im, al, bufx, bufy)`

* `elhead.x`: prints a STDOUT header. `el_head ()`
* `elimname.x`: adds image name to the table header. `el_imname (tp, im, is)`

* `elinrow.x`: reads in one row from the input `ellipse` data and converts it into section coordinate.
    - `el_inrow (im, is, al, tp, row)`

* `elrtable.x`: reads the input `ellipse` table and initializes isophote structures with data taken from the first and last line.
    - `el_rtable (im, is, al, ellname, tpin, inrow)`

* `elgeometry.x`: routines to translate coordinates between image (physical), section, and display systems.
    - `el_s2p (im, sec_coord, axis)` converts section coordinates to physical
    - `el_p2s (im, phys_coor, axis)` converts physical coordinates to section
    - Similarly: `el_s2d()`, `el_g2s()`, `el_s2g()`
    - `el_shrink()`, `el_unshrink()` adjust coordinate between physical and section space.

* `elshow.x`: shows the isophote parameters. `el_show (is, al, file)`

* `elzero.x`: deals with the special case of object center (SMA == 0).
    - `el_zero (im, sec, is, al, list)`

* `extnstr.x`: deals with strings.
    - `extnstr (str, n, substr)` extracts n-th sub-string from a list of dictionary.
    - `strnidx (ch, n, str)` returns the index of the n-th occurrence of a character in a string.

* `elodqf.x`: opens a data quality file (DQF). `el_odqf (im, sec, dqf)`
    - These are for HST observations and are used to mask out pixels. **Note**: not useful any more.

##### Display and interactive fitting

* `elgopen.x`: open appropriate graphic device
    - `el_gopen (device)`
* `eldisplay.x`: routines that display images.
    - `el_display (im, dp, is, sec, command, zoom, cx, cy)` is a general-purpose image display function.
    - `el_opend (im, dev)` opens image display structure.
    - `el_dimage (im, dp, is, sec, command, zoom, cx, cy)` displays the image itself.
    - `el_dscale` computes image display scaling parameters.
    - `el_mark (img, dp, sec, gp)` marks flagged pixels.

* `elintfit.x`: interactively fits one elliptical isophote.
    - `el_intfit (im, sec, dp, is, al, file, splot, sgraphics, inter, plot, list)`

* `elmask.x`: flags regions and plot it. 
    - `el_mask()`, `el_unmask()`, `el_limmsk()`

* `elplot.x`: plots the current ellipse stored in the isophote structure.
    - `el_plot (im, gp, dp, is)`
* `elppix.x`: plots a single pixel. `el_ppix (im, dp, gp, x, y)`
* `elpsample.x`: plots an azimuthal intensity sample. `el_psample (im, is, al, dev, x, y)`

* `elgraphics.x`: routines for plotting on the image.
* `elcursor.x`: `el_cursor(im, dp, x, y, key, command, sz_comm)` handels details of image curor input.
* `elcolon.x`: `el_colon()` processes cursor colon commands.

* `ipcolon.x`: processes cursor color commands.
* `ipdisplay.x`: displays image for the `isoexam` task.
* `ipdoit.x`: displays image and ellipses with cursor interaction.

##### Misc

* `elcurhelp.key`: cursor options for interactive isophote fitting
* `ipcurhelp.key`: cursor options for interactive isophote plotting

## Digest the `photutils.isophote` module

* `geometry.py`: provides a container class to store parameters for the geometry of an ellipse.
    - `_area(sma, eps, phi, r)`: gets the elliptical sector area.
    - `IN_MASK` and `OUT_MASK` used for finding the center of the galaxy in the beginning.
        * This can be replaced with something better.
    - `EllipseGeometry()` class: stores parameters for the geometry of an ellipse.
        * `__init__(self, x0, y0, sma, eps, pa, astep=0.1, linear_growth=False)`
        * `find_center(self, image, threshold=0.1, verbose=True)`: finds the center of a galaxy.
        * `radius(self, angle)`: convert polar angle into polar radius. 
        * `initialize_sector_geometry(self, phi)`: initializes geometry attributes associated with an elliptical sector at the given polar angle `phi`, including the four vertices that define the elliptical section, the sector area, and the angular width of the sector.
        * `bounding_ellipses(self)`: computes the semimajor axis of the two ellipses that bound the annulus where integrations take place.
        * `polar_angle_sector_limits(self)`: returns the two polar angles that bound the sector.
        * `to_polar(self, x, y)`: convert x, y coordinate to polar coordinate. There is a scalar (`_to_polar_scalar(self, x, y)`) and a vector (`_to_polar_vectorized(self, x, y)`) version.
        * `update_sma(self, step)`: calculates an updated value for the semimajor axis, given the current value and the step value.
        * `reset_sma(self. step)`: changes the direction of semimajor axis growth, from outwards to inwards.

* `harmonics.py`: provides tools for computing and fitting harmonic functions.
    * 没必要自己独立作为一个module, 计算部分可以放进数学包，优化部分可以和所有的优化工具放在一起.
    - `_least_squares_fit(optimize_func, parameters)`: use `scipy.optimize.leastsq` to optimize a function.
    - `first_and_second_harmonic_function(phi, c)`: simple equation to convert the harmonic coefficients into a correction value at given polar angle
    - `fit_first_and_second_harmonics(phi, intensities)`: fits the first and second harmonic function values to a set of (angle, intensity) pairs.
        * Defines a simple function to optimize over: the differences between the corrected value calculated by `first_and_second_harmonics()` and the images data.
        * Then uses `_least_squares_fit()` to optimize over it.
    - `fit_upper_harmonic(phi, intensities, order)`: similar, but for higher order harmonics.
        * For order 3 or 4, including both `an` and `bn` components.

* `integrator.py`:
    * 完全可以同时给出用不同方法得到的intensity，而且可以考虑得到沿着elliptical path的更多信息
    - `INTEGRATORS`: a dict that containts all integrators, connecting the class to a name.
    - `_Integrator` Class: supports different kinds of pixel integration methods.
        * `__init__(image, geometry, angles, radii, intensities)`.
            - `geometry`: `EllipseGeometry` instance, shape of the current ellipse.
            - `angles`, `radii`, `intensities`: list of output values extracted along the elliptical path.
        * `integrate(self, radius, phi)`: 最基本的method
        * `_reset()`: resets `_angles`, `_radii`, `_intensities`
        * `_store_results(self, phi, radius, sample)`: append values to the list; sample is for the `_intensities`.
        * `_get_polar_angle_step()`: returns the polar angle step used to walk over the elliptical path.
        * `get_sector_area()`: returns the area of elliptical sectors where the integration takes place.
        * `is_area()`: check if the integrator is an area integrator.
    - `_NearestNeighborIntegrator()` Class: inherits from `_Integrator`
        * Need to convert between polar and image coordinates. 
        * Just takes the pixel value. We never use this. 
    - `_BiLinearIntegrator()` Class: inherits from `_Integrator`.
        * Also needs to convert between polar and image coordinates.
        * Needs to check image boundaries and check if image is masked.
        * After that it is just a simple bi-linear interpolation.
    - `_AreaIntegrator()` Class: inherits from `_Integrator`
        * An auxiliary bilinear integrator to be used when sector areas contain a too small number of valid pixels
        * This is a little more complicated, need to study more, but probably is not necessary.
    - `_MeanIntegrator()` and `_MedianIntegrator()` Classes: inherits from `_Integrator`.
        * Pretty simple, has a `accumulate` method to gather pixel values. Then use a `compute_sample_value()` method to get the mean or median. This can be simplified.
    
* `sample.py`: provides a class to sample data along an elliptical path
    - `EllipseSample` Class: sample image data along an elliptical path
        * `__init__(self, image, sma, x0=None, y0=None, astep=0.1, eps=0.2, position_angle=0., sclip=3., nclip=0, linear_growth=False, integrmode='bilinear', geometry=None)`
            - If `geometry` is from somewhere else, change the `sma` too. This is not very good choice, should have a more consistent way to define ellipse.
            - This just defines a `EllipseGeometry` object.
            - Prepare `values`, `mean`, `gradient`, `gradient_error`, `gradient_relative_error`, `sector_area` attributes.
            - Also have `total_points` and `actual_points` attributes to report the total number of pairs angle-radius that were attempted and the actual number of samples used.
        * `_extract()`: just a wrapper; calls the `_extract(phi_min=0.05)` function to sample the elliptical path over the image array.
            - Defines the `integrator` and the initial polar radius and angle.
            - Walk along the elliptical path to get the intensities, radii, and angle list; also get the `sector_area`.
            - Perform sigma clipping as designed and compile the results into a 2-D array.
        * `_sigma_clip(angles, radii, intensities)`: sigma clipping function; calls `_iter_sigma_clip(angle, radii, intensities)` to do the job. 
            - Can use the `scipy` version for speed.
        * `update()`: update the `EllipseSample` instance after `extract`
            - Get the gradient and the error of the gradient; and check if the gradient is meaningful.
        * `_get_gradient(step)`: get the gradient and gradient error.
            - Define a new `EllipseSample` at a larger radius and calculate gradient. 
            - Does not feel like the most efficient approach.
        * `coordinates()`: convert radii and angles to (x, y). 
            - This definitely needs to become a common function.
    - `CentralEllipseSample(EllipseSample)`: special case for the central pixel.
        * `gradient` and `gradient_error` are `None`; just use the central pixel value.
    
- `isophote.py`:  provides classes to store the results of isophote fits.
    - Now depends on `astropy.table` and `units`. Can remove these dependences.
    - `Isophote` Class: stores the results of a single isophote.
        * `__init__(self, sample, niter, valid, stop_code)`
            - And prepare all the attributes that are going into the output. Similar to `Ellipse` standard output.
        * Need to have the capability to sort based on `sma` since the fitting starts from an intermediate radius. Using a `__lt__` method.
        * Properties: `sma`, `eps`, `pa`, `x0`, `y0`
        * `_compute_fluxes`: computes integrated fluxes in elliptical and circular apertures. 
            - Returns `tflux_e`, `tflux_c`, `npix_e`, `npix_c`
            - The "integrations" are done manually; can be improved. 
        * `_compute_deviations(sample, n)`: computes deviations from a perfect ellipse.
            - Based on `fit_first_and_second_harmonics` to calculate harmonic coefficients and `first_and_second_harmonic_function` to get the model and residual.
            - Then use `fit_upper_harmonic` to get the high-order harmonic values.
            - returns the `a`, `b`, `a_err`, `b_err` values for a given order `n`.
        * `_compute_errors()`: computes parameter errors based on the diagonal of the covariance matrix of the four harmonic coefficients for harmonics n=1 and n=2.
            - Get the `x0_err`, `y0_err`, `pa_err`, and `eps_err`. 
            - The harmonic fitting process is redundant. Need to think more about this, but looks like there is room to improve.
        * `fix_geometry(isophote)`: simply fix the geometry. 
            - Now only use this when the isophote is problematic.
        * `sampled_coordinates()`: calls the `sample.coordinates` function to get the (x,y) coordinates.
        * `to_table()`: convert the output to `astropy.table.QTable`.
            - Also has a `__str__` shortcut.
    - `CentralPixel(Isophote)` Class: for the special central pixel.
        * Most of the values and properties are fixed at 0 or `None`
    - `IsophoteList(Isophote)` Class: the container for a list of `isophote` objects.
        * The list itself is stored in `._list`.
        * Has standard `__len__`, `__delitem__`, `__setitem__`, `__getitem__`, `__iter__`, `__add__` methods.
        * Has `sort()`, `insert(index, value)`, `append(value)`, `etend(value)` methods.
        * `get_closest(sma)`: return the closest isophote to a radii.
        * `_collect_as_array(name)` and `_collect_as_list(name)`: collect certain attribute from the isophote object, and keep them in a list or array. The results are presented as properties of the object.
    - `_isophote_list_to_table(isophote_list)`: converts the isophote object list to a `astropy` table.

- `fitter.py`: provides a class to fit ellipses.
    - A bunch of DEFAULT parameters are defined here.
    - `_CORRECTORS = [_PositionCorrector0(), _PositionCorrector1(), _AngleCorrector(), _EllipticityCorrector()]`: in every iteration, the central position, PA, and ellipticity are changing; and everytime they change, you need a new sample.
        * `_PositionCorrector(_ParameterCorrector)`: change the center to `geometry.x0 + dx` and `geometry.y0 + dy`, redefine a new `EllipseSample`.
        * `_PositionCorrector0(_PositionCorrector)`: calculate the `dx` and `dy` based on `harmonic`, `eps`, and `gradient`.
        * `_PositionCorrector1(_PositionCorrector)`: similar...
        * `_AngleCorrector(_ParameterCorrector)`: provides a new PA to define new sample.
        * `_EllipticityCorrector(_ParameterCorrector)`: provides a new ellipticity value to define a new sample.
    - `CentralEllipseFitter(EllipseFitter)`: special class for central pixel. Just returns the `CentralPixel` value.
    - `EllipseFitter` Class: to fit ellipses.
        * `__init__(sample)`: takes an `isophote.EllipseSample` instance.
        * ` _check_conditions(sample, maxgerr, going_inwards, lexceed)`: 
            * `lexceed`: this flag signals that limiting gradient error (`maxgerr`) wasn't exceeded yet.
            * Check the center, the gradient error, the geometry; returns a flag `proceed` to tell the fitter to proceed or not.
        * `fit(conver=DEFAULT_CONVERGENCE, minit=DEFAULT_MINIT, maxit=DEFAULT_MAXIT, fflag=DEFAULT_FFLAG, maxgerr=DEFAULT_MAXGERR, going_inwards=False)`: core function.
            * Use the current `self._sample`
            * Also tracks the sample that caused the minimu harmonic amplitude(in absolute value). Use it for isophote if `maxit` is exceeded or `fflag` are set.
            * Now start the iteration in range of `maxit`:
            * 1. `sample.update()`: compute its gradient and associated values.
            * 2. `sample.extract()`: get values for angles, radii, and intensity.
            * 3. `fit_first_and_second_harmonics`: not allowed to fail at this step; get `coeffs`.
            * 4. Find the largest harmonic amplitude value and see if it decreases compared to `minimum_amplitude_value`; if so, update the value; `largest_harmonic_index` holds the index of the parameter that drives the large amplitude.
            * 5. `first_and_second_harmonic_function` to compute the `model`; and get the residual; check if the result is converged.; it may not have converged yet, but the sample contains too many invalid data points: return.
            * 6. Choose a `corrector` using the `largest_harmonic_index`; correct and update the sample.
            * 7. Check the condition using `_check_conditions`; if not `proceed`, return the current sample.

- `ellipse.py`: provides a class to fit elliptical isophotes.
    - This is the "user interface", it depends on `fitter`, `geometry`, `integrator`, `isophote`, `sample`.
    - `Ellipse` Class: to fit elliptical isophotes to a galaxy image.
        * `__init__(image, geometry=None, threshold=0.1)`: geometry is stored in `_geometry`.
        * `set_threshold(self, threshold)`: `threshold` is only used for central position.
        * `fit_image()`: key function to fit multiple isophotes to image. It takes almost all the contral parameters as input.

- `model.py`: profiles tools for building a model elliptical galaxy image from a list of isophotes.