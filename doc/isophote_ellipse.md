# Analysis of the `isophote.ellipse` IRAF code

----

## Key Information

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

* `elcopy.x`: routines that handle isophote structure.
    - `el_copy (is, iso)` copies data from one isophote structure to another.
    - `el_alloc (nharm)` allocates structure for holding isophote parameters.
    - `el_free (is)` frees structure that holds isophote parameters.

##### IO and other helpers

* `eldump.x`: writes isophote data to the main output table.
    - `el_dump (is, al, im, outname, tp, column, extracols, rownumber)`
* `eldsample.x`: dumps the azimutal sample of intensity as a binary table.
    - `el_dsample (file, is, im, al, bufx, bufy)`

* `elgeometry.x`: routines to translate coordinates between image (physical), section, and display systems.
    - `el_s2p (im, sec_coord, axis)` converts section coordinates to physical
    - `el_p2s (im, phys_coor, axis)` converts physical coordinates to section
    - Similarly: `el_s2d()`, `el_g2s()`, `el_s2g()`
    - `el_shrink()`, `el_unshrink()` adjust coordinate between physical and section space.

##### Math
* `elbilinear.x`: bi-linear interpolation
    - `el_bilnear (image, pixel_pointer, sma, ip, jp, fx, fy)`
        - `ip, jp`: pixel address; `fx, fy`: fractional pixel position on pixel array.
        - Depends on `el_subpix()`
* `elclip.x`: sigma-clipping function
    - `el_clip (data, usclip, lsclip, nclip, npoint, ndata, mean, sigma)`

##### Display and interactive fitting

* `eldisplay.x`: routines that display images.
    - `el_display (im, dp, is, sec, command, zoom, cx, cy)` is a general-purpose image display function.
    - `el_opend (im, dev)` opens image display structure.
    - `el_dimage (im, dp, is, sec, command, zoom, cx, cy)` displays the image itself.
    - `el_dscale` computes image display scaling parameters.
    - `el_mark (img, dp, sec, gp)` marks flagged pixels.
* `elcursor.x`: `el_cursor(im, dp, x, y, key, command, sz_comm)` handels details of image curor input.
* `elcolon.x`: `el_colon()` processes cursor colon commands.