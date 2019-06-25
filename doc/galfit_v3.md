# Analysis of the Galfit v3.0.7 Source Code

----

## File Contents

* Only considers the C source code and the header file.

### Top-Level

* `Makefile`: how to compile `galfit`. Looks pretty straightforward.
* `galfit.c`: top-level code for making the executable file `galfit`
    - Dependes on `stdio.h`, `stdlib.h`, and `curses.h` (definitions for screen handling and optimisation functions)
    - Uses `nrutil.h`, `structs.h`, and `debug.h`
    - Besides the `main()` function, there is just one `blank_img(image, inpars)` function.

### Headers

* `const.h`: just `PI`
* `debug.h`: just `DEBUG` and `CKIMG`
* `mymath.h`: just one `NINT(f)` function that add/minus 0.5 to `f` based on whether `f>0.`.
* `modelfunc.h`: definest the functions for each model component: `sersic`, `exponential`, `gaussian`, `moffat`, `ferrer`, `king`, `edgedisk`, `nuker`.
    - **Note**: some of these, like the `ferrer` and `nuker` are not really necessary. But should consider including something new from `imfit` or `libprofit`, or `galsim`.
* `nrutil.h`: the C version of the Numerical Recipes utility. Just some basic numerical functions.
    - **Note**: `galfit` makes use some of the algorithms from the book "Numerical Recipes in C". This header is for using these algorithms. Some of these can be improved or updated.
* `optimization.h`: defines functions used for optimization.
* `sersic.h`: defines the look-up table for Sersic index and the `kappa` parameter.
* `struct.h`: defines all the structures used in `galfit`.
    - Depends on `fitsio.h` from `cfitsio`.
    - Including `inpars{}`, `image{}`, `sampars{}`, `derivs{}`, `fitpars{}`, `cons{}`, `convpars{}`, `bicubic_grid{}`, `trunc_val{}`, `trunc_links{}`, `runflags{}`, `rotpars{}`. Comments are very helpful.

### Algorithms

* `assign_err.c`: scale the errorbars by `sqrt(chi2nu)` when `chi2nu != 1`.
    - `void assign_err (struct fitpars *fpar, double *sig, float chi2nu)`

* `calc_angles.c`: calculates the parameters for coordinate rotation.
    - `void calc_rotpars (float a[], struct fitpars *r, struct rotpars *rp)`

* `chi2calc.c`: calculates the chi2 of the model.
    - `double chi2calc (int nfree, float **data, float **sig, float **model, long naxes[], int *ndof)`

* `convolve.c`: convolves a model with the PSF. The input PSF must already be Fourier transformed using the `fftpsf.c`.
    - Complex_naxes[] stores the size of the FFT image, and the size is the number of real+complex pairs. So the physical storage in the x-direction is twice as large as the number in complex_naxes[2].
    - `void convolve (struct image *model, double **psfft, long complex_naxes[])`
    - `void complmultipl(double ar1[], double ar2[], unsigned long n)`
    - `void mkimg (float **array, long naxes[], char *outname, char *tname)`
    - `void copy_convregion (struct convpars *cpar, float **model, struct image convolved_img, int os)`
        - copy the convolved region back to the mini-size image

#### Components

#### High-level features

* `constraints.c`: deals with parameter constraints.
    - `void constraints (int mfit, float atry[], float a[], int ia[], struct cons *constr, float *aorig)`
    - **Note**: some of these are not very often used and can be ignored first.

* `boxiness.c`: extracts the diskyness/boxyness that the user wants to fit.
    - `void boxy_par (char *string, struct fitpars *fpar)`

* `bending.c`: extracts the Banana-fana parameter values that the user wants to fit.
    - `void bending_pars (char *string, struct fitpars *fpar)`

#### From the Numerical Recipe in C

* `bcucof.c`: constructs two-dimensional bicubic
* `bcuint.c`: two-dimensional bicubic interpolation
* `bessi1.c`: modified Bessel function I1
* `bessk0.c`: modified Bessel function K0
* `bessk1.c`: modified Bessel function K1
* `beta.c`: based on the Numerical Recipes.
    - Calculates the mathematical Beta function, which is ultimately used to compute the integrated magnitude of the generalized elliptical Sersic, Gaussian, Expdisk, etc. functions.
    - `float beta(float z, float w)`
    - `float ratio (float c)`: returns `PI * c / (4. * beta(1./c, 1 + 1./c))`

### Data IO and Other Helpers

* `assign_pars.c`: assigns parameters to structures.
    - -1 means this parameter is not being used
    - `int assign_pars (struct inpars *input, struct image *data, struct image *sigma, struct image *psf, struct image *badpix, struct fitpars *fpar, struct convpars *cpar, struct cons *constr)`

* `badpixlist.c`: reads in the mask file and convert it into a bad pixel list.
    - `badpixlist (struct image *badpix)`
    - **Note**: needs to investigate more.

* `check_constraints.c`: looks for strong parameter couplings and sets the fitting flag to 2 for all except the first component.
    - `void check_constraints (struct cons *constr, float a[], int ia[])`

* `checkprob.c`: checks to see if there may be problems with the fitting parameters that may cause convergence issues.
    - `int checkprob (char *objtype, int parnum, float val, float limit)`
    - `int check_ioffset (int ia, int prob)`

* `clearbuff.c`: clears the keyboard buffer. Used after a `scanf` statement to get rid of extra characters the user accidentally entered.

* `clearmodel.c`: clears the model image and fitting parameters.
    - `void clear_fpar_images (struct derivs *df, struct fitpars *fpar)`
    - `void clearmodel (struct image *model, struct derivs *df, struct fitpars *fpar)`

* `convregion.c`: figures out how big of a convolution box we have for each object.
    - `void convregion (struct derivs *df, float a[], struct image *d, struct convpars *cpar)`

* `copymat.c`: `void copymat (float **fmat, double **dmat, long naxes[], int op)`
* `copy_pars.c`: `void copy_pars (float a[], int ia[], struct fitpars *fpar, int dir)`
* `copy_struct.c`: copies the parameter structs.
    - `void copy_struct (struct fitpars *to, struct fitpars *from)`
* `count_pars.c`: counts the number of parameters.
    - `void count_pars (struct fitpars *fpar, int *nfree, int *nfixed, int *unused)`