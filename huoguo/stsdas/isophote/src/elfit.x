# Copyright restrictions apply - see stsdas$copyright.stsdas 

include	<imhdr.h>
include <imio.h>
include	<error.h>

include "ellipse.h"

# EL_FIT -- Fits one elliptical "isophote" over the image. 
# This is where the details of the fitting algorithm are handled.
# The stopping condition code is output in STOP(al). On entry,
# STOP(al) is used to set non-iterative mode, in this way
# allowing information to be passed from one isophote to the
# next. The extracted intensity sample may be written to table.

procedure el_fit (im, sec, is, al, file, sgraphics, splot, list)

pointer	im				# IMIO pointer
pointer	sec				# in-memory pixel structure
pointer	is				# isophote structure
pointer	al				# algorithm control structure
char	file[ARB]			# file name with int. sample
char    sgraphics[ARB]                  # graphics device for sample ploting.
bool	splot				# plot azimuthal samples ?
bool	list				# list at STDOUT ?

pointer bufx,  bufy			# extracted data buffers
pointer	bufx1, bufy1
pointer bufx2, bufy2
pointer bufyd, bufyd2
pointer	is2
bool	fixpar[4]
int	nbuffer				# buffer size
int	nbuffer1, nbuffer2
int	intmode				# mode for gradient computation
int	err				# error code from harmonic fit
int	gndata				# data points in gradient 2nd sample
int	i
int	largestc			# largest harmonic coefficient
real	coef[MAX_HARM], cerror[MAX_HARM]# harmonic coefficients
real	hrms				# rms of harmonic fit
real	mean1, sigma1, biga, aux
real	ea, eb				# errors in major and minor axis
real	gradient			# radial gradient
real	grad_radius			# incremental radius for grad comp.
real	dx, dy				# change in ellipse position.
bool	go				# go to next iteration ?

pointer	el_alloc()
errchk	malloc, realloc, mfree

begin
	# Alloc scratch space.
	call malloc (bufx,   SZ_BUFFER, TY_REAL)
	call malloc (bufy,   SZ_BUFFER, TY_REAL)
	call malloc (bufyd,  SZ_BUFFER, TY_REAL)
	call malloc (bufx1,  SZ_BUFFER, TY_REAL)
	call malloc (bufy1,  SZ_BUFFER, TY_REAL)
	call malloc (bufx2,  SZ_BUFFER, TY_REAL)
	call malloc (bufy2,  SZ_BUFFER, TY_REAL)
	call malloc (bufyd2, SZ_BUFFER, TY_REAL)
	is2 = el_alloc (NHARM(is))
	nbuffer  = SZ_BUFFER
	nbuffer1 = SZ_BUFFER

	# Initialize iteration control variables.
	go = true
	NITER(al) = 1
	biga = 1. / EPSILON

	# Main loop.
	while (go) {

	    # Get intensity sample from image array.
	    call el_get (im, sec, XC(is), YC(is), A(is), EPS(is), TETA(is), 
			bufx, bufy, nbuffer, NPOINT(is), NDATA(is), 
			MEAN(is), SIGMA(is), ASTEP(al), LINEAR(al), 
			INTEGR(al), USCLIP(al), LSCLIP(al), NCLIP(al), 
                        SAREA(is))

	    # Save raw intensity sample, for output at the end.
	    call realloc (bufyd, nbuffer, TY_REAL)
	    call amovr (Memr[bufy], Memr[bufyd], nbuffer)

	    # Fit 1st and 2nd harmonics.
	    if (STOP(al) != ST_NONITERATE) {
	        coef[1] = 1.
	        coef[2] = 1.
	        coef[3] = 2.
	        coef[4] = 2.
	        call el_harmonics (Memr[bufx], Memr[bufy], NPOINT(is), coef, 
	                           cerror, 4, hrms, err)
	        if (err == ERR) {
                    call el_free (is2)
	            call mfree (bufyd2, TY_REAL)
	            call mfree (bufy2,  TY_REAL)
	            call mfree (bufx2,  TY_REAL)
	            call mfree (bufy1,  TY_REAL)
	            call mfree (bufx1,  TY_REAL)
	            call mfree (bufyd,  TY_REAL)
	            call mfree (bufy,   TY_REAL)
	            call mfree (bufx,   TY_REAL)
	            STOP(al)  = ST_SINGULAR
	            VALID(is) = false
	            call eprintf ("Error in 1st and 2nd harmonic fit. ")
	            call eprintf ("Please, verify output table.\n ")
	            return
	        }
	    }

	    # Compute radial gradient:

	    # Enforce nearest-neighbor interpolation, only if sma
	    # is larger than 20 pixels and former gradient error
	    # is reasonable (less than 20%).
	    if ((A(is) <= 20.)           ||
	        (RESLOPE(is) > 0.2)      ||
	        (IS_INDEF(RESLOPE(is))))
	        intmode = INTEGR(al)
	    else
	        intmode = INT_NEIGHBOR

	    grad_radius = 0.1
	    call el_get (im, sec, XC(is), YC(is), (1. + grad_radius)*A(is), 
	                 EPS(is), TETA(is), bufx1, bufy1, nbuffer1, i, 
	                 gndata, mean1, sigma1, ASTEP(al), LINEAR(al), 
			 intmode, USCLIP(al), LSCLIP(al), NCLIP(al), aux)

	    gradient = (mean1 - MEAN(is)) / A(is) / grad_radius

	    # If no meaningful gradient, try another sample, this time
	    # using full sampling and slightly larger radius. Meaningful 
	    # gradient means something shallower, but still close to within 
	    # a factor 3 from previous gradient estimate. If no previous
	    # estimate is available, guess it.
	    if (!VALID(is))
#	        SLOPE(is) = - 0.1 * MEAN(is)
	        SLOPE(is) = - 0.05           # IB 9/22/97
	    if (gradient >= SLOPE(is)/3.) {
	        grad_radius = 0.2
	        call el_get (im, sec, XC(is), YC(is), (1. + grad_radius)*A(is), 
	                     EPS(is), TETA(is), bufx1, bufy1, nbuffer1, i, 
	                     gndata, mean1, sigma1, ASTEP(al), 
	                     LINEAR(al), INTEGR(al), USCLIP(al), LSCLIP(al), 
	                     NCLIP(al), aux)
	        gradient = (mean1 - MEAN(is)) / A(is) / grad_radius
	    }


	    # If still no meaningful gradient can be measured, try with 
	    # previous one, slightly shallower. A factor 0.8 is not too
	    # far from what is expected from geometrical sampling steps
	    # of 10-20% and a deVaucouleurs law or an exponential disk 
	    # (at least at its inner parts, r <~ 5 req). Gradient error 
	    # is meaningless in this case.
	    if (gradient >= SLOPE(is)/3.) {
	        SLOPE(is)   = 0.8 * SLOPE(is)
	        ESLOPE(is)  = INDEFR
	        RESLOPE(is) = INDEFR
	    } else {
	        SLOPE(is)   = gradient
	        if (!IS_INDEF(SIGMA(is))) {
	            ESLOPE(is)  = sqrt (SIGMA(is)**2 / NDATA(is) + 
	                          sigma1**2 / gndata) / A(is) / grad_radius
	            RESLOPE(is) = ESLOPE(is) / abs(gradient)
	        } else {
	            ESLOPE(is)  = INDEFR
	            RESLOPE(is) = INDEFR
	        }
	    }

	    # Compute ellipse parameter errors by direct projection 
	    # of coefficient errors. These showed to be the error
	    # formulae that best convey the true errors measured
	    # by Monte Carlo experiments (see reference in ellipse
	    # help page).
	    ea        = abs (cerror[2] / SLOPE(is))
	    eb        = abs (cerror[1] * (1.-EPS(is)) / SLOPE(is))
	    EX(is)    = sqrt ((ea * cos (TETA(is)))**2 + 
                              (eb * sin (TETA(is)))**2)
	    EY(is)    = sqrt ((ea * sin (TETA(is)))**2 + 
                              (eb * cos (TETA(is)))**2)
	    EEPS(is)  = abs (2. * cerror[4] *
			    (1. - EPS(is)) / A(is) / SLOPE(is))
	    if (abs (EPS(is)) > EPSILON) {
	        ETETA(is) = abs (2. * cerror[3] *
			    (1. - EPS(is)) / A(is) / SLOPE(is) / 
			    (1. - (1. - EPS(is))**2))
	    } else {
	        ETETA(is) = 0.
	    }

#  THESE ARE ABORTED ATTEMPTS TO FIND AN EMPIRICAL ERROR FORMULATION THAT
#  TRULY DEPICTS THE ACTUAL ERROR FOUND IN MONTE CARLO EXPERIMENTS.

#	    # Compute ellipse parameter errors by direct projection 
#	    # of coefficient errors combined with gradient error:
#	    if (!IS_INDEF(ESLOPE(is))) {
#	        aux = ESLOPE(is) * A(is) * grad_radius
#	        ea  = abs (sqrt(cerror[2]**2 + aux**2) / SLOPE(is))
#	        eb  = abs (sqrt(cerror[1]**2 + aux**2) * 
#	                  (1.-EPS(is)) / SLOPE(is))
#	        EX(is) = sqrt ((ea * cos (TETA(is)))**2 + 
#                              (eb * sin (TETA(is)))**2)
#	        EY(is) = sqrt ((ea * sin (TETA(is)))**2 + 
#                              (eb * cos (TETA(is)))**2)
#
#	        EEPS(is) = abs (2. * sqrt(cerror[4]**2 + aux**2) *
#		        	(1. - EPS(is)) / A(is) / SLOPE(is))
#
#	        if (abs (EPS(is)) > EPSILON) {
#	            ETETA(is) = abs (2. * sqrt(cerror[3]**2 + aux**2) *
#			        (1. - EPS(is)) / A(is) / SLOPE(is) / 
#			        (1. - (1. - EPS(is))**2))
#	        } else
#	            ETETA(is) = 0.
#	    } else {
#	        EEPS(is)  = INDEFR
#	        ETETA(is) = INDEFR
#	        EX(is)    = INDEFR
#	        EY(is)    = INDEFR
#	    }

	    # Compute ellipse parameter errors by direct projection 
	    # of coefficient errors with empirical correction:
#	    if (!IS_INDEF(ESLOPE(is))) {
#	        aux = abs (ESLOPE(is) / SLOPE(is))
#	        if (aux > 0.9)
#	            EEPS(is) = EEPS(is) + (aux - 0.9)**2 * 0.4
#	        if ((aux > 0.9) && (EPS(is) > 0.1))
#	                ETETA(is) = ETETA(is) + (aux - 0.9)**2 * 70./180.* PI
#	        if ((aux > 0.9) && (EPS(is) <= 0.1))
#	                ETETA(is) = ETETA(is) + (aux - 0.9)**2 * 70./180.* PI
#	        if (aux > 0.9) {
#	            EX(is) = EX(is) + (aux - 0.9)**2 * 40.
#	            EY(is) = EY(is) + (aux - 0.9)**2 * 40.
#	        }
#	    } else {
#	        EEPS(is)  = INDEFR
#	        ETETA(is) = INDEFR
#	        EX(is)    = INDEFR
#	        EY(is)    = INDEFR
#	    }


	    # Error propagation (coefficient and slope only):
#	    if (IS_INDEF(ESLOPE(is))) {
#	        ea  = abs (cerror[2] / coef[2])
#	        eb  = abs (cerror[1] / coef[1] * (1.-EPS(is)))
#	    } else {
#	        ea  = abs (sqrt((cerror[2] / coef[2]) ** 2 +
#	                         (ESLOPE(is) / SLOPE(is)) ** 2))
#	        eb  = abs (sqrt((cerror[1] / coef[1]) ** 2 +
#	                         (ESLOPE(is) / SLOPE(is)) ** 2) * (1.-EPS(is))) 
#	    }
#	    EX(is)    = sqrt ((ea * cos (TETA(is)))**2 + 
#                             (eb * sin (TETA(is)))**2)
#	    EY(is)    = sqrt ((ea * sin (TETA(is)))**2 + 
#                             (eb * cos (TETA(is)))**2)
#
#	    EEPS(is) = 2. * (1. - EPS(is)) * EPS(is) / A(is)
#	    if (IS_INDEF(ESLOPE(is)))
#	        EEPS(is) = EEPS(is) * cerror[4] / coef[4]
#	    else 
#	        EEPS(is) = EEPS(is) * 
#	                   sqrt ((cerror[4] / coef[4]) ** 2 +
#	                         (ESLOPE(is) / SLOPE(is)) ** 2)
#	    EEPS(is) = abs(EEPS(is))
#
#	    if (abs (EPS(is)) > EPSILON) {
#	        ETETA(is) = 2. * (1. - EPS(is)) * EPS(is) / A(is) /
#			    ((1. - EPS(is))**2 - 1.)
#	        if (IS_INDEF(ESLOPE(is)))
#	            ETETA(is) = ETETA(is) * cerror[3] / coef[3]
#	        else 
#	            ETETA(is) = ETETA(is) * 
#	                        sqrt ((cerror[3] / coef[3]) ** 2 +
#	                              (ESLOPE(is) / SLOPE(is)) ** 2)
#	        ETETA(is) = abs(ETETA(is))
#	    } else {
#	        ETETA(is) = 0.
#	    }

	    # Examine stopping criteria.
	    # First, get largest harmonic amplitude. Ignore if fixed.
	    ABIG(is) = -1./EPSILON
	    fixpar[1] = FIXX(al)	# these map the geometric
	    fixpar[2] = FIXY(al)	# parameters on the harmonic
	    fixpar[3] = FIXT(al)	# coefficient vector.
	    fixpar[4] = FIXE(al)
	    do i = 1, 4 {
	        if ((abs (coef[i]) > ABIG(is)) && (!fixpar[i])) {
	            ABIG(is) = abs (coef[i])
	            largestc = i
	        }
	    }

	    # This is the main convergency criterion. The largest
	    # harmonic amplitude is compared with the residual rms
	    # around the harmonic fit. Conparison takes into account
	    # user-supplied factor (usually a few percent) and also 
	    # the smoothing effect due to area averaging.
	    if (STOP(al) != ST_NONITERATE) {
	        if (abs (coef[largestc]) <= 
	            (CPAR(al) * sqrt(SAREA(is)) * hrms)) {
	            STOP(al) = ST_OK
	            go = false
	        }
	    } else 
	        go = false

	    # Examine other stopping conditions.
	    if (real(NDATA(is)) < 
	       (real(NPOINT(is)) * (1 - FBAD(al)))) { # less than FBAD
	        STOP(al) = ST_INDEF                   # sampled points 
	        go = false                            # with useful data.
	    } else if (NITER(al) >= NMAX(al)) {       # exceeded max. allowed
	        STOP(al) = ST_MAXIT		      # number of iterations.
	        go = false
	    }
	    if (STOP(al) == ST_OK) {		   # ensures minimum number
	        if (NITER(al) < NMIN(al))	   # of iterations.
	            go = true 
	    }

	    # If this is the isophote whith smallest amplitudes so far,
	    # store its parameters and associated data.
	    if (ABIG(is) < biga) {
	        call el_copy (is, is2)
	        biga     = ABIG(is)
	        nbuffer2 = nbuffer
	        call realloc (bufx2,  nbuffer2, TY_REAL)
	        call realloc (bufy2,  nbuffer2, TY_REAL)
	        call realloc (bufyd2, nbuffer2, TY_REAL)
	        call amovr (Memr[bufx],  Memr[bufx2],  nbuffer2)
	        call amovr (Memr[bufy],  Memr[bufy2],  nbuffer2)
	        call amovr (Memr[bufyd], Memr[bufyd2], nbuffer2)
	    }

	    # If everything fixed, no need to iterate.
	    if (FIXX(al) && FIXY(al) && FIXT(al) && FIXE(al))
	        go = false

	    # Iterate once more, if iterations still enabled.
	    if (go && (STOP(al) != ST_NONITERATE)) {
	        NITER(al) = NITER(al) + 1
	        # Apply correction to coeff. with largest amplitude.
	        switch (largestc) {
	        case 1:
	            aux = -coef[1] * (1. - EPS(is)) / SLOPE(is)
	            dx = -aux * sin (TETA(is))
	            dy =  aux * cos (TETA(is))
	            XC(is) = XC(is) + dx
	            YC(is) = YC(is) + dy
	        case 2:
	            aux = -coef[2] / SLOPE(is) 
	            dx = aux * cos (TETA(is))
	            dy = aux * sin (TETA(is))
	            XC(is) = XC(is) + dx
	            YC(is) = YC(is) + dy
	        case 3:
	            dx = 0.0
	            dy = 0.0
	            aux = coef[3] * 2. * (1. - EPS(is)) / 
			  A(is) / SLOPE(is) / ((1. - EPS(is))**2 - 1.) 
	            TETA(is) = TETA(is) + aux
	            # Restore angle to valid range (0 - PI).
	            while (TETA(is) > DPI)
	                TETA(is) = TETA(is) - DPI
	            if (TETA(is) > PI)
	                TETA(is) = TETA(is) - PI
	            while (TETA(is) < -DPI)
	                TETA(is) = TETA(is) + DPI
	            if (TETA(is) < -PI)
	                TETA(is) = TETA(is) + PI
	            if (TETA(is) < 0.)
	                TETA(is) = TETA(is) + PI
	        case 4:
	            dx = 0.0
	            dy = 0.0
	            aux = coef[4] * 2. * (1. - EPS(is)) / A(is) / SLOPE(is)
	            EPS(is) = EPS(is) - aux
	        }
	    }

	    # If center wandered more than allowed, put it back
	    # in place and signal the end of iterative mode.
	    if (!IS_INDEFR(WANDER(al))) {
	        if ((abs(dx) > WANDER(al)) || (abs(dy) > WANDER(al))) {
	            XC(is) = XC(is) - dx
	            YC(is) = YC(is) - dy
	            STOP(al) = ST_NONITERATE
	            go = false
	       }
	    }

	    # If ellipse diverged, or if maxrit was reached,
	    # signal the end of iterative mode for now on.
	    if (STOP(al) == ST_OK) {
	        if ((abs (EPS(is)) > (MAX_EPS)) ||
	            (XC(is) < 1.) || (XC(is) > IM_LEN (im, 1)) || 
	            (YC(is) < 1.) || (YC(is) > IM_LEN (im, 2)) ||
		    (A(is) >= MAXRIT(al))) {
	            STOP(al) = ST_NONITERATE
	            go = false
	        }
	    }

	    # See if eps == 0 (round isophote) was crossed.
	    if (EPS(is) < 0.) {
	        EPS(is) = -EPS(is)
	        if (TETA(is) < 0.) 
	            TETA(is) = TETA(is) + PI2
	        else
	            TETA(is) = TETA(is) - PI2
	    }

	    # If ellipse is an exact circle, computations diverge.
	    # Make it slightly flat.
	    if (EPS(is) == 0.0)
	        EPS(is) = 0.05
	}

	# If stopped abnormally, get data from ellipse with 
	# smallest harmonic amplitudes.
	if (STOP(al) != ST_OK) {
	    call el_copy (is2, is)
	    call amovr (Memr[bufx2],  Memr[bufx],  nbuffer2)
	    call amovr (Memr[bufy2],  Memr[bufy],  nbuffer2)
	    call amovr (Memr[bufyd2], Memr[bufyd], nbuffer2)
	}

	# Fit higher order harmonics. Third and fourth harmonics
	# are actual amplitudes in intensity but projected by gradient 
	# and normalized by semi-major axis.
	coef[1] = 3.				# 3rd harmonic
	coef[2] = 3.
	call el_harmonics (Memr[bufx], Memr[bufy], NPOINT(is), coef, cerror, 
			2, hrms, err)
	if (err == ERR) {
	    STOP(al) = ST_SINGULAR
	    call eprintf ("Error in 3rd harmonic fit. ")
	    call eprintf ("Please, verify output table.\n ")

	}
	A3(HH(is))  = coef[1]   / A(is) / abs(SLOPE(is)) 
	B3(HH(is))  = coef[2]   / A(is) / abs(SLOPE(is))
	if ((cerror[1] > 0.) && (coef[1] != 0.) &&
	    (cerror[2] > 0.) && (coef[2] != 0.)) {
	    if (!IS_INDEF(RESLOPE(is))) {
	        A3(EHH(is)) = abs(A3(HH(is))) * sqrt((cerror[1] / coef[1])**2 + 
	                      RESLOPE(is)**2)
	        B3(EHH(is)) = abs(B3(HH(is))) * sqrt((cerror[2] / coef[2])**2 + 
	                      RESLOPE(is)**2)
	    } else {
	        # Adopt reasonably large gradient error (~80%).
	        A3(EHH(is)) = abs(A3(HH(is))) * sqrt((cerror[1]/ coef[1])**2 + 
	                      0.64)
	        B3(EHH(is)) = abs(B3(HH(is))) * sqrt((cerror[2]/ coef[2])**2 + 
	                      0.64)
	    }
	} else {
	    A3(EHH(is)) = INDEFR
	    B3(EHH(is)) = INDEFR
	}

	coef[1] = 4.					# 4th harmonic
	coef[2] = 4.
	call el_harmonics (Memr[bufx], Memr[bufy], NPOINT(is), coef, cerror, 
			2, hrms, err)
	if (err == ERR) {
	    STOP(al) = ST_SINGULAR
	    call eprintf ("Error in 4th harmonic fit. ")
	    call eprintf ("Please, verify output table.\n ")
	}
	A4(HH(is))  = coef[1]   / A(is) / abs(SLOPE(is))
	B4(HH(is))  = coef[2]   / A(is) / abs(SLOPE(is))
	if ((cerror[1] > 0.) && (coef[1] != 0.) &&
	    (cerror[2] > 0.) && (coef[2] != 0.)) {
	    if (!IS_INDEF(RESLOPE(is))) {
	        A4(EHH(is)) = abs(A4(HH(is))) * sqrt((cerror[1]/ coef[1])**2 + 
	                      RESLOPE(is)**2)
	        B4(EHH(is)) = abs(B4(HH(is))) * sqrt((cerror[2]/ coef[2])**2 + 
	                      RESLOPE(is)**2)
	    } else {
	        # Adopt reasonably large gradient error (~80%).
	        A4(EHH(is)) = abs(A4(HH(is))) * sqrt((cerror[1]/ coef[1])**2 + 
	                      0.64)
	        B4(EHH(is)) = abs(B4(HH(is))) * sqrt((cerror[2]/ coef[2])**2 + 
	                      0.64)
	    }
	} else {
	    A4(EHH(is)) = INDEFR
	    B4(EHH(is)) = INDEFR
	}

	# Fit optional harmonics. They are NOT projected !
	if (NHARM(is) > 0) {
	    do i = 1, NHARM(is), 2 {
	        coef[1] = AI(HARM_NUMBERS(is),i)
	        coef[2] = BI(HARM_NUMBERS(is),i)
	        call el_harmonics (Memr[bufx], Memr[bufy], NPOINT(is), 
				   coef, cerror, 2, hrms, err)
	        if (err == ERR) {
	            STOP(al) = ST_SINGULAR
	            call eprintf ("Error in upper harmonic fit. ")
	            call eprintf ("Please, verify output table.\n ")
	        }
	        AI(HARM(is),i)  = coef[1]
	        BI(HARM(is),i)  = coef[2]
	        AI(EHARM(is),i) = cerror[1]
	        BI(EHARM(is),i) = cerror[2]
	    }
	}

	# If ellipse was held constant, discard parameter errors.
	if (STOP(al) == ST_NONITERATE) {
	    EX(is)     = INDEFR
	    EY(is)     = INDEFR
	    EEPS(is)   = INDEFR
	    ETETA(is)  = INDEFR
	}

	# If everything fixed, discard largest amplitude data.
	if (FIXX(al) && FIXY(al) && FIXT(al) && FIXE(al))
	    ABIG(is) = INDEFR

	# Dump intensity sample.
	call el_dsample (file, is, im, al, bufx, bufyd) 

	# Plot intensity sample. 
	if (splot)
	    call el_psample (im, is, al, sgraphics, Memr[bufx], Memr[bufyd])

	# Output some info to STDOUT.
	if (list) {
	    call el_unshrink (im, is, al)

	    call printf ("%7.2f %8.2f(%6.2f) %5.3f(%5.3f) %6.2f(%5.2f)")
	        call pargr (A(is))
	        call pargr (MEAN(is))
	        call pargr (SIGMA(is))
	        call pargr (EPS(is))
	        call pargr (EEPS(is))
	        call pargr (TETA(is))
	        call pargr (ETETA(is))


	    call printf (" %5.3f")
	        call pargr (RESLOPE(is))
	    call printf ("%5d %3d  %4d   %2d\n")
	        call pargi (NDATA(is))
	        call pargi (NPOINT(is) - NDATA(is))
	        call pargi (NITER(al))
	        call pargi (STOP(al))
	    call flush (STDOUT)

	    call el_shrink (im, is, al)
	}

	# Signal valid data now in isophote structure.
	VALID(is) = true

	# Reclaim scratch space
        call el_free (is2)
	call mfree (bufyd2, TY_REAL)
	call mfree (bufy2,  TY_REAL)
	call mfree (bufx2,  TY_REAL)
	call mfree (bufy1,  TY_REAL)
	call mfree (bufx1,  TY_REAL)
	call mfree (bufyd,  TY_REAL)
	call mfree (bufy,   TY_REAL)
	call mfree (bufx,   TY_REAL)
end






