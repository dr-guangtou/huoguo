# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<imhdr.h>
include	"ellipse.h"

define	PHI_MAX		0.2		# limits for sector angular whidth
define	PHI_MIN		0.05

# EL_GET -- Samples the image array over an elliptical path and returns two 
# buffers: one with the sampled azimutal angles, and the other with the sampled 
# intensities at the corresponding angles, with the mean subtracted. The mean, 
# standard deviation and number of sampled points are also returned. Samples 
# taken outside the image boundaries are given the value INDEF. Four sampling 
# modes are supported: bi-linear interpolation, mean or median over elliptical 
# annulus sectors, and nearest-neighbor. This last one is intended to be used 
# only when computing image gradient. Bi-linear mode integrates at the 
# sub-pixel level when the ellipse semi-major axis is smaller than 20 pixels.
# Sigma-clipping on extracted sample buffer can be optionally used. Output 
# buffers must be allocated at the calling program, with an initial size given 
# by nbuffer. nbuffer is also used as an incremental value for buffer size. 
# Isophote and algorithm parameters are passed explicitly in the parameter
# list instead of using structures, to allow sampling with any mode and at 
# any position in image, independent from current isophote and control 
# parameters.

procedure el_get (im, sec, x0, y0, a, eps, teta, bufx, bufy, nbuffer,
		  npoint, ndata, mean, sigma, astep, linear, integrmode, 
	          usclip, lsclip, nclip, aarea)
pointer	im				# i: IMIO pointer
pointer sec				# i: pointer to in-memory pixel struct
real	x0, y0				# i: center of ellipse on frame
real	a				# i: semi-major axis (in pixels)
real	eps				# i: ellipticity
real	teta				# i: position angle (in radians)
pointer bufx, bufy			# io: extracted data buffer pointers
int	nbuffer				# io: buffer size
int	npoint				# o: no. of acquired samples
int	ndata				# o: no. of valid samples (!= INDEFR)
real	mean, sigma			# o: sample statistics
real	astep				# i: distance between ellipses
int	integrmode			# i: sampling mode
real	usclip				# i: upper sigma-clip criterion
real	lsclip				# i: lower sigma-clip criterion
int	nclip				# i: number of s-clip iterations
real	aarea				# o: average sector area
bool	linear				# i: linear grow of a ?

pointer	medbuf				# buffer for sector median
int	i, j, i1, j1, i2, j2
int	bufsiz
int	npix
real	r, phi				# polar coord. of sector center
real	dphi				# angular whidth of sector
real	phistep				# step in polar angle
real	sarea				# actual sector area
real	rp, phip			# pixel polar coordinates
real	r1, r2, r3, r4, phi1, phi2	# polar coord. of sector vertices
real	x[4], y[4]			# image coordinates of sector vertices
real	a1, a2, rp1, rp2
real	fx, fy
real	area				# rigorous sector area
real	sa1, sa2, sa3, sa4
real	pixel, sample
real	aux
double	s, s2

real	el_sarea(), el_bilinear()
int     el_comparer()		# function to be called to compare elements
extern	el_comparer

errchk	el_getsec, realloc

begin
	# Initialize ellipse scanning.
	r = a
	if (linear) {				# limiting annulus ellipses.
	    a1 = a - astep / 2. 
	    a2 = a + astep / 2.
	} else {
	    a1 = a * (1. - ((1. - 1./astep) / 2.)) 
	    a2 = a * (1. + (astep - 1.) / 2.)
	}
	r3     = a2				# for building first sector.
	r4     = a1
	aux    = min ((a2 - a1), 3.)
	sarea  = (a2 - a1) * aux
	dphi   = max (min ((aux / a), PHI_MAX), PHI_MIN)
	phi    = dphi / 2.
	phi2   = phi - dphi / 2.
	aux    = 1. - eps
	r3     = a2 * aux / sqrt ((aux * cos (phi2))**2 + (sin (phi2))**2)
	r4     = a1 * aux / sqrt ((aux * cos (phi2))**2 + (sin (phi2))**2)
	npoint = 0
	s      = 0.d0
	s2     = 0.d0
	aarea  = 0.0
	bufsiz = nbuffer

	# Sample ellipse over image array.
	while (phi < DPI) {

	    switch (integrmode) {

	    case INT_NEIGHBOR:
	        area = 1.
	        # Step in angle is coarser in nearest-neighbor mode.
	        phistep = 2. / r
	        # Get image coordinates of (r, phi) pixel.
	        i = int (r * cos (phi + teta) + x0)
	        j = int (r * sin (phi + teta) + y0)
	        # If outside image boundaries, invalidate data point.
	        if ((i > 0) && (i < IM_LEN (im, 1)) &&
	            (j > 0) && (j < IM_LEN (im, 2))) {
	            call el_getpix (im, sec, i, j)
	            if (!IS_INDEFR (Memr[SUBRASTER(sec)]))
	                sample = Memr[SUBRASTER(sec)] 
	            else 
	                sample = INDEFR
	        } else
	            sample = INDEFR

	    case INT_LINEAR:
	        area    = 2.
	        phistep = 1. / r
	        # Get image coordinates of (r, phi) pixel.
	        x[1] = r * cos (phi + teta) + x0
	        y[1] = r * sin (phi + teta) + y0
	        i    = x[1]
	        j    = y[1]
	        # If outside image boundaries, invalidate data point.
	        if ((i > 0) && (i < IM_LEN (im, 1)-1) &&
	            (j > 0) && (j < IM_LEN (im, 2)-1)) {
	            fx = x[1] - real(i)
	            fy = y[1] - real(j)
	            sample = el_bilinear (im, sec, a, i, j, fx, fy)
	        } else
	            sample = INDEFR

	    case INT_MEAN, INT_MED:
	        # Get image coordinates of the four corners of the
	        # subraster which contains the elliptical sector. The
	        # sector has a whidth in phi such that it comes out
	        # whith roughly constant area along the ellipse.
	        phi1 = phi2
	        r1   = r4
	        r2   = r3
	        phi2 = phi + dphi / 2.
	        aux  = 1. - eps
	        r3   = a2 * aux / sqrt ((aux * cos (phi2))**2 + (sin (phi2))**2)
	        r4   = a1 * aux / sqrt ((aux * cos (phi2))**2 + (sin (phi2))**2)
	        x[1] = r1 * cos (phi1 + teta) + x0
	        y[1] = r1 * sin (phi1 + teta) + y0
	        x[2] = r2 * cos (phi1 + teta) + x0
	        y[2] = r2 * sin (phi1 + teta) + y0
	        x[3] = r3 * cos (phi2 + teta) + x0
	        y[3] = r3 * sin (phi2 + teta) + y0
	        x[4] = r4 * cos (phi2 + teta) + x0
	        y[4] = r4 * sin (phi2 + teta) + y0
		call el_qsortr (x, 4, el_comparer)
		call el_qsortr (y, 4, el_comparer)
	        i1 = x[1] - 1.
	        j1 = y[1] - 1.
	        i2 = x[4] + 1.
	        j2 = y[4] + 1.

	        # Compute sector area.
	        sa1  = el_sarea (a1, eps, phi1, r1)
	        sa2  = el_sarea (a2, eps, phi1, r2)
	        sa3  = el_sarea (a2, eps, phi2, r3)
	        sa4  = el_sarea (a1, eps, phi2, r4)
	        area = abs ((sa3 - sa2) - (sa4 - sa1))

	        # Compute step to next sector and its angular span.
	        dphi = max (min ((sarea / (r3 - r4) / r4), PHI_MAX), PHI_MIN)
	        phistep = dphi / 2. + phi2 - phi

	        # Read subraster. If outside image boundaries, 
	        # invalidate data point.
	        if ((i1 > 0) && (i1 < IM_LEN (im, 1)) &&
	            (j1 > 0) && (j1 < IM_LEN (im, 2)) &&
	            (i2 > 0) && (i2 < IM_LEN (im, 1)) &&
	            (j2 > 0) && (j2 < IM_LEN (im, 2))) {
	            call el_getsec (im, sec, i1, i2, j1, j2)
	            # Create buffer for median computation.
	            if (integrmode == INT_MED)
	                call malloc (medbuf, (i2-i1+1)*(j2-j1+1) , TY_REAL)
	            # Scan subraster, compute mean or median intensity.
	            sample = 0.
	            npix   = 0
	            do j = j1, j2 {
	                do i = i1, i2 {
	                    # Check if polar coordinates of each subraster
	                    # pixel put it inside elliptical sector.
	                    call el_polar (float(i), float(j), x0, y0, teta,
					   rp, phip)
	                    if ((phip < phi2) && (phip >= phi1)) {
	                        aux = (1. - eps) / sqrt (((1. - eps) * 
				      cos (phip))**2 + (sin (phip))**2)
	                        rp1 = a1 * aux 
	                        rp2 = a2 * aux 
	                        if ((rp < rp2) && (rp >= rp1)) {
                                    pixel = Memr[SUBRASTER(sec) + 
						 (j-j1) * (i2 - i1 + 1) +
						 i - i1]
	                            # Add valid pixel to sample.
	                            if (!IS_INDEFR (pixel)) {
	                                switch (integrmode) {
	                                case INT_MED:
	                                    Memr[medbuf+npix] = pixel
	                                    npix = npix + 1
	                                case INT_MEAN:
	                                    sample = sample + pixel 
	                                    npix = npix + 1
	                                }
	                            }
	                        }
	                    }
	                }
	            }
	            # If 6 or less pixels were sampled, get the
	            # bi-linear interpolated value instead.
	            if (npix <= 6) {
	                r       = (r1 + r2 + r3 + r4) / 4.
	                area    = 2.
	                x[1]    = r * cos (phi + teta) + x0
	                y[1]    = r * sin (phi + teta) + y0
	                i       = x[1]
	                j       = y[1]
	                fx      = x[1] - real(i)
	                fy      = y[1] - real(j)
	                sample  = el_bilinear (im, sec, a, i, j, fx, fy)
	                if (integrmode == INT_MED)
	                    call mfree (medbuf, TY_REAL)
	            } else {
	                # Compute mean or median.
	                switch (integrmode) {
	                case INT_MED:
	                    call el_qsortr (Memr[medbuf], npix, el_comparer)
	                    switch (mod (npix,2)) {
	                    case 0:
	                        sample = (Memr[medbuf + npix/2 - 1] +
	                                  Memr[medbuf + npix/2]) / 2.
	                    case 1:
	                        sample = Memr[medbuf + npix/2]
	                    }
	                    call mfree (medbuf, TY_REAL)
	                case INT_MEAN:
	                    sample = sample / float (npix)
	                }
	            }
	        } else
	            sample = INDEFR
	    }

	    # Add to output sample vectors.
	    Memr[bufx + npoint] = phi
	    Memr[bufy + npoint] = sample
	    npoint = npoint + 1
	    if (npoint > nbuffer) {
		nbuffer = nbuffer + bufsiz
	        call realloc (bufx, nbuffer, TY_REAL)
	        call realloc (bufy, nbuffer, TY_REAL)
	    }

	    # Update total area.
	    aarea = aarea + area

	    # Step over ellipse and calculate new radius vector of
	    # sector center. Step in phi must have a maximum value,
	    # otherwise subsampling may occur.
	    phi = phi + min (phistep, 0.5)
	    r   = a * (1. - eps) / sqrt (((1. - eps) * cos (phi))**2 + 
	          (sin (phi))**2)
	}

	# Compute sample statistics. VOPS operators can't be
	# used because of possible INDEF-valued points.
	s     = 0.D0
	s2    = 0.D0
	ndata = 0
	do i = 0, npoint-1 {
	    if (!IS_INDEFR (Memr[bufy+i])) {
	        s  = s  + Memr[bufy+i]
	        s2 = s2 + Memr[bufy+i] ** 2
	        ndata = ndata + 1
	    }
	}
	switch (ndata) {
	case 0:
	    mean  = 0.
	    sigma = INDEFR
	case 1:
	    mean  = s
	    sigma = INDEFR
	default:
	    mean  = real (s / double (ndata))
	    sigma = real((s2 - s * mean) / double (ndata - 1))
	    if (sigma > 0.)
	        sigma = sqrt(sigma)
	    else
	        sigma = 0.
	}

	# Clip off deviant points, only if it makes sense.
	# Statistics are updated accordingly.
	if ((nclip > 0) && (npoint > 15))
	    call el_clip (Memr[bufy], usclip, lsclip, nclip, npoint, ndata, 
                          mean, sigma)

	# Subtract mean from y data.
	do i = 0, npoint-1 {
	    if (!IS_INDEFR (Memr[bufy+i]))
	        Memr[bufy+i] = Memr[bufy+i] - mean
	}

	# Compute average sector area. This disregards the fact
	# that there might be clipped-off data points, since the 
	# sum was computed before clipping. For the average
	# computation, however, this is not a drawback.
	switch (integrmode) {
	case INT_MEAN,INT_MED:
	    aarea = aarea / float (npoint)
	case INT_LINEAR:
	    aarea = 2.
	case INT_NEIGHBOR:
	    aarea = 1.
	}
end


# EL_SAREA -- Compute area of elliptical sector.

real procedure el_sarea (a, eps, phi, r)

real	a, eps, phi, r
real	aux, saux

begin
	aux  = r * cos (phi) / a
	saux = aux / abs(aux)
	if (abs (aux) >= 1.)
	    aux = saux
	return (abs (a ** 2 * (1.-eps)/2. * acos (aux)))
end

