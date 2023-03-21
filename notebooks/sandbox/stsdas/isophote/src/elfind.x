# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imhdr.h>
include "ellipse.h"

define	RWINDOW		10	# re-centering window size.
define	INNER_RADIUS	4.	# concentric sample radii.
define	OUTER_RADIUS	8.


# EL_FIND  --  Atempts to find object center in frame. 

procedure el_find (im, is, al, sec, recenter, list, thresh)

pointer	im			# IMIO
pointer	is			# isophote data
pointer	al			# algorithm control
pointer	sec			# pixel/mask data
bool	recenter		# recenter x0,y0 ?
bool	list			# write msg at STDOUT ?
real	thresh			# k-sigma threshold

pointer	bufx, bufy
int	nbuffer
int	i,j, i1,j1, i2,j2
int	ic, jc			# re-centered coordinates
real	mean1, mean2		# averages along inner and outer circles.
real	std1, std2		# std. dev. along inner and outer circles.
real	fom			# figure of merit for detection
real	aux

real	el_p2s(), el_s2p()

begin
	# First, locate probable object center. It is either pointed
	# to by valid XC,YC pair, or sits in frame center.
	if ((!IS_INDEFR (XC(is))) && (!IS_INDEFR (YC(is)))) {

	    # Center coordinates are defined. Must check if they
	    # point to somewhere inside the frame. If not, set
	    # then to frame center.
	    if (PHYSICAL(is)) {
	        XC(is) = el_p2s (im, XC(is), 1)
	        YC(is) = el_p2s (im, YC(is), 2)
	    }
            if ((XC(is) < 1.) || (XC(is) > IM_LEN(im,1)))
	        XC(is) = IM_LEN(im,1) / 2
            if ((YC(is) < 1.) || (YC(is) > IM_LEN(im,2)))
	        YC(is) = IM_LEN(im,2) / 2
	} else {

	    # If center coordinates not defined, assume object is 
	    # in frame center. 
	    XC(is) = IM_LEN(im,1) / 2
	    YC(is) = IM_LEN(im,2) / 2
	}

	# Check to see if valid object is there.
	call malloc (bufx, SZ_BUFFER, TY_REAL)
	call malloc (bufy, SZ_BUFFER, TY_REAL)
	nbuffer = SZ_BUFFER

	# Limits of re-centering window.
	i1 = max (1, int(XC(is)) - RWINDOW / 2)
	j1 = max (1, int(YC(is)) - RWINDOW / 2)
	i2 = min (IM_LEN(im,1), int(XC(is)) + RWINDOW / 2)
	j2 = min (IM_LEN(im,2), int(YC(is)) + RWINDOW / 2)
	ic = 0
	jc = 0
	fom = 0.

	if (list) {
	    call printf ("Running object locator... ")
	    call flush (STDOUT)
	}

	# Scan window.
	do j = j1, j2 {
	    do i = i1, i2 {

	        # Extract two concentric circular samples. 
	        call el_get (im, sec, real(i), real(j), INNER_RADIUS, 0.0, 0.0, 
	                     bufx, bufy, nbuffer, NPOINT(is), NDATA(is), 
	                     mean1, std1, ASTEP(al), LINEAR(al), 
	                     INT_LINEAR, 4., 4.0, 0, SAREA(is))
	        call el_get (im, sec, real(i), real(j), OUTER_RADIUS, 0.0, 0.0, 
	                     bufx, bufy, nbuffer, NPOINT(is), NDATA(is), 
	                     mean2, std2, ASTEP(al), LINEAR(al), 
	                     INT_LINEAR, 4., 4.0, 0, SAREA(is))

	        # Figure of merit measures if there is reasonable
	        # signal at position i,j.
	        if (IS_INDEF(std1)) std1 = 0.0
	        if (IS_INDEF(std2)) std2 = 0.0
	        aux = std1 * std1 + std2 * std2



	        if (aux > 0.0) {
	            aux = (mean1 - mean2) / sqrt(aux)
	            if (aux > fom) {
	                fom = aux
	                ic  = i
	                jc  = j
	            }
	        }
	    }
	}
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)

	if (list)
	    call printf ("Done.\n")

	# If valid object, re-center if asked for. Otherwise, no-detection
	# is signaled by setting center coordinates to INDEF.
	if (fom > thresh) {
	    if (recenter) {
	        XC(is) = real (ic)
	        YC(is) = real (jc)
	    }
	} else {
	    XC(is) = INDEFR
	    YC(is) = INDEFR
	}

	# Restore coordinates to physical system if needed.
	if ((!IS_INDEFR (XC(is))) && (!IS_INDEFR (YC(is)))) {
	    if (PHYSICAL(is)) {
	        XC(is) = el_s2p (im, XC(is), 1)
	        YC(is) = el_s2p (im, YC(is), 2)
	    }
	}
end

