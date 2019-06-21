# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imhdr.h>
include	"ellipse.h"

# EL_INTEGR -- Compute integrated flux inside current ellipse, as well
# as inside circle defined by semi-major axis. Pixels in a square section
# enclosing circle are scanned; the distance of each pixel to the isophote
# center is compared both with the semi-major axis lenght and with the
# lenght of the ellipse radius vector, and integrals are updated if the
# pixel distance is smaller. 

procedure el_integr (im, sec, is)

pointer	im				# IMIO pointer
pointer	sec				# in-memory pixel structure
pointer	is				# isophote structure

int	i,j
int	imin,jmin, imax,jmax
real	dp, de, steta, cteta, sphi, cphi
real	dx, dy
double	tfe, tfc
bool	flag

begin
	steta = sin (-TETA(is))
	cteta = cos (-TETA(is))
	tfe = 0.D0		# Integrated intensities.
	tfc = 0.D0
	TFEA(is) = 0		# Number of valid pixels 
	TFCA(is) = 0		# in each integral.

	# Compute limits of square array that encloses circle.
	imin = max (1,            int (XC(is) - A(is) - 0.5) - 1)
	jmin = max (1,            int (YC(is) - A(is) - 0.5) - 1)
	imax = min (IM_LEN(im,1), int (XC(is) + A(is) + 0.5) + 1)
	jmax = min (IM_LEN(im,2), int (YC(is) + A(is) + 0.5) + 1)

	# Scan pixel array.
	do j = jmin, jmax {
	    do i = imin, imax {

	        # Distance of pixel to isophote center.
	        dx = real(i - int(XC(is) + 0.5))
	        dy = real(j - int(YC(is) + 0.5))
	        dp = dx*dx + dy*dy
	        if (dp > 0.)
	            dp = sqrt (dp)

	        # Update circle photometry.
	        flag = false
	        if (dp <= A(is)) {
	            call el_getpix (im, sec, i, j)
	            if (!IS_INDEFR (Memr[SUBRASTER(sec)])) {
	                tfc = tfc + double(Memr[SUBRASTER(sec)]) 
	                TFCA(is) = TFCA(is) + 1
	                flag = true
	            }
	        }	        

	        # Ellipse distance to isophote center (radius vector).
	        de = 0.
	        if (dp > 0.) {
	            sphi = dy / dp
	            cphi = dx / dp
	            de = A(is) * (1. - EPS(is))  / 
	                 sqrt (((cphi*cteta - sphi*steta) * (1. - EPS(is)))**2 +
	                       (sphi*cteta + cphi*steta)**2)
	        }

	        # Update ellipse photometry.
	        if (dp <= de) {
	            if (flag) {
	                tfe = tfe + double(Memr[SUBRASTER(sec)]) 
	                TFEA(is) = TFEA(is) + 1
	            } else {
	                call el_getpix (im, sec, i, j)
	                if (!IS_INDEFR (Memr[SUBRASTER(sec)])) {
	                    tfe = tfe + double(Memr[SUBRASTER(sec)]) 
	                    TFEA(is) = TFEA(is) + 1
	                }
	            }
	        }
	    }
	}

	TFE(is) = real (tfe)
	TFC(is) = real (tfc)
end
