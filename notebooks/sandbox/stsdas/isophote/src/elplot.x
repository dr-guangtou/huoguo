# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<imhdr.h>
include	"ellipse.h"
include	"eldisplay.h"

# EL_PLOT  --  Plot current ellipse stored in isophote structure.

procedure el_plot (im, gp, dp, is)

pointer	im				# IMIO pointer
pointer	gp				# GIO pointer
pointer dp                              # display control parameters
pointer is				# isophote structure

real	phi, r, x, y

begin
	if (IS_INDEF(XC(is)) || IS_INDEF(XC(is)))
	    return

	phi = 0.
	r   = A(is)

	# Move pen to beginning
	x   = (r * cos (TETA(is)) + XC(is))
	y   = (r * sin (TETA(is)) + YC(is))
	call el_gamove (im, dp, gp, x, y)

	# Scan phi angle, leaving some superposition
	# at the end.
	while (phi <= (1.1 * DPI)) {

	    # Get image coordinates of (r, phi) pixel
	    x = (r * cos (phi + TETA(is)) + XC(is))
	    y = (r * sin (phi + TETA(is)) + YC(is))

	    # Plot
	    call el_gadraw (im, dp, gp, x, y)

	    # Step over ellipse.
	    phi = phi + 2./ r
	    r = A(is) * (1. - EPS(is)) / sqrt (((1. - EPS(is)) * cos (phi))**2 + 
		(sin (phi))**2)
	}
end
