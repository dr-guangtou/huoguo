# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imhdr.h>
include <imset.h>
include <tbset.h>
include	<ctype.h>

include "ellipse.h"
include	"elcolnames.h"

define	PI		3.141592654
define	DPI		(2. * PI)

# T_MAP -- Generates elliptical contours in the form of x,y pair
# lists. Input is a SDAS table whith isophote parameters, output
# is a set of ASCII files suitable for input in task SGRAPH.
# This task was designed to work togheter whith the script tasks
# isomap.cl and isoimap.cl.

procedure t_map ()

pointer	tp, cptr			# table pointers
char	image[SZ_FNAME]			# image name
char	table[SZ_FNAME]			# table whith isophotal data
char	output[SZ_FNAME]		# output files
int	normalize			# absolute or viewport coordinates ?

pointer	im				# IMIO pointer
char	str[SZ_FNAME]
char	section[SZ_FNAME]
int	i, row, maxrow
int	out				# output file pointer
long	x1,x2,y1,y2,stepx,stepy		# section specs
long	ox1,ox2,oy1,oy2,ostepx,ostepy	# original section specs
real	r, fi				# polar coordinates
real	x, y				# image coordinates
real	a0				# semi-major axis length
real	x0, y0				# center of isophote on frame
real	eps0, teta0			# elipticity, pos. angle
real	large, small			# section dimensions
real	offx, offy

pointer	tbtopn(), immap()
int	tbpsta(), open(), strlen(), clgeti()
int	checkdim()

begin
	call clgstr ("image",  image,  SZ_FNAME)
	call clgstr ("table",  table,  SZ_FNAME)
	call clgstr ("output", output, SZ_FNAME)
	normalize = clgeti ("normalize")	

	# Get image section dimensions
	large = 1.
	small = 1.
	offx  = 0.
	offy  = 0.
	if (normalize == 1) {
	    im = immap (image, READ_ONLY, 0)

	    # change IM_NDIM(im) to checkdim(), JC Hsu 11/29/94
	    if (checkdim(im) != 2) {
	        call imunmap (im)
	        call error (0, "Input image section is not two-dimensional")
	    }
	    large = real (IM_LEN(im, 1))
	    small = real (IM_LEN(im, 2))
	    if (small > large) {
	        large = real (IM_LEN(im, 2))
	        small = real (IM_LEN(im, 1))
	    }
	    if (IM_LEN(im, 1) < IM_LEN(im, 2))
	        offx = (1. - small/large) / 2.
	    else
	        offy = (1. - small/large) / 2.
	    call imunmap (im)
	}

	# Treatment of image section specification.
	call imgsection (image, section, SZ_FNAME)
	if (strlen(section) > 0) {
	    i = 2
	    while (IS_WHITE(section[i]))
	        i = i + 1
	    call im_decode_subscript (section, i, x1, x2, stepx)
	    call im_decode_subscript (section, i, y1, y2, stepy)
	} else {
	    x1    = long(1)
	    y1    = long(1)
	    stepx = long(1)
	    stepy = long(1)
	}

	# Open isophote table.
	tp = tbtopn (table, READ_ONLY, 0)
	maxrow = tbpsta (tp, TBL_NROWS)
	call tbhgtt (tp, "IMAGE", str, SZ_FNAME)

	# Treatment of original image section specification.
	call imgsection (str, section, SZ_FNAME)
	if (strlen(section) > 0) {
	    i = 2
	    while (IS_WHITE(section[i]))
	        i = i + 1
	    call im_decode_subscript (section, i, ox1, ox2, ostepx)
	    call im_decode_subscript (section, i, oy1, oy2, ostepy)
	} else {
	    ox1    = long(1)
	    oy1    = long(1)
	    ostepx = long(1)
	    ostepy = long(1)
	}

	# Scan table.
	do row = 1, maxrow {

	    # Get row data.
	    call tbcfnd (tp, ES_CA, cptr, 1)
	    call tbegtr (tp, cptr, row, a0)
	    call tbcfnd (tp, ES_CEPS, cptr, 1)
	    call tbegtr (tp, cptr, row, eps0)
	    call tbcfnd (tp, ES_CTETA, cptr, 1)
	    call tbegtr (tp, cptr, row, teta0)
	    teta0 = teta0 / 180. * PI + PI2
	    call tbcfnd (tp, ES_CX0, cptr, 1)
	    call tbegtr (tp, cptr, row, x0)
	    call tbcfnd (tp, ES_CY0, cptr, 1)
	    call tbegtr (tp, cptr, row, y0)

	    # Generate curve.

	    call sprintf (str, SZ_FNAME, "%s.%03d")
	        call pargstr (output)
	        call pargi (row)
	    out = open (str, NEW_FILE, TEXT_FILE)

	    fi = 0.
	    r  = a0

	    if (a0 > 0.5) {
		
		# add 0.01 to DPI so the ellipse is nearer to a closed loop.
		#   JC Hsu 6/15/94
	        while (fi <= DPI+0.01) {

	            # Get image coordinates of (r, fi) pixel
	            x = r * cos (fi + teta0) + x0 + 1.
	            y = r * sin (fi + teta0) + y0 + 1.

	            # Convert them to original image coordinates,
	            # without windowing or subsampling.
	            x = (x - 1.) * float(ostepx) + float (ox1)
	            y = (y - 1.) * float(ostepy) + float (oy1)

	            # Print
	            call fprintf (out, "%g %g\n")
	                call pargr ((x - float(x1)) / float(stepx) / large + offx)
	                call pargr ((y - float(y1)) / float(stepy) / large + offy)

	            # Step over ellipse.
		    # reduce the step from 2. to 0.3 to get better resolution.
		    #   JC Hsu 6/15/94
	            fi = fi + 0.3 / r
	            r = a0 * (1. - eps0) / sqrt (((1. - eps0) * cos (fi))**2 + 
		        (sin (fi))**2)
	        }
	    } else {
	        # Print fake ellipse.
	        call fprintf (out, "%g %g\n")
	            call pargr ((x0 - float(x1)) / float(stepx) / large + offx)
	            call pargr ((y0 - float(y1)) / float(stepy) / large + offy)
	    }
	    call close (out)
	}

	call tbtclo (tp)
end
