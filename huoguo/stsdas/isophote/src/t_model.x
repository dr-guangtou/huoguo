# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imhdr.h>
include <imset.h>
include <tbset.h>

include	"ellipse.h"
include	"elcolnames.h"

define	PI		3.141592654
define	DPI		(2. * PI)
define	EPSILON		1.e-8
define	MARGIN		10
define	PIXEL		Memr[$1+($3-jmin)*ni+($2-imin)]
define	SPIXEL		Memr[$1+$3*ni+$2]

# T_MODEL -- Builds model galaxy image from isophote table.
# This table must be equally spaced in semi-major axis length,
# whith step < 1. pixel, in order to get a proper area coverage of
# the output image. The table must be sorted in order of increasing
# semi-major axis length. This task was designed to work together with
# the script task bmodel.cl, which takes care of table sorting and
# interpolation.
#
# The algorithm scans the input table, and, for each ellipse in there, 
# fills up the output image array with the corresponding isophotal
# intensity. Pixels in the target array are in general only partially 
# covered by the isophote "pixel"; the algorithm takes care of this
# partial pixel overlap, and keeps track of how much intensity was
# added to each pixel by storing the partial area information in an
# auxiliary array. The information in this array is used later to
# normalize the pixel intensities.

procedure t_model ()

pointer	im, parim, bufi, bufo		# image pointers
pointer	area				# fractional pixel information
pointer	tp, cptr			# table pointers
char	newname[SZ_PATHNAME]		# new image name
char	parname[SZ_PATHNAME]		# parent image name
char	table[SZ_PATHNAME]		# table whith isophotal data
real	backgr				# background
bool	highar				# add higher harmonics ?

int	dimx, dimy			# image size
int	i, j, row, maxrow
int	ic, jc				# galaxy center coordinates
int	imin, imax, jmin, jmax		# coord. of subraster whith model
int	ni, nj				# subraster size
real	r, phi				# polar coordinates
real	a0				# semi-major axis length
real	x0, y0				# center of isophote on frame
real	eps0, teta0			# elipticity, pos. angle
real	slope				# local intensity slope
real	a3, b3, a4, b4			# higher-order harmonics
real	int0				# isophotal intensity
real	intc				# center intensity
real	x, y, fx, fy, aux
long	cpu, clock			# time variables
bool	verbose

pointer	immap(), impl2r(), imps2r(), imgs2r()
pointer	tbtopn()
int	tbpsta(), strlen()
real	clgetr()
bool	clgetb()
long	cputime(), clktime()

begin
	call clgstr ("table",    table,   SZ_PATHNAME)
	call clgstr ("output",   newname, SZ_PATHNAME)
	call clgstr ("parent",   parname, SZ_PATHNAME)
	backgr   = clgetr ("backgr")
	highar   = clgetb ("highar")
	verbose  = clgetb ("verbose")
	cpu   = cputime (0)
	clock = clktime (0)

	ic = 0
	jc = 0

	# Open isophote table.
	tp = tbtopn (table, READ_ONLY, 0)
	maxrow = tbpsta (tp, TBL_NROWS)

	# Create new image; fill it whith background value. If no parent 
	# image was specified, assume its name is in table header.
	if (verbose) {
	    call printf ("Creating image and filling background value...\n")
	    call flush (STDOUT)
	}
	if (strlen(parname) == 0)
	    call tbhfkr (tp, "IMAGE", i, parname, j)
	parim = immap (parname, READ_ONLY, 0)
	im    = immap (newname, NEW_COPY, parim)
	dimx  = IM_LEN (im, 1)
	dimy  = IM_LEN (im, 2)
	IM_PIXTYPE(im) = TY_REAL

	bufo = impl2r (im, 1)
	do i = 1, dimy {
	    call amovkr (backgr, Memr[bufo], dimx)
	    bufo = impl2r (im, i)
	}
	call imunmap (im)
	call imunmap (parim)

	if (verbose) {
	    call printf ("Building model...\n")
	    call flush (STDOUT)
	}

	# Reopen image.
	im = immap (newname, READ_WRITE, 28800)
	call imsetr (im, IM_ADVICE, real(RANDOM))
	call imseti (im, IM_NBNDRYPIX, 3)
	call imseti (im, IM_TYBNDRY, BT_REFLECT)

	# Read last (largest radius) isophote data.
	call tbcfnd (tp, ES_CA, cptr, 1)
	call tbegtr (tp, cptr, maxrow, a0)
	call tbcfnd (tp, ES_CEPS, cptr, 1)
	call tbegtr (tp, cptr, maxrow, eps0)
	call tbcfnd (tp, ES_CTETA, cptr, 1)
	call tbegtr (tp, cptr, maxrow, teta0)
	teta0 = teta0 / 180. * PI + PI2
	call tbcfnd (tp, ES_CX0, cptr, 1)
	call tbegtr (tp, cptr, maxrow, x0)
	call tbcfnd (tp, ES_CY0, cptr, 1)
	call tbegtr (tp, cptr, maxrow, y0)

	# Set subraster limits.
	phi = 0.
	r   = a0
	imin = IM_LEN(im, 1)
	jmin = IM_LEN(im, 2)
	imax = 1
	jmax = 1
	while (phi < DPI) {
	    # Get image coordinates of (r, phi) pixel
	    i = int (r * cos (phi + teta0) + x0 + 0.5)
	    j = int (r * sin (phi + teta0) + y0 + 0.5)
	    # Set limits.
	    if ((i > imax) && (i <= IM_LEN(im,1)))
	        imax = i
	    if ((j > jmax) && (j <= IM_LEN(im,2)))
	        jmax = j
	    if ((i < imin) && (i > 0))
	        imin = i
	    if ((j < jmin) && (j > 0))
	        jmin = j
	    # Step over ellipse.
	    phi = phi + 1. / r
	    r = a0 * (1. - eps0) / sqrt (((1. - eps0) * cos (phi))**2 + 
	        (sin (phi))**2)
	}
	imin = imin - MARGIN
	if (imin < 1)
	    imin = 1
	jmin = jmin - MARGIN
	if (jmin < 1)
	    jmin = 1
	imax = imax + MARGIN
	if (imax > IM_LEN(im,1))
	    imax = IM_LEN(im,1)
	jmax = jmax + MARGIN
	if (jmax > IM_LEN(im,2))
	    jmax = IM_LEN(im,2)
	ni = imax - imin + 1
	nj = jmax - jmin + 1

	# Read subraster and allocate area for fractional pixel information.
	bufi = imgs2r (im, imin, imax, jmin, jmax)
	bufo = imps2r (im, imin, imax, jmin, jmax)
	do j = jmin, jmax {
	    do i = imin, imax {
	        PIXEL(bufo,i,j) = PIXEL(bufi,i,j)
	    }
	}
	call calloc (area, ni * nj, TY_REAL)

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
	    call tbcfnd (tp, ES_CINT, cptr, 1)
	    call tbegtr (tp, cptr, row, int0)
	    if (highar) {
	        call tbcfnd (tp, ES_CSLOPE, cptr, 1)
	        call tbegtr (tp, cptr, row, slope)
	        call tbcfnd (tp, ES_CA3, cptr, 1)
	        call tbegtr (tp, cptr, row, a3)
	        call tbcfnd (tp, ES_CB3, cptr, 1)
	        call tbegtr (tp, cptr, row, b3)
	        call tbcfnd (tp, ES_CA4, cptr, 1)
	        call tbegtr (tp, cptr, row, a4)
	        call tbcfnd (tp, ES_CB4, cptr, 1)
	        call tbegtr (tp, cptr, row, b4)
	        # Return deviations from ellipticity to 
	        # their original intensity amplitude meaning.
	        a3 = - a3 * slope * a0
	        b3 = - b3 * slope * a0
	        a4 = - a4 * slope * a0
	        b4 = - b4 * slope * a0
	    }

	    # Print verbose info.
	    if (verbose) {
	        call printf ("\r                  ")
	        call printf ("\ra0 = %g")
	            call pargr (a0)
	        call flush (STDOUT)
	    }

	    # Generate isophote.
	    phi = 0.
	    r   = a0

	    if (r > 0.) {
	        while (phi <= DPI) {

	            # Get image coordinates of (r, phi) pixel
	            x = r * cos (phi + teta0) + x0
	            y = r * sin (phi + teta0) + y0
	            i = int (x)
	            j = int (y)

	            # If outside subraster boundaries, ignore.
	            if ((i > imin) && (i < imax) &&
	                (j > jmin) && (j < jmax)) {

	                # Get fractional deviations of isophote pixel
		        # relative to the target array
	                fx = x - real (i)	        
	                fy = y - real (j)	        

	                # Add up the isophote contribution to the nearest
	                # pixels in the target array
	                PIXEL(bufo,i,j)     = PIXEL(bufo,i,j) + 
		 		              int0 * (1. - fx) * (1. - fy)
	                PIXEL(bufo,i+1,j)   = PIXEL(bufo,i+1,j) + 
					      int0 *       fx  * (1. - fy)
	                PIXEL(bufo,i,j+1)   = PIXEL(bufo,i,j+1) + 
					      int0 * (1. - fx) *       fy
	                PIXEL(bufo,i+1,j+1) = PIXEL(bufo,i+1,j+1) + 
					      int0 *       fx  *       fy

	                # Add up the fractional area contribution to the nearest
	                # pixels in the target array area map
	                PIXEL(area,i,j)     = PIXEL(area,i,j) + 
					      (1. - fx) * (1. - fy)
	                PIXEL(area,i+1,j)   = PIXEL(area,i+1,j) +
					            fx  * (1. - fy)
	                PIXEL(area,i,j+1)   = PIXEL(area,i,j+1) +
					      (1. - fx) *       fy
	                PIXEL(area,i+1,j+1) = PIXEL(area,i+1,j+1) +
					            fx  *       fy

	                # Add higher harmonics.
	                if (highar) { 
	                    aux = a3 * sin (3. * phi) + b3 * cos (3. * phi) +
	                          a4 * sin (4. * phi) + b4 * cos (4. * phi) / 4.
	                    PIXEL(bufo,i,j)     = PIXEL(bufo,i,j)     + aux
	                    PIXEL(bufo,i+1,j)   = PIXEL(bufo,i+1,j)   + aux
	                    PIXEL(bufo,i,j+1)   = PIXEL(bufo,i,j+1)   + aux 
	                    PIXEL(bufo,i+1,j+1) = PIXEL(bufo,i+1,j+1) + aux 
	                }

	                PIXEL(bufi,i,j)     = PIXEL(bufo,i,j)
	                PIXEL(bufi,i+1,j)   = PIXEL(bufo,i+1,j)
	                PIXEL(bufi,i,j+1)   = PIXEL(bufo,i,j+1)
	                PIXEL(bufi,i+1,j+1) = PIXEL(bufo,i+1,j+1)
	            }

	            # Step over ellipse.
	            phi = phi + 0.75 / r
	            r = a0 * (1. - eps0) / sqrt (((1. - eps0) * cos (phi))**2 + 
		        (sin (phi))**2)
	        }
	    } else {

	         # SMA = 0. Store center coordinates and
	         # intensity to set center pixel later.
	         ic   = int (x0 + 0.5)
	         jc   = int (y0 + 0.5)
	         intc = int0
	    }
	}

	call tbtclo (tp)

	# Normalize for partial pixel areas
	do j = jmin, jmax {
	    do i = imin, imax {
	        if (PIXEL(area,i,j) > 0.)
	            PIXEL(bufo,i,j) = (PIXEL(bufo,i,j) - backgr) /
	                               PIXEL(area,i,j)
	        else
	            PIXEL(bufo,i,j) = backgr
	    }
	}

	# Set center pixel.
	if ((ic > 0) && (jc > 0))
	    PIXEL(bufo,ic,jc) = intc

	call mfree (area, TY_REAL)
	call imunmap (im)

	if (verbose) {
	    call printf ("\n%7.2f  CPU seconds.\n%7.2f  minutes elapsed.\n")
	        call pargr (real (cputime (cpu)) / 1000.)
	        call pargr (real (clktime (clock)) / 60.)
	}
end
