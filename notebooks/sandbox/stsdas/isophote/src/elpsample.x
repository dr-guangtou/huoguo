# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<gset.h>
include	"ellipse.h"

# EL_PSAMPLE  --  Plot azimuthal intensity sample. 

procedure el_psample (im, is, al, dev, x, y) 

pointer	im				# IMIO
pointer	is				# isophote structure
pointer	al				# algorithm structure
char	dev[ARB]			# device
real	x[ARB], y[ARB]			# angle and residual vectors

pointer	xs, ys
pointer	str
int	i
real	x1, x2, y1,y2			# WCS

pointer	gp, gopen()
real	gstatr()

begin
	call malloc (str, SZ_LINE, TY_CHAR)

	# Process plot vectors: add mean intensity to 
	# residuals, set angle origin and sort angle.
	call malloc (xs, NPOINT(is), TY_REAL)
	call malloc (ys, NPOINT(is), TY_REAL)
	do i = 1, NPOINT(is) {
	    if (AANGLE(al))
	        Memr[xs+i-1] = x[i] + TETA(is) - PI2
	    else
	        Memr[xs+i-1] = x[i]
	    Memr[xs+i-1] = Memr[xs+i-1] / PI * 180.
	    if (!IS_INDEF(y[i]))
	        Memr[ys+i-1] = y[i] + MEAN(is)
	    else
	        Memr[ys+i-1] = INDEFR
	}

	gp = gopen (dev, NEW_FILE, STDGRAPH)

	# Scale graphics.
	call gascale (gp, Memr[xs], NPOINT(is), 1) 
	call gascale (gp, Memr[ys], NPOINT(is), 2) 
	if (AANGLE(al))
	    call glabax  (gp, "", "PA (degrees)", "Intens.")
	else
	    call glabax  (gp, "", "PAe (degrees)", "Intens.")

	# Plot sample.
	call gpline (gp, Memr[xs], Memr[ys], NPOINT(is))

	# Write id data.
	call ggwind  (gp, x1, x2, y1, y2)
	call gsetr (gp, G_TXSIZE, 2. * gstatr (gp, G_TXSIZE))
	y2 = 0.9 * (y2 - y1) + y1
	call el_unshrink (im, is, al)
	call sprintf (Memc[str], SZ_LINE, "SMA = %7.2f")
	    call pargr (A(is))
	call el_shrink (im, is, al)
	call gtext (gp, ((Memr[xs] + Memr[xs+NPOINT(is)-1])/2.), y2, 
	            Memc[str], "")
	call mfree (ys, TY_REAL)
	call mfree (xs, TY_REAL)
	call gflush (gp)
	call gclose (gp)
	call mfree (str, TY_CHAR)
end
