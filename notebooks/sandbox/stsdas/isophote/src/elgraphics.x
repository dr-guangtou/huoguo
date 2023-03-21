# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <gio.h>
include	"eldisplay.h"
include "ellipse.h"

# Routines for ploting at image/graphics display, which emulate
# GIO routines but disguise differences between standard graphics 
# and the image server. All input coordinates are in section system.



procedure el_gamove (im, dp, gp, x, y)

pointer	im				# IMIO pointer
pointer	dp				# display structure
pointer	gp				# GIO pointer
real	x,y				# section coordinates

real	el_s2d(), el_s2g()

begin
	switch (DDEV(dp)) {
	case G_DEV:
	    call gamove (gp, el_s2g (im, dp, x, 1), el_s2g (im, dp, y, 2))
	case I_DEV:
	    call gamove (gp, el_s2d (im, dp, x, 1), el_s2d (im, dp, y, 2))
	}
end




procedure el_gadraw (im, dp, gp, x, y)

pointer	im				# IMIO pointer
pointer	dp				# display structure
pointer	gp				# GIO pointer
real	x,y				# section coordinates

real	el_s2d(), el_s2g()

begin
	switch (DDEV(dp)) {
	case G_DEV:
	    call gadraw (gp, el_s2g (im, dp, x, 1), el_s2g (im, dp, y, 2))
	case I_DEV:
	    call gadraw (gp, el_s2d (im, dp, x, 1), el_s2d (im, dp, y, 2))
	}
end




procedure el_gfill (im, dp, gp, x, y, npts, style)

pointer	im				# IMIO pointer
pointer	dp				# display structure
pointer	gp				# GIO pointer
real	x[ARB],y[ARB]			# polygon vertices
int	npts				# number of vertices
int	style				# GIO fill pattern

pointer	bx, by
int	i

real	el_s2d(), el_s2g()

begin
	call malloc (bx, npts, TY_REAL)
	call malloc (by, npts, TY_REAL)

	do i = 1, npts {
	    switch (DDEV(dp)) {
	    case G_DEV:
	        Memr[bx+i-1] = el_s2g (im, dp, x[i], 1)
	        Memr[by+i-1] = el_s2g (im, dp, y[i], 2)
	    case I_DEV:
	        Memr[bx+i-1] = el_s2d (im, dp, x[i], 1)
	        Memr[by+i-1] = el_s2d (im, dp, y[i], 2)
	    }
	}

	call gfill (gp, Memr[bx], Memr[by], npts, style)

	call mfree (bx, TY_REAL)
	call mfree (by, TY_REAL)
end


