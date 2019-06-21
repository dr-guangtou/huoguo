# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <gset.h>
include "eldisplay.h"

# EL_PPIX  --  Plot (enhance) a single pixel. Graphics device
#              must be opened and closed by caller.

procedure el_ppix (im, dp, gp, x, y)

pointer	im				# IMIO
pointer	dp				# display
pointer	gp				# GIO
real	x, y				# pixel coordinates

real	px[5], py[5]			# vertices
real	ax, ay

begin
	# Rounding correction.
	ax = real (int (x + 0.5))
	ay = real (int (y + 0.5))

	# Fill up vertex vector.
	px[1] = ax - 0.5
	py[1] = ay - 0.5
	
	px[2] = ax - 0.5
	py[2] = ay + 0.5
	
	px[3] = ax + 0.5
	py[3] = ay + 0.5
	
	px[4] = ax + 0.5
	py[4] = ay - 0.5
	
	px[5] = ax - 0.5
	py[5] = ay - 0.5
	
	# Do it.
	call el_gfill (im, dp, gp, px, py, 5, GF_HOLLOW)
end
