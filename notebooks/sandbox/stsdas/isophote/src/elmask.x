# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<imhdr.h>
include	"ellipse.h"
include	"eldisplay.h"


# EL_MASK  --  Flag region and plot it. Both square regions denoted by
# central coordinate pair and PMASK size, and rectangular regions denoted
# by two coordinate pairs, are allowed.

procedure el_mask (im, sec, al, dp, cx, cy, cxa, cya)

pointer	im		# IMIO pointer
pointer	sec		# in-memory pixel structure
pointer	al		# algorithm control structure
pointer	dp		# displayed section structure
real	cx,  cy		# first corner / center coordinates
real	cxa, cya	# second corner coordinates

int	i1, i2, j1, j2
int	i, j
pointer	gp, el_gopen()

begin
	# Compute region limits. 
	call el_limmsk (im, al, cx, cy, cxa, cya, i1, i2, j1, j2)

	gp = el_gopen (DDEV(dp))

	# Do it.
	do j = j1, j2 {
	    do i = i1, i2 {
	        call el_pflag (im, sec, i, j, true)
	        call el_ppix (im, dp, gp, real(i), real(j))
	    }
	}

	call gflush (gp)
	call gclose (gp)
end



# EL_UNMASK  --  Un-flag region. Do not update display.

procedure el_unmask (im, sec, al, cx, cy, cxa, cya)

pointer	im		# IMIO pointer
pointer	sec		# in-memory pixel structure
pointer	al		# algorithm control structure
real	cx,  cy		# first corner / center coordinates
real	cxa, cya	# second corner coordinates

int	i1, i2, j1, j2
int	i, j

begin
	# Compute region limits. 
	call el_limmsk (im, al, cx, cy, cxa, cya, i1, i2, j1, j2)

	# Do it.
	do j = j1, j2 {
	    do i = i1, i2 {
	        call el_punf (im, sec, i, j)
	    }
	}
end



# EL_LIMMSK  --  Compute mask limits.

procedure el_limmsk (im, al, cx, cy, cxa, cya, i1, i2, j1, j2)

pointer	im		# i: IMIO pointer
pointer	al		# i: algorithm control structure
real	cx,  cy		# i: first corner / center coordinates
real	cxa, cya	# i: second corner coordinates
int	i1,i2, j1,j2	# o: mask limits

int	i

begin
	if (REGION(al)) {
	    i1 = int (cx)
	    i2 = int (cxa)
	    j1 = int (cy)
	    j2 = int (cya)
	    if (i1 > i2) {
	        i  = i2
	        i2 = i1
	        i1 = i
	    }
	    if (j1 > j2) {
	        i  = j2
	        j2 = j1
	        j1 = i
	    }
	} else {
	    i = PMASK(al) / 2 + 1
	    i1 = max (int (cx + 0.5) - i + 1, 1)
	    i2 = min (int (cx + 0.5) + i - 1, IM_LEN(im,1))
	    j1 = max (int (cy + 0.5) - i + 1, 1)
	    j2 = min (int (cy + 0.5) + i - 1, IM_LEN(im,2))
	}
end
