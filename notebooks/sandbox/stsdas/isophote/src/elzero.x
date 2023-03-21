# Copyright restrictions apply - see stsdas$copyright.stsdas 

include "ellipse.h"

# EL_ZERO -- Special case of zero semi-major axis.

procedure el_zero (im, sec, is, al, list)

pointer	im				# IMIO pointer
pointer	sec				# in-memory pixel structure
pointer	is				# isphote structure
pointer	al				# algorithm control structure
bool	list				# list at STDOUT ?

int	ic, jc, i
real	mean, sma

begin
	# Central pixel coordinates.
	ic = int (XC(is) + 0.5)
	jc = int (YC(is) + 0.5)

	# Estimate central gradient. First, undo
	# SMA update performed by main loop.
        if (LINEAR(al))
            sma = A(is) + ASTEP(al)
        else
            sma = A(is) * ASTEP(al)
	mean = MEAN(is)
	call el_getpix (im, sec, ic, jc)
	MEAN(is) = Memr[SUBRASTER(sec)]
	if (!IS_INDEF(MEAN(is)) && (sma > 0.0)) {
	    SLOPE(is) = (mean - MEAN(is)) / sma
	    NDATA(is) = 1
	} else {
	    SLOPE(is) = INDEFR
	    NDATA(is) = 0
	}

	A(is)       = 0.
	NPOINT(is)  = 1

	# Set everything else to INDEF.
	EPS(is)     = INDEFR
	TETA(is)    = INDEFR
	EX(is)      = INDEFR
	EY(is)      = INDEFR
	EEPS(is)    = INDEFR
	ETETA(is)   = INDEFR
	SIGMA(is)   = INDEFR
	ESLOPE(is)  = INDEFR
	RESLOPE(is) = INDEFR
	SAREA(is)   = INDEFR
	ABIG(is)    = INDEFR
	TFE(is)     = INDEFR
	TFC(is)     = INDEFR
	TFEA(is)    = INDEFI
	TFCA(is)    = INDEFI
	STOP(al)    = INDEFI
	NITER(al)   = INDEFI
	do i = 0, 3 {
	    Memr[HH(is)+i]  = INDEFR
	    Memr[EHH(is)+i] = INDEFR
	}
	if (NHARM(is) > 0) {
	    do i = 1, NHARM(is), 2 {
	        AI(HARM(is),i)  = INDEFR
	        AI(EHARM(is),i) = INDEFR
	        BI(HARM(is),i)  = INDEFR
	        BI(EHARM(is),i) = INDEFR
	    }
	}
	# Output basic info to STDOUT.
	if (list) {
	    call printf ("%7.2f %8.2f\n")
	        call pargr (A(is))
	        call pargr (MEAN(is))
	}
end






