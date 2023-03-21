# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imio.h>
include	"ellipse.h"

# EL_DUMP -- Write isophote data to main output table.

procedure el_dump (is, al, im, outname, tp, column, extracols, rownumber, m0, 
                   refer, backgr)

pointer	is				# isophote data
pointer	al				# algorithm data
pointer	im				# IMIO
char	outname[ARB]			# output file name
pointer	tp, column[ARB]			# pointers for SDAS tables
pointer	extracols			# pointers for opt. harm. columns
int	rownumber			# SDAS table row number
real	m0				# magnitude of reference
real	refer				# reference intensity in one pixel
real	backgr				# "background" value in one pixel

int	i, j
real	ebar, mag, tmage, tmagc 
real	emagl, emags
real	pixvar, isoerr
real	aux

int	strlen()

begin
	rownumber = rownumber + 1

	# Compute isophote magnitude and errors.
	aux = (MEAN(is) - backgr) / (refer)
	if (!((aux) > 0.))
	    aux = 1. / (refer)
	mag = m0 - 2.5 * log10 (aux)

	if (A(is) > 0.) {
	    if (!IS_INDEF(SIGMA(is))) {
	        ebar = sqrt(SIGMA(is)**2 / NDATA(is))
	        aux = ((MEAN(is) - ebar) - backgr) / (refer)
	        if ((aux) > 0.)
	            emagl = -(mag - m0 + 2.5 * log10 (aux))
	        else
	            emagl = INDEFR
	        aux = ((MEAN(is) + ebar) - backgr) / (refer)
	        if ((aux) > 0.)
	           emags = mag - m0 + 2.5 * log10 (aux)
	        else
	            emags = INDEFR
	    } else {
	        emagl = INDEFR
	        emags = INDEFR
	    }
	} else {
	    emagl = INDEFR
	    emags = INDEFR
	}

	# Compute integrated magnitudes.
	if (A(is) > 0.) {
	    aux = (TFE(is) - backgr*TFEA(is)) / (refer)
	    if (!((aux) > 0.))
	        aux = 1. / (refer)
	    tmage = m0 - 2.5 * log10 (aux)
	    aux = (TFC(is) - backgr*TFCA(is)) / (refer)
	    if (!((aux) > 0.))
	        aux = 1. / (refer)
	    tmagc = m0 - 2.5 * log10 (aux)
	} else {
	    tmage = INDEFR
	    tmagc = INDEFR
	}

	# Compute pixel variance and isophote error.
	if (A(is) > 0.) {
	    pixvar = SIGMA(is) * sqrt (SAREA(is))
	    isoerr = SIGMA(is) / sqrt (float(NDATA(is)))
	} else {
	    pixvar = INDEFR
	    isoerr = INDEFR
	}

	# Translate to physical coord. system.
	call el_unshrink (im, is, al)

	# Output
	if (strlen (outname) > 0) {
	    call tbeptr(tp, column[1], rownumber, A(is))
	    call tbeptr(tp, column[2], rownumber, MEAN(is))
	    call tbeptr(tp, column[3], rownumber, isoerr)
	    call tbeptr(tp, column[4], rownumber, pixvar)
	    call tbeptr(tp, column[5], rownumber, SIGMA(is))
	    call tbeptr(tp, column[6], rownumber, EPS(is))
	    call tbeptr(tp, column[7], rownumber, EEPS(is))
	    call tbeptr (tp, column[8], rownumber, TETA(is))
	    call tbeptr (tp, column[9], rownumber, ETETA(is))
	    call tbeptr (tp, column[10], rownumber, XC(is))
	    call tbeptr (tp, column[11], rownumber, EX(is))
	    call tbeptr (tp, column[12], rownumber, YC(is))
	    call tbeptr (tp, column[13], rownumber, EY(is))
	    call tbeptr (tp, column[14], rownumber, SLOPE(is))
	    call tbeptr (tp, column[15], rownumber, ESLOPE(is))
	    call tbeptr (tp, column[16], rownumber, RESLOPE(is))
	    call tbeptr (tp, column[17], rownumber, A(is)**0.25)
	    call tbeptr (tp, column[18], rownumber, mag)
	    call tbeptr (tp, column[19], rownumber, emagl)
	    call tbeptr (tp, column[20], rownumber, emags)
	    call tbeptr (tp, column[21], rownumber, TFE(is))
	    call tbeptr (tp, column[22], rownumber, TFC(is))
	    call tbeptr (tp, column[23], rownumber, tmage)
	    call tbeptr (tp, column[24], rownumber, tmagc)
	    call tbepti (tp, column[25], rownumber, TFEA(is))
	    call tbepti (tp, column[26], rownumber, TFCA(is))
	    do i = 0, 6, 2 {
	        call tbeptr (tp, column[i+27], rownumber, Memr[HH(is)+i/2])
	        call tbeptr (tp, column[i+28], rownumber, Memr[EHH(is)+i/2])
	    }
	    call tbepti (tp, column[35], rownumber, NDATA(is))
	    call tbepti (tp, column[36], rownumber, NPOINT(is) - NDATA(is))
	    call tbepti (tp, column[37], rownumber, NITER(al))
	    call tbepti (tp, column[38], rownumber, STOP(al))
	    call tbeptr (tp, column[39], rownumber, ABIG(is))
	    call tbeptr (tp, column[40], rownumber, SAREA(is))

	    # Optional harmonics
	    if (NHARM(is) > 0) {
	        j = 0
	        do i = 1, NHARM(is), 2 {
	            call tbeptr (tp, Memi[extracols+j], rownumber, 
	                         AI(HARM(is),i))
	            j = j + 1
	            call tbeptr (tp, Memi[extracols+j], rownumber, 
	                         AI(EHARM(is),i))
	            j = j + 1
	            call tbeptr (tp, Memi[extracols+j], rownumber, 
	                         BI(HARM(is),i))
	            j = j + 1
	            call tbeptr (tp, Memi[extracols+j], rownumber, 
	                         BI(EHARM(is),i))
	            j = j + 1
	        }
	    }
	}

	# Translate back to section coord. system.
	call el_shrink (im, is, al)
end
