# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<imio.h>
include <tbset.h>
include	"ellipse.h"
include	"elcolnames.h"

# EL_DSAMPLE -- Dump intensity azimutal sample as binary table.

procedure el_dsample (file, is, im, al, bufx, bufy) 

char	file[ARB]			# output table name
pointer	is				# isophote pointer
pointer	im				# IMIO pointer
pointer	al				# algorithm pointer
pointer	bufx, bufy			# data buffer pointers

pointer	tp, col1, col2			# binary table pointers

pointer	str, k[1], index
int	i, nrows
real	angle

int	tbtopn(), tbpsta(), strlen()
real	el_s2p()

errchk	tbtopen, tbtclo, tbcdef

begin
	call malloc (str, SZ_LINE, TY_CHAR)

	if (strlen (file) > 0) {

	    # Create table.
            tp = tbtopn (file, NEW_FILE, 0)

	    if (AANGLE(al))
	        call tbcdef (tp,col1,"PA","degrees","%10.6f",TY_REAL, 1, 1)
	    else
	        call tbcdef (tp,col1,"PAe","degrees","%10.6f",TY_REAL, 1, 1)
	    call tbcdef (tp, col2, "INTENS", "",     "",       TY_REAL, 1, 1)
	    call tbtcre (tp)

	    # Add image name to header.
	    call el_imname (tp, im, is)

	    # Add ellipse geometry to header.
	    call el_unshrink (im, is, al)

	    call sprintf (Memc[str], SZ_LINE, "%g")
	        call pargr (A(is))
	    call tbhanp (tp, ES_CA, TY_REAL, Memc[str], i)
	    call sprintf (Memc[str], SZ_LINE, "%g")
	        call pargr (EPS(is))
	    call tbhanp (tp, ES_CEPS, TY_REAL, Memc[str], i)
	    call sprintf (Memc[str], SZ_LINE, "%g")
	        call pargr (TETA(is)/PI*180.+90.)
	    call tbhanp (tp, ES_CTETA, TY_REAL, Memc[str], i)
	    call sprintf (Memc[str], SZ_LINE, "%g")
	        call pargr (el_s2p(im, XC(is), 1))
	    call tbhanp (tp, ES_CX0, TY_REAL, Memc[str], i)
	    call sprintf (Memc[str], SZ_LINE, "%g")
	        call pargr (el_s2p(im, YC(is), 2))
	    call tbhanp (tp, ES_CY0, TY_REAL, Memc[str], i)

	    call el_shrink (im, is, al)

	    # Dump sample.
	    do i = 0, NPOINT(is) - 1 {
	        # Angle.
	        if (AANGLE(al)) {
	            angle = Memr[bufx+i] + TETA(is) - PI2
	            if (angle > DPI)
	                angle = angle - DPI
	            if (angle < 0.)
	                angle = angle + DPI
	        }else
	            angle = Memr[bufx+i]

	        angle = angle / PI * 180.
	        call tbeptr (tp, col1, i+1, angle)

	        # Intensity.
	        if (!IS_INDEF(Memr[bufy+i]))
	            call tbeptr (tp, col2, i+1, Memr[bufy+i] + MEAN(is))
	        else
	            call tbeptr (tp, col2, i+1, INDEFR)
	    }

	    # Sort output table in ascending angle (thanks to Phil Hodge !)
	    nrows = tbpsta (tp, TBL_NROWS)
	    call malloc (index, nrows, TY_INT)
	    do i = 1, nrows
	        Memi[index+i-1] = i
	    call tbtflu (tp)			   # flush the buffer
	    if (AANGLE(al))
	        call tbcfnd (tp, "PA", k, 1)	   # find column ES_CA
	    else
	        call tbcfnd (tp, "PAe", k, 1)	   # find column ES_CA
	    call tbtsrt (tp, 1, k, false, nrows, Memi[index]) # sort the index
	    call reorder (tp, nrows, Memi[index])  # reorder the rows
	    call mfree (index, TY_INT)

            call tbtclo (tp)
	}

	call mfree (str, TY_CHAR)
end
