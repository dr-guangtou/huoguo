# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <tbset.h>

include	"ellipse.h"
include	"elcolnames.h"


# EL_INROW  --  Read one row with ellipse data from input table and
#               convert to section coordinate system.

procedure el_inrow (im, is, al, tp, row)

pointer	im				# IMIO pointer
pointer	is				# isophote structure
pointer	al				# algorithm control structure
pointer	tp, cptr			# input table pointers
int	row 				# current row number

begin
	call tbcfnd (tp, ES_CA, cptr, 1)
	call tbegtr (tp, cptr, row, A(is))
	call tbcfnd (tp, ES_CEPS, cptr, 1)
	call tbegtr (tp, cptr, row, EPS(is))
	call tbcfnd (tp, ES_CTETA, cptr, 1)
	call tbegtr (tp, cptr, row, TETA(is))
	call tbcfnd (tp, ES_CX0, cptr, 1)
	call tbegtr (tp, cptr, row, XC(is))
	call tbcfnd (tp, ES_CY0, cptr, 1)
	call tbegtr (tp, cptr, row, YC(is))

	call el_shrink (im, is, al)

	# Bump row counter backwards. This assumes that
	# initialization routine positions at last line.
	row = row - 1
end
