# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 

include	<tbset.h>
include	"ellipse.h"

# EL_RTABLE  --  Open input ellipse table and initialize isophote
#                structure with data taken from first and last line. 

procedure el_rtable (im, is, al, ellname, tpin, inrow)

pointer	im				# IMIO pointer
pointer	is				# isophote structure
pointer	al				# algorithm control structure
char	ellname[SZ_PATHNAME]		# table name
pointer	tpin				# table pointer (output)
int	inrow				# current row (output)

pointer	tbtopn()
int	tbpsta()

begin
	# Open input table.
	tpin = tbtopn (ellname, READ_ONLY, 0)

	# Set minimum semi-major axis.
	inrow = 1
	call el_inrow (im, is, al, tpin, inrow)
	MINA(al) = A(is)

	# Initializes with last row's data. This assumes
	# that table will be scanned backwards. Don't 
	# forget to bump inrow outwards since each call
	# to el_inrow bumps it inwards automatically.
	inrow = tbpsta (tpin, TBL_NROWS)
	call el_inrow (im, is, al, tpin, inrow)
	inrow = inrow + 1
end
