# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imhdr.h>
include	"ellipse.h"

# EL_GETPIX -- Gets a single pixel from image. If memory-intensive, gets
# pixel from array pointed by PIXARRAY(sec). Otherwhise, uses IMIO call.
# Out-of-bound conditions must me checked prior to calling this procedure.

procedure el_getpix (im, sec, ip, jp)

pointer	im				# IMIO pointer
pointer sec				# pointer to in-memory pixel struct
int	ip, jp				# pixel coordinates, in section units

long	offset

pointer	imgs2r(), imgs2i()

begin
	# Pixel taken form in-memory pixel array.
	if (PIXARRAY(sec) != NULL) {
	    call realloc (SUBRASTER(sec), 1, TY_REAL)
	    offset = long ((jp - 1) * (IM_LEN(im,1)))
	    offset = offset + long (ip - 1)
	    Memr[SUBRASTER(sec)] = Memr[PIXARRAY(sec)+offset]

	# Pixel taken directly from file. Both DQF and pixel mask,
	# if existent, must be queried.
	} else {
	    SUBRASTER(sec) = imgs2r(im, ip, ip, jp, jp)
	    if (DQF(sec) != NULL) {
	        if (Memi[imgs2i(DQF(sec), ip,ip, jp,jp)] != 0)
	            Memr[SUBRASTER(sec)] = INDEFR
	    }
	    if (MASK(sec) != NULL) {
	        if (Memi[imgs2i(MASK(sec), ip,ip, jp,jp)] != 0)
	            Memr[SUBRASTER(sec)] = INDEFR
	    }
	}
end
