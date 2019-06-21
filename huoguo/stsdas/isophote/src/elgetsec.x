# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imhdr.h>
include	"ellipse.h"

# EL_GETSEC -- Gets an image subraster. If in memory-intensive mode,
# subraster is taken from array pointed by PIXARRAY(sec). Otherwhise, 
# it is taken directly from the pixel file by IMIO. Out-of-bound 
# conditions must me checked prior to calling this procedure.


procedure el_getsec (im, sec, i1, i2, j1, j2)

pointer	im				# IMIO pointer
pointer sec				# pointer to in-memory pixel struct
int	i1, i2, j1, j2			# subraster corners, in image units

pointer	mask1, mask2
int	i, j, k
int	ras_size			# subraster size
long	offset, index

pointer	imgs2r(), imgs2i()

begin
	ras_size = (i2 - i1 + 1) * (j2 -j1 + 1)

	# Pixels taken form in-memory pixel array. Buffer
	# might be re-allocated if size is not sufficient.
	if (PIXARRAY(sec) != NULL) {
	    call realloc (SUBRASTER(sec), ras_size, TY_REAL)
	    k      = 0
	    index  = 0
	    offset = long ((j1 - 1) * (IM_LEN(im,1)))
	    offset = offset + long (i1 - 1)
	    do j = j1, j2 {
	        do i = i1, i2 {
	            Memr[SUBRASTER(sec)+k] = Memr[PIXARRAY(sec)+offset+index]
	            index = index + 1
	            k     = k + 1
	        }
	        index = index + long ((IM_LEN(im,1) - 1) - (i2-i1))
	    }

	# Pixels taken directly from file. IMIO manages the buffer.
	# Both DQF and pixel mask, if existent, must be queried.
	} else {
	    SUBRASTER(sec) = imgs2r (im, i1, i2, j1, j2)
	    if (DQF(sec) != NULL) {
	        mask1 = imgs2i (DQF(sec), i1, i2, j1, j2)
	        do i = 0, ras_size-1 {
	            if (Memi[mask1+i] != 0)
	                Memr[SUBRASTER(sec)+i] = INDEFR
	        }
	    }
	    if (MASK(sec) != NULL) {
	        mask2 = imgs2i (MASK(sec), i1, i2, j1, j2)
	        do i = 0, ras_size-1 {
	            if (Memi[mask2+i] != 0)
	                Memr[SUBRASTER(sec)+i] = INDEFR
	        }
	    }
	}
end
