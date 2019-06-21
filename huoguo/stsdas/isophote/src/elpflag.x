# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imhdr.h>
include	"ellipse.h"

# EL_PFLAG  --  Flag a single pixel.

procedure el_pflag (im, sec, i, j, interact)

pointer im				# IMIO pointer
pointer sec				# pointer to in-memory pixel struct
int	i, j				# pixel coords, in section units
bool	interact			# signals interactive request

pointer	pix, imps2i()

begin
	if (!interact) {
	    # If not requested by an interactive action, assumes that
	    # caller already tested presence of in-memory pixel array.
	    # This entry point is used at task initialization only, when
	    # first reading DQF / mask.
	    call el_pmem (im, sec, i, j, INDEFR)

	} else {
	    # Interactive request. First, flag in-memory pixel.
	    if (PIXARRAY(sec) != NULL)
	        call el_pmem (im, sec, i, j, INDEFR)
	    # Now, update pixel mask. If pixel mask does not exist,
	    # create it.
	    if (MASK(sec) == NULL)
	        call el_opm (im, sec, READ_WRITE)
	    pix = imps2i (MASK(sec), i,i, j,j)
	    Memi[pix] = 1
	}
end



# EL_PUNF  --  Un-flag a single pixel.

procedure el_punf (im, sec, i, j)

pointer im				# IMIO pointer
pointer sec				# pointer to in-memory pixel struct
int	i, j				# pixel coords, in section units

pointer	imgs2i(), imps2i()
pointer	imgs2r()

begin
	# Pixel flagged in DQF: ignore request.
	if (DQF(sec) != NULL) {
	    if (Memi[imgs2i (DQF(sec), i,i, j,j)] != 0)
	        return
	}

	# Un-flag in-memory pixel, reading it afresh from file.
	if (PIXARRAY(sec) != NULL)
	    call el_pmem (im, sec, i, j, Memr[imgs2r(im, i,i, j,j )])

	# Update pixel mask.
	if (MASK(sec) != NULL)
	    Memi[imps2i (MASK(sec), i,i, j,j)] = 0
end



# EL_PMEM  --  Set a pixel in memory array.

procedure el_pmem (im, sec, i, j, value)

pointer im, sec
int	i, j
real	value

long	offset

begin
	offset = long ((j - 1) * IM_LEN(im,1))
	offset = offset + i - 1
	Memr[PIXARRAY(sec)+offset] = value
end


