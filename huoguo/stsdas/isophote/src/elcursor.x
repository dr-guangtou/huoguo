# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include "eldisplay.h"

# EL_CURSOR  --  Handles details of image/graphics cursor input.
# Graphics cursor returns coordinates in plotted section units,
# image cursor returns then in physical units; both have to be
# translated to input image section system.
# Routine maps 'n' key as EOF, to signal next isophote.

int procedure el_cursor (im, dp, x, y, key, command, sz_comm)

pointer	im			# i: IMIO pointer
pointer	dp			# i: displayed section info
real	x,y			# o: cursor coordinates
int	key			# o: key pressed
char	command[sz_comm]	# o: colon command string
int	sz_comm			# i: size of command string

int	i, wcs, clgcur()
real	el_p2s(), el_g2s()

begin
	switch (DDEV(dp)) {

	case G_DEV:
	    i = clgcur (EL_GCUR, x, y, wcs, key, command, sz_comm)
	    x = el_g2s (im, dp, x, 1)
	    y = el_g2s (im, dp, y, 2)

	case I_DEV:
	    i = clgcur (EL_IMCUR, x, y, wcs, key, command, sz_comm)
	    x = el_p2s (im, x, 1)
	    y = el_p2s (im, y, 2)
	}

	if (key == 'n')
	    i = EOF

	return (i)
end
