# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 

include	"ellipse.h"

# EL_BILINEAR  --  Bi-linear interpolation. If semi-major axis is
# smaller than 20 pixels, use sub-pixel integration.

real procedure el_bilinear (im, sec, sma, ip, jp, fx, fy)

pointer	im			# IMIO pointer
pointer sec			# pointer to in-memory pixel struct
real	sma			# semi-major axis
int	ip, jp			# pixel address
real	fx, fy			# fractional pixel position on pixel array 

real	qx, qy
real	sample

real	el_subpix()

begin
	qx = 1. - fx
	qy = 1. - fy
	call el_getsec (im, sec, ip, ip+1, jp, jp+1)
	if (!IS_INDEFR (Memr[SUBRASTER(sec)  ])  && 
	    !IS_INDEFR (Memr[SUBRASTER(sec)+1])  &&
	    !IS_INDEFR (Memr[SUBRASTER(sec)+2])  && 
	    !IS_INDEFR (Memr[SUBRASTER(sec)+3])) {
	    if (sma > 20.) {
	        sample = Memr[SUBRASTER(sec)]   * qx * qy + 
			 Memr[SUBRASTER(sec)+1] * fx * qy +
	                 Memr[SUBRASTER(sec)+2] * qx * fy + 
	                 Memr[SUBRASTER(sec)+3] * fx * fy
	    } else {
	        sample = el_subpix (Memr[SUBRASTER(sec)],
	                            Memr[SUBRASTER(sec)+1],
	                            Memr[SUBRASTER(sec)+2],
	                            Memr[SUBRASTER(sec)+3],
	                            fx, fy)
	    }
	} else
	    sample = INDEFR

	return (sample)
end

