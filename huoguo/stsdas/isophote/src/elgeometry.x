include	<imio.h>
include "ellipse.h"
include "eldisplay.h"


# Routines to translate coordinates and shrink/unshrink geometry 
# in between display, graphics, section and physical coordinates 
# (sort of "poor man's WCS").
#
#
# Physical:   1           to   IM_PHYSLEN
# Section:    IM_VOFF+1   to   IM_VOFF+IM_VSTEP*IM_LEN
# Graphics:   D1-IM_VOFF  to   D2-IM_VOFF
# Display:    0.0         to   1.0





##########################################################################
#                                                                        #
#                   COORDINATE  TRANSLATION                              #
#                                                                        #
##########################################################################


# EL_S2P  --  Section to physical.

real procedure el_s2p (im, sec_coor, axis)

pointer	im			# IMIO pointer
real	sec_coor		# coordinate in section system
int	axis			# image axis

begin
	return (real(IM_VOFF(im,axis)) + real(IM_VSTEP(im, axis)) * sec_coor)
end




# EL_P2S  --  Physical to section.

real procedure el_p2s (im, phys_coor, axis)

pointer	im			# IMIO pointer
real	phys_coor		# coordinate in physical system
int	axis			# image axis

begin
	return ((phys_coor - real(IM_VOFF(im,axis))) / real(IM_VSTEP(im, axis)))
end




# EL_S2D  --  Section to display.

real procedure el_s2d (im, dp, sec_coor, axis)

pointer	im			# IMIO pointer
pointer	dp			# displayed section info
real	sec_coor		# coordinate in section system
int	axis			# image axis

real	el_s2p()

begin
	return ((el_s2p (im, sec_coor, axis) - 
	        real(DOFF(dp,axis))) / DSCALE(dp,axis))
end




# EL_G2S  --  Graphics to section.

real procedure el_g2s (im, dp, graph_coor, axis)

pointer	im			# IMIO pointer
pointer	dp			# displayed section info
real	graph_coor		# coordinate in plotted system
int	axis			# image axis

real	el_p2s()

begin
	return (el_p2s (im,
	graph_coor * IM_VSTEP(im,axis) + real(D1(dp,axis)),
	axis))
end




# EL_S2G  --  Section to graphics.

real procedure el_s2g (im, dp, sec_coor, axis)

pointer	im			# IMIO pointer
pointer	dp			# displayed section info
real	sec_coor		# coordinate in section system
int	axis			# image axis

real	el_s2p()

begin
	return ((el_s2p (im, sec_coor, axis) - real(D1(dp,axis)) + 
	        real (IM_VSTEP(im,axis))) / real (IM_VSTEP(im,axis)))
end






##########################################################################
#                                                                        #
#                       SHRINK / UNSHRINK                                #
#                                                                        #
##########################################################################

# These routines assume that section step is the same in both axes.



# EL_SHRINK  --  From physical space to section space.

procedure el_shrink (im, is, al)

pointer	im			# IMIO
pointer	is			# isophote structure
pointer	al			# algorithm structure

real	scale
real	el_p2s()

begin
	# PA has units and origin always translated,
	# independent of PHYSICAL flag.
	if (!IS_INDEF(TETA(is)))
	    TETA(is) = (TETA(is)/180.* PI) + PI2
	if (!IS_INDEF(ETETA(is)))
	    ETETA(is) = ETETA(is) * PI / 180.

	# Shrink geometry.
	if (PHYSICAL(is)) {
	    # Center coordinates.
            if ((!IS_INDEFR (XC(is))) && (!IS_INDEFR (YC(is)))) {
	        XC(is) = el_p2s (im, XC(is), 1)
	        YC(is) = el_p2s (im, YC(is), 2)
	    }
	    # SMA and related parameters.
	    scale = real(IM_VSTEP(im,1))
            if (!IS_INDEFR(A(is)))
	        A(is) = A(is) / scale
	    if (al != NULL) {
                if (!IS_INDEFR(MINA(al)))
	            MINA(al) = MINA(al) / scale
                if (!IS_INDEFR(MAXA(al)))
	            MAXA(al) = MAXA(al) / scale
                if (!IS_INDEFR(MAXRIT(al)))
	            MAXRIT(al) = MAXRIT(al) / scale
                if (LINEAR(al))
	            ASTEP(al) = ASTEP(al) / scale
	    }
	}
end




# EL_UNSHRINK  --  From section space to physical space.

procedure el_unshrink (im, is, al)

pointer	im			# IMIO
pointer	is			# isophote structure
pointer	al			# algorithm structure

real	scale
real	el_s2p()

begin
	# Unshrink geometry.
	if (PHYSICAL(is)) {
	    # Center coordinates.
            if ((!IS_INDEFR (XC(is))) && (!IS_INDEFR (YC(is)))) {
	        XC(is) = el_s2p (im, XC(is), 1)
	        YC(is) = el_s2p (im, YC(is), 2)
	    }
	    # SMA and related parameters.
	    scale = real(IM_VSTEP(im,1))
            if (!IS_INDEFR(A(is)))
	        A(is) = A(is) * scale
	    if (al != NULL) {
                if (!IS_INDEFR(MINA(al)))
	            MINA(al) = MINA(al) * scale
                if (!IS_INDEFR(MAXA(al)))
	            MAXA(al) = MAXA(al) * scale
                if (!IS_INDEFR(MAXRIT(al)))
	            MAXRIT(al) = MAXRIT(al) * scale
                if (LINEAR(al))
	            ASTEP(al) = ASTEP(al) * scale
	    }
	}

	# PA has units and origin always translated,
	# independent of PHYSICAL flag.
	if (!IS_INDEF(TETA(is)))
	    TETA(is) = (TETA(is) - PI2) / PI * 180.
	if (!IS_INDEF(ETETA(is)))
	    ETETA(is) = ETETA(is) / PI * 180.
end

