# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	"ellipse.h"

# EL_POLAR -- Given x,y coordinates on image grid, returns radius 
# and polar angle on ellipse coordinate system. Takes care of 
# different definitions for theta and phi:
#
#		-PI < teta < PI
#		  0 < phi  < DPI

procedure el_polar (x, y, x0, y0, teta, radius, angle)

real	x, y, radius, angle
real	x0, y0, teta		# ellipse center and position angle

real x1, y1, teta1

begin	
	x1 = x - x0
	y1 = y - y0
	radius = x1**2 + y1**2
	if (radius > 0.0) {
	    radius = sqrt (radius)
	    angle  = asin (abs (y1) / radius)
	} else {
	    radius = 0.0
	    angle  = 1.0
	}
	if ((x1 >= 0.) && (y1 < 0.))
	    angle = DPI - angle
	else if ((x1 < 0.) && (y 1>= 0.))
	    angle = PI - angle
	else if ((x1 < 0.) && (y1 < 0.))
	    angle = PI + angle
	else {}
	teta1 = teta
	if (teta < 0.)
	    teta1 = teta + DPI
	angle = angle - teta1
	if (angle < 0. )
	    angle = angle + DPI
end
