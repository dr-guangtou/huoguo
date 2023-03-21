# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
define	NCELL	8		# sqrt(number of cells) in target pixel

# EL_SUBPIX -- Integrates a bi-linear interpolation on a square 2 x 2 grid,
# using sub-pixel resolution.

real procedure el_subpix (z1, z2, z3, z4, fx, fy)

real	z1, z2, z3, z4		# pixel contents 
real	fx, fy			# fractional position of target pixel

real	za, zb, z, sum
real	x, y, a1, a2, a3
real	correction
int	i, j

begin	
	sum = 0.
	a1  = z2 - z1
	a2  = z4 - z3
	a3  = 1./ NCELL
	correction = 0.5 + a3 / 2.
	do j = 1, NCELL {
	    y = j * a3 + fy - correction
	    do i = 1, NCELL {
	        x = i * a3 + fx - correction
	        za = a1 * x + z1
	        zb = a2 * x + z3
	        z  = (zb - za) * y + za
	        sum = sum + z
	    }
	}
	sum = sum / NCELL**2
	return (sum)
end
