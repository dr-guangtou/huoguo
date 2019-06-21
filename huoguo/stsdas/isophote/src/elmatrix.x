# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
define	FIK		Memi[ik+$1-1]
define	FJK		Memi[jk+$1-1]

# EL_MATINV -- Invert a symmetric matrix and calculate the determinant.
# The inverse is returned in place. This SPP version keeps the 
# non-structured style of the original FORTRAN II code.
#
# Original: Bevington, p. 302-303

procedure el_matinv (array, norder, det)

double	array[norder,ARB]	# io: double precision matrix
int	norder			# i:  degree of matrix
real	det			# o:  determinant

pointer	sp, ik, jk
int	i, j, k, l
double	amax, save

errchk	smark, salloc, sfree

begin
	call smark (sp)
	call salloc (ik, norder, TY_INT)
	call salloc (jk, norder, TY_INT)
	det = 1.

	do k = 1, norder {

	    # Find largest element in rest of matrix
	    amax = 0.d0
21
	    do i = k, norder {
	        do j = k, norder {
	            if (abs (amax) < abs (array[i,j])) {
	                amax = array[i,j]
	                FIK(k) = i
	                FJK(k) = j
	            }
	        }
	    }

	    # Interchange rows and columns to put amax in array(k,k)
	    if (abs (amax) < 1.d-10) {
	        det = 0.
	        call sfree (sp)
	        return
	    }
	    i = FIK(k)
	    if (i < k)
	        goto 21
	    if (i > k) {
	        do j = 1, norder {
	            save       = array[k,j]
	            array[k,j] = array[i,j]
	            array[i,j] = -save
	        }
	    }
	    j = FJK(k)
	    if (j < k)
	        goto 21
	    if (j > k) {
	        do i = 1, norder {
	            save       = array[i,k]
	            array[i,k] = array[i,j]
	            array[i,j] = -save
	        }
	    }

	    # Accumulate elements of inverse matrix
	    do i = 1, norder {
	        if (i != k)
	            array[i,k] = -array[i,k] / amax
	    }
	    do i = 1, norder {
	        do j = 1, norder {
	            if ((i != k) && (j != k)) 
	                array[i,j] = array[i,j] + array[i,k] * array[k,j]
	       }
	    }
	    do j = 1, norder {
	        if (j != k)
	            array[k,j] = array[k,j] / amax
	    }
	    array[k,k] = 1.d0 / amax
	    det = det * amax
	}

	# Restore ordering of matrix
	do l = 1, norder {
	    k = norder -l + 1
	    j = FIK(k)
	    if (j > k) {
	        do i = 1, norder {
	            save       =  array[i,k]
	            array[i,k] = -array[i,j]
	            array[i,j] =  save
	        }
	    } 
	    i = FJK(k)
	    if (i > k) {
	        do j = 1, norder {
	            save       =  array[k,j]
	            array[k,j] = -array[i,j]
	            array[i,j] =  save
	        }
	    }
	}
	call sfree (sp)
end


# EL_MATMUL -- Multiplies a square matrix by a columnar matrix.
# The result is returned in the columnar input matrix.

procedure el_matmul (array, column, norder)

double	array[norder,ARB]
double	column[ARB]
int	norder

pointer	sp, result
int	i, j

begin
	call smark (sp)
	call salloc (result, norder, TY_DOUBLE)

	do j = 0, norder - 1
	    Memd[result+j] = 0.d0

	do j = 1, norder {
	    do i = 1, norder {
	        Memd[result+j-1] = Memd[result+j-1] +  array[i,j] * column[i]
	    }
	}

	call amovd (Memd[result], column, norder)
	call sfree (sp)
end
