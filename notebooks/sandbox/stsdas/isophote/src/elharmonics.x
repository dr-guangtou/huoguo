# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	"ellipse.h"

define	FARRAY		Memd[array+$2*norder+$1]	# coeff. matrix
define	FRSIDE		Memd[rside+$1]			# right side
define	FBF		Memd[bf+$1]			# basis functions
define	FHN		Memr[hn+$1]			# harmonics

# EL_HARMONICS -- Decompose data vectors in harmonic components.
# Solution of a complete system whith cross terms is implemented,
# to allow complete freedom on the spacing of data points.
# The harmonic number(s) must be passed from the calling
# program in array coef. Data values flagged with INDEFR are
# discarded. The error parameter is returned with the value ERR 
# in the case of a singular system or if ndata < norder+1.
# The fitted harmonics are removed from the input bufy vector.

procedure el_harmonics (bufx, bufy, ndata, coef, cerror, norder, rms, err)

real	bufx[ARB]		#i:  input data vectors
real	bufy[ARB]		#i
int	ndata			#i:  size of data vectors
real	coef[ARB]		#io: coefficient vector and harmonic numb.
real	cerror[ARB]		#o:  coeff. error vector
int	norder			#i:  2 * number of harmonics
real	rms			#o:  rms of fit
int	err			#o:  error code

pointer	sp, array, rside	# coefficient matrix and right side matrix
pointer	bf			# basis functions
pointer hn			# harmonic numbers
int	i, j, k, n
real	det, yc
double	drms

errchk	smark, salloc, sfree

begin
	if (ndata < (norder + 1)) {
	    err = ERR
	    return
	}

	# Create scratch space
	call smark (sp)
	call salloc (array, norder * norder, TY_DOUBLE)
	call salloc (rside, norder, TY_DOUBLE)
	call salloc (bf,    norder, TY_DOUBLE)
	call salloc (hn,    norder, TY_REAL)
	n = norder - 1

	# Save harmonic numbers
	do j = 0, n
	    FHN(j) = coef[j+1]

	# Zero sums
	do j = 0, n {
	    FRSIDE(j) = 0.d0
	    do i = 0, n {
	        FARRAY(i,j) = 0.d0
	    }
	}

	# Scan data vectors, building cross sums
	do i = 1, ndata {
	    # Ignore INDEF data
	    if (!IS_INDEFR (bufy[i])) {

	        # Build basis functions: sin (harm * data), cos (harm * data)
	        do j = 0, n
	            FBF(j) = double (sin (FHN(j) * bufx[i] + 
					  PI2 * float(mod(j,2))))

	        # Add to upper triangle of coefficient matrix
	        do j = 0, n {
	            do k = 0, n
	                FARRAY(k,j) = FARRAY(k,j) + FBF(k) * FBF(j)
	        }

	        # Add to right side vector
	        do j = 0, n
	            FRSIDE(j) = FRSIDE(j) + FBF(j) * bufy[i]
	    }
	}

	# Build lower triangle of coefficient matrix
	do j = 0, n {
	    do i = 0, n {
	        if (i != j)
	            FARRAY(i,j) = FARRAY(j,i)
	    }
	}

	# Solve system	    
	if (ndata > 5) {
	    call el_matinv (FARRAY(0,0), norder, det)
	    if (abs (det) < EPSILON) {
	        err = ERR
	        return
	    }
	    call el_matmul (FARRAY(0,0), FRSIDE(0), norder)
	} else {
	    do i = 0, n {
	        FRSIDE(i)   = 0.d0
	        FARRAY(i,i) = 0.d0
	    }
	}

	# Obtain r.m.s. of fit
	drms = 0.d0
	k    = 0
	do i = 1, ndata {
	    yc = 0.
	    do j = 0, n
	        yc = yc + real (FRSIDE(j)) *
			  sin (FHN(j) * bufx[i] + PI2 * float(mod(j,2)))
	    if (!IS_INDEFR (bufy[i])) {
	        bufy[i] = bufy[i] - yc
	        drms = drms + double(bufy[i]) * double(bufy[i])
	        k = k + 1
	    } 
	}
	i = k - norder - 1
	if (i > 0) 
	    drms = drms / double (i)
	else
	    drms = 0.d0
	rms = real (sqrt (drms))

	# Get solution
	do j = 0, n {
	    coef[j+1]   = real (FRSIDE(j))
	    cerror[j+1] = real (sqrt (drms * FARRAY(j,j)))
	}

	call sfree (sp)
	err = OK
end
