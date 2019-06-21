# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 

# EL_CLIP --  Sigma-clipping. Input vector has its deviant points
# marked by setting then to INDEF. Sample statistics are updated.

procedure el_clip (data, usclip, lsclip, nclip, npoint, ndata, mean, sigma)

real	data[ARB]		# io: data vector to be clipped
real	usclip			# i:  upper sigma-clip criterion
real	lsclip			# i:  lower sigma-clip criterion
int	nclip			# i:  number of clip iterations
int	npoint			# i:  data vector dimension
int	ndata			# o:  number of valid data points
real	mean, sigma		# o:  mean and sigma after clipping 

int	i, iter
real	low, up
double	sum, sumsq

begin
	# Iterate.
	do iter = 1, nclip {

	    # Sums.
            sum   = 0.0D0
            sumsq = 0.0D0
            ndata = 0
            do i = 1, npoint {
	        if (!IS_INDEF(data[i])) {
	            sum   = sum   + double(data[i])
	            sumsq = sumsq + double(data[i]) ** 2
	            ndata = ndata + 1
	        }
	    }

	    # Mean and stddev.
            switch (ndata) {
            case 0:
                mean  = INDEFR
                sigma = INDEFR
            case 1:
                mean  = real(sum)
                sigma = INDEFR
            default:
                mean  = real((sum) / double(ndata))
                sigma = real((sumsq - mean * sum) / double(ndata - 1))
                if (sigma > 0.0)  
                    sigma = sqrt(sigma)
                else
                    sigma = 0.0
            }

	    # Flag deviant data points.
	    if (!IS_INDEF(mean) && !IS_INDEF(sigma)) {
	        low   = mean - lsclip * sigma
	        up    = mean + usclip * sigma
	        do i = 1, npoint {
	            if (!IS_INDEF(data[i])) {
	                if ((data[i] > up) || (data[i] < low)) 
	                    data[i] = INDEFR
	            }
	        }
	    }
	}

	# Final count of good points.
        ndata = 0
        do i = 1, npoint {
	    if (!IS_INDEF(data[i]))
	        ndata = ndata + 1
	}

	if (IS_INDEF(mean))
	    mean = 0.0
end
