# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
define	LOGPTR 	20			# log2(maxpts)  (1e6)

# EL_QSORTR -- Reorders the elements of the real X array.

procedure el_qsortr (x, nelem, el_comparer)

real	x[ARB]			# array to be sorted
int	nelem			# number of elements in array
int     el_comparer()		# function to be called to compare elements
extern	el_comparer 		# function to be called to compare elements
#--

int	i, j, k, lv[LOGPTR], p, uv[LOGPTR]
real	temp, pivot
define	swap {temp=$1;$1=$2;$2=temp}


begin
	lv[1] = 1
	uv[1] = nelem
	p = 1

	while (p > 0) {
	    if (lv[p] >= uv[p])			# only one elem in this subset
		p = p - 1			# pop stack
	    else {
		# Dummy loop to trigger the optimizer.
		do p = p, ARB {
		    i = lv[p] - 1
		    j = uv[p]

		    # Select as the pivot the element at the center of the
		    # subfile, to avoid quadratic behavior on an already
		    # sorted list.

		    k = (lv[p] + uv[p]) / 2
		    swap (x[j], x[k])
		    pivot = x[j]			# pivot line

		    while (i < j) {
			for (i=i+1;  el_comparer (x[i], pivot) < 0;  i=i+1)
			    ;
			for (j=j-1;  j > i;  j=j-1)
			    if (el_comparer (x[j], pivot) <= 0)
				break
			if (i < j)			# out of order pair
			    swap (x[i], x[j])	# interchange elements
		    }

		    j = uv[p]			# move pivot to position i
		    swap (x[i], x[j])		# interchange elements

		    if (i-lv[p] < uv[p] - i) {	# stack so shorter done first
			lv[p+1] = lv[p]
			uv[p+1] = i - 1
			lv[p] = i + 1
		    } else {
			lv[p+1] = i + 1
			uv[p+1] = uv[p]
			uv[p] = i - 1
		    }

		    break
		}

		p = p + 1			# push onto stack
	    }
	}
end


# EL_COMPARER -- Real comparison procedure for tqsort

int procedure el_comparer (a, b)

real	a, b

begin
	if (a > b)
	    return ( 1 )
	else if (a < b)
	    return ( -1)
	else
	    return ( 0 )
end

