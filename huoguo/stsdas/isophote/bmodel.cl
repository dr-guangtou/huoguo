# BMODEL -- Builds a galaxy model image.
#
# Results from the isophotal analysis generated by task ELLIPSE
# are used to build a model image. The original SDAS table created 
# by ELLIPSE is sorted in place by TSORT, and is next interpolated 
# in a fine grid in the semi-major axis length by task TREBIN. The 
# resulting temporary table is feed as input to the hidden task
# MODEL, which creates the image.
#
#						Ivo Busko   11/89
#
#
#  11/95  Added "fulltable" and optional parent image name (IB)
# 26-May-2004  Phil Hodge  Remove code to save/restore parameters.


procedure bmodel (table, output)

char	table  = ""       {prompt="isophote table"}
char	output = ""       {prompt="output (model) image"}
char	parent = ""       {prompt="parent image"}
bool	fulltable = yes   {prompt="use full range of `SMA' from table ?"}
real	minsma = 1.       {min=0., prompt="minimum modelling SMA"}
real	maxsma = 1.       {min=0., prompt="maximum modelling SMA"}
real	backgr = 0.       {prompt="background value"}
char	interp = "spline" {prompt="interpolation algorithm",enum="nearest|linear|poly3|spline"}
bool	highar  = no      {prompt="add higher harmonics ?"}
bool	verbose = no      {prompt="print info ?"}

begin
	file	out, par, tab, tfile, tfile2
	real	min, max
	bool	verb, fta

	# Check for the presence of pre-requisite packages.
	if (!deftask("tsort")) {
	    print ("Package 'ttools' must be loaded first !\n")
	    bye
	}

	# Read task parameters.
	out = output
	par = parent
	tab = table
	fta = fulltable
	min = minsma
	max = maxsma
	verb = verbose
	if (!fta) {
	    if (((max - min) <= 0.) || (max <= 0.) || (min < 0.)) {
	        print ("Error in maxsma, minsma.\n")
	        bye
	    }
	}

	# Do it.
	if (verb)
	    print ("Creating rebinned temporary table...\n")
	tfile = mktemp ("tmp$treb")
	tcopy (tab, tfile, verbose-)
	tsort (tfile, "SMA")
	# Get minsma and maxsma from table
	if (fta) {
	    tinfo (tfile, ttout=no)
	    tabpar (tfile, "SMA", 1)
	    min = real (tabpar.value)
	    tabpar (tfile, "SMA", tinfo.nrows)
	    max = real (tabpar.value)
	}
	# First row is sma=0., copy parameters
	# from second row.
	tabpar (tfile, "SMA", 1)
	if (real(tabpar.value) <= 0.) {
	    tabpar (tfile, "ELLIP", 2)
	    partab (tabpar.value, tfile, "ELLIP", 1)
	    tabpar (tfile, "PA", 2)
	    partab (tabpar.value, tfile, "PA", 1)	    
	}
	tfile2 = mktemp ("tmp$treb")
	tcopy (tfile, tfile2, verbose-)
	tcalc (tfile2, "PA", "if PA<0. then PA+180. else PA", datatype="real")
	delete (tfile // ".tab", verify = no)
	trebin (tfile2, tfile, "SMA", min, max, 0.50, xtable="",
		function=interp, extrapolate=yes,
		value=INDEF, padvalue=INDEF, verbose=no)

	model.backgr  = backgr
	model.highar  = highar
	model.verbose = verb
	model (tfile, out, par)

	# Delete temporary files.
	delete (tfile // ".tab", verify = no)
	delete (tfile2 // ".tab", verify = no)
end
