# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imset.h>
include <imhdr.h>
include <imio.h>
include <ctype.h>
include <tbset.h>
include <error.h>

include "ellipse.h"
include "eldisplay.h"
include "elcolnames.h"

# ELLIPSE -- Fits elliptical isophotes to images. Main program handles
# initialization, higher level I/O and main loop which scans in semi-major 
# axis space. Details of isophote fitting and interactive cursor control are 
# handled in subroutines.
#
# Current version uses clcmdw routine to send image display commands to the CL.
#
#
# 20 Oct 89  -  I.Busko - Task created.
#  6 Nov 89  -          - Added new integration modes and output of
#                         intensity samples.
# 20 Dez 89  -          - In-memory pixel subraster, faster gradient
#                         computation.
# 20 Feb 90  -          - Optional harmonic computation.
# 26 Mar 90  -          - Clipping of extreme data points
# 29 Mar 90  -          - Linear/geometric growing radius, fixed
#                         ellipse parameters
#  6 Aug 90  -          - Added structure for algorithm control
#  1 Sep 90  -          - Interactive image/graphic output
#  6 Sep 90  -          - Sub-pixel integration
# 21 Sep 90  -          - Bad pixel list
# 01 Jul 94  - J-C Hsu  - restrict ellip0 to >= 1.e-6 (in elfit.x)
#    Nov 95  - IB       - Major upgrade to 2nd version
# 23 Jul 96  -          - Limited ellipse center variations.
# 13 Sep 96  -          - Read ellipses from input table and don't fit.
# 20 Oct 97  -          - Variable threshold in object locator.
# 08 Jan 98  -          - Ignore some checks when reading from table.
# 13 Jan 98  -          - Upper/lower sigma clip.
# 04 Mar 98  -          - New logic for stopping when reading from table.
# 24 Jul 98  -          - First ellipse fit ignores gradient error result.
#  3 Feb 04  - P Hodge  - Remove "real" from the definition of el_bclean.

procedure t_ellipse ()

char    inname[SZ_PATHNAME]             # input image name
char    outname[SZ_PATHNAME]            # output table name
char    ellname[SZ_PATHNAME]            # input table with ellipses
char    sample[SZ_PATHNAME]             # intensity sample files
char    intmode[SZ_LINE]                # integration mode
char    harmonics[SZ_LINE]              # optional harmonics
char    dqfile[SZ_PATHNAME]             # data quality file name / extension
real    mag0                            # magnitude of reference
real    refer                           # reference intensity
real    backgr                          # zero (DC) bias level
bool    hcenter, heps, hteta            # hold parameters fixed ?
char    graphics[SZ_FNAME]              # Graphics device for isophote ploting.
char    sgraphics[SZ_FNAME]             # Graphics device for sample ploting.
bool    interactive                     # Interactive?
bool    splot                           # Simultaneous plot ?
bool	list				# List at STDOUT ?
bool	samplot				# Plot azimuthal samples ?
bool	memory				# memory-intensive ?

pointer is                              # isophote structure pointer
pointer im                              # image file pointer
pointer sec                             # in-memory pixel structure
pointer dp                              # display control parameters
pointer al                              # algorithm control structure
pointer tp, column[NUM_COLS]            # SDAS table, column pointers
pointer extracols
pointer	index, nrows, k[1]
pointer tpin                            # input table with ellipses
pointer	pp				# pset pointer
pointer	dqm_line			# DQF/mask line
char    pfile[SZ_PATHNAME]
char    str[SZ_PATHNAME]
int     pfindex                         # suffix of sample table name
int     rownumber                       # SDAS table row number
int     i, j
int     isignal                         # interactive signal
int	minit				# minimum # iterations
int	dev				# graphics device number
int	inrow				# current row in input ellipse table
int	krow				# row counter for input ellipse table
int	nurows				# # of rows in input ellipse table
long    cpu, clock                      # time variables
real    x1, y1, eps1, teta1             # first fitted isophote at A0
real    nexta                           # next radius to be fitted
real	maxsma				# maximum semi-major axis lenght
real	olksig				# object locator threshold
real    aux1, aux2
bool	grow				# current state of SMA incr/decr
bool	go_inwards			# signals change in SMA state
bool	lexceed				# last isoph. already exceeded limit
bool    firstellipse			# signals first isophote
bool	centercoords			# x0,y0 input by pset or prompt ?
bool	recenter			# re-center x0,y0 by finding routine ?
bool	xylearn				# update pset with x0,y0 ?
bool	extellip			# external ellipse input ?

pointer clopset()
pointer	immap(), imgs2r(), imgl2i()
pointer	 el_opend()
int      strlen(), strdic()
int     el_intfit(), nscan()
int	checkdim(), clgeti(), tbpsta()
long    cputime(), clktime()
real    clgetr()
bool	clgetb(), streq()

errchk  immap, malloc, el_intfit

string  maginfo	"Magnitudes will be computed with refer = 1. and zerolevel = 0.\n"

begin
        # Read main I/O parameters. 
        call clgstr ("input",  inname,  SZ_PATHNAME)
	if (strlen (inname) == 0)
	    call error (0, "No input image.")
        call clgstr ("output", outname, SZ_PATHNAME)
	if (strlen (outname) == 0)
	    call error (1, "No output table.")

        # Create structure for algorithm control.
        call malloc (al, LEN_ALGCONTROL, TY_STRUCT)

        # Create structure for isophote parameters.
        call malloc (is, LEN_ISTRUCT, TY_STRUCT)
        call malloc (HH(is),  MAX_HARM, TY_REAL)
        call malloc (EHH(is), MAX_HARM, TY_REAL)
	# Mark it as empty.
	VALID(is) = false

        # Read task parameters from main parameter file: 
	# some I/O, and interactive controls.
        call clgstr ("dqf", dqfile, SZ_PATHNAME)
	call el_bclean (dqfile, SZ_PATHNAME)
        call clgstr ("inellip", ellname, SZ_PATHNAME)
	call el_bclean (ellname, SZ_PATHNAME)
	if (strlen (ellname) == 0)
	    extellip = false
	else
	    extellip = true
        memory = clgetb ("memory")
        interactive = clgetb ("interactive")
        if (interactive)
	    splot = true
        else
            splot = false
        call clgstr ("device", graphics, SZ_FNAME)
	dev = strdic (graphics, graphics, SZ_LINE, GDEV_NAMES)
        PMASK(al)  = clgeti ("masksz")
        REGION(al) = clgetb ("region")
        list = clgetb ("verbose")

	# Read parameters from geometry pset.
	centercoords = false
        XC(is) = clgetr ("x0")
        YC(is) = clgetr ("y0")
	if (!IS_INDEF(XC(is)) && !IS_INDEF(YC(is)))
	    centercoords = true
        EPS(is)  = clgetr ("ellip0")
        TETA(is) = clgetr ("pa0")
	# Semi-major axis controls.
        A(is)      = clgetr ("sma0")
        MINA(al)   = clgetr ("minsma")
        MAXA(al)   = clgetr ("maxsma")
        ASTEP(al)  = clgetr ("step")
        LINEAR(al) = clgetb ("linear")
        MAXRIT(al) = clgetr ("maxrit")
        if (!LINEAR(al))
            ASTEP(al) = 1. + ASTEP(al)
        recenter = clgetb ("recenter")
        xylearn  = clgetb ("xylearn")
        PHYSICAL(is) = clgetb ("physical")

        # Read parameters from control pset.
        minit      = clgeti ("minit")
        NMAX(al)   = clgeti ("maxit")
        hcenter    = clgetb ("hcenter")
        heps       = clgetb ("hellip")
        hteta      = clgetb ("hpa")
        LSLOPE(al) = clgetr ("maxgerr")
        olksig     = clgetr ("olthresh")
        SOFT(al)   = clgetb ("soft")
        CPAR(al)   = clgetr ("conver")
        WANDER(al) = clgetr ("wander")

        # Read parameters from samples pset.
        call clgstr ("integrmode", intmode, SZ_LINE)
        USCLIP(al) = clgetr("usclip")
        LSCLIP(al) = clgetr("lsclip")
        NCLIP(al) = clgeti("nclip")
        FBAD(al)  = clgetr("fflag")
        call clgnone ("sdevice", sgraphics, SZ_FNAME)
        call clgnone ("tsample", sample, SZ_PATHNAME)
        call clgnone ("harmonics", harmonics, SZ_LINE)
	if (strlen (sgraphics) > 0)
	    samplot = true
	else
	    samplot = false
        AANGLE(al) = clgetb ("absangle")

	# Turn off sample ploting in case of graphics device conflict.
	i = strdic (sgraphics, str, SZ_PATHNAME, GDEVICES)
	if (i == dev) {
	    samplot = false
	    call eprintf ("WARNING: Sample plotting turned off.\n")
	}

	# Read parameters from magnitude pset.
        mag0   = clgetr ("mag0")
        refer  = clgetr ("refer")
        backgr = clgetr ("zerolevel")

        # Set fix flags
        FIXX(al) = false
        FIXY(al) = false
        FIXE(al) = false
        FIXT(al) = false
        if (hcenter) {
            FIXX(al) = true
            FIXY(al) = true
        }
        if (hteta)
            FIXT(al) = true
        if (heps)
            FIXE(al) = true

        # Map image.
        im = immap (inname, READ_ONLY, 0)
        call imseti (im, IM_TYBNDRY, BT_NEAREST) # boundary extension
        call imseti (im, IM_NBNDRYPIX, 10)
	# Change IM_NDIM(im) to checkdim(im), JC Hsu 11/29/94
        if (checkdim(im) != 2) {
            call imunmap (im)
            call error (3, "Input image is not two-dimensional")
        }
        call imsetr (im, IM_ADVICE, real(RANDOM))

	# Structure to hold pixel/dqf/mask info.
        call malloc (sec, LEN_PIXSTRUCT, TY_STRUCT)
	SUBRASTER(sec) = NULL

	# If memory-intensive, read pixel array. 
	if (memory)
	    PIXARRAY(sec)  = imgs2r (im, 1, IM_LEN(im,1), 1, IM_LEN(im,2))
	else
	    PIXARRAY(sec)  = NULL

	# Open DQF.
	call el_odqf (im, sec, dqfile)

        # Read dq file and flag pixels, only if memory-intensive.
        if ((PIXARRAY(sec) != NULL) && (DQF(sec) != NULL)) {
	    do j = 1, IM_LEN(DQF(sec),2) {
	        dqm_line = imgl2i (DQF(sec), j)
	        do i = 1, IM_LEN(DQF(sec),1) {
	            if (Memi[dqm_line+i-1] != 0)
	                call el_pflag (im, sec, i, j, false)
	        }
	    }
	}

	# Open pixel mask.
	call el_opm (im, sec, READ_ONLY)

        # Read pixel mask and flag pixels, only if memory-intensive.
        if ((PIXARRAY(sec) != NULL) && (MASK(sec) != NULL)) {
	    do j = 1, IM_LEN(MASK(sec),2) {
	        dqm_line = imgl2i (MASK(sec), j)
	        do i = 1, IM_LEN(MASK(sec),1) {
	            if (Memi[dqm_line+i-1] != 0)
	                call el_pflag (im, sec, i, j, false)
	        }
	    }
	}

	# Atempts to find object center.
	if (!extellip && olksig > 0.0)
	    call el_find (im, is, al, sec, recenter, list, olksig)

        # If still not successful, try once more to get center of initial 
	# isophote. If in interactive mode, prints warning message; user 
	# is supposed to define initial isophote center with cursor. If in 
	# non-interactive mode, and when learn-enabled, prompts at STDIN, 
	# so user is prompted even when task is run with mode=h. In this
	# case, updates pset too.
	if (!extellip) {
            if ((IS_INDEFR (XC(is))) || (IS_INDEFR (YC(is)))) {
                if (interactive) {
                    call printf ("\n\t\t Can not find object. Please, use\n")
                    call printf (
                          "\t\t `x' cursor key to mark galaxy center, and\n")
                    call printf ("\t\t next `f' key to fit.\n\n")
                } else {
	            if (xylearn) {
                        call printf (
                              "Can not find object. Please, enter initial\n")
                        call printf ("guess for starting isophote center:\n")
                        call printf ("X0: ")
	                call flush (STDOUT)
	                call scan ()
	                call gargr (XC(is))
                        call printf ("Y0: ")
	                call flush (STDOUT)
	                call scan ()
	                call gargr (YC(is))
	                pp = clopset ("geompar") 
	                call clppsetr (pp, "x0", XC(is))
	                call clppsetr (pp, "y0", YC(is))
	                call clcpset (pp)
	                centercoords = true
	            } else
	                call error (2, "Can not find object center.")
	        }
            }
        }

	# Initialize from input ellipse table.
	if (extellip) {
	    call el_rtable (im, is, al, ellname, tpin, inrow)
	    nurows = inrow
            krow = 1
	}

	# If not reading from table, translate geometry 
	# parameters to section system.
	if (!extellip) 
	    call el_shrink (im, is, al)

        # Check against some obvious inconsistencies.
        if (refer <= backgr) {
            call eprintf ("WARNING: reference <= zero (bias) level. \n")
            call eprintf (maginfo) 
            refer = 1.
            backgr  = 0.
        }
        if ((!IS_INDEFR (XC(is))) && (!IS_INDEFR (YC(is)))) {
            if ((XC(is) < 1.) || (XC(is) > IM_LEN(im,1)) ||
                (YC(is) < 1.) || (YC(is) > IM_LEN(im,2))) {
                call imunmap (im)
                call error (4, "Ellipse center is outside image.")
            }
        }
	if (IM_VSTEP(im,1) != IM_VSTEP(im,2))
            call error (9, "Different step in each axis is not allowed.")
	if (!extellip) {
	    if (IS_INDEF(MAXA(al)))
	        maxsma = max (IM_LEN(im, 1), IM_LEN(im, 2))
	    else
	        maxsma = MAXA(al)
            if (MINA(al) > maxsma)
                call error (5, "Incompatible minimum and maximum sma.")
            if ((A(is) < MINA(al)) || (A(is) > maxsma))
                call error (6, "Incompatible starting sma.")
            if ((ASTEP(al) < EPSILON) || (IS_INDEFR (ASTEP(al))))
                call error (7, "Error in sma step specification.")
	}

	# Test special case of minimum SMA = zero.
        if (MINA(al) < MIN_SMA) {
	    MINA(al) = MIN_SMA
	    ZEROA(al) = true
	} else
	    ZEROA(al) = false
        if (A(is) < MIN_SMA)
	    A(is) = MINA(al)

	# Open display and load image.
	if (interactive) {
	    dp = el_opend (im, dev)
	    call el_display (im, dp, is, sec, "", DZ_RESET, INDEFR, INDEFR)
	}

        # Get optional harmonic numbers.
        call sscan (harmonics)
        NHARM(is) = 0
        HARM_NUMBERS(is) = NULL
        call gargr (aux1)                       # this code is just to
        if (nscan() > 0) {                      # count the number of
            NHARM(is) = 1                       # optional harmonics.
            call gargr (aux2)
            while (aux1 != aux2) {
                NHARM(is) = NHARM(is) + 1
                aux1 = aux2
                call gargr (aux2)
            }
            NHARM(is) = 2 * NHARM(is)
            call malloc (HARM_NUMBERS(is), NHARM(is), TY_REAL)
            call malloc (HARM(is),         NHARM(is), TY_REAL)
            call malloc (EHARM(is),        NHARM(is), TY_REAL)
            call malloc (extracols, 2 * NHARM(is), TY_INT)
            call sscan (harmonics)              # here begins actual input.
            do i = 1, NHARM(is), 2 {
                call gargr (AI(HARM_NUMBERS(is),i))
                BI(HARM_NUMBERS(is),i) = AI(HARM_NUMBERS(is),i)
            }
        }

        # Open output table and initialize output lines.
        call el_opens (is, im, outname, tp, column, extracols, list)
        rownumber = 0

        # Initialize variables. 
	A0(is)      = A(is)
	RESLOPE(is) = 0.			# enable neighbor interp.
        if (streq (intmode, "bi-linear"))
            INTMODE(al) = INT_LINEAR
        else if (streq (intmode, "mean"))
            INTMODE(al) = INT_MEAN
        else if (streq (intmode, "median"))
            INTMODE(al) = INT_MED
        if (LINEAR(al))                      # radius of next isophote.
            nexta = A(is) + ASTEP(al)
        else
            nexta = A(is) * ASTEP(al)
        grow         = true                  # enable radius growing.
        firstellipse = true
	lexceed      = true
        call strcpy ("", pfile, SZ_PATHNAME)
        pfindex = 1
        cpu     = cputime (0)
        clock   = clktime (0)
        isignal = GO_OK
	# Enable iterative mode, or disable it 
	# if external ellipse input exists.
	if (!extellip) 
            STOP(al) = ST_OK
	else
	    STOP(al) = ST_NONITERATE

        # Main loop: fit isophote and bump SMA.

        while ((A(is) >= MINA(al)) && (isignal == GO_OK)) {

	    # If existent, read ellipse data from input table.
	    if (extellip) {
                if (krow > nurows)
                    break
	        call el_inrow (im, is, al, tpin, inrow)
                krow = krow + 1
	        if (A(is) == 0.0) 
	            break
            }

            # Reset stop flag, if iterations not disabled.
            if (STOP(al) != ST_NONITERATE)
                STOP(al) = ST_OK

            # Enforce bi-linear interpolation, in case next isophote
            # is closer than 2 pixels.
            if (abs (nexta - A(is)) < 2.)
                INTEGR(al) = INT_LINEAR
            else
                INTEGR(al) = INTMODE(al)

            # Creates table name for intensity sample. It's composed of
            # root 'sample' plus suffix given by pfindex.
            if (strlen (sample) > 0) {
                call sprintf (pfile, SZ_PATHNAME, "%s_%03d")
                    call pargstr (sample)
                    call pargi (pfindex)
                pfindex = pfindex + 1
            }

	    # First isophote is fitted with augmented minimum 
	    # number of iterations, to maximize chances of finding
	    # correct geometric parameters (in case first guess
	    # is too bad).
            if ((NMIN(al) != 1) && (firstellipse))
	        NMIN(al) = max (8, (minit * 2))
	    else 
	        NMIN(al) = minit

            # Fit isophote.
	    isignal = el_intfit (im, sec, dp, is, al, pfile, samplot,
	                         sgraphics, interactive, splot, list)

	    # If no meaningful solution can be found
	    # at initial guess, can't help but abort.
	    if (firstellipse  && !VALID(is)) 
                call error (8, "No meaningful solution.")

            # Store data from first isophote, to be used when
            # changing from growing to shrinking sma. Do not need
	    # to store full isophote data, just ellipse geometry.
	    # Also, if center coordinates were not input from pset or
	    # from prompting, store then on pset to be eventually 
	    # used next time task is executed (only if learn-enabled).
            if (firstellipse            && 
	        VALID(is)               && 
	        (isignal != GO_DISCARD) && 
	        (isignal != GO_QUIT))   {
                x1     = XC(is)
                y1     = YC(is)
                eps1   = EPS(is)
                teta1  = TETA(is)
	        A0(is) = A(is)
	        if (!centercoords && xylearn) {
	            call el_unshrink (im, is, al)
	            pp = clopset ("geompar")
	            call clppsetr (pp, "x0", XC(is))
	            call clppsetr (pp, "y0", YC(is))
	            call clcpset (pp)
	            call el_shrink (im, is, al)
	        }
            }

	    # Take appropriate course of action depending on
	    # interactive signal returned by fitting routine.
	    switch (isignal) {

	    case GO_OK, GO_INWARDS:
	        # Compute integrated fluxes and write results to output file.
	        if (VALID(is)) {
	            call el_integr (im, sec, is)
	            call el_dump (is, al, im, outname, tp, column, extracols,
	                          rownumber, mag0, refer, backgr)
	        }

            case GO_DISCARD:
                pfindex = pfindex - 1   # to overwrite last sample table.
	        isignal = GO_OK         # reset for next isophote.

            case GO_GIVEUP:
                pfindex = pfindex - 1   # to overwrite last sample table.
	        isignal = GO_INWARDS    # change grow sense.

	    case GO_QUIT:
	        if (VALID(is)) {
	            call el_integr (im, sec, is)
	            call el_dump (is, al, im, outname, tp, column, extracols,
	                          rownumber, mag0, refer, backgr)
	        }
	        break   # to avoid next semi-major axis bumping. 
	    }

            # Grow/shrink semi-major axis.
            if (grow) {
                if (LINEAR(al)) {
                    A(is) = A(is) + ASTEP(al)
                    nexta = A(is) + ASTEP(al)
                } else {
                    A(is) = A(is) * ASTEP(al)
                    nexta = A(is) * ASTEP(al)
                }
            } else {
                if (LINEAR(al)) {
                    A(is) = A(is) - ASTEP(al)
                    nexta = A(is) - ASTEP(al)
                } else {
                    A(is) = A(is) / ASTEP(al)
                    nexta = A(is) / ASTEP(al)
                }
            }

	    # Here begins logic to detect end of SMA growing.

	    go_inwards = false                  # default: no change

	    # "Hard" criteria, which apply no matter what happened
	    # whith gradient error and maximum SMA: either sampling 
	    # is being done mainly off-boundaries, or `i' cursor key 
	    # was pressed.
	    if ((STOP(al) == ST_INDEF)    ||
                (isignal  == GO_INWARDS))
	        go_inwards = true

	    # Enter here only if SMA still growing and go-inwards signal
            # is not set yet. But do not do anything if this is the first
            # ellipse. From the first ellipse the fit must always go
            # outwards no matter what.
	    if (grow && !go_inwards && !firstellipse) {
	        if (IS_INDEF(MAXA(al))) {
	            # Maximum SMA was not defined; gradient error
	            # triggers change of SMA growing state. 
	            if (((RESLOPE(is) > LSLOPE(al)) || 
	                IS_INDEF(RESLOPE(is)))) {
	                if (SOFT(al)) {
	                    if (lexceed)                      # soft stop.
	                        go_inwards = true
	                } else                                # hard stop.
	                    go_inwards = true
	            }
	        } else {
	            # Maximum SMA was defined; gradient error
	            # triggers entry to no-iterate mode; SMA
	            # can still grows up to MAXA(al).
	            if (((RESLOPE(is) > LSLOPE(al)) || 
	                 IS_INDEF(RESLOPE(is)))) {
	                if (SOFT(al)) {
	                    if (lexceed)                      # soft stop.
	                         STOP(al) = ST_NONITERATE
	                } else                                # hard stop.
	                    STOP(al) = ST_NONITERATE
	            }
	            # Tests now maximum SMA condition.
	            if ((A(is) >  MAXA(al)))
	                go_inwards = true
	        }
	    }

	    # Soft stop handling. This sets lexceed to signal that current
	    # isophote already exceeded gradient error criterion.
	    # Two consecutive isophotes have to exceed criterion in
	    # order to trigger change of SMA growing state.
	    if (((RESLOPE(is) > LSLOPE(al)) || IS_INDEF(RESLOPE(is))))
	        lexceed = true
	    else
	        lexceed = false

	    # Set SMA to go inwards.
	    if (go_inwards) {
                isignal = GO_OK
                if (grow) {
                    grow = false
	            if (!extellip)
                        STOP(al) = ST_OK
                    if (LINEAR(al)) {
                        A(is) = A0(is) - ASTEP(al)
                        nexta = A(is)  - ASTEP(al)
                    } else {
                        A(is) = A0(is) / ASTEP(al)
                        nexta = A(is)  / ASTEP(al)
                    }
                    XC(is)   = x1       # Recover starting isophote data.
                    YC(is)   = y1
                    EPS(is)  = eps1
                    TETA(is) = teta1
                }
            }

	    # First ellipse processed; lower flag.
            firstellipse = false
        }

	# Special case of zero SMA. Add one row to
	# output table with central brightness data.
        if ((isignal != GO_QUIT) && (isignal != GO_DISCARD) && (ZEROA(al))) {
            call el_zero (im, sec, is, al, list)
	    call el_dump (is, al, im, outname, tp, column, extracols,rownumber, 
	                  mag0, refer, backgr)
	}

	# Sort output table in ascending order of semi-major 
	# axis (thanks to Phil Hodge !).
	nrows = tbpsta (tp, TBL_NROWS)
	call malloc (index, nrows, TY_INT)
	do i = 1, nrows
	    Memi[index+i-1] = i
	call tbtflu (tp)				   # flush the buffer
	call tbcfnd (tp, ES_CA, k, 1)			   # find column ES_CA
	call tbtsrt (tp, 1, k, false, nrows, Memi[index])  # sort the index
	call reorder (tp, nrows, Memi[index])		   # reorder the rows
	call mfree (index, TY_INT)

        # Finished !
	if (list && VALID(is)) {
	    call printf ("%7.2f  CPU seconds.\n%7.2f  minutes elapsed.\n")
                call pargr (real (cputime (cpu)) / 1000.)
                call pargr (real (clktime (clock)) / 60.)
	}

        # Close files and reclaim space.
        call tbtclo (tp)
	if (DQF(sec) != NULL)
	    call imunmap (DQF(sec))
	if (MASK(sec) != NULL)
	    call imunmap (MASK(sec))
        call imunmap (im)
	if (extellip)
            call tbtclo (tpin)
	if (interactive)
	    call mfree (dp, TY_STRUCT)
	if (PIXARRAY(sec) != NULL)
            call mfree (SUBRASTER(sec), TY_REAL)
        call mfree (sec, TY_STRUCT)
        call mfree (HH(is),  TY_REAL)
        call mfree (EHH(is), TY_REAL)
        if (HARM_NUMBERS(is) != NULL) {
            call mfree (HARM_NUMBERS(is), TY_REAL)
            call mfree (HARM(is),         TY_REAL)
            call mfree (EHARM(is),        TY_REAL)
            call mfree (extracols, TY_INT)
        }
        call mfree (is,  TY_STRUCT)
        call mfree (al,  TY_STRUCT)
end




# EL_BCLEAN  --  Clean string from leading blanks.

procedure el_bclean (str, size)

char	str[ARB]
int	size

char	tempst[SZ_FNAME]
int	i, j

begin
	j = 1
	i = 1
	while (IS_WHITE(str[i]))
	    i = i + 1
	while ((str[i]!=EOS) && (i<=size)) {
	    tempst[j] = str[i]
	    i = i + 1
	    j = j + 1
	}
	tempst[j] = EOS
	call strcpy (tempst, str, size)
end
