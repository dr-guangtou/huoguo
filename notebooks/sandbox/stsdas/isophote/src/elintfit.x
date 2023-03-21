# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	"ellipse.h"
include "eldisplay.h"

# EL_INTFIT -- Interactive fit of one elliptical isophote.
# Calls el_fit, which takes care of fitting algorithm details.
# Handles the image/graphics cursor processing. Routine el_cursor 
# maps `n' key as EOF, ending cursor loop. Returned value is 
# GO_-type interactive signal.

int procedure el_intfit (im, sec, dp, is, al, file, splot, sgraphics, 
	                 inter, plot, list)

pointer	im				# IMIO pointer
pointer	sec				# in-memory pixel structure
pointer dp                              # display control parameters
pointer	is				# isophote structure
pointer	al				# algorithm control structure
char	file[ARB]			# file whith intens. sample
bool	inter				# interactive ?
bool	splot				# plot azimuthal samples ?
char    sgraphics[ARB]                  # graphics device for sample ploting.
bool	plot				# plot ellipses ?
bool	list				# list at STDOUT ?

pointer	command				# cursor readback string
pointer	strdisp				# display command enhacements string
int	key				# cursor key
int	i
real	cx,cy, cxa,cya			# cursor coordinates
real	newa
bool	replot
pointer	gp

int	el_cursor()
pointer	el_gopen()

begin
	call calloc (command, SZ_LINE, TY_CHAR)
	call calloc (strdisp, SZ_LINE, TY_CHAR)
	Memc[command] = EOS
	Memc[strdisp] = EOS

	if (!inter && !plot) {
	    # Fit isophote without ploting and interaction.
	    call el_fit (im, sec, is, al, file, sgraphics, splot, list)

	} else if (!inter && plot) {
	    # Fit, plot and go.
	    call el_fit (im, sec, is, al, file, sgraphics, splot, list)
	    gp = el_gopen (DDEV(dp))
	    call el_plot (im, gp, dp, is)
	    call gflush (gp)
	    call gclose (gp)

	} else if (inter && plot) {
	    if (!IS_INDEF(XC(is)) && !IS_INDEF(YC(is)))
	        key = 'f'     # Center is defined, go directly to fit routine.
	    else 
	        key = ','     # Center not defined, read cursor first.

	    # Main interactive loop.
	    repeat {

	        switch (key) {

	        case '?':	# Cursor help.
	            switch (DDEV(dp)) {
	            case I_DEV:
	                call pagefile (EL_CURHELP, EL_CURPROMPT)
	            case G_DEV:
	                gp = el_gopen (DDEV(dp))
	                call gpagefile (gp, EL_CURHELP, EL_CURPROMPT)
#	                call gflush (gp)
	                call gclose (gp)
	            }

	        case ':':	# Access control and isophote parameters.
	            call el_unshrink (im, is, al)
	            call el_colon (is, al, DDEV(dp),
	                           Memc[command], Memc[strdisp])
	            call el_shrink (im, is, al)

	        case 'a':	# Redefine ellipse sma and angle.
	            if (!IS_INDEF(XC(is)) && !IS_INDEF(YC(is))) {
	                newa = (cx - XC(is))**2 + (cy - YC(is))**2
	                if (newa > 0.)
	                    newa = sqrt (newa)
	                if (newa < MINA(al))
	                    call eprintf ("WARNING: Smaller than minsma.\n")
	                if (!IS_INDEF(MAXA(al))) {
	                    if (newa > MAXA(al))
	                        call eprintf ("WARNING: Larger than maxsma.\n")
	                }
	                A(is) = newa
	                TETA(is) = asin ((cy - YC(is)) / newa)
	                if ((cx - XC(is)) < 0.)
	                    TETA(is) = PI - TETA(is)
	                if (TETA(is) > PI)
	                    TETA(is) = TETA(is) - PI
	                if (TETA(is) < 0.)
	                    TETA(is) = TETA(is) + PI
	                STOP(al) = ST_CHNG
	            }

	        case 'c':	     # Continue in non-interactive mode,
	                             # with disabled ploting.
	            if (!IS_INDEF(XC(is)) && !IS_INDEF(YC(is))) {
	                inter = false    
	                plot  = false
	                call mfree (strdisp, TY_CHAR)
	                call mfree (command, TY_CHAR)
	                return (GO_OK)
	            }

	        case 'd':	    # Discard.
	            if (VALID(is)) {
	                call printf ("Current isophote discarded.\n")
	                call mfree (strdisp, TY_CHAR)
	                call mfree (command, TY_CHAR)
	                return (GO_DISCARD)
	            }

	        case 'f':	# Fit current isophote.
	            if (!IS_INDEF(XC(is)) && !IS_INDEF(YC(is))) {
	                call el_fit (im, sec, is, al, file, sgraphics, 
	                             splot, list)
	                if (plot)
	                    replot = true
	                else
	                    replot = false
	            }

	        case 'g':	    # Give up: discard, go inwards.
	            if (VALID(is)) {
	                call printf ("Current isophote discarded.\n")
	                call mfree (strdisp, TY_CHAR)
	                call mfree (command, TY_CHAR)
	                return (GO_GIVEUP)
	            }

	        case 'h':	   # Continue in non-interactive mode,
	                           # but plot ellipses anyway.
	            if (!IS_INDEF(XC(is)) && !IS_INDEF(YC(is))) {
	                inter = false  
	                call mfree (strdisp, TY_CHAR)
	                call mfree (command, TY_CHAR)
	                return (GO_OK)
	            }

	        case 'i':	           	# Go inwards.
	            if (!IS_INDEF(XC(is)) && !IS_INDEF(YC(is))) {
	                call mfree (strdisp, TY_CHAR)
	                call mfree (command, TY_CHAR)
	                return (GO_INWARDS)
	            }

	        case 'm':	                # Mask pixels.
	            if (REGION(al)) {
	                cxa = cx
	                cya = cy
	                while ((cxa == cx) && (cya == cy)) {
                            call printf (" again:\n")
	                    call flush(STDOUT) 
	                    i = el_cursor (im, dp, cxa, cya, key, 
	                                   Memc[command], SZ_LINE)
	                }
	                call el_mask (im, sec, al, dp, cx, cy, cxa, cya)
	            } else {
	                if (PMASK(al) > 0) {
	                    call el_mask (im, sec, al, dp, cx, cy, cx, cy)
	                }
	            }
	            replot = false

	        case 'o':                  	# Reset zoom.
	            call el_display (im, dp, is, sec, Memc[strdisp], 
	                             DZ_RESET, cx, cy)

	        case 'p':	                # Plot current ellipse.
	            if (!IS_INDEF(XC(is)) && !IS_INDEF(YC(is)))
	                replot = true

	        case 'q':	                # Quit program.
	            call mfree (strdisp, TY_CHAR)
	            call mfree (command, TY_CHAR)
	            return (GO_QUIT)

	        case 'r':                       # Roam, Re-display image.
	            call el_display (im, dp, is, sec, Memc[strdisp], 
	                             DZ_CENTER, cx, cy)

	        case 'u':	                # Un-mask pixels.
	            if (REGION(al)) {
	                cxa = cx
	                cya = cy
	                while ((cxa == cx) && (cya == cy)) {
                            call printf (" again:\n")
	                    call flush(STDOUT) 
	                    i = el_cursor (im, dp, cxa, cya, key, 
	                                   Memc[command], SZ_LINE)
	                }
	                call el_unmask (im, sec, al, cx, cy, cxa, cya)
	            } else {
	                if (PMASK(al) > 0) {
	                    call el_unmask (im, sec, al, cx, cy, cx, cy)
	                }
	            }
	            replot = false

	        case 'x':	                # Define ellipse center.
	            XC(is) = cx
	            YC(is) = cy
	            STOP(al) = ST_CHNG

	        case 'z':                  	# Zoom display.
	            call el_display (im, dp, is, sec, Memc[strdisp], 
	                             DZ_ZOOM, cx, cy)

	        default:
	            replot = false
	        }

	        if (replot) {
	            gp = el_gopen (DDEV(dp))
	            call el_plot (im, gp, dp, is)
	            call gflush (gp)
	            call gclose (gp)
	            replot = false
	        }

	    } until (el_cursor (im, dp, cx, cy, key, Memc[command], SZ_LINE)
		     == EOF)
	}

	call mfree (strdisp, TY_CHAR)
	call mfree (command, TY_CHAR)

	return (GO_OK)
end

