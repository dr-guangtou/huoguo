# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	"ellipse.h"
include "eldisplay.h"

# IP_DOIT  --  Display image and ellipses with cursor interaction.

procedure ip_doit (im, sec, is, dp, tp)

pointer	im				# IMIO pointer
pointer	sec				# in-memory pixel structure
pointer	is				# isophote structure
pointer dp                              # display control parameters
pointer	tp				# table pointer

pointer	command				# cursor readback string
pointer	strdisp				# display command enhacements string
int	key				# cursor key
real	cx,cy				# cursor coordinates
pointer	gp

int	el_cursor()
pointer	el_gopen()

begin
	call calloc (command, SZ_LINE, TY_CHAR)
	call calloc (strdisp, SZ_LINE, TY_CHAR)
	Memc[command] = EOS
	Memc[strdisp] = EOS

	# Force display loading.
	key = 'o'

	# Main interactive loop.
	repeat {

	    switch (key) {

	    case '?':	# Cursor help.
	        switch (DDEV(dp)) {
	        case I_DEV:
	            call pagefile (IP_CURHELP, IP_CURPROMPT)
	        case G_DEV:
	            gp = el_gopen (DDEV(dp))
	            call gpagefile (gp, IP_CURHELP, IP_CURPROMPT)
	            call gflush (gp)
	            call gclose (gp)
	        }

	    case ':':	# Colon command.
	        call ip_colon (DDEV(dp), Memc[command], Memc[strdisp])

	    case 'o':                  	# Reset zoom.
	        call ip_display (im, dp, is, sec, tp, Memc[strdisp], DZ_RESET, 
                                 cx, cy)

	    case 'q':	                # Quit program.
	        call mfree (strdisp, TY_CHAR)
	        call mfree (command, TY_CHAR)
	        return

	    case 'r':                       # Roam, Re-display image.
	        call ip_display (im, dp, is, sec, tp, Memc[strdisp], 
	                         DZ_CENTER, cx, cy)

	    case 'z':                  	# Zoom display.
	        call ip_display (im, dp, is, sec, tp, Memc[strdisp], 
	                             DZ_ZOOM, cx, cy)

	    default:
	    }

	} until (el_cursor (im, dp, cx, cy, key, Memc[command], SZ_LINE)==EOF)

	call mfree (strdisp, TY_CHAR)
	call mfree (command, TY_CHAR)
end

