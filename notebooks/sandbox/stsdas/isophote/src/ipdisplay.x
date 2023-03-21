include <imio.h>
include <imhdr.h>
include <tbset.h>
include	"ellipse.h"
include	"eldisplay.h"

# IP_DISPLAY  --  Image display for 'isoexam' task.
#
# This is a variant of the general-purpose el_display routine.

procedure ip_display (im, dp, is, sec, tp, command, zoom, cx, cy)

pointer	im			# i: IMIO pointer
pointer	dp			# i: displayed section info
pointer	is			# i: isophote structure
pointer	sec			# i: pixel structure
pointer	tp			# i: table pointer
char	command[ARB]		# i: command line parameters
int	zoom			# i: zoom control
real	cx, cy			# i: cursor coordinates

pointer	gp			# GIO pointer
int	row			# current row (ellipse) in table

pointer	el_gopen()
int	tbpsta()

begin
	# Begin by displaying image.
	call el_dimage (im, dp, is, sec, command, zoom, cx, cy)

	# Open graphic overlay/display.
	gp = el_gopen (DDEV(dp))

	# Mark flagged pixels.
	call el_mark (im, dp, sec, gp)

	# Plot all ellipses in table, beginning by
	# last one (largest SMA).
	row = tbpsta (tp, TBL_NROWS)
	while (row > 0) {

	    # Read current table row, bump row number backwards.
	    call el_inrow (im, is, NULL, tp, row)

	    # Plot current ellipse.
	    if (A(is) > 0.0)
	        call el_plot (im, gp, dp, is)
	}

	# Done !
	call gflush (gp)
	call gclose (gp)
end

