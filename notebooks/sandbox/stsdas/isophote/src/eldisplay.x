include <imio.h>
include <imhdr.h>
include	"ellipse.h"
include	"eldisplay.h"

# EL_DISPLAY  --  General-purpose image display.
#
# Image section spec is stripped off from image name, and internal
# section spec is appended to it. This is used to implement the 
# zoom/roam capability. 
#
# If flagged pixels exist, they are marked.

procedure el_display (im, dp, is, sec, command, zoom, cx, cy)

pointer	im			# i: IMIO pointer
pointer	dp			# i: displayed section info
pointer	is			# i: isophote structure
pointer	sec			# i: in-memory pixel structure
char	command[ARB]		# i: command line parameters
int	zoom			# i: zoom control
real	cx, cy			# i: cursor coordinates

pointer	gp			# GIO pointer

pointer	el_gopen()

begin
	# Begin by displaying image. This must be the first routine
	# to be called because it updates the display control structure
	# to reflect the current section beign loaded into display.
	call el_dimage (im, dp, is, sec, command, zoom, cx, cy)

	# Open graphic overlay/display.
	gp = el_gopen (DDEV(dp))

	# Now plot current ellipse.
	if (VALID(is))
	    call el_plot (im, gp, dp, is)

	# Mark flagged pixels.
	call el_mark (im, dp, sec, gp)

	# Done !
	call gflush (gp)
	call gclose (gp)
end





# EL_DIMAGE  --  Display the image itself.
#
# Non-standard clcmdw routine is used to access either `tv.display' 
# or `plot.contour' tasks.

procedure el_dimage (im, dp, is, sec, command, zoom, cx, cy)

pointer	im			# i: IMIO pointer
pointer	dp			# i: displayed section info
pointer	is			# i: isophote structure
pointer	sec			# i: in-memory pixel structure
char	command[ARB]		# i: command line parameters
int	zoom			# i: zoom control
real	cx, cy			# i: cursor coordinates

pointer	str
pointer	name
int	x1, x2, y1, y2
real	pcx, pcy		# cursor coordinates in physical units.
real	pl			# physical axis lenght of displayed section.

int	strlen()
real	el_s2p()

begin
	call malloc (str,  SZ_COMMAND, TY_CHAR)
	call malloc (name, SZ_FNAME,   TY_CHAR)

	# Process zoom/roam/reset.
	switch (zoom) {
	case DZ_RESET:
	    x1 = int(IM_VOFF(im,1)) + IM_VSTEP(im,1)
	    x2 = int(IM_VOFF(im,1)) + int(IM_LEN(im,1)) * IM_VSTEP(im,1)
	    y1 = int(IM_VOFF(im,2)) + IM_VSTEP(im,2)
	    y2 = int(IM_VOFF(im,2)) + int(IM_LEN(im,2)) * IM_VSTEP(im,2)

	case DZ_CENTER,DZ_ZOOM:
	    # center of displayed section.
	    pcx = el_s2p (im, cx, 1)
	    pcy = el_s2p (im, cy, 2)
	    # size of displayed section.
	    pl = real(max (D2(dp,1)-D1(dp,1), D2(dp,2)-D1(dp,2)))
	    if (zoom == DZ_ZOOM)
	        pl  = pl / (1. + DZ_FACTOR)
	    # everything must be inside input image section.
	    x1 = max (int(IM_VOFF(im,1)) + IM_VSTEP(im,1), 
	              int (pcx - pl / 2.))
	    y1 = max (int(IM_VOFF(im,2)) + IM_VSTEP(im,2), 
	              int (pcy - pl / 2.))
	    x2 = min (int(IM_VOFF(im,1)) + int(IM_LEN(im,1)) * IM_VSTEP(im,1), 
	              int (pcx + pl / 2.))
	    y2 = min (int(IM_VOFF(im,2)) + int(IM_LEN(im,2)) * IM_VSTEP(im,2), 
	              int (pcy + pl / 2.))
	}

	# Section beign displayed.
	D1(dp,1) = x1
	D2(dp,1) = x2
	D1(dp,2) = y1
	D2(dp,2) = y2
	call el_dscale (im, dp)

	# Build section spec.
	call imgimage (IM_NAME(im), Memc[name], SZ_FNAME)
	call sprintf (Memc[str], SZ_COMMAND, "%s[%d:%d:%d,%d:%d:%d]")
	    call pargstr (Memc[name])
	    call pargi (x1)
	    call pargi (x2)
	    call pargi (IM_VSTEP(im,1))
	    call pargi (y1)
	    call pargi (y2)
	    call pargi (IM_VSTEP(im,2))
	call strcpy (Memc[str], Memc[name], SZ_COMMAND)

	# Assemble command line. It includes task name, image name
	# with explicit section, image buffer (`display' task only),
	# user-supplied parameter sub-string, and default parameters.
	call strcpy  ("", Memc[str], SZ_COMMAND)
	switch (DDEV(dp)) {

	case G_DEV:
	    call strcat (CCOMM1, Memc[str], SZ_COMMAND)
	    call strcat (Memc[name], Memc[str], SZ_COMMAND)
	    if (strlen (command) > 0)
	        call strcat (command, Memc[str], SZ_COMMAND)
	    call strcat (CCOMM2, Memc[str], SZ_COMMAND)

	case I_DEV:
	    call strcat (DCOMM1, Memc[str], SZ_COMMAND)
	    call strcat (Memc[name], Memc[str], SZ_COMMAND)
	    call strcat (DCOMM2, Memc[str], SZ_COMMAND)
	    if (strlen (command) > 0)
	        call strcat (command, Memc[str], SZ_COMMAND)
	    call strcat (DCOMM3, Memc[str], SZ_COMMAND)
	}

	# Display.
	call clcmdw (Memc[str])

	call mfree (name, TY_CHAR)
	call mfree (str,  TY_CHAR)
end





# EL_MARK  --  Mark flagged pixels, only if DQF/mask activated.

procedure el_mark (im, dp, sec, gp)

pointer	im			# IMIO pointer
pointer	dp			# displayed section info
pointer	sec			# in-memory pixel structure
pointer	gp			# GIO pointer

int	x1, x2, y1, y2
int	i, j
real	el_p2s()

begin
	if ( (DQF(sec) != NULL) || (MASK(sec) != NULL)) {
	    y1 = int (el_p2s (im, real (D1(dp,2)), 2))
	    y2 = int (el_p2s (im, real (D2(dp,2)), 2))
	    x1 = int (el_p2s (im, real (D1(dp,1)), 1))
	    x2 = int (el_p2s (im, real (D2(dp,1)), 1))
	    call printf ("Loading display with pixel mask...\n")
	    call flush (STDOUT)
	    do j = y1, y2 {
	        do i = x1, x2 {
	            call el_getpix (im, sec, i, j)
	            if (IS_INDEF(Memr[SUBRASTER(sec)]))
	                call el_ppix (im, dp, gp, real(i), real(j))
	        }
	    }
	}
end





# EL_OPEND  --  Open image display structure.

pointer procedure el_opend (im, dev)

pointer	im			# i: IMIO pointer
int	dev			# i: device number

pointer	dp

begin
	# Create structure for display parameters.
	# This has to be freed by caller.
	call malloc (dp, LEN_DSTRUCT, TY_STRUCT)

	# Initialize display to full input section.
	D1(dp,1) = IM_VOFF(im,1) + IM_VSTEP(im,1)
	D1(dp,2) = IM_VOFF(im,2) + IM_VSTEP(im,2)
	D2(dp,1) = IM_VOFF(im,1) + IM_LEN(im,1) * IM_VSTEP(im,1)
	D2(dp,2) = IM_VOFF(im,2) + IM_LEN(im,2) * IM_VSTEP(im,2)

	# Compute scaling parameters.
	call el_dscale (im, dp)

	# Graphics overlay.
	DDEV(dp) = dev

	return (dp)
end




# EL_DSCALE  --  Compute image display scaling parameters.
#
# Scale and offset variables are used to plot graphics in image 
# device viewport, which is set to 0.0 - 1.0 range. Definitions 
# here are dependent on the way task `tv.display' implements its 
# `fill' mode. Also, they assume a square display frame buffer.

procedure el_dscale (im, dp)

pointer	im			# i: IMIO pointer
pointer	dp			# i: displayed section info

pointer	tty, ttygdes()
int	px, py, ttygeti()
real	fact


begin
	# Size in physical pixels of each displayed axis. 
	px = (D2(dp,1) - D1(dp,1) + 1) / IM_VSTEP(im,1)
	py = (D2(dp,2) - D1(dp,2) + 1) / IM_VSTEP(im,2)

	# Scale is defined by largest axis.
	fact = real (IM_VSTEP(im,1)) / real(IM_VSTEP(im,2))
	if (px > py) { 
	    DSCALE(dp,1) = real (D2(dp,1) - D1(dp,1) + 1)
	    DSCALE(dp,2) = real (D2(dp,1) - D1(dp,1) + 1) / fact
	} else if (py > px) {
	    DSCALE(dp,2) = real (D2(dp,2) - D1(dp,2) + 1)
	    DSCALE(dp,1) = real (D2(dp,2) - D1(dp,2) + 1) * fact
	} else {
	    DSCALE(dp,1) = real (D2(dp,1) - D1(dp,1) + 1)
	    DSCALE(dp,2) = DSCALE(dp,1)
	}

	# Offsets for display buffer totaly filled up by square section.
	DOFF(dp,1) = real(D1(dp,1)) - 0.5
	DOFF(dp,2) = real(D1(dp,2)) - 0.5

	# Get device resolution from graphcap entry. Correct 
	# offset by one device pixel when plotting.
	if (DDEV(dp) != G_DEV) {
            tty = ttygdes ("stdimage")
	    DOFF(dp,1) = DOFF(dp,1) - DSCALE(dp,1) * 
	                             (1. / real(ttygeti (tty,"xr")))
	    DOFF(dp,2) = DOFF(dp,2) - DSCALE(dp,2) * 
	                             (1. / real(ttygeti (tty,"yr")))
	    call ttycdes (tty)
	}

	# When displayed section is not square, smaller
	# axis' offset has to be modified to take into
	# account the dark band that results when using
	# fill mode in a square frame buffer.
	if (px < py) {
	    DOFF(dp,1) = DOFF(dp,1) - 
	                 (DSCALE(dp,1) - real (D2(dp,1)-D1(dp,1)+1)) / 2.
	} else if (py < px) { 
	    DOFF(dp,2) = DOFF(dp,2) - 
	                 (DSCALE(dp,2) - real (D2(dp,2)-D1(dp,2)+1)) / 2.
	}
end
