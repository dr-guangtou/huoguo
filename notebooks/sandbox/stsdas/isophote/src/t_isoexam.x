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

# T_ISOEXAM  --  Read isophote table and plot contents on image display.
#
# This is basicaly a stripped-down version of the code in t_ellipse.x,
# retaining only the functionality associated with the interactive display
# functions. However, this code has to deal simultaneously with *all* ellipses
# in the input table, while t_ellipse.x is capable of handling only one
# ellipse at a time.
#
# Current version uses clcmdw routine to send image display commands to the CL.
#
#
# 16 Sep 96  -  I.Busko         - Task created.

procedure t_isoexam ()

char    inname[SZ_PATHNAME]             # input image name
char    ellname[SZ_PATHNAME]            # input table with ellipses
char    dqfile[SZ_PATHNAME]             # data quality file name / extension
char    graphics[SZ_FNAME]              # Graphics device for isophote ploting.
pointer is                              # isophote structure pointer
pointer sec                             # pixel array structure
pointer dp                              # display control parameters
pointer	im				# IMIO pointer
pointer	dqm_line			# DQF/mask line
pointer tp				# table pointer
int     i, j
int	dev				# graphics device number

pointer	immap(), el_opend()
pointer	tbtopn(), imgl2i(), imgs2r()
int      strlen(), strdic()

begin
        # Read main I/O parameters. 
	call clgstr ("table", ellname, SZ_PATHNAME)
	if (strlen (ellname) == 0)
	    call error (0, "No table.")
        call clgstr ("dqf", dqfile, SZ_PATHNAME)
	call el_bclean (dqfile, SZ_PATHNAME)
        call clgstr ("device", graphics, SZ_FNAME)
	dev = strdic (graphics, graphics, SZ_LINE, GDEV_NAMES)

	# Open isophote table.
	tp = tbtopn (ellname, READ_ONLY, 0)

	# Open image, getting name from table header.
	call tbhfkr (tp, "IMAGE", i, inname, j)
	im = immap (inname, READ_ONLY, 0)
        call imsetr (im, IM_ADVICE, real(RANDOM))

        # Create structure for isophote parameters.
        call malloc (is, LEN_ISTRUCT, TY_STRUCT)
	VALID(is) = false

	# Structure to hold pixel/dqf/mask info.
        call malloc (sec, LEN_PIXSTRUCT, TY_STRUCT)
	SUBRASTER(sec) = NULL

	# Read pixel array. If we use disk-intensive mode instead, display
	# updates become *too* slow, particularly with large sections.
	PIXARRAY(sec)  = imgs2r (im, 1, IM_LEN(im,1), 1, IM_LEN(im,2))

	# Open DQF.
	call el_odqf (im, sec, dqfile)

        # Read dq file and flag pixels.
        if (DQF(sec) != NULL) {
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

        # Read pixel mask and flag pixels.
        if (MASK(sec) != NULL) {
	    do j = 1, IM_LEN(MASK(sec),2) {
	        dqm_line = imgl2i (MASK(sec), j)
	        do i = 1, IM_LEN(MASK(sec),1) {
	            if (Memi[dqm_line+i-1] != 0)
	                call el_pflag (im, sec, i, j, false)
	        }
	    }
	}

	# Open display.
	dp = el_opend (im, dev)

	# Do it.
	call ip_doit (im, sec, is, dp, tp)

        # Close everything.
	call mfree (dp,  TY_STRUCT)
        call mfree (sec, TY_STRUCT)
        call mfree (is,  TY_STRUCT)
        call imunmap (im)
        call tbtclo (tp)
end


