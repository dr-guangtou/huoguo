# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<error.h>
include	<gio.h>
include	<gset.h>
include	"ellipse.h"
include	"eldisplay.h"

# List of colon commands.
define	CMDS "|color|dispars"

define	CCOLOR		 1	# Set/list graphics color
define	CDISPARS	 2	# Display parameters


# IP_COLON  --  Processes cursor colon commands. 
#
# This is a stripped-down version of el_colon, to be
# used with the 'isoexam' display task.

procedure ip_colon (dev, cmdstr, dcmd)

int	dev			# image/graphics device
char	cmdstr[ARB]		# string after colon
char	dcmd[ARB]		# display command

pointer	cmd
int	ncmd, i

int	nscan(), strdic()

begin
	call malloc (cmd, SZ_LINE, TY_CHAR)

	# Use formated scan to parse the command string.
	# The first word is the command and it may be minimum match
	# abbreviated with the list of commands.

	call sscan (cmdstr)
	call gargwrd (Memc[cmd], SZ_LINE)
	ncmd = strdic (Memc[cmd], Memc[cmd], SZ_LINE, CMDS)

	switch (ncmd) {

	case CCOLOR:
	    call gargwrd (Memc[cmd], SZ_LINE)
	    if (nscan() == 1) {
		call printf ("color = %s\n")
	            call extnstr (GDEV_NAMES, dev, Memc[cmd])
		    call pargstr (Memc[cmd])
	    } else {
	        i = strdic (Memc[cmd], Memc[cmd], SZ_LINE, GDEV_NAMES)
	        if (i < 2)
	            call eprintf ("Non-recognizable color.\n")
	        else
	            dev = i
	    }

	case CDISPARS:
	    call gargstr (dcmd, SZ_LINE)

	default:
	    call eprintf ("Non-recognizable colon command.\n")
	}

	call mfree (cmd, TY_CHAR)
end
