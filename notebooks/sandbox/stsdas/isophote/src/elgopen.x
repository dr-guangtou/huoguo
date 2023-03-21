# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <gset.h>
include	"eldisplay.h"

# EL_GOPEN  --  Open appropriate graphics device.

pointer procedure el_gopen (device)

int	device

pointer	str, gp, gopen()

begin
	call malloc (str, SZ_LINE, TY_CHAR)
	call extnstr (GDEVICES, device, Memc[str])

	switch (device) {
	case G_DEV:
	    gp = gopen (Memc[str], APPEND, STDGRAPH)
	case I_DEV:
#	    gp = gopen (Memc[str], NEW_FILE+AW_DEFER, STDIMAGE)
            gp = gopen (Memc[str], APPEND, STDIMAGE)
	}

	call mfree (str, TY_CHAR)
	return (gp)
end
