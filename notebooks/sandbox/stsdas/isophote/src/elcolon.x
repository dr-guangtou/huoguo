# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<error.h>
include	<gio.h>
include	<gset.h>
include	"ellipse.h"
include	"eldisplay.h"

# List of colon commands.
define	CMDS "|show|maxit|minit|integrmode|minsma|maxsma|step|maxrit        \
	      |maxgerr|soft|cpar|usclip|lsclip|nclip|fflag|linear|hcenter \
	      |hellip|hpa|x0|y0|ellip|pa|masksz|region|color|dispars|wander"

define	CSHOW		1	# Show values of isophote parameters
define	CMAXIT		2	# Set/list max. no. of iterations
define	CMINIT		3	# Set/list min. no. of iterations
define	CINTEGR		4	# Set/list integration mode
define	CMINA		5	# Set/list minimum radius
define	CMAXA		6	# Set/list maximum radius
define	CASTEP		7	# Set/list radius step
define	CMAXRIT		8	# Set/list max. radius for iterative fit
define	CLSLOPE		9	# Set/list maximum acceptable gradient error
define	CSSTOP		10	# Set/list soft stop flag
define	CCPAR		11	# Set/list convergency sensitivity
define	CUSCLIP		12	# Set/list sigma-clipping criterion
define	CLSCLIP		13	# Set/list sigma-clipping criterion
define	CNCLIP		14	# Set/list number of clipping iterations
define	CFBAD		15	# Set/list bad data fraction
define	CLINEAR		16	# Set/list radius growing mode
define	CHCENTER	17	# Set/list fixed center flag
define	CHEPS		18	# Set/list fixed ellipticity flag
define	CHTETA		19	# Set/list fixed pos. angle flag
define	CCX0		20	# Set/list ellipse parameters
define	CCY0		21
define	CCEPS		22
define	CCTETA		23
define	CMASK		24	# Set/list pixel mask size
define	CREGION		25	# Set/list regionmask flag
define	CCOLOR		26	# Set/list graphics color
define	CDISPARS	27	# Display parameters
define	CWANDER		28	# Maximum center wandering


# EL_COLON -- Processes cursor colon commands. Most parameters are
# silently range-enforced on input; constants must match the ranges 
# defined in psets.

procedure el_colon (is, al, dev, cmdstr, dcmd)

pointer	is			# isophote pointer
pointer	al			# algorithm control structure
int	dev			# image/graphics device
char	cmdstr[ARB]		# string after colon
char	dcmd[ARB]		# display command

pointer	cmd
int	ival, ncmd, i
real	rval
bool	bval
pointer	gp1

int	nscan(), strdic()
pointer	el_gopen()

begin
	call malloc (cmd, SZ_LINE, TY_CHAR)

	# Use formated scan to parse the command string.
	# The first word is the command and it may be minimum match
	# abbreviated with the list of commands.

	call sscan (cmdstr)
	call gargwrd (Memc[cmd], SZ_LINE)
	ncmd = strdic (Memc[cmd], Memc[cmd], SZ_LINE, CMDS)

	switch (ncmd) {

	case CSHOW:
	    call gargwrd (Memc[cmd], SZ_LINE)
	    if (nscan() == 1) {
	        switch (dev) {
	        case I_DEV:
		    call el_show (is, al, "STDOUT")
	        case G_DEV:
	            gp1 = el_gopen (dev)
	            bval = (and (GP_GFLAGS(gp1), GF_WSACTIVE) != 0)
	            if (bval)
	                call gdeactivate (gp1, 0)
		    call el_show (is, al, "STDOUT")
	            if (bval)
	                call greactivate (gp1, AW_PAUSE)
	            call gflush (gp1)
	            call gclose (gp1)
	        }
	    } else {
		iferr (call el_show (is, al, Memc[cmd]))
		    call erract (EA_WARN)
	    }

	case CMAXIT:
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("maxit = %d\n")
		    call pargi (NMAX(al))
	    } else {
		NMAX(al) = max (ival, 2)
	    }

	case CMINIT:
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("minit = %d\n")
		    call pargi (NMIN(al))
	    } else {
		NMIN(al) = max (ival, 1)
	    }

	case CINTEGR:
	    call gargwrd (Memc[cmd], SZ_LINE)
	    if (nscan() == 1) {
		call printf ("integrmode = %s\n")
	            call extnstr (IMODES, INTMODE(al), Memc[cmd])
		    call pargstr (Memc[cmd])
	    } else {
	        i = strdic (Memc[cmd], Memc[cmd], SZ_LINE, IMODES)
	        if (i != 0)
	            INTMODE(al) = i
	        else 
	            call eprintf ("Non-recognizable integration mode.\n")
	    }

	case CMINA:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("minsma = %g\n")
		    call pargr (MINA(al))
	    } else
		MINA(al) = max (rval, 0.)

	case CMAXA:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("maxsma= %g\n")
		    call pargr (MAXA(al))
	    } else
		MAXA(al) = max (rval, 1.)

	case CASTEP:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("step = %g\n")
	            if (LINEAR(al))
		        call pargr (ASTEP(al))
	            else
		        call pargr (ASTEP(al)-1.)
	    } else {
	            if (LINEAR(al))
		        ASTEP(al) = max (rval, 1.E-1)
	            else
		        ASTEP(al) = max (rval+1., 1.01)
	    }

	case CMAXRIT:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("maxrit = %g\n")
		    call pargr (MAXRIT(al))
	    } else
		MAXRIT(al) = max (rval, 0.)

	case CLSLOPE:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("maxgerr = %g\n")
		    call pargr (LSLOPE(al))
	    } else
		LSLOPE(al) = max (rval, 0.)

	case CSSTOP:
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("soft = %b\n")
		    call pargb (SOFT(al))
	    } else
	        SOFT(al) = bval

	case CCPAR:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("cpar = %g\n")
		    call pargr (CPAR(al))
	    } else
		CPAR(al) = max (rval, 0.)

	case CUSCLIP:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("usclip = %g\n")
		    call pargr (USCLIP(al))
	    } else
		USCLIP(al) = max (rval, 0.)

	case CLSCLIP:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("lsclip = %g\n")
		    call pargr (LSCLIP(al))
	    } else
		LSCLIP(al) = max (rval, 0.)

	case CNCLIP:
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("nclip = %d\n")
		    call pargi (NCLIP(al))
	    } else {
		NCLIP(al) = max (ival, 0)
	    }

	case CFBAD:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("fflag = %g\n")
		    call pargr (FBAD(al))
	    } else
		FBAD(al) = min (max (rval, 0.), 1.)

	case CCX0:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("x0 = %g (%g)\n")
	            call pargr (XC(is))
		    call pargr (EX(is))
	    } else {
		XC(is) = max (rval, 1.)
	        STOP(al) = ST_CHNG
	    }

	case CCY0:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("y0 = %g (%g)\n")
	            call pargr (YC(is))
		    call pargr (EY(is))
	    } else {
		YC(is) = max (rval, 1.)
	        STOP(al) = ST_CHNG
	    }

	case CCEPS:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("ellip = %g (%g)\n")
		    call pargr (EPS(is))
		    call pargr (EEPS(is))
	    } else {
		EPS(is) = min (max (rval, 0.), 1.)
	        STOP(al) = ST_CHNG
	    }

	case CCTETA:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("pa = %g (%g)\n")
		    call pargr (TETA(is))
		    call pargr (ETETA(is))
	    } else {
	        rval = min (max (rval, -90.), 90.)
	        rval = (rval)
		TETA(is) = rval
	        STOP(al) = ST_CHNG
	    }

	case CWANDER:
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("wander = %g\n")
		    call pargr (WANDER(al))
	    } else
		WANDER(al) = max (rval, 0.)

	case CLINEAR:
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("linear = %b\n")
		    call pargb (LINEAR(al))
	    } else
	        LINEAR(al) = bval

	case CHCENTER:
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("hcenter = %b\n")
		    call pargb (FIXX(al))
	    } else {
	        FIXX(al) = bval
	        FIXY(al) = bval
	    }

	case CHEPS:
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("hellip = %b\n")
		    call pargb (FIXE(al))
	    } else
	        FIXE(al) = bval

	case CHTETA:
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("hpa = %b\n")
		    call pargb (FIXT(al))
	    } else
	        FIXT(al) = bval

	case CMASK:
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("masksz = %d\n")
		    call pargi (PMASK(al))
	    } else
		PMASK(al) = max (ival, 1)

	case CREGION:
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("region = %b\n")
		    call pargb (REGION(al))
	    } else
	        REGION(al) = bval

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
