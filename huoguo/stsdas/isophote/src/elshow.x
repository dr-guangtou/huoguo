# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	"ellipse.h"

# EL_SHOW -- Show the isophote parameters. This routine assumes
#            that caller worries about translation to/from physical
#            coordinate system.

procedure el_show (is, al, file)

pointer	is			# isophote pointer
pointer	al			# algorithm pointer
char	file[ARB]		# output file

pointer	fd, sp, str

int	open()
long	clktime()

errchk	open()

begin
	fd = open (file, APPEND, TEXT_FILE)
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)

	call cnvtime (clktime(0), Memc[str], SZ_LINE)
	call fprintf (fd, "\n%s\n")
	    call pargstr (Memc[str])

	call fprintf (fd, "semi-major axis      = %g\n")
	    call pargr (A(is))
	call fprintf (fd, "mean intensity       = %g  (%g)\n")
	    call pargr (MEAN(is))
	    call pargr (SIGMA(is))
	call fprintf (fd, "ellipticity          = %g  (%g)\n")
	    call pargr (EPS(is))
	    call pargr (EEPS(is))
	call fprintf (fd, "position angle       = %g  (%g)\n")
	    call pargr (TETA(is))
	    call pargr (ETETA(is))
	call fprintf (fd, "X center             = %g  (%g)\n")
	    call pargr (XC(is))
	    call pargr (EX(is))
	call fprintf (fd, "Y center             = %g  (%g)\n")
	    call pargr (YC(is))
	    call pargr (EY(is))
	call fprintf (fd, "gradient             = %g  (%g)\n")
	    call pargr (SLOPE(is))
	    call pargr (ESLOPE(is))
	call fprintf (fd, "grad. relative error = %g\n")
	    call pargr (RESLOPE(is))
	call fprintf (fd, "Stopping code        = %d\n")
	    call pargi (STOP(al))
	call fprintf (fd, "Total no. of samples = %d,  useful samples = %d\n")
	    call pargi (NPOINT(is))
	    call pargi (NDATA(is))

	call sfree (sp)
	call close (fd)
end
