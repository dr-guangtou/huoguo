# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imio.h>
include	"ellipse.h"
include	"elcolnames.h"
	
# EL_OPENS -- Opens main output table.

procedure el_opens (is, im, outname, tp, column, extracols, list) 

pointer	is				# isophote structure pointer
pointer	im				# IMIO pointer
char	outname[ARB]			# output table name
pointer	tp, column[ARB]			# pointers for SDAS tables
pointer	extracols
bool	list				# list at STDOUT ?

char	str[STRSIZE]
char	form1[STRSIZE], form2[STRSIZE]
char	unit[STRSIZE]
int	i, j

int	strlen(), tbtopn()

errchk	tbtopen, tbtclo, tbcdef

begin
	call strcpy (ES_FORM1, form1, STRSIZE)
	call strcpy (ES_FORM2, form2, STRSIZE)
	call strcpy (unit, unit, STRSIZE)

	# Print header at STDOUT
	if (list)
	    call el_head ()

	if (strlen (outname) > 0) {
            tp = tbtopn (outname, NEW_FILE, 0)

	    # Standard columns
	    call tbcdef (tp, column[1], ES_CA, ES_UA, ES_FA, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[2], ES_CINT, unit, ES_FINT, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[3], ES_CEINT, unit, ES_FEINT, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[4], ES_CPXV, unit, form1, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[5], ES_CRMS, unit, form1, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[6], ES_CEPS, unit, ES_FEPS, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[7], ES_CEEPS, unit, ES_FEEPS, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[8], ES_CTETA, ES_UTETA, ES_FTETA, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[9], ES_CETETA, ES_UETETA, ES_FETETA,
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[10], ES_CX0, ES_UX0, ES_FX0, TY_REAL, 1, 1)
	    call tbcdef (tp, column[11], ES_CEX, ES_UEX, ES_FEX, TY_REAL, 1, 1)
	    call tbcdef (tp, column[12], ES_CY0, ES_UY0, ES_FY0, TY_REAL, 1, 1)
	    call tbcdef (tp, column[13], ES_CEY, ES_UEY, ES_FEY, TY_REAL, 1, 1)
	    call tbcdef (tp, column[14], ES_CSLOPE, unit, ES_FSLOPE, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[15], ES_CESLOPE, unit, ES_FESLOPE, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[16], ES_CRESLOPE, unit, ES_FRESLOPE, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[17], ES_CR4, ES_UR4, ES_FR4, TY_REAL, 1, 1)
	    call tbcdef (tp, column[18], ES_CMAG, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[19], ES_CEML, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[20], ES_CEMU, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[21], ES_CTFE, unit, ES_FTFE, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[22], ES_CTFC, unit, ES_FTFC, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[23], ES_CTME, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[24], ES_CTMC, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[25], ES_CNPE, unit, ES_FNPE, TY_INT, 1, 1)
	    call tbcdef (tp, column[26], ES_CNPC, unit, ES_FNPC, TY_INT, 1, 1)
	    call tbcdef (tp, column[27], ES_CA3, unit, form1, TY_REAL, 1, 1)
	    call tbcdef (tp, column[28], ES_CEA3, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[29], ES_CB3, unit, form1, TY_REAL, 1, 1)
	    call tbcdef (tp, column[30], ES_CEB3, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[31], ES_CA4, unit, form1, TY_REAL, 1, 1)
	    call tbcdef (tp, column[32], ES_CEA4, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[33], ES_CB4, unit, form1, TY_REAL, 1, 1)
	    call tbcdef (tp, column[34], ES_CEB4, unit, form2, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[35], ES_CNDATA, unit, ES_FNDATA, 
	                 TY_INT, 1, 1)
	    call tbcdef (tp, column[36], ES_CNFLAG, unit, ES_FNFLAG, 
	                 TY_INT, 1, 1)
	    call tbcdef (tp, column[37], ES_CNITER, unit, ES_FNITER, 
	                 TY_INT, 1, 1)
	    call tbcdef (tp, column[38], ES_CSTOP, unit, ES_FSTOP, 
	                 TY_INT, 1, 1)
	    call tbcdef (tp, column[39], ES_CBIGA, unit, form1, 
	                 TY_REAL, 1, 1)
	    call tbcdef (tp, column[40], ES_CSAREA, ES_USAREA, ES_FSAREA, 
	                 TY_REAL, 1, 1)

	    # Extra columns for optional harmonics
	    if (NHARM(is) > 0) {
	        j = 0
	        do i = 1, NHARM(is), 2 {
	            call sprintf (str, STRSIZE, ES_CAI)
	                call pargi (int(AI(HARM_NUMBERS(is),i)))
	            call tbcdef (tp, Memi[extracols+j], str, unit, form1, 
				 TY_REAL, 1, 1)
	            j = j + 1
	            call sprintf (str, STRSIZE, ES_CEAI)
	                call pargi (int(AI(HARM_NUMBERS(is),i)))
	            call tbcdef (tp, Memi[extracols+j], str, unit, form1, 
				 TY_REAL, 1, 1)
	            j = j + 1
	            call sprintf (str, STRSIZE, ES_CBI)
	                call pargi (int(BI(HARM_NUMBERS(is),i)))
	            call tbcdef (tp, Memi[extracols+j], str, unit, form1, 
				 TY_REAL, 1, 1)
	            j = j + 1
	            call sprintf (str, STRSIZE, ES_CEBI)
	                call pargi (int(BI(HARM_NUMBERS(is),i)))
	            call tbcdef (tp, Memi[extracols+j], str, unit, form1, 
				 TY_REAL, 1, 1)
	            j = j + 1
	        }
	    }

	    call tbtcre (tp)

	    # Add input image name to table header.
	    call el_imname (tp, im, is)
	}
end
