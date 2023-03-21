# Column names, units and formats in isophote table.

define	NUM_COLS	40		# minimum number of 
					# columns in SDAS table
define	STRSIZE		10		# format string size

define	ES_FORM1	"%9.3g"		# often used formats
define	ES_FORM2	"%7.3g"
define	ES_EUNIT	""		# empty unit string

define	ES_CA		"SMA"		# semi-major axis length
define	ES_UA		"pixel"
define	ES_FA		"%7.2f"

define	ES_CINT		"INTENS"	# mean isophotal intensity
define	ES_FINT		"%10.3g"

define	ES_CEINT	"INT_ERR"	# isophotal intensity error
define	ES_FEINT	"%10.3g"

define	ES_CRMS		"RMS"		# scatter about mean intensity

define	ES_CPXV		"PIX_VAR"	# pixel variance

define	ES_CEPS		"ELLIP"		# ellipticity
define	ES_FEPS		"%6.4f"

define	ES_CEEPS	"ELLIP_ERR"	# ellipticity error
define	ES_FEEPS	"%6.4f"

define	ES_CTETA	"PA"		# position angle
define	ES_UTETA	"degrees"
define	ES_FTETA	"%6.2f"

define	ES_CETETA	"PA_ERR"	# error in position angle
define	ES_UETETA	"degrees"
define	ES_FETETA	"%6.2f"

define	ES_CX0		"X0"		# ellipse center
define	ES_UX0		"pixel"
define	ES_FX0		"%7.2f"

define	ES_CEX		"X0_ERR"	# center error
define	ES_UEX		"pixel"
define	ES_FEX		"%6.2f"

define	ES_CY0		"Y0"		# ellipse center
define	ES_UY0		"pixel"
define	ES_FY0		"%7.2f"

define	ES_CEY		"Y0_ERR"	# center error
define	ES_UEY		"pixel"
define	ES_FEY		"%6.2f"

define	ES_CSLOPE	"GRAD"		# local intensity gradient
define	ES_FSLOPE	"%8.3g"

define	ES_CESLOPE	"GRAD_ERR"	# error in gradient
define	ES_FESLOPE	"%6.3g"

define	ES_CRESLOPE	"GRAD_R_ERR"	# relative error in gradient
define	ES_FRESLOPE	"%6.3g"

define	ES_CR4		"RSMA"		# (semi-major axis length) ** 1/4
define	ES_UR4		"pixel**1/4"
define	ES_FR4		"%7.5f"

define	ES_CMAG		"MAG"		# isophotal magnitude

define	ES_CEML		"MAG_LERR"	# lower mag error bar

define	ES_CEMU		"MAG_UERR"	# upper mag error bar

define	ES_CTFE		"TFLUX_E"	# integrated flux in
define	ES_FTFE		"%12.5g"

define	ES_CTFC		"TFLUX_C"	# integrated flux in
define	ES_FTFC		"%12.5g"

define	ES_CTME		"TMAG_E"	# integrated magnitude in isoph.

define	ES_CTMC		"TMAG_C"	# integrated magnitude in circle

define	ES_CNPE		"NPIX_E"	# no. of valid pixels in
define	ES_FNPE		"%6d"		# ellipse integral

define	ES_CNPC		"NPIX_C"	# no. of valid pixels in
define	ES_FNPC		"%6d"		# circle integral

define	ES_CA3		"A3"		# harmonic amplitude

define	ES_CEA3		"A3_ERR"	# error in harmonic amplitude

define	ES_CB3		"B3"		# harmonic amplitude

define	ES_CEB3		"B3_ERR"	# error in harmonic amplitude

define	ES_CA4		"A4"		# harmonic amplitude

define	ES_CEA4		"A4_ERR"	# error in harmonic amplitude

define	ES_CB4		"B4"		# harmonic amplitude

define	ES_CEB4		"B4_ERR"	# error in harmonic amplitude

define	ES_CNDATA	"NDATA"		# no. of points
define	ES_FNDATA	"%5d"

define	ES_CNFLAG	"NFLAG"		# no. of flagged points
define	ES_FNFLAG	"%5d"

define	ES_CNITER	"NITER"		# no. of iterations
define	ES_FNITER	"%3d"

define	ES_CSTOP	"STOP"		# stopping condition code
define	ES_FSTOP	"%2d"

define	ES_CBIGA	"A_BIG"		# largest harm. amplitude

define	ES_CSAREA	"SAREA"		# sector area 
define	ES_USAREA	"pixel"
define	ES_FSAREA	"%5.1f"

define	ES_CAI		"AI%d"		# optional harm. amplitudes

define	ES_CEAI		"AI%d_ERR"	# error in optional harm. amplitudes

define	ES_CBI		"BI%d"		# optional harm. amplitudes

define	ES_CEBI		"BI%d_ERR"	# error in optional harm. amplitudes
