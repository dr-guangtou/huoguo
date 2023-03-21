# Structure for holding isophote parameters.

define	LEN_ISTRUCT	31

define	A		Memr[$1+0]	# semi-major axis
define	A0		Memr[$1+1]	# semi-major axis at startup
define	MEAN		Memr[$1+2]	# mean intensity
define	SIGMA		Memr[$1+3]	# intensity scatter
define	EPS		Memr[$1+4]	# ellipticity
define	EEPS		Memr[$1+5]	# error
define	TETA		Memr[$1+6]	# position angle
define	ETETA		Memr[$1+7]	# error
define	XC		Memr[$1+8]	# coordinates of center
define	EX		Memr[$1+9]	# and errors
define	YC		Memr[$1+10]
define	EY		Memr[$1+11]
define	SLOPE		Memr[$1+12]	# local gradient
define	ESLOPE		Memr[$1+13]	# gradient error
define	RESLOPE		Memr[$1+14]	# gradient relative error
define	TFE		Memr[$1+15]	# total flux inside isophote
define	TFC		Memr[$1+16]	# total flux inside circle
define	SAREA		Memr[$1+17]	# average sector area in pixels
define	ABIG		Memr[$1+18]	# largest harmonic amplitude
define	TFEA		Memi[$1+19]	# no. of valid pixels inside ellipse
define	TFCA		Memi[$1+20]	# no. of valid pixels inside circle
define	NPOINT		Memi[$1+21]	# total no. of points in sample
define	NDATA		Memi[$1+22]	# no. of useful points in sample
define	HH		Memi[$1+23]	# arrays with higher order 
define	EHH		Memi[$1+24]	# harmonics and errors
define	NHARM		Memi[$1+25]	# no. of opt. harmonic amplitudes
define	HARM_NUMBERS	Memi[$1+26]	# array whith optional harmonic numbers
define	HARM		Memi[$1+27]	# arrays whith optional harmonic 
define	EHARM		Memi[$1+28]	# amplitudes and errors
define	PHYSICAL	Memb[$1+29]	# physical coordinate I/O ?
define	VALID		Memb[$1+30]	# valid data in structure ?

# Harmonic array elements.

define	A3		Memr[$1+0]	# non-optional harmonics
define	B3		Memr[$1+1]
define	A4		Memr[$1+2]
define	B4		Memr[$1+3]
define	MAX_HARM	4		# max. no. of non-optional harmonics
define	AI		Memr[$1+$2-1]	# optional harmonics
define	BI		Memr[$1+$2]

# Structure for holding algorithm control and state parameters.

define	LEN_ALGCONTROL	28

define	INTMODE		Memi[$1+0]	# integration mode chosen by user
define	INTEGR		Memi[$1+1]	# current integration mode
define	NMAX		Memi[$1+2]	# maximum no. of iterations
define	NMIN		Memi[$1+3]	# minimum no. of iterations
define	NITER		Memi[$1+4]	# actual no. of iterations
define	STOP		Memi[$1+5]	# stopping condition code
define	NCLIP		Memi[$1+6]	# number of sigma-clipping iterations
define	PMASK		Memi[$1+7]	# square pixel mask size
define	MINA		Memr[$1+8]	# minimum radius to fit
define	MAXA		Memr[$1+9]	# maximum radius to fit
define	MAXRIT		Memr[$1+10]	# maximum radius for iterative fit
define	ASTEP		Memr[$1+11]	# step for next isophote
define	LSLOPE		Memr[$1+12]	# max. acceptable slope relative error
define	CPAR		Memr[$1+13]	# convergency sensitivity
define	USCLIP		Memr[$1+14]	# upper sigma-clip criterion
define	LSCLIP		Memr[$1+15]	# lower sigma-clip criterion
define	FBAD		Memr[$1+16]	# acceptable fraction of bad points
define	WANDER		Memr[$1+17]	# maximum change in ellipse center
define	LINEAR		Memb[$1+18]     # linear/geometric radius growing
define	ZEROA		Memb[$1+19]	# zero radius enabled
define	FIXX		Memb[$1+20]	# flags for holding parameters fixed
define	FIXY		Memb[$1+21]
define	FIXE		Memb[$1+22]
define	FIXT		Memb[$1+23]
define	SOFT		Memb[$1+24]	# soft stop flag
define	AANGLE		Memb[$1+25]	# sample angles in image coord. system
define	REGION		Memb[$1+26]	# region masking mode
#			empty space

# Structure for holding pixel subraster and pixel mask info.

define	LEN_PIXSTRUCT	5

define	PIXARRAY	Memi[$1+0]	# pixel array (Mem) pointer
define	SUBRASTER	Memi[$1+1]	# subraster (Mem) pointer
define	DQF		Memi[$1+2]	# data quality file (IMIO) pointer
define	MASK		Memi[$1+3]	# masked pixel list (IMIO) pointer
#			empty space

# Stopping condition codes.
define	ST_CHNG		-1		# there are changed parameters
define	ST_OK		0		# converged
define	ST_INDEF	1		# > NPOINT*FBAD are INDEF
define	ST_MAXIT	2		# exceeded max. iterations
define	ST_SINGULAR	3		# singular matrix in harm. fit
define	ST_NONITERATE	4		# stop iterations

# Interactive signals.
define	GO_OK		0		# continue
define	GO_QUIT		1		# quit program
define	GO_INWARDS	2		# signals to change direction
define	GO_DISCARD	3		# discard current isophote
define	GO_GIVEUP	4		# discard, go inwards

# Integration and averaging modes.
define	IMODES		"|bi-linear|mean|median"
define	INT_LINEAR	1		# bi-linear interpolation
define	INT_MEAN	2		# mean over elliptical sector
define	INT_MED		3		# median over elliptical sector
define	INT_NEIGHBOR	4		# nearest neighbor

# General constants.
define	MAX_EPS		0.95		# maximum allowed ellipticity
define	MIN_SMA		0.5		# minimum allowed semi-major axis
define	PI		3.141592654
define	PI2		(PI / 2.)
define	DPI		(2 * PI)
define	EPSILON		1.e-10
define	SZ_BUFFER	300		# initial buffer size and buffer
					# grow step
