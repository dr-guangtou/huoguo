# Cursor help file.

define	EL_CURHELP	"isophote$src/elcurhelp.key"
define	EL_CURPROMPT	"isophote fitting cursor commands"
define	IP_CURHELP	"isophote$src/ipcurhelp.key"
define	IP_CURPROMPT	"isophote display cursor commands"


# Graphics devices.

define	GDEVICES	"|stdgraph|imdwhite|imdred|imdgreen|imdblue|imdyellow"
define	GDEV_NAMES	"|stdgraph|white|red|green|blue|yellow"

define	G_DEV		1
define	WHITE		2
define	RED		3
define	GREEN		4
define	BLUE		5
define	YELLOW		6

define	I_DEV		WHITE,RED,GREEN,BLUE,YELLOW


# Cursors.

define	EL_GCUR		"gcommands"
define	EL_IMCUR	"icommands"


# Display command sub-strings.

define	SZ_COMMAND	3 * SZ_LINE
define	DCOMM1  	"display "
define	DCOMM2  	" 1 "
define	DCOMM3		" erase+ border+ fill+ xcenter=0.5 ycenter=0.5 \
			xsize=1. ysize=1. xmag=1. ymag=1. mode=h >& dev$null"
define	CCOMM1  	"contour "
define	CCOMM2  	" append- mode=h >& dev$null"


# Display zoom modes and factor.

define	DZ_CENTER	1	# re-center with same zoom factor.
define	DZ_ZOOM		2	# zoom by DZ_FACTOR.
define	DZ_RESET	3	# redisplay full image section.

define	DZ_FACTOR	0.5	# 50% increase in visualized pixel size.


# Structure for holding displayed image section info.

define	LEN_DSTRUCT	10

define	DDEV		Memi[$1]          # device number
define	D1		Memi[$1+$2]       # coordinates of displayed section,
define	D2		Memi[$1+$2+2]     # in physical reference system.
define	DOFF		Memr[$1+$2+4]     # display offset.
define	DSCALE		Memr[$1+$2+6]     # display scale factor.
#			empty space
