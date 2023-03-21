/* General constants */
#define PI 3.14159265358979323846
#define DPI (2. * PI)
#define PI2	(PI / 2.)

/* Constants related to the geometry of the isophote */
/* Minimum and maximum of ellipticity */
#define MAX_EPS	0.95
#define MIN_EPS	0.05
/* Minimum allowed semi-major axis length */
#define MIN_SMA 0.5

/* Stopping condition codes. */
#define ST_CHNG       -1      /* There are changed parameters */
#define ST_OK          0      /* Converged */
#define ST_INDEF       1      /* NPOINT * FBAD are INDEF */
#define ST_MAXIT       2      /* Exceeded max. iterations */
#define ST_SINGULAR    3      /* Singular matrix in harm. fit */
#define ST_NONITERATE  4      /* Stop iterations */

#define	EPSILON		1.e-10
/* Initial buffer size and buffer growth size */
#define	SZ_BUFFER	300