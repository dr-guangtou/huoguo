# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	"ellipse.h"


########################################################################
#                                                                      #
#            Routines to handle isophote structure                     #
#                                                                      #
########################################################################


# EL_COPY -- Copy data from one isophote structure to another. 

procedure el_copy (is, iso)

pointer	is, iso
int	i

begin
	A(iso)       = A(is)
	XC(iso)      = XC(is)
	YC(iso)      = YC(is)
	EPS(iso)     = EPS(is)
	TETA(iso)    = TETA(is)
	EX(iso)      = EX(is)
	EY(iso)      = EY(is)
	EEPS(iso)    = EEPS(is)
	ETETA(iso)   = ETETA(is)
	NPOINT(iso)  = NPOINT(is)
	NDATA(iso)   = NDATA(is)
	MEAN(iso)    = MEAN(is)
	SIGMA(iso)   = SIGMA(is)
	SLOPE(iso)   = SLOPE(is)
	ESLOPE(iso)  = ESLOPE(is)
	RESLOPE(iso) = RESLOPE(is)
	SAREA(iso)   = SAREA(is)
	ABIG(iso)    = ABIG(is)
	VALID(iso)   = VALID(is)

	A3(HH(iso))  = A3(HH(is))     # This code depends somewhat on the
	B3(HH(iso))  = B3(HH(is))     # number of non-optional harmonics,
	A4(HH(iso))  = A4(HH(is))     # MAX_HARM, which is set currently
	B4(HH(iso))  = B4(HH(is))     # to 4. Definitions in ellipse.h also
	A3(EHH(iso)) = A3(EHH(is))    # have to comply to this number.
	B3(EHH(iso)) = B3(EHH(is))
	A4(EHH(iso)) = A4(EHH(is))
	B4(EHH(iso)) = B4(EHH(is))

	NHARM(iso)   = NHARM(is)
	if (NHARM(iso) > 0) {
            do i = 1, NHARM(iso), 2 {
                AI(HARM_NUMBERS(iso),i) = AI(HARM_NUMBERS(is),i)
                BI(HARM_NUMBERS(iso),i) = BI(HARM_NUMBERS(is),i)
	        AI(HARM(iso),i)  = AI(HARM(is),i)
	        BI(HARM(iso),i)  = BI(HARM(is),i)
	        AI(EHARM(iso),i) = AI(EHARM(is),i)
	        BI(EHARM(iso),i) = BI(EHARM(is),i)
           }
	}
end



# EL_ALLOC  --  Allocs structure for holding isophote parameters.

pointer procedure el_alloc (nharm)

int	nharm		# number of optional harmonics

pointer	is

begin
	call malloc (is, LEN_ISTRUCT, TY_STRUCT)
        call malloc (HH(is), MAX_HARM, TY_REAL)
        call malloc (EHH(is),MAX_HARM, TY_REAL)
        if (nharm > 0) {
            call malloc (HARM_NUMBERS(is), nharm, TY_REAL)
            call malloc (HARM(is),         nharm, TY_REAL)
            call malloc (EHARM(is),        nharm, TY_REAL)
	}
	return (is)
end



# EL_FREE  --  Frees structure that holds isophote parameters.

procedure el_free (is)

pointer	is		# structure pointer

begin
        call mfree (HH(is), TY_REAL)
        call mfree (EHH(is), TY_REAL)
        if (NHARM(is) > 0) {
            call mfree (HARM_NUMBERS(is), TY_REAL)
            call mfree (HARM(is),         TY_REAL)
            call mfree (EHARM(is),        TY_REAL)
	}
	call mfree (is, TY_STRUCT)
end
