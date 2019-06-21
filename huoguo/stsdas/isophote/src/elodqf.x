# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<imio.h>
include	"ellipse.h"

	
# EL_ODQF  --  Opens DQF.

procedure el_odqf (im, sec, dqf)

pointer	im				# input image IMIO pointer
pointer	sec				# pixel/dqf/mask structure pointer
char	dqf[ARB]			# dqf name or extension

pointer	imname				# input image name
pointer	dqfname				# DQF name
pointer	section				# input image section spec
pointer	str
int	ip, i
bool	flag

pointer	immap()
int	strlen(), imaccess()

begin
	call malloc (dqfname, SZ_PATHNAME, TY_CHAR)
	call malloc (section, SZ_PATHNAME, TY_CHAR)
	call malloc (str,     SZ_PATHNAME, TY_CHAR)

	# Put aside section spec, to be used later. 
	call imgsection (IM_NAME(im), Memc[section], SZ_PATHNAME)

	# Prepare input image name for parsing.
	call malloc (imname,  SZ_PATHNAME, TY_CHAR)
	call strcpy (IM_NAME(im), Memc[imname], SZ_PATHNAME)

	# To be a valid DQF extension, it must begin with a dot and be 3
	# characters long. Everything else is assumed to be a plain
	# file name.

	flag = false
	switch (strlen (dqf)) {

	# No DQF spec; do nothing.
	case 0:
	    DQF(sec) = NULL
	    call mfree (imname,  TY_CHAR)
	    call mfree (str,     TY_CHAR)
	    call mfree (section, TY_CHAR)
	    call mfree (dqfname, TY_CHAR)
	    return

	# Extension only.
	case 4:
	    if (dqf[1] == '.') {

	        # Parse input image name, looking for a 3-character extension.
	        ip = 1
	        Memc[dqfname] = Memc[imname]
	        flag = true
	        while (Memc[imname+ip] != EOS) {

	            # A dot was found, if an EOS or an '[' is 4 characters
	            # away, assume this is the input image extension. 
	            if (Memc[imname+ip] == '.')  {
	                if ((Memc[imname+ip+4] == '[') ||
	                    (Memc[imname+ip+4] == EOS)) {
	                    # Extension found ! Replace input image's
	                    # by the one supplied in dqf parameter.
	                    do i = 0, 3 {
	                        Memc[dqfname+ip+i] = dqf[i+1]
	                    }
	                    ip = ip + 4
	                    flag = false
	                } else {
	                    Memc[dqfname+ip] = '.'
	                    ip = ip + 1
	                }
	            } else {
	                Memc[dqfname+ip] = Memc[imname+ip]
	                ip = ip + 1
	            }
	        }
	        Memc[dqfname+ip] = EOS
	    } else {
	        # 4-character file name.
	        call strcpy (dqf, Memc[dqfname], 4)
	        flag = false
	    }

	# Full name already supplied.
	default:
	    call strcpy (dqf, Memc[dqfname], strlen(dqf))
	    flag = false

	}

	# Strip off any section spec eventually in dqf name.
	call imgimage (Memc[dqfname], Memc[str], strlen(Memc[dqfname]))

	# Append extension, if not done before.
	if (flag)
	    call strcat (dqf, Memc[str], SZ_PATHNAME)

	# Now, complete by appendig original image section spec.
	if (strlen(Memc[section]) > 0)
	    call strcat (Memc[section], Memc[str], SZ_PATHNAME)

	# Map DQF, if existant.
	if (imaccess (Memc[str], READ_ONLY) == YES)
	    DQF(sec) = immap (Memc[str], READ_ONLY, 0)
	else
	    DQF(sec) = NULL

	# Done.
	call mfree (imname,  TY_CHAR)
	call mfree (str,     TY_CHAR)
	call mfree (section, TY_CHAR)
	call mfree (dqfname, TY_CHAR)
end

