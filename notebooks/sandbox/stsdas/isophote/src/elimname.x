# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include <imio.h>
include	"ellipse.h"

	
# EL_IMNAME  --  Add image name to table header. Section spec will
# depend on PHYSICAL flag. If set, specs for the appropriate axes
# will be replaced by '*', otherwise full spec will be used.

procedure el_imname (tp, im, is) 

pointer	tp				# table pointer
pointer	im				# IMIO pointer
pointer	is				# isophote structure pointer

pointer	name
pointer	section
pointer	axis_spec
int	ip, jp, axis_count, sec_count

int	strlen()

begin
	call malloc (name,      SZ_PATHNAME, TY_CHAR)
	call malloc (section,   SZ_PATHNAME, TY_CHAR)
	call malloc (axis_spec, SZ_FNAME,    TY_CHAR)

	# Isolate section spec from image name.
	call imgimage (IM_NAME(im), Memc[name], SZ_PATHNAME)
	call imgsection (IM_NAME(im), Memc[section], SZ_PATHNAME)

	# Enter here only if section spec was actually
	# appended to image name. Section string is scanned
	# and the two mapped axes specs are either replaced
	# by '*' or left as original.
	if (strlen (Memc[section]) > 0) {
	    sec_count = 1
	    ip = 1
	    call strcat ("[", Memc[name], SZ_PATHNAME)

	    # Process each axis in physical image.
	    do axis_count = 1, IM_NPHYSDIM(im) {

	        # Get section spec for current axis.
	        jp = 0
	        while ((Memc[section+ip] != ',') && 
	               (Memc[section+ip] != ']')) {
	            Memc[axis_spec+jp] = Memc[section+ip]
	            ip = ip + 1
	            jp = jp + 1
	        }
	        Memc[axis_spec+jp] = EOS
	        ip = ip + 1

	        # Append appropriate spec to name string.
	        if ((axis_count == IM_VMAP(im,sec_count)) &&
	            (sec_count <= 2)) {
	            if (PHYSICAL(is))
	                call strcat ("*,", Memc[name], SZ_PATHNAME)
	            else {
	                call strcat (Memc[axis_spec], Memc[name], SZ_PATHNAME)
	                call strcat (",", Memc[name], SZ_PATHNAME)
	            }
	            sec_count = sec_count + 1
	        } else {
	            call strcat (Memc[axis_spec], Memc[name], SZ_PATHNAME)
	            call strcat (",", Memc[name], SZ_PATHNAME)
	        }
	    }

	    # Close section spec.
	    ip = 0
	    while (Memc[name+ip] != EOS)
	        ip = ip + 1
	    Memc[name+ip]   = ']'
	    Memc[name+ip+1] = EOS
	}

	# Write table header keyword.
	call tbhadt (tp, "IMAGE", Memc[name])

	call mfree (axis_spec, TY_CHAR)
	call mfree (section,   TY_CHAR)
	call mfree (name,      TY_CHAR)
end

