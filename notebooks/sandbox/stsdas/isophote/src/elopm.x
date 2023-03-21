# Copyright restrictions apply - see stsdas$copyright.stsdas 
# 
include	<imio.h>
include <imhdr.h>
include	"ellipse.h"

# EL_OPM  --  Opens pixel mask. This mask has the same name as the input 
# image, however with extension ".pl". Input image extension, if existent, 
# is assumed to be 3 characters long. `mode' controls what to return: if
# READ_ONLY, routine either maps a pre-existing pixel mask file for input
# only, or returns a null pointer. If READ_WRITE, routine either maps a 
# pre-existing file for both input and output, or creates an empty pixel 
# mask and maps it. When creating a new pixel mask, any section spec 
# associated with the input image is ignored. This will generate pixel masks 
# that are always 2-dim sections with no cluster/group attributes, with 
# physical size equal to the physical size of the input 2-dim full image
# plane.

procedure el_opm (im, sec, mode)

pointer	im				# input image IMIO pointer
pointer	sec				# pixel/dqf/mask structure pointer
int	mode				# operation mode

pointer	pmname
pointer	section
pointer	fullname
pointer	line
int	i, len

pointer	immap(), impl2i()
int	imaccess(), strlen()

begin
	# Assemble pixel mask name: get cluster name stripped off eventual
	# section spec, and look for a 3-character name extension. Replace
	# it by (or just append) ".pl" extension.
	call malloc (pmname, SZ_PATHNAME, TY_CHAR)
	call imgcluster (IM_NAME(im), Memc[pmname], SZ_PATHNAME)
	len = strlen(Memc[pmname])
	if (Memc[pmname+len-4] == '.')
	    Memc[pmname+len-4] = EOS
	call strcat (".pl", Memc[pmname], SZ_PATHNAME)

	# Assemble section spec, from input image section data.
	call malloc (section, SZ_PATHNAME, TY_CHAR)
	call sprintf (Memc[section], SZ_PATHNAME, "[%d:%d:%d,%d:%d:%d]")
	    call pargl (IM_VOFF(im,1)+IM_VSTEP(im,1))
	    call pargl (IM_VOFF(im,1)+IM_LEN(im,1)*IM_VSTEP(im,1))
	    call pargi (IM_VSTEP(im,1))
	    call pargl (IM_VOFF(im,2)+IM_VSTEP(im,2))
	    call pargi (IM_VOFF(im,2)+IM_LEN(im,2)*IM_VSTEP(im,2))
	    call pargi (IM_VSTEP(im,2))

	# Assemble full name (pixel mask + section).
	call malloc (fullname, SZ_PATHNAME, TY_CHAR)
	call strcpy (Memc[pmname], Memc[fullname], SZ_PATHNAME)
	call strcat (Memc[section], Memc[fullname], SZ_PATHNAME)

	switch (mode) {

	case READ_ONLY:
	    # Attempts to open pixel mask with appropriate section; if
	    # non-existent, return.
	    if (imaccess (Memc[fullname], READ_ONLY) == YES)
	        MASK(sec) = immap (Memc[fullname], READ_WRITE, 0)
	    else
	        MASK(sec) = NULL

	case READ_WRITE:
	    # Attempts to open pixel mask with appropriate section; if
	    # non-existent, create pixel mask with physical size equal
	    # to input image's, and map it according to input image section.
	    if (imaccess (Memc[fullname], READ_WRITE) == YES)
	        MASK(sec) = immap (Memc[fullname], READ_WRITE, 0)
	    else {
	        # Open brand new mask with physical size equal to input image.
	        MASK(sec) = immap (Memc[pmname], NEW_IMAGE, 0)
	        IM_LEN(MASK(sec),1) = IM_PHYSLEN(im,1)
	        IM_LEN(MASK(sec),2) = IM_PHYSLEN(im,2)
	        # Ensure mask is cleared.
	        do i = 1, IM_LEN(MASK(sec),2) {
	            line = impl2i (MASK(sec), i)
	            call amovki (0, Memi[line], IM_LEN(MASK(sec),1))
	        }
	        call imunmap (MASK(sec))
	        # Re-open mask, now with appropriate section spec.
	        MASK(sec) = immap (Memc[fullname], READ_WRITE, 0)
	    }
	}

	# Done.
	call mfree (fullname, TY_CHAR)
	call mfree (section,  TY_CHAR)
	call mfree (pmname,   TY_CHAR)
end

