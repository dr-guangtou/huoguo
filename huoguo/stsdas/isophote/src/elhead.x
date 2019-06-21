# Copyright restrictions apply - see stsdas$copyright.stsdas 
	
# EL_HEAD -- Print STDOUT header.

procedure el_head () 

begin
	call printf ("#\n")
	call printf ("# Semi-    Isophote      Ellipticity     ")
	call printf ("Position   Grad.  Data Flag Iter. Stop\n")
	call printf ("# major      mean                         ")
	call printf ("Angle      rel.                  code\n")
	call printf ("# axis     intensity                                ")
	call printf ("error\n")
	call printf ("#(pixel)                                 ")
	call printf ("(degree)\n#\n")
	call flush (STDOUT)
end
