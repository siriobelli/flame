PRO flame_simple_stack, fuel=fuel

	; read in the filenames
	readcol, fuel.science_filelist, filenames, format='A'

	; identify the A and B and X positions
	diagnostics = *fuel.diagnostics

	; select all A frames
	w_A = where(diagnostics.offset_pos eq 'A', /null)
	
	; select all B frames
	w_B = where(diagnostics.offset_pos eq 'B', /null)

	; select all X frames
	w_X = where(diagnostics.offset_pos eq 'X', /null)



	; stack all A, B, and X frames
	;*************************************

	; read in first file
	frame0 = readfits(filenames[0], header)

	if w_A ne !NULL then begin
		cube_A = dblarr( (size(frame0))[1], (size(frame0))[2], n_elements(w_A) )
		for i=0,n_elements(w_A)-1 do cube_A[*,*,i] = readfits(filenames[w_A[i]])
		if n_elements(w_A) eq 1 then stack_A = cube_A $
		else stack_A = median(cube_A, dimension=3)
	endif

	if w_B ne !NULL then begin
		cube_B = dblarr( (size(frame0))[1], (size(frame0))[2], n_elements(w_B) )
		for i=0,n_elements(w_B)-1 do cube_B[*,*,i] = readfits(filenames[w_B[i]])
		if n_elements(w_B) eq 1 then stack_B = cube_B $
		else stack_B = median(cube_B, dimension=3)
	endif

	if w_X ne !NULL then begin
		cube_X = dblarr( (size(frame0))[1], (size(frame0))[2], n_elements(w_X) )
		for i=0,n_elements(w_X)-1 do cube_X[*,*,i] = readfits(filenames[w_X[i]])
		if n_elements(w_X) eq 1 then stack_x = cube_x $
		else stack_X = median(cube_X, dimension=3)
	endif



	; combine A, B, and X stacks
	;*************************************

	if w_A NE !NULL and w_B ne !NULL then begin
		writefits, fuel.intermediate_dir + 'quickstack_A-B.fits', stack_A - stack_B, header
		print, 'I wrote ', fuel.intermediate_dir + 'quickstack_A-B.fits'
	endif

	if w_A NE !NULL and w_X ne !NULL then begin
		writefits, fuel.intermediate_dir + 'quickstack_A-X.fits', stack_A - stack_X, header
		print, 'I wrote ', fuel.intermediate_dir + 'quickstack_A-X.fits'
	endif

	if w_B NE !NULL and w_X ne !NULL then begin
		writefits, fuel.intermediate_dir + 'quickstack_B-X.fits', stack_B - stack_X, header
		print, 'I wrote ', fuel.intermediate_dir + 'quickstack_B-X.fits'
	endif



END