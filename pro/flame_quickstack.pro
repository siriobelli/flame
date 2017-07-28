PRO flame_quickstack, fuel

;
; Takes the A and B positions from the fuel.diagnostics structure,
; stacks all the frames in each position, and outputs the difference A-B.
; If AB subtraction is not set, then simply stack all the frames
; Useful for a quick look at the data.
; If some frames are missing the star trace, then the X position is also considered
;

	flame_util_module_start, fuel, 'flame_quickstack'


	; identify the A and B and X positions
	;*************************************

	; get the diagnostics
	diagnostics = fuel.diagnostics

	; select all A frames
	w_A = where(diagnostics.offset_pos eq 'A', /null)

	; select all B frames
	w_B = where(diagnostics.offset_pos eq 'B', /null)

	; select all X frames
	w_X = where(diagnostics.offset_pos eq 'X', /null)



	; stack all A, B, and X frames
	;*************************************

	; read in first file
	frame0 = readfits(fuel.util.science.raw_files[0], header)

	; total stack
	cube_tot = dblarr( (size(frame0))[1], (size(frame0))[2], fuel.util.science.n_frames )
	for i=0, fuel.util.science.n_frames-1 do cube_tot[*,*,i] = readfits(fuel.util.science.raw_files[i])
	if fuel.util.science.n_frames eq 1 then stack_tot = cube_tot $
		else stack_tot = median(cube_tot, dimension=3)

	; stack all A frames
	if w_A ne !NULL then begin
		cube_A = dblarr( (size(frame0))[1], (size(frame0))[2], n_elements(w_A) )
		for i=0,n_elements(w_A)-1 do cube_A[*,*,i] = readfits(fuel.util.science.raw_files[w_A[i]])
		if n_elements(w_A) eq 1 then stack_A = cube_A $
			else stack_A = median(cube_A, dimension=3)
	endif

	; stack all B frames
	if w_B ne !NULL then begin
		cube_B = dblarr( (size(frame0))[1], (size(frame0))[2], n_elements(w_B) )
		for i=0,n_elements(w_B)-1 do cube_B[*,*,i] = readfits(fuel.util.science.raw_files[w_B[i]])
		if n_elements(w_B) eq 1 then stack_B = cube_B $
			else stack_B = median(cube_B, dimension=3)
	endif

	; stack all X frames
	if w_X ne !NULL then begin
		cube_X = dblarr( (size(frame0))[1], (size(frame0))[2], n_elements(w_X) )
		for i=0,n_elements(w_X)-1 do cube_X[*,*,i] = readfits(fuel.util.science.raw_files[w_X[i]])
		if n_elements(w_X) eq 1 then stack_x = cube_x $
			else stack_X = median(cube_X, dimension=3)
	endif



	; combine A, B, and X stacks
	;*************************************

	; if there is no A-B subtraction, then simply output the total stack
	if fuel.input.AB_subtraction eq 0 then begin

		writefits, fuel.util.intermediate_dir + 'quickstack.fits', stack_tot, header
		print, 'I wrote ', fuel.util.intermediate_dir + 'quickstack.fits'

	; if doing A-B subtraction, then find theright combination
	endif else begin

		if w_A NE !NULL and w_B ne !NULL then begin
			writefits, fuel.util.intermediate_dir + 'quickstack_A-B.fits', stack_A - stack_B, header
			print, 'I wrote ', fuel.util.intermediate_dir + 'quickstack_A-B.fits'
		endif

		if w_A NE !NULL and w_X ne !NULL then begin
			writefits, fuel.util.intermediate_dir + 'quickstack_A-X.fits', stack_A - stack_X, header
			print, 'I wrote ', fuel.util.intermediate_dir + 'quickstack_A-X.fits'
		endif

		if w_B NE !NULL and w_X ne !NULL then begin
			writefits, fuel.util.intermediate_dir + 'quickstack_B-X.fits', stack_B - stack_X, header
			print, 'I wrote ', fuel.util.intermediate_dir + 'quickstack_B-X.fits'
		endif

	endelse


  flame_util_module_end, fuel

END
