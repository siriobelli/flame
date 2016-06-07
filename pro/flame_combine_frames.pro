;
;
; TO-DO:
;
; - frames are combined by shifting A-B-A by an integer number of pixels - can you do better?
; - Ideally, need to de-shift to the nearest integer during the rectification step
; - dithering length is calculated from first A frame and first B frame - we need a better definition
; - handle better the case in which A or B are actually sky
;



FUNCTION flame_combine_stack, filenames=filenames, sigma_clip=sigma_clip, rejected_im=rejected_im
;
; read in FITS files and mean-stack them doing a sigma clipping
; optionally, outputs an image with the number of masked pixels
;

	; number of frames 
	N_frames = n_elements(filenames)

	; read first frame
	im_0 = readfits(filenames[0], header)

	; make big cube containing all frames
	im_cube = dblarr( N_frames, (size(im_0))[1], (size(im_0))[2] )

	; read all frames 
	for i_frame=0,N_frames-1 do im_cube[i_frame,*,*] = readfits(filenames[i_frame])

	; calculate median at each pixel of the image
	median_im = median(im_cube, dimension=1)

	; calculate standard deviation at each pixel of the image
	sigma_im = stddev(im_cube, dimension=1, /nan)

	; make a cube with the corresponding median at each (pixel,frame) position
	median_cube = im_cube
	median_cube[*]=0.0
	for i_frame=0,N_frames-1 do median_cube[i_frame, *, *] = median_im

	; make a cube with the corresponding sigma at each (pixel,frame) position
	sigma_cube = im_cube
	sigma_cube[*]=0.0
	for i_frame=0,N_frames-1 do sigma_cube[i_frame, *, *] = sigma_im

	; for each pixel in the cube, calculate the deviation from the median in units of sigma
	deviation = (im_cube - median_cube)/sigma_cube 

	; and now make a mask that indicates those pixels that are more than sigma_clip off from the median
	mask_cube = fix(im_cube)
	mask_cube[*] = 0
	if keyword_set(sigma_clip) then $
		mask_cube[where( abs(deviation) GT sigma_clip, /null )] = 1

	; turn masked pixels into NaNs
	im_cube[where(mask_cube, /null)] = !values.d_nan

	; make an image with the number of rejected pixels
	rejected_im = total(mask_cube, 1)

	; finally stack the frames
	return,	mean(im_cube, dimension=1, /nan)

END


; -------------------------------------------------------------------------------------------


PRO flame_combine_frames_oneslit, slit=slit, fuel=fuel

	; prefix for output file names
	filename_prefix = fuel.output_dir + 'slit' + string(slit.number, format='(I02)') + '-' + slit.name

	; input filenames for this slit
	filenames = *slit.filenames

	; identify the A and B and X positions
	diagnostics = *fuel.diagnostics

	; select all A frames
	w_A = where(diagnostics.offset_pos eq 'A', /null)
	
	; select all B frames
	w_B = where(diagnostics.offset_pos eq 'B', /null)

	; select all X frames
	w_X = where(diagnostics.offset_pos eq 'X', /null)


	; combine all frames for sky spectrum
	;*************************************

	; rectified but non-sky-subtracted files
	sky_filenames = flame_util_replace_string(filenames, '.fits', '_rectified.fits')

	; read in the header with the wavelength calibration
	header = headfits(sky_filenames[0])

	; stack and get the sky spectrum
	stack_sky = flame_combine_stack(filenames=sky_filenames, sigma_clip=3.0)

	; write out the sky spectrum
	writefits, filename_prefix + '_stack_sky.fits', stack_sky, header


	; stack all A, B, and X frames
	;*************************************

	if w_A ne !NULL then begin
		stack_A_filenames = flame_util_replace_string(filenames[w_A], '.fits', '_skysub_rectified.fits')
		stack_A = flame_combine_stack(filenames=stack_A_filenames, sigma_clip=3.0, rejected_im=rejected_im_A)
		writefits, filename_prefix + '_stack_A.fits', stack_A, header
	endif

	if w_B ne !NULL then begin
		stack_B_filenames = flame_util_replace_string(filenames[w_B], '.fits', '_skysub_rectified.fits')
		stack_B = flame_combine_stack(filenames=stack_B_filenames, sigma_clip=3.0, rejected_im=rejected_im_B)
		writefits, filename_prefix + '_stack_B.fits', stack_B, header
	endif

	if w_X ne !NULL then begin
		stack_X_filenames = flame_util_replace_string(filenames[w_X], '.fits', '_skysub_rectified.fits')
		stack_X = flame_combine_stack(filenames=stack_X_filenames, sigma_clip=3.0, rejected_im=rejected_im_X)
		writefits, filename_prefix + '_stack_X.fits', stack_X, header
	endif


	; combine A, B, and X stacks
	;*************************************

	if w_A NE !NULL and w_B ne !NULL then begin
		writefits, filename_prefix + '_stack_A-B.fits', stack_A - stack_B, header
		writefits, filename_prefix + '_rejectedpixels_A-B.fits', rejected_im_A + rejected_im_B, header
	endif

	if w_A NE !NULL and w_X ne !NULL then begin
		writefits, filename_prefix + '_stack_A-X.fits', stack_A - stack_X, header
		writefits, filename_prefix + '_rejectedpixels_A-X.fits', rejected_im_A + rejected_im_X, header
	endif

	if w_B NE !NULL and w_X ne !NULL then begin
		writefits, filename_prefix + '_stack_B-X.fits', stack_B - stack_X, header
		writefits, filename_prefix + '_rejectedpixels_B-X.fits', rejected_im_B + rejected_im_X, header
	endif


	; combine A and B into negative-positive-negative
	;*************************************

	if w_A ne !NULL and w_B ne !NULL then begin

	; A-B stack
	stack_AB = stack_A - stack_B

	; invert the A-B stack
	stack_BA = -stack_AB

	; find the dithering length, which for now we take from the first A and the first B frames
	; (A is always the one on top)
	dithering_length = round(diagnostics[w_A[0]].position - diagnostics[w_B[0]].position)

	; zero padding on top
	padding = dblarr( (size(stack_A))[1], dithering_length )
	padding[*] = 0.0
	AB_cleansum_padded = [ [stack_AB], [padding] ]
	BA_cleansum_padded = [ [stack_BA], [padding] ]

	BA_cleansum_shifted = shift(BA_cleansum_padded, 0, dithering_length)

	AB_combined = mean( [ [[AB_cleansum_padded]], [[BA_cleansum_shifted]] ], dimension=3, /nan )

	; ; last detail: the beginning and the end of the spectrum, outside the observed range, are NaN.
	; ; but the zero padding added zeroes there. Let's delete them

	; ; first find the x coordinates of the NaN regions, defined as columns where >90% of pixels (pre-padding) are NaNs
	; number_NaNs = total( ~finite(AB_cleansum), 2)
	; w_nans = where(number_NaNs GT 0.9 * (size(AB_cleansum))[2], /null )

	; ; now let's set all those columns equal to NaN
	; if w_nans NE !NULL then $
	; 	for i=0, n_elements(w_nans)-1 do AB_combined[w_nans[i],*] = !values.d_NaN

	; output final result
	writefits, filename_prefix + '_ABcombined.fits', AB_combined, header

	endif

END


; ---------------------------------------------------------------------------------------------------------------------------



PRO flame_combine_frames, fuel=fuel
 

	print, ' '
	print, 'Combine frames'
	print, '****************'

	; extract the slits 
	slits = *fuel.slits

 	; loop through all slits
	for i_slit=0, n_elements(slits)-1 do begin

		this_slit = slits[i_slit]
		print, 'Combining slit ', this_slit.number, ' - ', this_slit.name

		flame_combine_frames_oneslit, slit=this_slit, fuel=fuel

	endfor


END



