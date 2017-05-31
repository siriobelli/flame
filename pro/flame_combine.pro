;
;
; TO-DO:
;
; - handle better the case in which A or B are actually sky
;



FUNCTION flame_combine_stack, filenames=filenames, sigma_clip=sigma_clip, rejected_im=rejected_im, sigma_im=sigma_im
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
	for i_frame=0,N_frames-1 do im_cube[i_frame,*,*] = mrdfits(filenames[i_frame], 0, /silent)

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

	; check if there are more than one extensions in the FITS file (i.e. sigma images)
	fits_info, filenames[0], N_ext=N_ext, /silent

	; if sigma images exist, then calculate the combined sigma image
	if N_ext GE 2 then begin

			; read in the sigma images
			im_cube_sigma = dblarr( N_frames, (size(im_0))[1], (size(im_0))[2] )
			for i_frame=0,N_frames-1 do im_cube_sigma[i_frame,*,*] = mrdfits(filenames[i_frame], 1, /silent)

			; make an image with the number of non-rejected pixels
			nonrejected_im = N_frames - rejected_im

			; make the final sigma image by adding the uncertainties in quadrature
			sigma_im = sqrt( total(im_cube_sigma^2, 3, /nan)  ) / float(nonrejected_im)

	endif

	; finally stack the frames
	return,	mean(im_cube, dimension=1, /nan)

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_combine_oneslit, slit=slit, fuel=fuel

	; parameter for sigma-clipping when combining the frames
	sigma_clip = fuel.util.combine_sigma_clip

	; prefix for output file names
	filename_prefix = fuel.input.output_dir + 'slit' + string(slit.number, format='(I02)') + '-' + slit.name

	; input filenames for this slit
	filenames = slit.cutouts.filename

	; identify the A and B and X positions
	diagnostics = fuel.diagnostics

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
	stack_sky = flame_combine_stack(filenames=sky_filenames, sigma_clip=sigma_clip, sigma_im=stack_sky_sigma)

	; write out the sky spectrum
	writefits, filename_prefix + '_stack_sky.fits', stack_sky, header
	writefits, filename_prefix + '_stack_sky.fits', stack_sky_sigma, /append


	; stack all A, B, and X frames
	;*************************************

	if w_A ne !NULL then begin

		stack_A_filenames = flame_util_replace_string(filenames[w_A], '.fits', '_rectified.fits')
		stack_A = flame_combine_stack(filenames=stack_A_filenames, sigma_clip=sigma_clip, rejected_im=rejected_im_A, sigma_im=stack_A_sigma)
		writefits, filename_prefix + '_stack_A.fits', stack_A, header
		writefits, filename_prefix + '_stack_A.fits', stack_A_sigma, /append

		stack_A_skysub_filenames = flame_util_replace_string(filenames[w_A], '.fits', '_skysub_rectified.fits')
		stack_A_skysub = flame_combine_stack(filenames=stack_A_skysub_filenames, sigma_clip=sigma_clip, rejected_im=rejected_im_A, sigma_im=stack_A_skysub_sigma)
		writefits, filename_prefix + '_stack_A_skysub.fits', stack_A_skysub, header
		writefits, filename_prefix + '_stack_A_skysub.fits', stack_A_skysub_sigma, /append

	endif

	if w_B ne !NULL then begin

		stack_B_filenames = flame_util_replace_string(filenames[w_B], '.fits', '_rectified.fits')
		stack_B = flame_combine_stack(filenames=stack_B_filenames, sigma_clip=sigma_clip, rejected_im=rejected_im_B, sigma_im=stack_B_sigma)
		writefits, filename_prefix + '_stack_B.fits', stack_B, header
		writefits, filename_prefix + '_stack_B.fits', stack_B_sigma, /append

		stack_B_skysub_filenames = flame_util_replace_string(filenames[w_B], '.fits', '_skysub_rectified.fits')
		stack_B_skysub = flame_combine_stack(filenames=stack_B_skysub_filenames, sigma_clip=sigma_clip, rejected_im=rejected_im_B, sigma_im=stack_B_skysub_sigma)
		writefits, filename_prefix + '_stack_B_skysub.fits', stack_B_skysub, header
		writefits, filename_prefix + '_stack_B_skysub.fits', stack_B_skysub_sigma, /append

	endif

	if w_X ne !NULL then begin

		stack_X_filenames = flame_util_replace_string(filenames[w_X], '.fits', '_rectified.fits')
		stack_X = flame_combine_stack(filenames=stack_X_filenames, sigma_clip=sigma_clip, rejected_im=rejected_im_X, sigma_im=stack_X_sigma)
		writefits, filename_prefix + '_stack_X.fits', stack_X, header
		writefits, filename_prefix + '_stack_X.fits', stack_X_sigma, /append

		stack_X_skysub_filenames = flame_util_replace_string(filenames[w_X], '.fits', '_skysub_rectified.fits')
		stack_X_skysub = flame_combine_stack(filenames=stack_X_skysub_filenames, sigma_clip=sigma_clip, rejected_im=rejected_im_X, sigma_im=stack_X_skysub_sigma)
		writefits, filename_prefix + '_stack_X_skysub.fits', stack_X_skysub, header
		writefits, filename_prefix + '_stack_X_skysub.fits', stack_X_skysub_sigma, /append

	endif


	; combine A, B, and X stacks
	;*************************************

	if w_A NE !NULL and w_B ne !NULL then begin

		writefits, filename_prefix + '_stack_A-B.fits', stack_A - stack_B, header
		writefits, filename_prefix + '_stack_A-B.fits', sqrt(stack_A_sigma^2 + stack_B_sigma^2), /append

		writefits, filename_prefix + '_rejectedpixels_A-B.fits', rejected_im_A + rejected_im_B, header

		writefits, filename_prefix + '_stack_A-B_skysub.fits', stack_A_skysub - stack_B_skysub, header
		writefits, filename_prefix + '_stack_A-B_skysub.fits', sqrt(stack_A_skysub_sigma^2 + stack_B_skysub_sigma^2), /append

	endif

	if w_A NE !NULL and w_X ne !NULL then begin

		writefits, filename_prefix + '_stack_A-X.fits', stack_A - stack_X, header
		writefits, filename_prefix + '_stack_A-X.fits', sqrt(stack_A_sigma^2 + stack_X_sigma^2), /append

		writefits, filename_prefix + '_rejectedpixels_A-X.fits', rejected_im_A + rejected_im_X, header

		writefits, filename_prefix + '_stack_A-X_skysub.fits', stack_A_skysub - stack_X_skysub, header
		writefits, filename_prefix + '_stack_A-X_skysub.fits', sqrt(stack_A_skysub_sigma^2 + stack_X_skysub_sigma^2), /append

	endif

	if w_B NE !NULL and w_X ne !NULL then begin

		writefits, filename_prefix + '_stack_B-X.fits', stack_B - stack_X, header
		writefits, filename_prefix + '_stack_B-X.fits', sqrt(stack_B_sigma^2 + stack_X_sigma^2), /append

		writefits, filename_prefix + '_rejectedpixels_B-X.fits', rejected_im_B + rejected_im_X, header

		writefits, filename_prefix + '_stack_B-X_skysub.fits', stack_B_skysub - stack_X_skysub, header
		writefits, filename_prefix + '_stack_B-X_skysub.fits', sqrt(stack_B_skysub_sigma^2 + stack_X_skysub_sigma^2), /append

	endif

	; if only one among the A and B positions is available, then we are done
	if w_A eq !NULL or w_B eq !NULL then return


	; combine A and B into negative-positive-negative
	;*************************************

	; A-B stack
	stack_AB = stack_A_skysub - stack_B_skysub
	stack_AB_sigma = sqrt(stack_A_skysub_sigma^2 + stack_B_skysub_sigma^2)

	; NOTE: if A is not the "top" position, then this is actually B-A.
	; This ensures that the central trace of AB_combined is always positive
	if diagnostics[w_A[0]].position LT diagnostics[w_B[0]].position then $
		stack_AB = stack_B_skysub - stack_A_skysub

	; invert the A-B stack
	stack_BA = -stack_AB
	stack_BA_sigma = stack_AB_sigma

	; find the dithering length
	; (keep in mind that the rectification step already shifted each frame to the floor() of the reference position)
	dithering_length = abs( floor(diagnostics[w_A[0]].position) - floor(diagnostics[w_B[0]].position) )

	; if the dithering was not along the slit, then we are done
	if dithering_length GE (size(stack_A))[2] then return

	; zero padding on top
	padding = dblarr( (size(stack_A))[1], dithering_length )
	padding[*] = 0.0
	AB_cleansum_padded = [ [stack_AB], [padding] ]
	AB_cleansum_padded_sigma = [ [stack_AB_sigma], [padding] ]
	BA_cleansum_padded = [ [stack_BA], [padding] ]
	BA_cleansum_padded_sigma = [ [stack_BA_sigma], [padding] ]

	; shift by the dithering length
	BA_cleansum_shifted = shift(BA_cleansum_padded, 0, dithering_length)
	BA_cleansum_shifted_sigma = shift(BA_cleansum_padded_sigma, 0, dithering_length)

	; combine positive and negative
	AB_combined = mean( [ [[AB_cleansum_padded]], [[BA_cleansum_shifted]] ], dimension=3, /nan )
	AB_combined_sigma = 0.5 * sqrt( AB_cleansum_padded_sigma^2 + BA_cleansum_shifted_sigma^2 )

	; output final result
	writefits, filename_prefix + '_ABcombined.fits', AB_combined, header
	writefits, filename_prefix + '_ABcombined.fits', AB_combined_sigma, /append

	; also output the SNR map
	writefits, filename_prefix + '_ABcombined_SNR.fits', AB_combined/AB_combined_sigma, header


END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_combine_multislit, fuel=fuel
;
; if the dithering length matches the distance between two slits, then it mean-stack
; that these are the A and B positions for the same object, and we need to combine them
;

	; identify the A and B positions
	diagnostics = fuel.diagnostics
	w_A = where(diagnostics.offset_pos eq 'A', /null)
	w_B = where(diagnostics.offset_pos eq 'B', /null)

	; if either the A or B positions do not exist, then we are done
	if w_A eq !NULL or w_B eq !NULL then return

	; dithering length (by definition; see flame_combine_oneslit)
	dithering_length = abs( floor(diagnostics[w_A[0]].position) - floor(diagnostics[w_B[0]].position) )

	; number of pixel along the spatial position
	Nx = n_elements(*fuel.slits[0].rough_wavecal)

	; calculate the vertical coordinate of the geometric center of each slit
	slit_center = fltarr(n_elements(fuel.slits))
	for i_slit=0, n_elements(fuel.slits)-1 do slit_center[i_slit] = $
		poly( 0.5*Nx, fuel.slits[i_slit].bottom_poly) + 0.5*fuel.slits[i_slit].height

	; loop through the slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		; if the dithering length is smaller than the slit height, then it is an on-slit dithering
		if dithering_length LE fuel.slits[i_slit].height then continue
		; otherwise, check whether we have a good match among the slits

		; this is the distance of this slit to every other slit
		distance = slit_center-slit_center[i_slit]

		; this is the difference between the distance and the dithering length. in units of the slit height
		; (only match with slits that are above this one, to avoid double counting)
		delta = (distance - dithering_length) / fuel.slits.height

		; select the slit with delta closest to zero
		mindelta = min(abs(delta), j_slit)

		; if the dithering length falls within the central 50% of the slit, then we have a match
		if abs(mindelta) LT 0.5 then begin

			print, ''
			print, 'Combining slit ' + strtrim(fuel.slits[i_slit].number, 2) + ' - ' + fuel.slits[i_slit].name + $
				' with slit ' + strtrim(fuel.slits[j_slit].number, 2) + ' - ' + fuel.slits[j_slit].name

			; prefix for file names
			filename_prefix_i = fuel.input.output_dir + 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + $
			 	'-' + fuel.slits[i_slit].name
			filename_prefix_j = fuel.input.output_dir + 'slit' + string(fuel.slits[j_slit].number, format='(I02)') + $
			 	'-' + fuel.slits[j_slit].name

			; calculate the signs so that the stacked A-B has positive signal
			if floor(diagnostics[w_A[0]].position) GT floor(diagnostics[w_B[0]].position) then $
				signs = [-1, 1] else $
				signs = [1, -1]

			; combine tha A-B stacks
			flame_util_combine_slits, [filename_prefix_i + '_stack_A-B.fits', filename_prefix_j + '_stack_A-B.fits'], $
			 	output = fuel.input.output_dir + 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + '+slit' + $
				 string(fuel.slits[j_slit].number, format='(I02)') + '_stack_A-B.fits', signs = signs, $
				 sky_filenames=[filename_prefix_i + '_stack_sky.fits', filename_prefix_j + '_stack_sky.fits']

			; combine tha skysubtracted A-B stacks
 			flame_util_combine_slits, [filename_prefix_i + '_stack_A-B_skysub.fits', filename_prefix_j + '_stack_A-B_skysub.fits'], $
 			 	output = fuel.input.output_dir + 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + '+slit' + $
 				 string(fuel.slits[j_slit].number, format='(I02)') + '_stack_A-B_skysub.fits', signs = signs

		endif

	endfor


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_combine, fuel

		flame_util_module_start, fuel, 'flame_combine'


 	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		this_slit = fuel.slits[i_slit]
		print, 'Combining slit ' + strtrim(this_slit.number, 2) + ' - ' + this_slit.name

		flame_combine_oneslit, slit=this_slit, fuel=fuel

	endfor

	; if there is more than one slit, it may be necessary to combine two different slits together
	if n_elements(fuel.slits) GT 1 and fuel.input.AB_subtraction then flame_combine_multislit, fuel=fuel


	flame_util_module_end, fuel


  print, '-------------------------------------'
	print, '-------------------------------------'
	print, '-------------------------------------'

	; print total execution time
	print, ' '
	print, 'The data reduction took a total of ', $
		cgnumber_formatter((systime(/seconds) - fuel.util.start_time)/60.0, decimals=2), ' minutes.'
		print, ' '


END
