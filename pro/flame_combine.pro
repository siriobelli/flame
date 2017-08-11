;
;
; TO-DO:
;
; - handle better the case in which A or B are actually sky
;



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_combine_stack, filenames=filenames, output_filename=output_filename, $
	sigma_clip=sigma_clip
;
; read in FITS files and mean-stack them doing a sigma clipping
; write multi-HDU output FITS file:
; HDU 0: stacked spectrum
; HDU 1: error spectrum
; HDU 2: sigma map [i.e., standard deviation of values for each pixel]
; HDU 3: exptime map
;

	; number of frames
	N_frames = n_elements(filenames)

	; check that there is more than one frame
	if N_frames LE 1 then message, 'cannot stack only one frame'


	; construct the grid for the output image
	; ----------------------------------------------------------------------------

	; read header of first frame and get the grid of (lambda,gamma)
	header0 = headfits(filenames[0])
	lambda_min = sxpar(header0, 'CRVAL1')
	lambda_step = sxpar(header0, 'CDELT1')
	N_x = sxpar(header0, 'NAXIS1')
	gamma_min = sxpar(header0, 'CRVAL2')
	gamma_max = sxpar(header0, 'CRVAL2') + sxpar(header0, 'NAXIS2')

	; for each frame, read header and compare the grid
	for i_frame=1, N_frames-1 do begin

		header = headfits(filenames[i_frame])

		; check that they have the same lambda axis
		if (sxpar(header, 'CRVAL1')-lambda_min)/lambda_min GT 0.001 then message, 'wavelength axes not identical'
		if (sxpar(header, 'CDELT1')-lambda_step)/lambda_step GT 0.001 then message, 'wavelength axes not identical'

		; get the new extreme values for the gamma range
		if sxpar(header, 'CRVAL2') LT gamma_min then gamma_min = sxpar(header, 'CRVAL2')
		if sxpar(header, 'CRVAL2') + sxpar(header, 'NAXIS2') GT gamma_max then gamma_max = sxpar(header, 'CRVAL2') + sxpar(header, 'NAXIS2')

	endfor

	; check that final gamma range is reasonable
	if gamma_max - gamma_min GT 3000 then message, 'vertical dimension of the combined image is too large'

	; make big cube containing all images
	im_cube = dblarr( N_frames, N_x, gamma_max-gamma_min+1 )
	im_cube[*] = !values.d_nan

	; cube for the error spectra
	error_cube = dblarr( N_frames, N_x, gamma_max-gamma_min+1 )
	error_cube[*] = !values.d_nan

	; cube for the exptime
	exptime_cube = dblarr( N_frames, N_x, gamma_max-gamma_min+1 )
	exptime_cube[*] = !values.d_nan


	; read in all frames
	; ----------------------------------------------------------------------------

	; read in all frames
	for i_frame=0, N_frames-1 do begin

		; read in image and header
		im = mrdfits(filenames[i_frame], 0, header, /silent)

		; determine the spatial shift needed to align this with the cube
		y_start = sxpar(header, 'CRVAL2') - gamma_min
		y_end = y_start + sxpar(header, 'NAXIS2') - 1

		; insert the image at the right place in the cube
		im_cube[i_frame, *, y_start:y_end] = im

		; if present, read the error spectrum
		fits_info, filenames[i_frame], N_ext=N_ext, /silent
		if N_ext GE 1 then error_cube[i_frame, *, y_start:y_end] = $
		 	mrdfits(filenames[i_frame], 1, /silent)

		; make map of exposure time
		exptime = sxpar(header, 'EXPTIME')
		map_exptime = im*0.0 + exptime ; account for NaNs
		exptime_cube[i_frame, *, y_start:y_end] = map_exptime

	endfor


	; sigma-clipping and stacking
	; ----------------------------------------------------------------------------

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
	error_cube[where(mask_cube, /null)] = !values.d_nan
	exptime_cube[where(mask_cube, /null)] = !values.d_nan


	; make stack and output files
	; ----------------------------------------------------------------------------

	; stack the frames
	im_stack =	mean(im_cube, dimension=1, /nan)

	; make the error spectrum
	error_stack = sqrt( total(error_cube^2, 1, /nan)  ) / float( total(finite(mask_cube), 1))

	; make final map of exptime
	exptime_stack = total(exptime_cube, 1, /nan)

	; make header array with the correct grid
	sxaddpar, header0, 'EXTNAME', 'DATA'
	sxaddpar, header0, 'CRVAL2', gamma_min
	sxaddpar, header0, 'NAXIS2', gamma_max-gamma_min+1
	sxaddpar, header0, 'CRDELT2', 1.0

	; write out FITS file with stack
	writefits, output_filename, im_stack, header0

	; make header for extensions
	mkhdr, xten_hdr, error_stack, /image

	; add extension with error
	sxaddpar, xten_hdr, 'EXTNAME', 'NOISE'
	writefits, output_filename, error_stack, xten_hdr, /append

	; add extension with the pixel standard deviation
	sxaddpar, xten_hdr, 'EXTNAME', 'SIGMA'
	writefits, output_filename, sigma_im, xten_hdr, /append

	; add extension with exptime
	sxaddpar, xten_hdr, 'EXTNAME', 'EXPTIME'
	writefits, output_filename, exptime_stack, xten_hdr, /append

	return, im_stack


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************

FUNCTION flame_combine_AB, filename_top=filename_top, filename_bottom=filename_bottom, $
		dithering_length=dithering_length, error_image=error_image
	;
	; For on-source dithering, combine A and B traces together.
	; filename_top is the filename of the frame with the trace on top.
	; Optionally output error spectrum
	;

	; read in frames with error spectra too
	top = mrdfits(filename_top, 0, /silent)
	top_sigma = mrdfits(filename_top, 1, /silent)
	bottom = mrdfits(filename_bottom, 0, /silent)
	bottom_sigma = mrdfits(filename_bottom, 1, /silent)

	; A-B stack
	AB = top - bottom
	AB_sigma = sqrt(top_sigma^2 + bottom_sigma^2)

	; change NaNs into zeros, to properly account for the empty areas at the edges
	AB[where(~finite(AB), /null)] = 0.0

	; invert the A-B stack
	BA = -AB
	BA_sigma = AB_sigma

	; zero padding on top
	padding = dblarr( (size(top))[1], dithering_length )
	padding[*] = 0.0
	AB_padded = [ [AB], [padding] ]
	AB_padded_sigma = [ [AB_sigma], [padding] ]
	BA_padded = [ [BA], [padding] ]
	BA_padded_sigma = [ [BA_sigma], [padding] ]

	; shift by the dithering length
	BA_shifted = shift(BA_padded, 0, dithering_length)
	BA_shifted_sigma = shift(BA_padded_sigma, 0, dithering_length)

	; combine positive and negative
	AB_combined = mean( [ [[AB_padded]], [[BA_shifted]] ], dimension=3, /nan )
	error_image = sqrt( total( [ [[AB_padded_sigma^2]], [[BA_shifted_sigma^2]] ], 3, /nan ) ) / $
	 	double( total( finite([ [[AB_padded_sigma]], [[BA_shifted_sigma]] ]), 3 ) )

	; revert the true bad pixels (i.e., both A and B are NaN) into NaNs
	AB_combined[where(~finite(error_image), /null)] = !values.d_nan

	return, AB_combined

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_combine_oneslit, i_slit=i_slit, fuel=fuel

	; parameter for sigma-clipping when combining the frames
	sigma_clip = fuel.settings.combine_sigma_clip

	; prefix for output file names
	filename_prefix = fuel.util.output_dir + 'slit' + $
		string(fuel.slits[i_slit].number, format='(I02)') + '-' + fuel.slits[i_slit].name

	; input filenames for this slit
	filenames = fuel.slits[i_slit].cutouts.filename
	if fuel.slits[i_slit].cutouts[0].illcorr_applied then $
				filenames = flame_util_replace_string(filenames, '_corr', '_illcorr')

	; identify the A and B and X positions
	diagnostics = fuel.diagnostics

	; select all A frames
	w_A = where(diagnostics.offset_pos eq 'A', /null)

	; select all B frames
	w_B = where(diagnostics.offset_pos eq 'B', /null)

	; select all X frames
	w_X = where(diagnostics.offset_pos eq 'X', /null)

	; read in a header with the wavelength and gamma calibration
	header = headfits(flame_util_replace_string(filenames[0], '.fits', '_rectified.fits'))


	; combine sky spectra
	;*************************************

	; rectified but non-sky-subtracted files
	sky_filenames = flame_util_replace_string(filenames, '.fits', '_skymodel_rectified.fits')

	; stack and get the sky spectrum
	stack_sky = flame_combine_stack(filenames=sky_filenames, $
		output_filename = filename_prefix + '_sky.fits', sigma_clip=sigma_clip)


	; stack all A, B, and X frames
	;*************************************

	; go in reverse order of preference so that the output_file field at the end has
	; the preferred version (stack_A > stack_B > stack_X)

	if w_X ne !NULL then begin

		stack_X_filenames = flame_util_replace_string(filenames[w_X], '.fits', '_rectified.fits')
		stack_X = flame_combine_stack(filenames=stack_X_filenames, $
			output_filename=filename_prefix + '_X.fits', sigma_clip=sigma_clip)
		fuel.slits[i_slit].output_file = filename_prefix + '_X.fits'

		stack_X_skysub_filenames = flame_util_replace_string(filenames[w_X], '.fits', '_skysub_rectified.fits')
		stack_X_skysub = flame_combine_stack(filenames=stack_X_skysub_filenames, $
			output_filename=filename_prefix + '_skysub_X.fits', sigma_clip=sigma_clip)

	endif

	if w_B ne !NULL then begin

		stack_B_filenames = flame_util_replace_string(filenames[w_B], '.fits', '_rectified.fits')
		stack_B = flame_combine_stack(filenames=stack_B_filenames, $
		 	output_filename=filename_prefix + '_B.fits', sigma_clip=sigma_clip)
		fuel.slits[i_slit].output_file = filename_prefix + '_B.fits'

		stack_B_skysub_filenames = flame_util_replace_string(filenames[w_B], '.fits', '_skysub_rectified.fits')
		stack_B_skysub = flame_combine_stack(filenames=stack_B_skysub_filenames, $
		 	output_filename=filename_prefix + '_skysub_B.fits', sigma_clip=sigma_clip)

	endif

	if w_A ne !NULL then begin

		stack_A_filenames = flame_util_replace_string(filenames[w_A], '.fits', '_rectified.fits')
		stack_A = flame_combine_stack(filenames=stack_A_filenames, $
			output_filename=filename_prefix + '_A.fits', sigma_clip=sigma_clip)
		fuel.slits[i_slit].output_file = filename_prefix + '_A.fits'

		stack_A_skysub_filenames = flame_util_replace_string(filenames[w_A], '.fits', '_skysub_rectified.fits')
		stack_A_skysub = flame_combine_stack(filenames=stack_A_skysub_filenames, $
			output_filename=filename_prefix + '_skysub_A.fits', sigma_clip=sigma_clip)

	endif

stop

	; combine A, B, and X stacks
	;*************************************

	if w_B NE !NULL and w_X ne !NULL then begin

		writefits, filename_prefix + '_B-X.fits', stack_B - stack_X, header
		writefits, filename_prefix + '_B-X.fits', sqrt(stack_B_sigma^2 + stack_X_sigma^2), /append

		writefits, filename_prefix + '_rejectedpixels_B-X.fits', rejected_im_B + rejected_im_X, header

		writefits, filename_prefix + '_skysub_B-X.fits', stack_B_skysub - stack_X_skysub, header
		writefits, filename_prefix + '_skysub_B-X.fits', sqrt(stack_B_skysub_sigma^2 + stack_X_skysub_sigma^2), /append

		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_B-X.fits'

	endif

	if w_A NE !NULL and w_X ne !NULL then begin

		writefits, filename_prefix + '_A-X.fits', stack_A - stack_X, header
		writefits, filename_prefix + '_A-X.fits', sqrt(stack_A_sigma^2 + stack_X_sigma^2), /append

		writefits, filename_prefix + '_rejectedpixels_A-X.fits', rejected_im_A + rejected_im_X, header

		writefits, filename_prefix + '_skysub_A-X.fits', stack_A_skysub - stack_X_skysub, header
		writefits, filename_prefix + '_skysub_A-X.fits', sqrt(stack_A_skysub_sigma^2 + stack_X_skysub_sigma^2), /append

		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_A-X.fits'

	endif

	if w_A NE !NULL and w_B ne !NULL then begin

		writefits, filename_prefix + '_A-B.fits', stack_A - stack_B, header
		writefits, filename_prefix + '_A-B.fits', sqrt(stack_A_sigma^2 + stack_B_sigma^2), /append

		writefits, filename_prefix + '_rejectedpixels_A-B.fits', rejected_im_A + rejected_im_B, header

		writefits, filename_prefix + '_skysub_A-B.fits', stack_A_skysub - stack_B_skysub, header
		writefits, filename_prefix + '_skysub_A-B.fits', sqrt(stack_A_skysub_sigma^2 + stack_B_skysub_sigma^2), /append

		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_A-B.fits'

	endif

	; if only one among the A and B positions is available, then we are done
	if w_A eq !NULL or w_B eq !NULL then return


	; combine A and B into negative-positive-negative
	;*************************************

	; find the dithering length
	; (keep in mind that the rectification step already shifted each frame to the floor() of the reference position)
	dithering_length = abs( floor(diagnostics[w_A[0]].position) - floor(diagnostics[w_B[0]].position) )

	; if the dithering was not along the slit, then we are done
	if dithering_length GE (size(stack_A))[2] then return

	; combine A and B traces

	; make sure to pick the right offset position for the "top" trace
		if diagnostics[w_A[0]].position GT diagnostics[w_B[0]].position then begin
			AB = flame_combine_AB( filename_top=filename_prefix + '_A.fits', filename_bottom=filename_prefix + '_B.fits', $
				dithering_length=dithering_length, error_image=AB_sigma )
			AB_skysub = flame_combine_AB( filename_top=filename_prefix + '_skysub_A.fits', filename_bottom=filename_prefix + '_skysub_B.fits', $
				dithering_length=dithering_length, error_image=AB_skysub_sigma )
		endif else begin
			AB = flame_combine_AB( filename_top=filename_prefix + '_B.fits', filename_bottom=filename_prefix + '_A.fits', $
				dithering_length=dithering_length, error_image=AB_sigma )
			AB_skysub = flame_combine_AB( filename_top=filename_prefix + '_skysub_B.fits', filename_bottom=filename_prefix + '_skysub_A.fits', $
				dithering_length=dithering_length, error_image=AB_skysub_sigma )
		endelse

	; output final result
	writefits, filename_prefix + '_ABcombined.fits', AB, header
	writefits, filename_prefix + '_ABcombined.fits', AB_sigma, /append
	writefits, filename_prefix + '_skysub_ABcombined.fits', AB_skysub, header
	writefits, filename_prefix + '_skysub_ABcombined.fits', AB_skysub_sigma, /append

	; also output the SNR map
	writefits, filename_prefix + '_ABcombined_SNR.fits', AB/AB_sigma, header
	writefits, filename_prefix + '_skysub_ABcombined_SNR.fits', AB_skysub/AB_skysub_sigma, header

	; update the filename of the final output
	fuel.slits[i_slit].output_file = filename_prefix + '_skysub_ABcombined.fits'


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_combine_multislit, fuel=fuel
;
; if the dithering length matches the distance between two slits, then it means
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

	; number of pixels along the spatial position
	Nx = n_elements(*fuel.slits[0].rough_skylambda)

	; calculate the vertical coordinate of the geometric center of each slit
	slit_center = fltarr(n_elements(fuel.slits))
	for i_slit=0, n_elements(fuel.slits)-1 do slit_center[i_slit] = $
		poly( 0.5*Nx, fuel.slits[i_slit].bottom_poly) + 0.5*fuel.slits[i_slit].height

	; loop through the slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		if fuel.slits[i_slit].skip then continue

		; if the dithering length is clearly smaller than the slit height, then it is an on-slit dithering
		if dithering_length LE 0.85*fuel.slits[i_slit].height then continue
		; otherwise, check whether we have a good match among the slits

		; this is the distance of this slit to every other slit
		distance = slit_center-slit_center[i_slit]

		; this is the difference between the distance and the dithering length. in units of the slit height
		; (only match with slits that are above this one, to avoid double counting)
		delta = (distance - dithering_length) / fuel.slits.height

		; select the slit with delta closest to zero
		mindelta = min(abs(delta), j_slit)

		; check that the candidate slit is not skipped
		if fuel.slits[j_slit].skip then continue

		; if the dithering length falls within the central 50% of the slit, then we have a match
		if abs(mindelta) LT 0.5 then $
			; still need to check that these two slits are horizontally aligned
			; do this by comparing the wavelength solutions
			if fuel.slits[i_slit].outlambda_min eq fuel.slits[j_slit].outlambda_min then begin

			print, ''
			print, 'Combining slit ' + strtrim(fuel.slits[i_slit].number, 2) + ' - ' + fuel.slits[i_slit].name + $
				' with slit ' + strtrim(fuel.slits[j_slit].number, 2) + ' - ' + fuel.slits[j_slit].name

			; prefix for file names
			filename_prefix_i = fuel.util.output_dir + 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + $
			 	'-' + fuel.slits[i_slit].name
			filename_prefix_j = fuel.util.output_dir + 'slit' + string(fuel.slits[j_slit].number, format='(I02)') + $
			 	'-' + fuel.slits[j_slit].name

			; calculate the signs so that the stacked A-B has positive signal
			if floor(diagnostics[w_A[0]].position) GT floor(diagnostics[w_B[0]].position) then $
				signs = [-1, 1] else $
				signs = [1, -1]

			; combine the A-B stacks
			outname = fuel.util.output_dir + 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + '+slit' + $
				string(fuel.slits[j_slit].number, format='(I02)') + '_A-B.fits'
			flame_util_combine_slits, [filename_prefix_i + '_A-B.fits', filename_prefix_j + '_A-B.fits'], $
			 	output = outname, signs = signs, $
				 sky_filenames=[filename_prefix_i + '_sky.fits', filename_prefix_j + '_sky.fits']

			; combine the skysubtracted A-B stacks
			outname = fuel.util.output_dir + 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + '+slit' + $
			 string(fuel.slits[j_slit].number, format='(I02)') + '_skysub_A-B.fits'
 			flame_util_combine_slits, [filename_prefix_i + '_skysub_A-B.fits', filename_prefix_j + '_skysub_A-B.fits'], $
 			 	output = outname, signs = signs

			; update the file name of the final output
			fuel.slits[i_slit].output_file = outname

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

		if fuel.slits[i_slit].skip then continue

		print, 'Combining slit ' + strtrim(fuel.slits[i_slit].number, 2) + ' - ' + fuel.slits[i_slit].name

		; handle errors by ignoring that slit
		if fuel.settings.debugging eq 0 then begin
			catch, error_status
			if error_status ne 0 then begin
				print, ''
		    print, '**************************'
		    print, '***       WARNING      ***'
		    print, '**************************'
		    print, 'Error found. Skipping slit ' + strtrim(fuel.slits[i_slit].number,2), ' - ', fuel.slits[i_slit].name
				fuel.slits[i_slit].skip = 1
				catch, /cancel
				continue
			endif
		endif

		flame_combine_oneslit, i_slit=i_slit, fuel=fuel

	endfor

	; if there is more than one slit, it may be necessary to combine two different slits together
	if n_elements(where(fuel.slits.skip eq 0, /null)) GT 1 and fuel.input.AB_subtraction then flame_combine_multislit, fuel=fuel


	flame_util_module_end, fuel

END
