

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_combine_stack, fuel=fuel, filenames=filenames, sky_filenames=sky_filenames, diagnostics=diagnostics, $
	 output_filename=output_filename, noalign=noalign
;
; Read in FITS files and mean-stack them after a sigma clipping.
; Alignment is done using the gamma coordinate, which is the vertical wcs coordinate,
; unless the option /noalign is set.
; NB: the grid of the first filename is assumed for the output file
;
; Write multi-HDU output FITS file:
; HDU 0: stacked spectrum
; HDU 1: error spectrum
; HDU 2: sigma map [i.e., standard deviation of values for each pixel, including correction for correlated noise]
; HDU 3: model sky
; HDU 4: exptime map
; HDU 5: weight map
;

	; number of frames
	N_frames = n_elements(filenames)

	; check that there are frames
	if N_frames EQ 0 then message, 'no frames to stack'

	; if there is only one frame, then print a warning
	if N_frames EQ 1 then print, 'Warning: stacking only one frame'

	; check that, if present, the sky frames are in the correct number
	if n_elements(sky_filenames) NE 0 and n_elements(sky_filenames) NE N_frames then $
	 	message, 'sky_filenames: wrong number of frames'


	; construct the grid for the output image
	; ----------------------------------------------------------------------------

	; read header of first frame and get the grid of (lambda,gamma)
	header0 = headfits(filenames[0])
	lambda_min = sxpar(header0, 'CRVAL1')
	lambda_step = sxpar(header0, 'CDELT1')
	N_x = sxpar(header0, 'NAXIS1')
	N_y = sxpar(header0, 'NAXIS2')
	gamma_min = sxpar(header0, 'CRVAL2')
	gamma_max = gamma_min + N_y - 1

	; make big cube containing all images
	im_cube = dblarr( N_frames, N_x, N_y )
	im_cube[*] = !values.d_nan

	; cube for the error spectra
	error_cube = dblarr( N_frames, N_x, N_y )
	error_cube[*] = !values.d_nan

	; cube for the sky spectra
	sky_cube = dblarr( N_frames, N_x, N_y )
	sky_cube[*] = !values.d_nan

	; cube for the exptime
	exptime_cube = dblarr( N_frames, N_x, N_y )
	exptime_cube[*] = !values.d_nan


	; read in all frames
	; ----------------------------------------------------------------------------

	; read in all frames
	for i_frame=0, N_frames-1 do begin

		; read in image and header
		im = mrdfits(filenames[i_frame], 0, header, /silent)

		; check that the lambda axis is the same as for the first file
		if (sxpar(header, 'CRVAL1')-lambda_min)/lambda_min GT 0.001 then message, 'wavelength axes not identical'
		if (sxpar(header, 'CDELT1')-lambda_step)/lambda_step GT 0.001 then message, 'wavelength axes not identical'

		; read in the spatial grid
		N_y_i = sxpar(header, 'NAXIS2')
		gamma_min_i = sxpar(header, 'CRVAL2')
		gamma_max_i = gamma_min_i + N_y_i - 1

		if keyword_set(noalign) then begin

			; do not align, simply cut out the top if needed
			bot_i = 0
			bot_ref = 0
			top_i = min([N_y,N_y_i])-1
			top_ref = min([N_y,N_y_i])-1

		endif else begin

			; determine if the vertical dimensions are on the same grid
			shift = gamma_min_i - gamma_min
			if shift ne fix(shift) then message, 'Before stacking the files need to be resampled along the spatial direction!'

			; determine the starting and ending pixels for the proper alignment of the frame
			if gamma_min_i GE gamma_min then bot_i = 0 else bot_i = gamma_min-gamma_min_i
			if gamma_min_i GE gamma_min then bot_ref = gamma_min_i-gamma_min else bot_ref = 0
			if gamma_max_i GE gamma_max then top_i = gamma_max-gamma_min_i  else top_i = gamma_max_i-gamma_min_i
			if gamma_max_i GE gamma_max then top_ref = gamma_max-gamma_min else top_ref = gamma_max_i-gamma_min

		endelse

		; add the image to the cube
		im_cube[i_frame, *, bot_ref:top_ref] = im[*, bot_i:top_i]

		; if present, read the error spectrum
		fits_info, filenames[i_frame], N_ext=N_ext, /silent
		if N_ext GE 1 then begin
			err = mrdfits(filenames[i_frame], 1, /silent)
			error_cube[i_frame, *, bot_ref:top_ref] = err[*, bot_i:top_i]
		endif

		; if present, read the model sky
		if n_elements(sky_filenames) GE 1 then begin
			sky = mrdfits(sky_filenames[i_frame], 0, /silent)
			sky_cube[i_frame, *, bot_ref:top_ref] = sky[*, bot_i:top_i]
		endif

		; make map of exposure time
		exptime = sxpar(header, 'EXPTIME')
		map_exptime = im*0.0 + exptime ; account for NaNs
		exptime_cube[i_frame, *, bot_ref:top_ref] = map_exptime[*, bot_i:top_i]


	endfor


	; sigma-clipping
	; ----------------------------------------------------------------------------

	; calculate median at each pixel of the image
	median_im = median(im_cube, dimension=1)

	; make a cube with the corresponding median at each (pixel,frame) position
	median_cube = im_cube
	median_cube[*]=0.0
	for i_frame=0,N_frames-1 do median_cube[i_frame, *, *] = median_im

	; calculate standard deviation at each pixel of the image
	sigma_im = stddev(im_cube, dimension=1, /nan)

	; make a cube with the corresponding sigma at each (pixel,frame) position
	sigma_cube = im_cube
	sigma_cube[*]=0.0
	for i_frame=0,N_frames-1 do sigma_cube[i_frame, *, *] = sigma_im

	; for each pixel in the cube, calculate the deviation from the median in units of sigma
	deviation = (im_cube - median_cube)/sigma_cube

	; and now make a mask that indicates those pixels that are more than sigma_clip off from the median
	mask_cube = fix(im_cube)
	mask_cube[*] = 0
	if fuel.settings.combine_sigma_clip GT 0.0 then $
		mask_cube[where( abs(deviation) GT fuel.settings.combine_sigma_clip, /null )] = 1

	; turn masked pixels into NaNs
	im_cube[where(mask_cube, /null)] = !values.d_nan
	error_cube[where(mask_cube, /null)] = !values.d_nan
	sky_cube[where(mask_cube, /null)] = !values.d_nan
	exptime_cube[where(mask_cube, /null)] = !values.d_nan

	; make a cube that flags pixels that are actually good
	goodpix_cube = byte(im_cube*0.0)
	goodpix_cube = finite(im_cube)

	; for each pixel of the final image, how many frames did actually contribute?
	im_goodpix = total(goodpix_cube, 1)

	; set the weights for each frame
	case strlowcase(fuel.settings.frame_weights) of
		'none': weights = replicate(1.0, N_frames)
		'flux': weights = diagnostics.flux
		'seeing': weights = 1.0/diagnostics.seeing
		'peak': weights = diagnostics.flux / diagnostics.seeing
		else: message, 'fuel.settings.frame_weights: only valid options are: none, flux, seeing, peak'
	endcase

	; for convenience, scale the weights so that the maximum is one
	weights /= max(weights, /nan)
	w_noweight = where( ~finite(weights), /null )
	if w_noweight NE !NULL then begin
		print, 'WARNING: some of the frames have zero weight!'
		weights[w_noweight] = 0.0
	endif

	; make the weight cube
	weight_cube = im_cube
	weight_cube[*]=0.0
	for i_frame=0,N_frames-1 do weight_cube[i_frame, *, *] = weights[i_frame]

	; add the good pixel mask, so that weights are zero for bad pixels
	weight_cube *= finite(im_cube)


	; make stack and output files
	; ----------------------------------------------------------------------------

	; mean-stack the frames
	im_stack = total( weight_cube * im_cube, 1, /nan ) / total(weight_cube, 1, /nan)

	; make the error spectrum
	error_stack = sqrt( total( weight_cube^2 * error_cube^2, 1, /nan)  ) / total( weight_cube, 1, /nan)

	; mean-stack the sky model
	sky_stack = total( weight_cube * sky_cube, 1, /nan ) / total(weight_cube, 1, /nan)

	; make a clean sigma image (i.e., excluding rejected pixels)
	weight_tot_sq = total( weight_cube^2, 1, /nan )
	weight_tot = total( weight_cube, 1, /nan )
	im_stack_cube = im_cube
	im_stack_cube[*] = 0.0
	for i_frame=0,N_frames-1 do im_stack_cube[i_frame, *, *] = im_stack
	sigma_stack = sqrt( weight_tot_sq / ( weight_tot^2 - weight_tot_sq ) ) * $
		sqrt( total( weight_cube * (im_cube-im_stack_cube)^2, 1, /nan ) / weight_tot )

	; correct for the correlated noise
	; see Eq. 9 in Fruchter & Hook 2002, and also footnote 20 in Kriek et al. 2015
	sigma_stack *= 1.5

	; make final map of exptime
	exptime_stack = total(exptime_cube, 1, /nan)

	; require a minimum contribution in terms of number of frames for a pixel to be valid
	w_void = where(im_goodpix LE fuel.settings.combine_min_framefrac*N_frames)
	im_stack[w_void] = !values.d_nan
	error_stack[w_void] = !values.d_nan
	sigma_stack[w_void] = !values.d_nan
	sky_stack[w_void] = !values.d_nan
	exptime_stack[w_void] = 0.0

	; make header array with the correct grid
	sxaddpar, header0, 'EXTNAME', 'DATA'
	sxaddpar, header0, 'CRVAL2', gamma_min
	sxaddpar, header0, 'NAXIS2', gamma_max-gamma_min+1
	sxaddpar, header0, 'CRDELT2', 1.0
	; leave the YCUTOUT value of the first frame in the header

	; write out FITS file with stack
	writefits, output_filename, im_stack, header0

	; make header for extensions
	mkhdr, xten_hdr, error_stack, /image

	; add extension with error
	sxaddpar, xten_hdr, 'EXTNAME', 'NOISE'
	writefits, output_filename, error_stack, xten_hdr, /append

	; add extension with the pixel standard deviation
	sxaddpar, xten_hdr, 'EXTNAME', 'SIGMA'
	writefits, output_filename, sigma_stack, xten_hdr, /append

	; add extension with sky
	sxaddpar, xten_hdr, 'EXTNAME', 'SKY'
	writefits, output_filename, sky_stack, xten_hdr, /append

	; add extension with exptime
	sxaddpar, xten_hdr, 'EXTNAME', 'EXPTIME'
	writefits, output_filename, exptime_stack, xten_hdr, /append

	; add extension with total weight map
	sxaddpar, xten_hdr, 'EXTNAME', 'WEIGHT'
	writefits, output_filename, weight_tot, xten_hdr, /append

	print, output_filename, ' written'

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_combine_diff, filename1=filename1, filename2=filename2, $
		output_filename=output_filename, combined_filename=combined_filename
;
; make difference image
; the input files must be written by flame_combine_stack and have six extensions
; if combined_filename is specified, then the ABcombined image is also written,
; using the gamma coordinate to align
; NB: the alignment is done in the "observed" frame using the YCUTOUT keyword
; from the FITS header
;

	; read in the two files
	; ----------------------------------------------------------------------------

	im1 = mrdfits(filename1, 0, hdr1, /silent)
	im2 = mrdfits(filename2, 0, hdr2, /silent)

	err1 = mrdfits(filename1, 1, err_hdr, /silent)
	err2 = mrdfits(filename2, 1, /silent)

	sig1 = mrdfits(filename1, 2, sig_hdr, /silent)
	sig2 = mrdfits(filename2, 2, /silent)

	sky1 = mrdfits(filename1, 3, sky_hdr, /silent)
	sky2 = mrdfits(filename2, 3, /silent)

	exptime1 = mrdfits(filename1, 4, exptime_hdr, /silent)
	exptime2 = mrdfits(filename2, 4, /silent)

	weight1 = mrdfits(filename1, 5, weight_hdr, /silent)
	weight2 = mrdfits(filename2, 5, /silent)


	; align the two frames
	; ----------------------------------------------------------------------------

	; read the y coordinate of the (0,0) pixel for each frame
	ymin1 = round(sxpar(hdr1, 'YCUTOUT'))
	ymin2 = round(sxpar(hdr2, 'YCUTOUT'))

	; read the y coordinate of the (0,Ny) pixel for each frame
	ymax1 = ymin1 + sxpar(hdr1, 'NAXIS2') - 1
	ymax2 = ymin2 + sxpar(hdr2, 'NAXIS2') - 1

	; determine the starting and ending pixels for the proper alignment of the frame
	if ymin2 GE ymin1 then bot2 = 0 else bot2 = ymin1-ymin2
	if ymin2 GE ymin1 then bot1 = ymin2-ymin1 else bot1 = 0
	if ymax2 GE ymax1 then top2 = ymax1-ymin2  else top2 = ymax2-ymin2
	if ymax2 GE ymax1 then top1 = ymax1-ymin1 else top1 = ymax2-ymin1

	; read the gamma coordinate (to be used later on)
	gamma_min1 = sxpar(hdr1, 'CRVAL2')
	gamma_min2 = sxpar(hdr2, 'CRVAL2')


	; combine the data
	; ----------------------------------------------------------------------------

	; make the difference image
	imdiff = im1[*,bot1:top1] - im2[*,bot2:top2]

	; combine the noise
	errdiff = sqrt( err1[*,bot1:top1]^2 + err2[*,bot2:top2]^2 )

	; combine the sigma
	sigdiff = sqrt( sig1[*,bot1:top1]^2 + sig2[*,bot2:top2]^2 )

	; combine the sky
	skydiff = sky1[*,bot1:top1] + sky2[*,bot2:top2]

	; take only the exptime of A
	exptimediff1 = exptime1[*,bot1:top1]
	exptimediff2 = exptime2[*,bot2:top2]
	exptimediff = exptimediff1

	; take only the weight of A
	weightdiff1 = weight1[*,bot1:top1]
	weightdiff2 = weight2[*,bot2:top2]
	weightdiff = weightdiff1


	; output file
	; ----------------------------------------------------------------------------

	; make the header for the combined spectrum
	hdr_diff = hdr1
	if ymin1 ne 0 then sxaddpar, hdr_diff, 'YCUTOUT', sxpar(hdr_diff, 'YCUTOUT')+bot1
	if ymin1 ne 0 then sxaddpar, hdr_diff, 'CRVAL2', sxpar(hdr_diff, 'CRVAL2')+bot1

	writefits, output_filename, imdiff, hdr_diff
	writefits, output_filename, errdiff, err_hdr, /append
	writefits, output_filename, sigdiff, sig_hdr, /append
	writefits, output_filename, skydiff, sky_hdr, /append
	writefits, output_filename, exptimediff, exptime_hdr, /append
	writefits, output_filename, weightdiff, weight_hdr, /append
	print, output_filename, ' written'


	; double-combine frames
	; ----------------------------------------------------------------------------

	; if the keyword is not provided, then skip
	if ~keyword_set(combined_filename) then return

	; what is the nod amplitude?
	nod = (gamma_min1-ymin1) - (gamma_min2-ymin2)
	print, 'Nod amplitude: ' + strtrim(nod, 2) + ' rectified pixels'

	; height of the A-B image
	Ny = (size(imdiff))[2]

	; if the nod amplitude is zero then something is wrong
	if nod eq 0 then message, 'Nod amplitude cannot be zero'

	; if the nod amplitude is almost the same as the frame height, then we are done
	if abs(nod+0.0) GT 0.9*Ny then return

	; make an empty cube for the stacking
	cube_imdiff = dblarr( 2, (size(imdiff))[1], Ny + abs(nod) )
	cube_imdiff[*] = !values.d_nan
	cube_errdiff = cube_imdiff
	cube_sigdiff = cube_imdiff
	cube_skydiff = cube_imdiff
	cube_exptimediff = cube_imdiff
	cube_weightdiff = cube_imdiff

	; determine the amount of vertical shifting
	shiftA=0
	shiftB=0
	if nod GT 0 then shiftA+=nod else shiftB+=abs(nod)

	; add the difference image in the first layer of the cube
	cube_imdiff[0,*,0+shiftA:Ny-1+shiftA] = imdiff
	cube_errdiff[0,*,0+shiftA:Ny-1+shiftA] = errdiff
	cube_sigdiff[0,*,0+shiftA:Ny-1+shiftA] = sigdiff
	cube_skydiff[0,*,0+shiftA:Ny-1+shiftA] = skydiff
	cube_exptimediff[0,*,0+shiftA:Ny-1+shiftA] = exptimediff1
	cube_weightdiff[0,*,0+shiftA:Ny-1+shiftA] = weightdiff1

	; add the negative of the difference image in the second layer of the cube
	cube_imdiff[1,*,0+shiftB:Ny-1+shiftB] = -imdiff
	cube_errdiff[1,*,0+shiftB:Ny-1+shiftB] = errdiff
	cube_sigdiff[1,*,0+shiftB:Ny-1+shiftB] = sigdiff
	cube_skydiff[1,*,0+shiftB:Ny-1+shiftB] = skydiff
	cube_exptimediff[1,*,0+shiftB:Ny-1+shiftB] = exptimediff2
	cube_weightdiff[1,*,0+shiftB:Ny-1+shiftB] = weightdiff2

	; weights for stacking
	cube_weights = cube_weightdiff

	; finally stack the cubes
	dbl_imdiff = total(cube_weights*cube_imdiff, 1, /nan) / total(cube_weights, 1, /nan)
	dbl_errdiff = sqrt( total(cube_weights^2 * cube_errdiff^2, 1, /nan) ) / total(cube_weights, 1, /nan)
	dbl_sigdiff = sqrt( total(cube_weights^2 * cube_sigdiff^2, 1, /nan) ) / total(cube_weights, 1, /nan)
	dbl_skydiff = total(cube_weights*cube_skydiff, 1, /nan) / total(cube_weights, 1, /nan)
	dbl_exptimediff = total(cube_exptimediff, 1, /nan)
	dbl_weightdiff = total(cube_weights, 1, /nan)

	; set to NaNs pixels with exptime=0
	w_nan = where(dbl_exptimediff eq 0.0, /null)
	if w_nan NE !NULL then begin
		dbl_imdiff[w_nan] = !values.d_nan
		dbl_errdiff[w_nan] = !values.d_nan
		dbl_sigdiff[w_nan] = !values.d_nan
		dbl_skydiff[w_nan] = !values.d_nan
	endif

	; make the header for the double-combined spectrum
	hdr_dbl = hdr_diff
	if shiftA ne 0 then sxaddpar, hdr_dbl, 'YCUTOUT', sxpar(hdr_dbl, 'YCUTOUT')+shiftA
	if shiftA ne 0 then sxaddpar, hdr_dbl, 'CRVAL2', sxpar(hdr_dbl, 'CRVAL2')+shiftA

	; write out FITS file
	writefits, combined_filename, dbl_imdiff, hdr_dbl
	writefits, combined_filename, dbl_errdiff, err_hdr, /append
	writefits, combined_filename, dbl_sigdiff, sig_hdr, /append
	writefits, combined_filename, dbl_skydiff, sky_hdr, /append
	writefits, combined_filename, dbl_exptimediff, exptime_hdr, /append
	writefits, combined_filename, dbl_weightdiff, weight_hdr, /append
	print, combined_filename, ' written'


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
	flame_combine_stack, fuel=fuel, filenames=sky_filenames, diagnostics=fuel.diagnostics, $
		output_filename = filename_prefix + '_sky.fits', /noalign


	; stack all A, B, and X frames
	;*************************************

	; get all the rectified and skysub filenames
	filenames_rectified = flame_util_replace_string(filenames,  '.fits', '_rectified.fits')
	filenames_skysub_rectified = flame_util_replace_string(filenames,  '.fits', '_skysub_rectified.fits')

	; go in reverse order of preference so that the output_file field at the end has
	; the preferred version (stack_A > stack_B > stack_X)

	if w_X ne !NULL then begin
		flame_combine_stack, fuel=fuel, filenames=filenames_rectified[w_X], sky_filenames=sky_filenames[w_X], $
		 	diagnostics=fuel.diagnostics[w_X], output_filename=filename_prefix + '_X.fits'
		flame_combine_stack, fuel=fuel, filenames=filenames_skysub_rectified[w_X], sky_filenames=sky_filenames[w_X], $
		 	diagnostics=fuel.diagnostics[w_X], output_filename=filename_prefix + '_skysub_X.fits'
		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_X.fits'
	endif

	if w_B ne !NULL then begin
		flame_combine_stack, fuel=fuel, filenames=filenames_rectified[w_B], sky_filenames=sky_filenames[w_B], $
			diagnostics=fuel.diagnostics[w_B], output_filename=filename_prefix + '_B.fits'
		flame_combine_stack, fuel=fuel, filenames=filenames_skysub_rectified[w_B], sky_filenames=sky_filenames[w_B], $
			diagnostics=fuel.diagnostics[w_B], output_filename=filename_prefix + '_skysub_B.fits'
		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_B.fits'
	endif

	if w_A ne !NULL then begin
		flame_combine_stack, fuel=fuel, filenames=filenames_rectified[w_A], sky_filenames=sky_filenames[w_A], $
			diagnostics=fuel.diagnostics[w_A], output_filename=filename_prefix + '_A.fits'
		flame_combine_stack, fuel=fuel, filenames=filenames_skysub_rectified[w_A], sky_filenames=sky_filenames[w_A], $
			diagnostics=fuel.diagnostics[w_A], output_filename=filename_prefix + '_skysub_A.fits'
		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_A.fits'
	endif


	; make difference images
	;*************************************


	if w_B NE !NULL and w_X ne !NULL then begin
		flame_combine_diff, filename1=filename_prefix+'_B.fits', filename2=filename_prefix+'_X.fits', $
			output_filename=filename_prefix + '_B-X.fits'
		flame_combine_diff, filename1=filename_prefix+'_skysub_B.fits', filename2=filename_prefix+'_skysub_X.fits', $
			output_filename=filename_prefix + '_skysub_B-X.fits'
		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_B-X.fits'
	endif


	if w_A NE !NULL and w_X ne !NULL then begin
		flame_combine_diff, filename1=filename_prefix+'_A.fits', filename2=filename_prefix+'_X.fits', $
			output_filename=filename_prefix + '_A-X.fits'
		flame_combine_diff, filename1=filename_prefix+'_skysub_A.fits', filename2=filename_prefix+'_skysub_X.fits', $
			output_filename=filename_prefix + '_skysub_A-X.fits'
		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_A-X.fits'
	endif


	if w_A NE !NULL and w_B ne !NULL then begin

		flame_combine_diff, filename1=filename_prefix+'_A.fits', filename2=filename_prefix+'_B.fits', $
			output_filename=filename_prefix + '_A-B.fits', combined_filename=filename_prefix + '_ABcombined.fits'
		flame_combine_diff, filename1=filename_prefix+'_skysub_A.fits', filename2=filename_prefix+'_skysub_B.fits', $
			output_filename=filename_prefix + '_skysub_A-B.fits', combined_filename=filename_prefix + '_skysub_ABcombined.fits'
		fuel.slits[i_slit].output_file = filename_prefix + '_skysub_A-B.fits'

		; if written, then use the ABcombined file
		if file_test(filename_prefix + '_skysub_ABcombined.fits') then $
			fuel.slits[i_slit].output_file = filename_prefix + '_skysub_ABcombined.fits'

	endif


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
	dithering_length = floor(diagnostics[w_B[0]].position) - floor(diagnostics[w_A[0]].position)

	; number of pixels along the spatial position
	Nx = n_elements(*fuel.slits[0].rough_skylambda)

	; calculate the top and bottom edges of the slits, in the middle of the detector
	slit_bottom = fltarr(n_elements(fuel.slits))
	for i_slit=0, n_elements(fuel.slits)-1 do slit_bottom[i_slit] = $
		poly( 0.5*Nx, fuel.slits[i_slit].bottom_poly)
	slit_top = slit_bottom + fuel.slits.height

	; make plot that shows slit pairing
	cgplot, [0], /nodata, xra=[0, 2.5], $
	 	yra= [ min(slit_bottom)-1.2*abs(dithering_length), max(slit_top)+1.2*abs(dithering_length)], $
		ytit='pixel position along the vertical direction', charsize=1

	; for each slit show their name and position on the detector
	for i_slit=0, n_elements(fuel.slits)-1 do begin
		cgplot, [1,1], [slit_bottom[i_slit], slit_top[i_slit]], /overplot, thick=3
		cgplot, [0,3], slit_bottom[i_slit] + [0,0], /overplot, linestyle=2
		cgplot, [0,3], slit_top[i_slit] + [0,0], /overplot, linestyle=2
		cgtext, 0.9, 0.5*(slit_top+slit_bottom)[i_slit], 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + $
			' - ' + 'pos A', charsize=1, alignment=1
	endfor

	; now shift them by the dithering length
	for i_slit=0, n_elements(fuel.slits)-1 do begin
		cgplot, [1.2,1.2], [slit_bottom[i_slit], slit_top[i_slit]] + dithering_length, $
		 	/overplot, thick=3, color='red'
		cgtext, 1.3, 0.5*(slit_top+slit_bottom)[i_slit] + dithering_length, 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + $
			' - ' + 'pos B', charsize=1, alignment=0, color='red'
	endfor


	; calculate the quality of each possible overlap and store it in a matrix
	overlap_coefficient = fltarr(n_elements(fuel.slits), n_elements(fuel.slits))

	; loop through the slits and fill in the matrix
	for i_slit=0, n_elements(fuel.slits)-1 do $
		for j_slit=0, n_elements(fuel.slits)-1 do begin

			; calculate the top and bottom edges of the possible pairing
			top_edges = [slit_top[i_slit], slit_top[j_slit] + dithering_length]
			bottom_edges = [slit_bottom[i_slit], slit_bottom[j_slit] + dithering_length]

			; "overlap coefficient":
			; 1: perfect overlap
			; 0.01: tiny overlap
			; less than 0: no overlap
			; -1: infinitely far apart
			overlap_coefficient[i_slit,j_slit] = ( min(top_edges)-max(bottom_edges) ) / ( max(top_edges)-min(bottom_edges) )

	endfor

	; check whether there are any overlaps:
	if total( overlap_coefficient GT 0.0 ) EQ 0 then begin
		print, 'WARNING: no slits could be paired'
		return
	endif

	; we don't care about negative coefficients (meaning no overlap)
	overlap_coefficient = overlap_coefficient > 0.0

	print, ''
	print, 'Overlap matrix:'
	print, overlap_coefficient
	print, 'Possible overlaps found: ', n_elements( where(overlap_coefficient GT 0.0) )

	; let's keep track of the slits that have been paired
	slit_paired = bytarr(n_elements(fuel.slits))

	; go through all possible pairings, in decreasing order of overlap coefficient
	; stop when the pairing involves a slit that has already been paired
	for pairing_number=1, n_elements(where(overlap_coefficient GT 0.0)) do begin

		; find the highest overlap coefficient
		top_overlap = max(overlap_coefficient, ind2d)
	 	i_slit = (array_indices(overlap_coefficient, ind2d))[0]
		j_slit = (array_indices(overlap_coefficient, ind2d))[1]

		; remove this overlap
		overlap_coefficient[ind2d] = 0.0

		; check if these slits have already been paired - in that case we are done
		if slit_paired[i_slit] or slit_paired[j_slit] then break

		; mark these slits as paired
		slit_paired[i_slit] = 1
		slit_paired[j_slit] = 1

		print, ''
		print, 'Pairing #' + strtrim(pairing_number)
		print, 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + $
			' - ' +	'slit' + string(fuel.slits[j_slit].number, format='(I02)')
		print, 'Overlap: ' + cgnumber_formatter(top_overlap*100.0, decimals=2) + ' %'

		; check that both slits have actually been reduced
		if fuel.slits[i_slit].skip or fuel.slits[j_slit].skip then begin
			print, 'no combination; slits have been skipped'
			continue
		endif

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
		fuel.slits[j_slit].output_file = outname

	endfor


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_combine_mask, fuel
	;
	; combine the flux of all reduced slits into one large FITS file and do the same for the SNR map
	; NB: data will be resampled, this output should not be used for science
	;

	; select all the slits that were reduced
	w_slits = where(fuel.slits.skip eq 0, /null)
	if w_slits EQ !NULL then return
	slits = fuel.slits[w_slits]

	; sort the slits by y position on the mask
	slits = slits[sort(slits.approx_top + slits.approx_bottom)]

	; empty arrays with the lambda scale of all frames
	lambda_min = []
	lambda_step = []
	N_x = []

	; read all headers and get the wavelength scales
	for i_slit=0, n_elements(slits)-1 do begin
		header = headfits(slits[i_slit].output_file)
		lambda_min = [lambda_min, sxpar(header, 'CRVAL1')]
		lambda_step = [lambda_step, sxpar(header, 'CDELT1')]
		N_x = [N_x, sxpar(header, 'NAXIS1')]
	endfor

	; calculate the lambda_max for each slit
	lambda_max = lambda_min + (N_x-1)*lambda_step

	; get the wavelength properties of the output grid
	lambda_0 = min([lambda_min])
	lambda_1 = max([lambda_max])
	lambda_delta = median([lambda_step])

	; make wavelength grid for the output file
	new_lambda_axis = lambda_0 + lambda_delta * dindgen((lambda_1-lambda_0)/lambda_delta+1)

	; make the one-pixel row to insert between slits
	separation_row = dblarr(n_elements(new_lambda_axis)) + !values.d_nan

	; empty final maps
	flux_map = separation_row
	snr_map = separation_row

	; loop through the slits
	for i_slit=0, n_elements(slits)-1 do begin

		; read in the data and the error and calculate snr
		flux = mrdfits(slits[i_slit].output_file, 0, /silent)
		error = mrdfits(slits[i_slit].output_file, 1, /silent)
		snr = flux/error

		; make the original wavelength grid
		lambda_axis = lambda_min[i_slit] + lambda_step[i_slit] * dindgen( (size(flux))[1] )

		; set edges to NaNs to avoid crazy extrapolations
		flux[0:2,*] = !values.d_nan
		flux[-3:-1,*] = !values.d_nan
		snr[0:2,*] = !values.d_nan
		snr[-3:-1,*] = !values.d_nan

		; resample on new lambda grid
		flux_resampled = dblarr( n_elements(new_lambda_axis), (size(flux))[2] )
		snr_resampled = flux_resampled
		for i_row=0, (size(flux))[2]-1 do flux_resampled[*,i_row] = interpol(flux[*,i_row], lambda_axis, new_lambda_axis )
		for i_row=0, (size(snr))[2]-1 do snr_resampled[*,i_row] = interpol(snr[*,i_row], lambda_axis, new_lambda_axis )

		; add to the final map
		flux_map = [ [flux_map], [flux_resampled], [separation_row] ]
		snr_map = [ [snr_map], [snr_resampled], [separation_row] ]

	endfor

	; add wavelength calibration to the header
	SXADDPAR, header, 'CRPIX1', 1
	SXADDPAR, header, 'CRVAL1', lambda_0
	SXADDPAR, header, 'CDELT1', lambda_delta

	; output flux map
	output_filename = fuel.util.output_dir + 'combined_flux.fits'
	writefits, output_filename, flux_map, header
	print, output_filename, ' written'

	; output SNR map
	output_filename = fuel.util.output_dir + 'combined_SNR.fits'
	writefits, output_filename, snr_map, header
	print, output_filename, ' written'


END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_combine, fuel

		flame_util_module_start, fuel, 'flame_combine'


 	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		if fuel.slits[i_slit].skip then continue

		print, ''
		print, 'Combining slit ' + strtrim(fuel.slits[i_slit].number, 2) + ' - ' + fuel.slits[i_slit].name

		; handle errors by ignoring that slit
		if fuel.settings.stop_on_error eq 0 then begin
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

	; combine the SNR map of all slits that were reduced into one large FITS file
 	flame_combine_mask, fuel


	flame_util_module_end, fuel

END
