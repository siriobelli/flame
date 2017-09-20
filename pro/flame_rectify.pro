


PRO flame_rectify_one, fuel, filename, lambda_coeff=lambda_coeff, gamma_coeff=gamma_coeff, output_name=output_name, slit=slit

	print, 'rectifying ', filename

	; check whether there is an extension (containing the error spectrum)
	rdfits_struct, filename, struct, /silent, /header_only
	Next = n_tags(struct)

	; read in file to rectify
	im = mrdfits(filename, 0, header, /silent)
	if Next GT 1 then im_sigma = mrdfits(filename, 1, /silent)

	; read dimensions of the observed frame
	N_imx = (size(im))[1]
	N_imy = (size(im))[2]

	; get the parameters for the output grid
	lambda_0 = slit.outlambda_min
	delta_lambda = slit.outlambda_delta
	Nx = slit.outlambda_Npix

	; create 2D arrays containing the observed coordinates of each pixel
	x_2d = indgen(N_imx) # replicate(1, N_imy)
	y_2d = replicate(1, N_imx) # indgen(N_imy)

	; create 2D arrays containing the rectified coordinates of each pixel
	lambda_2d = flame_util_transform_coord(x_2d, y_2d, lambda_coeff )
	gamma_2d = flame_util_transform_coord(x_2d, y_2d, gamma_coeff )

	; define grid on the gamma axis - note that the grid points are integer numbers
	gamma_min = floor( min(gamma_2d, /nan) )
	gamma_max = floor( max(gamma_2d, /nan)+0.5 )
	Ny = gamma_max - gamma_min

	; calculate the *absolute* y coordinate (i.e. in the raw frame) corresponding to gamma=gamma_min at x=0
	; if gamma is not linear in y then throw an error
	if gamma_coeff[1,0] eq 1.0 and (size(gamma_coeff))[1] eq 2 then $
		abs_y_gamma_min = gamma_min - gamma_coeff[0,0] + slit.yrange_cutout[0] $
		else message, 'Non-linear vertical rectification not supported'

	; normalize the lambda values (otherwise triangulate does not work well; maybe because the scale of x and y is too different)
	lambdax_2d = (lambda_2d-lambda_0) / delta_lambda

	; make sure that the interpolation method is a valid string
	if total(strlowcase(fuel.settings.interpolation_method) eq strlowcase(['Linear', 'NearestNeighbor', $
		'NaturalNeighbor', 'Quintic']) ) NE 1 then begin
		print, 'WARNING: interpolation_method not valid: ' + fuel.settings.interpolation_method
		print, 'using NaturalNeighbor instead'
		fuel.settings.interpolation_method = 'NaturalNeighbor'
	endif

	; resample image onto new grid using griddata
	triangulate, lambdax_2d, gamma_2d, triangles
	new_im = griddata(lambdax_2d, gamma_2d, im, triangles=triangles, start=[0.0, gamma_min], delta=[1.0, 1.0], dimension=[Nx, Ny], $
		method=fuel.settings.interpolation_method, missing=!values.d_nan)
	if Next GT 1 then new_im_sigma = griddata(lambdax_2d, gamma_2d, im_sigma, triangles=triangles, start=[0.0, gamma_min], delta=[1.0, 1.0], dimension=[Nx, Ny], $
		method=fuel.settings.interpolation_method, missing=!values.d_nan)

	; add the wavelength calibration to the FITS header
	SXADDPAR, Header, 'CTYPE1', 'AWAV    '
	SXADDPAR, Header, 'CUNIT1', 'MICRON'
	SXADDPAR, Header, 'CRPIX1', 1
	SXADDPAR, Header, 'CRVAL1', lambda_0
	SXADDPAR, Header, 'CDELT1', delta_lambda

	; add the spatial position to the FITS header
	SXADDPAR, Header, 'CUNIT2', 'PIXEL'
	SXADDPAR, Header, 'CRPIX2', 1
	SXADDPAR, Header, 'CRVAL2', gamma_min
	SXADDPAR, Header, 'CDELT2', 1.0
	SXADDPAR, Header, 'YCUTOUT', abs_y_gamma_min, 'Y coordinate of the first pixel; added by FLAME'

	; delete WCS keywords
	SXDELPAR, Header, 'CTYPE2'
	SXDELPAR, Header, 'CD1_1'
	SXDELPAR, Header, 'CD1_2'
	SXDELPAR, Header, 'CD2_1'
	SXDELPAR, Header, 'CD2_2'

	; write rectified image
	writefits, output_name, new_im, header
	if Next GT 1 then writefits, output_name, new_im_sigma, /append

END


; ---------------------------------------------------------------------------------------------------------------------------


PRO flame_rectify, fuel

	flame_util_module_start, fuel, 'flame_rectify'


	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		if fuel.slits[i_slit].skip then continue

		this_slit = fuel.slits[i_slit]

		print, 'Rectifying slit ', this_slit.number, ' - ', this_slit.name

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

		for i_frame=0, n_elements(this_slit.cutouts)-1 do begin

			this_cutout = this_slit.cutouts[i_frame]

			filename = this_cutout.filename
			if this_cutout.illcorr_applied then $
				filename = flame_util_replace_string(this_cutout.filename, '_corr', '_illcorr')

			; rectify observed frame
			flame_rectify_one, fuel, filename, lambda_coeff=*this_cutout.lambda_coeff, gamma_coeff=*this_cutout.gamma_coeff, $
				output_name = flame_util_replace_string(filename, '.fits', '_rectified.fits'), slit=this_slit

			; rectify sky model
			flame_rectify_one, fuel, flame_util_replace_string(filename, '.fits', '_skymodel.fits'), lambda_coeff=*this_cutout.lambda_coeff, gamma_coeff=*this_cutout.gamma_coeff, $
				output_name = flame_util_replace_string(filename, '.fits', '_skymodel_rectified.fits'), slit=this_slit

			; rectify sky-subtracted frame
			flame_rectify_one, fuel, flame_util_replace_string(filename, '.fits', '_skysub.fits'), lambda_coeff=*this_cutout.lambda_coeff, gamma_coeff=*this_cutout.gamma_coeff, $
				output_name = flame_util_replace_string(filename, '.fits', '_skysub_rectified.fits'), slit=this_slit

		endfor

	endfor


  flame_util_module_end, fuel

END
