


PRO flame_illumcorr_applycorrection, cutout, illumcorr, gamma_0
	;
	; apply the illumination correction to a cutout and save new file
	;

	; read in slit
	im = mrdfits(cutout.filename, 0, hdr, /silent)
	im_sigma = mrdfits(cutout.filename, 1, /silent)

	; cutout dimensions
	N_pixel_x = (size(im))[1]
	N_pixel_y = (size(im))[2]

	; calculate the gamma coordinate for each observed pixel, row by row
	gamma_coordinate = im * 0.0
	for i_row=0, N_pixel_y-1 do gamma_coordinate[*,i_row] = $
		flame_util_transform_coord(dindgen(N_pixel_x), replicate(i_row, N_pixel_x), *cutout.gamma_coeff )

	; calculate the illumination correction at each pixel
	illumination_correction = poly(gamma_coordinate-gamma_0, illumcorr)

	; set the correction to NaN when outside the boundary
	illumination_correction[where(gamma_coordinate LT min(gamma_coordinate, /nan) OR $
		gamma_coordinate GT max(gamma_coordinate, /nan), /null)] = !values.d_NaN

	; set the correction to NaN if it's more than 25%
	illumination_correction[where(illumination_correction GT 1.25 or $
		illumination_correction LT 0.75, /null)] = !values.d_NaN

	; apply illumination correction
	im /= illumination_correction
	im_sigma /= illumination_correction

	; filename for the illumination-corrected cutout
	illcorr_filename = flame_util_replace_string(cutout.filename, '_corr', '_illcorr')

	; write out the illumination-corrected cutout
  writefits, illcorr_filename, im, hdr
	writefits, illcorr_filename, im_sigma, /append
	print, illcorr_filename + ' written'

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_illumcorr_getcorrection_emlines, cutout, gamma_0=gamma_0
	;
	; use the speclines to derive an illumination correction
	; along the spatial slit axis
	; returns the polynomial coefficients and the zero-point gamma_0
	; such that the illumination field is:
	; I(gamma) = poly(gamma-gamma_0, poly_coeff)
	;

	; speclines from sky or arcs
	speclines = *cutout.speclines

	; sort by wavelength
	sorted_lambdas = speclines[ sort(speclines.lambda) ].lambda

	; find unique wavelengths
	lambdas = sorted_lambdas[uniq(sorted_lambdas)]

	; calculate Gussian fluxes
	line_flux = sqrt(2.0*3.14) * speclines.peak * speclines.sigma

	; for each line, normalize flux to the median
	line_norm = line_flux * 0.0
	for i_line=0, n_elements(lambdas)-1 do begin
		w_thisline = where(speclines.lambda eq lambdas[i_line], /null)
		line_norm[w_thisline] = line_flux[w_thisline] / median(line_flux[w_thisline])
	endfor

	; calculate the gamma coordinate of each OHline detection
	line_lambda = flame_util_transform_coord(speclines.x, speclines.y, *cutout.lambda_coeff )
	line_gamma = flame_util_transform_coord(speclines.x, speclines.y, *cutout.gamma_coeff )

	; sort by gamma
	sorted_gamma = line_gamma[sort(line_gamma)]
	sorted_illum = line_norm[sort(line_gamma)]

	; do not consider measurements that are more than a factor of three off
	w_tofit = where(sorted_illum GT 0.33 and sorted_illum LT 3.0, /null)

	; set the boundaries for a meaningful correction
	gamma_min = sorted_gamma[3]
	gamma_max = sorted_gamma[-4]

	; median filter the illumination values
	sorted_gamma_filtered = median(sorted_gamma[w_tofit], 5)

	; fit polynomial to the illumination correction as a function of gamma
	; NB: to increase robustness, set zero-point for gamma so that we are working with small numbers
	poly_coeff = robust_poly_fit(sorted_gamma_filtered-gamma_min, sorted_illum[w_tofit], 8)

	; scatter plot of the illumination (show all lines)
	cgps_open, flame_util_replace_string(cutout.filename, '.fits', '_illumcorr.ps'), /nomatch
	cgplot, sorted_gamma[w_tofit], sorted_illum[w_tofit], psym=3, /ynozero, charsize=1.2, $
		xtitle='gamma coordinate', ytitle='Illumination'

	; overplot the smooth illumination
	x_axis = gamma_min + (gamma_max-gamma_min)*dindgen(200)/199.
	cgplot, x_axis, poly(x_axis-gamma_min, poly_coeff), $
		/overplot, color='red', thick=3

	if median(poly(x_axis-gamma_min, poly_coeff)) GT 2.0 or $
		median(poly(x_axis-gamma_min, poly_coeff)) LT 0.5 then message, 'Illumination correction failed'

	; overplot flat illumination
	margin = abs(gamma_max-gamma_min)
	cgplot, [gamma_min - margin, gamma_max + margin], [1,1], /overplot, linestyle=2, thick=3

	; overplot the limit to the correction (25%)
	cgplot, [gamma_min - margin, gamma_max + margin], 0.75+[0,0], /overplot, linestyle=2, thick=1
	cgplot, [gamma_min - margin, gamma_max + margin], 1.25+[0,0], /overplot, linestyle=2, thick=1
	cgps_close

	; output the zero-point for the gamma axis
	gamma_0 = gamma_min

	; return the coefficients describing the illumination correction as a function of gamma
	return, poly_coeff

END




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_illumcorr_getcorrection_flat, cutout, slit=slit, gamma_0=gamma_0
	;
	; use the illumflat to derive an illumination correction
	; along the spatial slit axis
	; returns the polynomial coefficients and the zero-point gamma_0
	; such that the illumination field is:
	; I(gamma) = poly(gamma-gamma_0, poly_coeff)
	;

	; read the illumination flat cutout
	illumflat = mrdfits(cutout.filename, 0, hdr, /silent)

	; read dimensions of the observed frame
	N_imx = (size(illumflat))[1]
	N_imy = (size(illumflat))[2]

	; create 2D arrays containing the observed coordinates of each pixel
	x_2d = indgen(N_imx) # replicate(1, N_imy)
	y_2d = replicate(1, N_imx) # indgen(N_imy)

	; create 2D arrays containing the rectified coordinates of each pixel
	lambda_2d = flame_util_transform_coord(x_2d, y_2d, *cutout.lambda_coeff )
	gamma_2d = flame_util_transform_coord(x_2d, y_2d, *cutout.gamma_coeff )

	; get the parameters for the output grid
	lambda_0 = slit.outlambda_min
	delta_lambda = slit.outlambda_delta
	Nx = slit.outlambda_Npix

	; normalize the lambda values (otherwise triangulate does not work well; maybe because the scale of x and y is too different)
	lambdax_2d = (lambda_2d-lambda_0) / delta_lambda

	; find the range of gamma values in the image
	gamma_min = round( min(gamma_2d, /nan) )
	gamma_max = round( max(gamma_2d, /nan) )
	Ny = gamma_max-gamma_min

	; resample image onto new grid using griddata
	triangulate, lambdax_2d, gamma_2d, triangles
	rectified_illumflat = griddata(lambdax_2d, gamma_2d, illumflat, triangles=triangles, start=[0.0, gamma_min], delta=[1.0, 1.0], dimension=[Nx, Ny], /linear, missing=!values.d_nan)

	; get the median spectrum of the flat
	spectrum = median(rectified_illumflat, dimension=2)

	; make normalization image
	norm = spectrum # replicate(1, Ny)

	; normalize flat
	normalized_flat = rectified_illumflat / norm

	; get the median profile as a function of gamma
	profile_flat = median(normalized_flat, dimension=1)

	; re-normalized profile
	profile_flat /= median(profile_flat)

	; do not consider measurements that are more than a factor of three off
	w_tofit = where(profile_flat GT 0.33 and profile_flat LT 3.0, /null)

	; make gamma grid
	gamma_grid = gamma_min + dindgen(Ny)

	; fit polynomial to the illumination correction as a function of gamma
	; NB: to increase robustness, set zero-point for gamma so that we are working with small numbers
	poly_coeff = robust_poly_fit(gamma_grid[w_tofit]-gamma_min, profile_flat[w_tofit], 8)

	; plot of the illumination
	cgps_open, flame_util_replace_string(cutout.filename, '.fits', '_illumcorr.ps'), /nomatch
	cgplot, gamma_grid[w_tofit], profile_flat[w_tofit], psym=3, /ynozero, charsize=1.2, $
		xtitle='gamma coordinate', ytitle='Illumination'

	; overplot the smooth illumination
	x_axis = gamma_min + (gamma_max-gamma_min)*dindgen(200)/199.
	cgplot, x_axis, poly(x_axis-gamma_min, poly_coeff), $
		/overplot, color='red', thick=3

	if median(poly(x_axis-gamma_min, poly_coeff)) GT 2.0 or $
		median(poly(x_axis-gamma_min, poly_coeff)) LT 0.5 then message, 'Illumination correction failed'

	; overplot flat illumination
	margin = abs(gamma_max-gamma_min)
	cgplot, [gamma_min - margin, gamma_max + margin], [1,1], /overplot, linestyle=2, thick=3

	; overplot the limit to the correction (25%)
	cgplot, [gamma_min - margin, gamma_max + margin], 0.75+[0,0], /overplot, linestyle=2, thick=1
	cgplot, [gamma_min - margin, gamma_max + margin], 1.25+[0,0], /overplot, linestyle=2, thick=1
	cgps_close

	; output the zero-point for the gamma axis
	gamma_0 = gamma_min

	; return the coefficients describing the illumination correction as a function of gamma
	return, poly_coeff

END




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_illumcorr, fuel

	flame_util_module_start, fuel, 'flame_illumcorr'


	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		; this slit
		this_slit = fuel.slits[i_slit]

		if this_slit.skip then continue

	  print, 'Illumination correction for slit ', strtrim(this_slit.number,2), ' - ', this_slit.name
		print, ' '

		; handle errors by ignoring that slit
		if fuel.settings.stop_on_error eq 0 then begin
			catch, error_status
			if error_status ne 0 then begin
				print, ''
		    print, '**************************'
		    print, '***       WARNING      ***'
		    print, '**************************'
		    print, 'Error found. Skipping slit ' + strtrim(this_slit.number,2), ' - ', this_slit.name
				fuel.slits[i_slit].skip = 1
				catch, /cancel
				continue
			endif
		endif


		; if needed, use illumflat to get illumination correction
		if fuel.util.illumflat.n_frames GT 0 then begin

      print, 'Calculating illumination from the flat'
			illumcorr = flame_illumcorr_getcorrection_flat( fuel.slits[i_slit].illumflat_cutout, slit=this_slit, gamma_0=gamma_0)

		; if needed, use arcs to get illumination correction
		endif else $
		if fuel.util.arc.n_frames GT 0 then begin

      print, 'Calculating illumination from the arcs'
			illumcorr = flame_illumcorr_getcorrection_emlines( fuel.slits[i_slit].arc_cutout, gamma_0=gamma_0)

    endif else print, 'Calculating illumination from the science frames'

		; loop through the frames
	  for i_frame=0, n_elements(fuel.slits[i_slit].cutouts)-1 do begin

	      ; calculate the illumination correction from the science frame (unless arcs are used)
				if fuel.util.arc.n_frames EQ 0 then $
					illumcorr = flame_illumcorr_getcorrection_emlines( fuel.slits[i_slit].cutouts[i_frame], gamma_0=gamma_0)

				; if we are not applying the illumination correction, then skip remaining part
				if ~fuel.settings.illumination_correction then continue

				; apply the illumination correction
				flame_illumcorr_applycorrection, fuel.slits[i_slit].cutouts[i_frame], illumcorr, gamma_0

				; update the flag
				fuel.slits[i_slit].cutouts[i_frame].illcorr_applied = 1

    endfor

  endfor


  flame_util_module_end, fuel

END
