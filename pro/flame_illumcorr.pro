


PRO flame_illumcorr_one, fuel=fuel, i_slit=i_slit, i_frame=i_frame
	;
	; use the speclines to derive and apply an illumination correction
	; along the spatial slit axis
	;

	cutout = fuel.slits[i_slit].cutouts[i_frame]
	speclines = *cutout.speclines

	; read in slit
	im = mrdfits(cutout.filename, 0, hdr, /silent)
	im_sigma = mrdfits(cutout.filename, 1, /silent)

	; cutout dimensions
	N_pixel_x = (size(im))[1]
	N_pixel_y = (size(im))[2]

	; sort by wavelength
	sorted_lambdas = speclines[ sort(speclines.lambda) ].lambda

	; find unique wavelengths
	lambdas = sorted_lambdas[uniq(sorted_lambdas)]

	; calculate Gussian fluxes
	OHflux = sqrt(2.0*3.14) * speclines.peak * speclines.sigma

	; for each line, normalize flux to the median
	OHnorm = OHflux * 0.0
	for i_line=0, n_elements(lambdas)-1 do begin
		w_thisline = where(speclines.lambda eq lambdas[i_line], /null)
		OHnorm[w_thisline] = OHflux[w_thisline] / median(OHflux[w_thisline])
	endfor

	; calculate the gamma coordinate of each OHline detection
	OHlambda = flame_util_transform_coord(speclines.x, speclines.y, *cutout.lambda_coeff )
	OHgamma = flame_util_transform_coord(speclines.x, speclines.y, *cutout.gamma_coeff )

	; sort by gamma
	sorted_gamma = OHgamma[sort(OHgamma)]
	sorted_illum = OHnorm[sort(OHgamma)]

	; do not consider measurements that are more than a factor of three off
	w_tofit = where(sorted_illum GT 0.33 and sorted_illum LT 3.0, /null)

	; set the boundaries for a meaningful correction
	gamma_min = sorted_gamma[3]
	gamma_max = sorted_gamma[-4]

	; fit polynomial to the illumination correction as a function of gamma
	; NB: to increase robustness, set zero-point for gamma so that we are working with small numbers
	poly_coeff = robust_poly_fit(sorted_gamma[w_tofit]-gamma_min, sorted_illum[w_tofit], 8)

	; scatter plot of the illumination (show all OH lines)
	cgplot, sorted_gamma[w_tofit], sorted_illum[w_tofit], psym=3, /ynozero, charsize=1.2, $
		xtitle='gamma coordinate', ytitle='Illumination'

	; overplot the smooth illumination
	x_axis = gamma_min + (gamma_max-gamma_min)*dindgen(200)/199.
	cgplot, x_axis, poly(x_axis-gamma_min, poly_coeff), $
		/overplot, color='red', thick=3

	if median(poly(x_axis-gamma_min, poly_coeff)) GT 2.0 or $
		median(poly(x_axis-gamma_min, poly_coeff)) LT 0.5 then message, 'Illumination correction failed'

	; overplot flat illumination
	cgplot, [gamma_min - 0.5*gamma_max , gamma_max*1.5], [1,1], /overplot, linestyle=2, thick=3

	; overplot the limit to the correction (25%)
	cgplot, [gamma_min - 0.5*gamma_max , gamma_max*1.5], 0.75+[0,0], /overplot, linestyle=2, thick=1
	cgplot, [gamma_min - 0.5*gamma_max , gamma_max*1.5], 1.25+[0,0], /overplot, linestyle=2, thick=1

	; if we do not have to apply the illumination correction, then we are done
	if ~fuel.settings.illumination_correction then return

	; calculate the gamma coordinate for each observed pixel, row by row
	gamma_coordinate = im * 0.0
	for i_row=0, N_pixel_y-1 do gamma_coordinate[*,i_row] = $
		flame_util_transform_coord(dindgen(N_pixel_x), replicate(i_row, N_pixel_x), *cutout.gamma_coeff )

	; calculate the illumination correction at each pixel
	illumination_correction = poly(gamma_coordinate-gamma_min, poly_coeff)

	; set the correction to NaN when outside the boundary
	illumination_correction[where(gamma_coordinate LT gamma_min OR $
		gamma_coordinate GT gamma_max, /null)] = !values.d_NaN

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

	; update the flag
	fuel.slits[i_slit].cutouts[i_frame].illcorr_applied = 1


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
		if fuel.settings.debugging eq 0 then begin
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


		if fuel.util.arc.n_frames GT 0 then begin

      print, 'Here you should derive the illumination correction from the arcs'

    endif else $
		  for i_frame=0, n_elements(fuel.slits[i_slit].cutouts)-1 do begin

				print, ''

	      ; calculate and apply the illumination correction
		    flame_illumcorr_one, fuel=fuel, i_slit=i_slit, i_frame=i_frame

      endfor

  endfor



  flame_util_module_end, fuel

END
