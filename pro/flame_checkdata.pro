
PRO flame_checkdata_refstar, fuel

	; check if the reference star has been specified
	if fuel.input.star_y_A eq 0.0 then return

	; x coordinate where the star trace is certainly visible
	star_x = mean(fuel.util.star_x_range)

	; identify the slit with the reference star
	i_ref = -1
	for i_slit=0, n_elements(fuel.slits)-1 do $
		if fuel.input.star_y_A GE poly(star_x, fuel.slits[i_slit].bottom_poly) $
			and fuel.input.star_y_A LE poly(star_x, fuel.slits[i_slit].bottom_poly) + $
			fuel.slits[i_slit].height then $
				i_ref = i_slit

	; if there is no slit with the reference star, then exit
	if i_ref eq -1 then begin
		print, 'Did not find the slit with the reference star.'
		return
	endif

	print, 'Reference star is in slit number ', strtrim(fuel.slits[i_ref].number, 2)
  cgPS_open, fuel.input.output_dir + 'reference_star.ps', /nomatch


	; load output spectrum of reference star
	; --------------------------------------

	; load output file
	ref_spec = mrdfits(fuel.slits[i_ref].output_file, 0, header, /silent)

	; get the wavelength calibration from the header
 	lambda_unit = strlowcase( strtrim(sxpar(header, 'CUNIT1'), 2) )
	lambda_axis = sxpar(header, 'CRVAL1') + sxpar(header,'CDELT1') * $
		( findgen(sxpar(header,'NAXIS1')) - sxpar(header,'CRPIX1') + 1d )


 	; median profile and fit
 	; -----------------------------

	; make the y-axis
	yaxis = dindgen( (size(ref_spec))[2] )

	; extract the spatial profile
	ref_profile = median(ref_spec, dimension=1)

	; identify the peak (should be easy for the reference star)
	est_peak = max(ref_profile, est_center)

	; fit a gaussian to the integrated profile
	wfit = where( finite(ref_profile) )
  fit_result = gaussfit(yaxis[wfit], ref_profile[wfit], ref_coeff, nterms=4, $
    estimates=[ est_peak, est_center, 3.0, 0.0], $
    chisq=chisq, sigma=coeff_err)

	; calculate seeing in arcsec
	median_seeing = 2.355 * ref_coeff[2] * fuel.instrument.pixel_scale
	print, ''
	print, 'The final effective seeing calculated from the reference star is ' + $
	 	cgnumber_formatter(median_seeing, decimals=2) + ' arcsec.'

	; plot the median profile
	cgplot, yaxis, ref_profile, psym=16, charsize=1, xtit='y pixel coordinate', $
		ytit='median flux', title='reference star: median profile (seeing = ' + $
		cgnumber_formatter(median_seeing, decimals=2) + ' arcsec)'

	; overplot Gaussian fit
	cgplot, dindgen(300)/299.0*yaxis[-1], $
		ref_coeff[0] * exp( -0.5* ( (dindgen(300)/299.0*yaxis[-1]-ref_coeff[1])/ref_coeff[2] )^2 ) + ref_coeff[3], $
		/overplot, color='red', thick=3


	; fit profile as a function of wavelength
	; ----------------------------------------

	; set the bin size, in pixels, along the wavelength direction
	binsize = 100
	starting_pixel = 0

	; total number of pixels along the wavelength direction
	N_pixel_x = (size(ref_spec))[1]

	; empty arrays for seeing and position measurement
	coord_x = []
	seeing = []
	center = []

	while starting_pixel LT N_pixel_x do begin

    ; extract the bin
    end_pixel = min([starting_pixel + binsize - 1, N_pixel_x-1])
    cutout_bin = ref_spec[starting_pixel : end_pixel, *]

    ; spatial profile
    profile = median(cutout_bin, dimension=1)

		; fit a Gaussian
		wfit = where( finite(profile) )
  	fit_result = gaussfit(yaxis[wfit], profile[wfit], coeff, nterms=4, $
    	estimates=ref_coeff, $
    	chisq=chisq, sigma=coeff_err)

		; save the result of the fit
		coord_x = [coord_x, 0.5*(starting_pixel+end_pixel)]
		seeing = [seeing, 2.355 * coeff[2] * fuel.instrument.pixel_scale]
		center = [center, coeff[1]]

    ; advance to next bin
    starting_pixel += binsize

  endwhile

	; plot seeing as a function of wavelength
	cgplot, coord_x, seeing, thick=4, charsize=1, $
		ytit='seeing (arcsec)', $
		yra=[min(seeing, /nan)-0.1, max(seeing, /nan)+0.1], $
		xra=[0, N_pixel_x], xsty=1+8, xthick=4, ythick=4, $
		position = [0.15, 0.55, 0.9, 0.9], xtickformat='(A1)'

	; top x-axis with pixel coordinate
	cgaxis, xaxis=1, xra = [ 0 , N_pixel_x ], xsty=1, charsize=1, $
		xtit='x-coordinate (pixel)'

	; show the median seeing
	cgplot, [0, N_pixel_x], [0,0]+median_seeing, /overplot, thick=3, linestyle=2

	; make center relative to initial position
	center -= center[0]

	; show centering as a function of wavelength
	cgplot, coord_x, center, thick=4, charsize=1, $
		xtit='x-coordinate (pixel)', ytit='center position (pixel)', $
		yra=[min(center, /nan)-1.0, max(center, /nan)+1.0], $
		xra=[0, N_pixel_x], xsty=1+4, xthick=4, ythick=4, $
		position = [0.15, 0.15, 0.9, 0.55], /noerase

	; bottom x-axis with wavelength
	cgaxis, xaxis=0, xra = [ lambda_axis[0] , lambda_axis[-1] ], xsty=1, charsize=1, $
		xtit='wavelength (' + lambda_unit + ')'

	; overplot zero line
	cgplot, [0, N_pixel_x], [0,0], /overplot, thick=3, linestyle=2


	; extract and plot 1D spectrum of reference star
	; ----------------------------------------------

	; extract spectrum from +/- 2 sigma around the center
	window_min = ( ref_coeff[1] - 2.0*ref_coeff[2] ) > 0
	window_max = ref_coeff[1] + 2.0*ref_coeff[2] < N_pixel_x-1

	; extract boxcar spectrum
	spectrum = total(ref_spec[ * , window_min:window_max ], 2, /nan)

	; show star spectrum
	cgplot, lambda_axis, smooth(spectrum, 17), charsize=1, $
		xtit='wavelength (' + lambda_unit + ')', ytit='flux', /ynozero, $
		title='boxcar extraction of the reference star spectrum'


	cgPS_close

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_checkdata_sky, fuel, i_slit=i_slit

	print, 'Checking slit number ' + strtrim(fuel.slits[i_slit].number, 2)


	; load the sky spectrum for this slit
	;-------------------------------------

	; filename of the output sky stack
	skystack_filename = fuel.input.output_dir + 'slit' + $
		string(fuel.slits[i_slit].number, format='(I02)') + '-' + fuel.slits[i_slit].name + '_stack_sky.fits'

	; load the sky spectrum
	sky_spec2d = mrdfits(skystack_filename, 0, header, /silent)

	; extract 1D spectrum
	sky_spec =  mean(sky_spec2d, dimension=2, /nan)

	; get the wavelength calibration from the header
 	lambda_unit = strlowcase( strtrim(sxpar(header, 'CUNIT1'), 2) )
	lambda_axis = sxpar(header, 'CRVAL1') + sxpar(header,'CDELT1') * $
		( findgen(sxpar(header,'NAXIS1')) - sxpar(header,'CRPIX1') + 1d )

	; for now we only support micron
	if lambda_unit ne 'micron' then message, lambda_unit + ' not supported!'


	; fit the OH lines
	;-------------------------------------

	; load line list
	readcol, fuel.util.linelist_filename, line_list, format='D', /silent

  ; calculate approximate sky line width
	linewidth_um = median(lambda_axis) / (2.36 * fuel.instrument.resolution_slit1arcsec)

	; identify the OH lines that are in this wavelength range
	w_lines = where(line_list GT min(lambda_axis, /nan) $
		AND line_list LT max(lambda_axis, /nan), /null )

	; make sure there are sky lines here
	if w_lines EQ !NULL then begin
    print, 'Warning: wavelength range does not contain sky lines'
    return
  endif

	; keep only the OH lines of interest
	line_list = line_list[w_lines]

	; arrays for the fit results
	sky_lambda_th = []
	sky_lambda_obs = []
	sky_sigma = []

	; fit a Gaussian to every sky line
	for i_line=0,n_elements(line_list)-1 do begin

		; select the region to fit
		w_fit = where( abs(lambda_axis-line_list[i_line]) LT 6.0*linewidth_um, /null )

		; check that the region is within the observed range
		if w_fit eq !NULL then continue

    ; check that there actually is signal and it's not just a bunch of NaNs
    if n_elements( where( finite(sky_spec[w_fit]), /null ) ) LE 5 then continue

		; error handling for the gaussian fitting
		catch, error_gaussfit
		if error_gaussfit ne 0 then begin
			print, 'GAUSSFIT ERROR STATUS: ' + strtrim(error_gaussfit,2)
			catch, /cancel
			continue
		endif

		; estimate parameters of the Gaussian
		est_peak = max( median( sky_spec[w_fit], 3) , /nan)
		est_center = line_list[i_line]
		est_sigma = linewidth_um
		est_cont = min( median( sky_spec[w_fit], 3) , /nan)

		; Gaussian fit
		!NULL = gaussfit( lambda_axis[w_fit], sky_spec[w_fit], gauss_param, nterms=4, $
			estimates=[est_peak, est_center, est_sigma, est_cont], sigma=gauss_err, chisq=chisq )

		; check that chi square makes sense
		if ~finite(chisq) then continue

		; check that the peak of the Gaussian is positive
		if gauss_param[0] LT 0.0 then continue

		; check that the SNR is high
		if gauss_param[0] LT 5.0*gauss_err[0] then continue

		; check that the center of the Guassian is in the observed range
		if gauss_param[1] LT min(lambda_axis[w_fit]) or gauss_param[1] GT max(lambda_axis[w_fit]) then continue

		; check that the Gaussian width makes sense
		if gauss_param[2] LT linewidth_um/10.0 or gauss_param[2] GT linewidth_um*10.0 then continue

		; save the results
		sky_lambda_th = [sky_lambda_th, line_list[i_line] ]
		sky_lambda_obs = [sky_lambda_obs, gauss_param[1] ]
		sky_sigma = [sky_sigma, gauss_param[2] ]

	endfor

	; calculate the wavelength residuals in angstrom
	residuals = 1d4 * (sky_lambda_th-sky_lambda_obs)

	; calculate the spectral resolution R
	spectral_R = sky_lambda_th/(sky_sigma*2.36)

	; plot the result of the fit
	;-------------------------------------

	; x axis range
	xra=[lambda_axis[0], lambda_axis[-1]]


	; panel 1: plot the spectrum
	cgplot, lambda_axis, sky_spec, charsize=0.8, xsty=1, xtit='', ytit='sky flux', $
		title = skystack_filename, position = [0.10, 0.70, 0.95, 0.95], $
		xtickformat="(A1)", xra=xra, /nodata

	; show the OH lines that were identified
	for i_line=0, n_elements(sky_lambda_th)-1 do $
		cgplot, sky_lambda_th[i_line] + [0,0], [-2,2]*max(abs(sky_spec)), /overplot, color='red'

	; show the spectrum on top, for clarity
	cgplot, lambda_axis, sky_spec, /overplot


	; panel 2: show the residuals
	cgplot, sky_lambda_th, residuals, /ynozero, xra=xra, $
		xsty=1, psym=16, color='red', symsize=0.7, $
		ytit='residuals (angstrom)', charsize=0.8, $
		/noerase, position = [0.10, 0.45, 0.95, 0.70], xtickformat="(A1)"
	cgplot, [xra[0], xra[1]], [0,0], /overplot, thick=3, linestyle=2


	; panel 3: plot the spectral resolution
	cgplot, sky_lambda_th, spectral_R, /ynozero, xra=xra, $
		xsty=1, psym=16, color='red', symsize=0.7, $
		xtit='wavelength (micron)', ytit='spectral resolution R', charsize=0.8, $
		/noerase, position = [0.10, 0.20, 0.95, 0.45]
	cgplot, [xra[0], xra[1]], [0,0]+median(spectral_R), /overplot, thick=3, linestyle=2


	; print some stats on wavelength calibration
	cgtext, 0.10, 0.11, 'wavelength calibration residuals: ', /normal, charsize=0.7
	cgtext, 0.10, 0.08, 'standard deviation = ' + $
		cgnumber_formatter( stddev(residuals, /nan), decimals=3) + ' ' + STRING("305B), /normal, charsize=0.7
	cgtext, 0.10, 0.06, 'root mean square = ' + $
		cgnumber_formatter( sqrt( mean(residuals^2, /nan)), decimals=3) + ' ' + STRING("305B), /normal, charsize=0.7
	cgtext, 0.10, 0.04, 'median absolute deviation = ' + $
		cgnumber_formatter( median(abs(residuals)), decimals=3) + ' ' + STRING("305B), /normal, charsize=0.7


	; print some stats on spectral resolution
	cgtext, 0.50, 0.11, 'spectral resolution: ', /normal, charsize=0.7
	cgtext, 0.50, 0.08, 'median R = ' + $
		cgnumber_formatter( median(spectral_R), decimals=0), /normal, charsize=0.7
	cgtext, 0.50, 0.06, 'stddev R = ' + $
		cgnumber_formatter( stddev(spectral_R, /nan), decimals=0), /normal, charsize=0.7


	; print some stats on velocity resolution
	cgtext, 0.75, 0.11, 'median velocity resolution: ', /normal, charsize=0.7
	cgtext, 0.75, 0.08, 'FWHM = ' + $
		cgnumber_formatter( median(3d5/spectral_R), decimals=1) + ' km/s', /normal, charsize=0.7
	cgtext, 0.75, 0.06, 'sigma = ' + $
		cgnumber_formatter( median(3d5/spectral_R)/2.36, decimals=1) + ' km/s', /normal, charsize=0.7


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_checkdata_speclines, fuel, i_slit=i_slit

return

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_checkdata, fuel

	flame_util_module_start, fuel, 'flame_checkdata'

	; calculate diagnostics from reference star
	flame_checkdata_refstar, fuel

	; calculate diagnostics for each slit
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		cgPS_open, fuel.input.output_dir + 'slit' + string(fuel.slits[i_slit].number, format='(I02)') + $
			'-' + fuel.slits[i_slit].name +  '_datacheck.ps', /nomatch

		; calculate diagnostics from the sky spectrum
		flame_checkdata_sky, fuel, i_slit=i_slit

		; show result of wavelength calibration on individual frames
		flame_checkdata_speclines, fuel, i_slit=i_slit

		cgPS_close

	endfor


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
