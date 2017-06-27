
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


PRO flame_checkdata, fuel

	flame_util_module_start, fuel, 'flame_checkdata'

	; calculate diagnostics from reference star
	flame_checkdata_refstar, fuel


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
