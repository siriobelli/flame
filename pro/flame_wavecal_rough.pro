;
; First pass at a rough wavelength calibration.
; For each slit, load only the first frame, extract the spectrum
; of the central few pixel rows, and cross-correlate that
; with a model spectrum of the sky until you find a reasonable match
;


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



FUNCTION flame_generate_lambda_axis, $
	N_pixels=N_pixels, lambda_central=lambda_central, delta_lambda=delta_lambda, delta_delta_lambda=delta_delta_lambda

	; returns a non-uniform wavelength axis centered on lambda_central, with pixel scale at the central pixel
	; given by delta_lambda, and its fractional variation is given by delta_delta_lambda as in:
	; lambda[x_c+1] - lambda[x_c] = delta_lambda
	; lambda[x_c+2] - lambda[x_c+1] = (1 + delta_delta_lambda) * ( lambda[x_c+1] - lambda[x_c] )

	if delta_delta_lambda eq 0. then $
		 lambda_axis = lambda_central + delta_lambda * (dindgen(N_pixels) - double(N_pixels)/2.) $
	 else $
		 lambda_axis = lambda_central + delta_lambda / alog(1. + delta_delta_lambda) * ( (1. + delta_delta_lambda)^(dindgen(N_pixels) - double(N_pixels)/2.) - 1. )

	return, lambda_axis

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_crosscorr, sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
	 approx_lambda_central=lambda_central, pix_scale_grid=pix_scale_grid, pix_scale_variation_grid=pix_scale_variation_grid, $
	 lambda_axis=lambda_axis_output, title=title

	;
	; Estimate the wavelength solution by comparing the observed sky spectrum with a model sky spectrum
	; Loop through a grid of pixel scale values and at each step cross-correlate the spectra to find the lambda shift
	; The wavelength solution is parameterized as an exponential function of the pixel coordinate.
	; The free parameters are central lambda, pixel scale at central pixel, and relative variation in
	; pixel scale between adjacent pixels
	; Return the wavelength axis.
	;

	; if the pixel scale variation grid is not specified, then fix pix_scale_variation to 0 and using a uniform lambda axis
	if ~keyword_set(pix_scale_variation_grid) then pix_scale_variation_grid = 0.

	; wavelength shift values to probe
	lag = -n_elements(sky)/2 + indgen(n_elements(sky))	; in pixels

	crosscorr_array = dblarr(n_elements(pix_scale_grid), n_elements(pix_scale_variation_grid))
	delta_array = intarr(n_elements(pix_scale_grid), n_elements(pix_scale_variation_grid))

	print, 'Finding wavelength solution via cross-correlation...'

	for i=0,n_elements(pix_scale_grid)-1 do for j=0,n_elements(pix_scale_variation_grid)-1 do begin

		print, strjoin(replicate(string(8B),10)) + strtrim( round( 100. * float(i)/float(n_elements(pix_scale_grid)-1) ) , 2) + '%' , format='(a,$)'

		; this is the new wavelength axis, with the chosen pixel scale and pixel scale variation
		lambda_axis = flame_generate_lambda_axis(N_pixels=n_elements(sky), lambda_central=lambda_central, delta_lambda=pix_scale_grid[i], delta_delta_lambda=pix_scale_variation_grid[j])

		; need to interpolate the model onto the new wavelength axis
		new_model_flux = interpol( model_flux, model_lambda, lambda_axis)

		; find the best wavelength shift
		crosscorr = c_correlate(sky, new_model_flux, lag)
		max_cc = max(crosscorr, ind, /nan)
		crosscorr_array[i,j] = max_cc
		delta = lag[ind]
		delta_array[i,j] = delta

	endfor

	; find the best point in the grid, the one with the highest value of cross-correlation
	max_cc = max(crosscorr_array, ind)

	; now find the best-fit parameters

	; wavelength shift in pixels
	best_delta = delta_array[ind]

	; pixel scale and its variation
	if size(crosscorr_array, /n_dimensions) eq 2 then begin
		ind_2d = array_indices(crosscorr_array, ind)
		best_pix_scale = pix_scale_grid[ ind_2d[0] ]
		best_pix_scale_variation = pix_scale_variation_grid[ ind_2d[1] ]
	endif else begin
		best_pix_scale = pix_scale_grid[ ind ]
		best_pix_scale_variation = pix_scale_variation_grid
	endelse

	print, ''
  	print, 'max crosscorr = ', max_cc
  	print, 'central wavelength = ', lambda_central
	print, 'delta = ', best_delta
  	print, 'pixel scale = ', best_pix_scale
	print, 'pixel scale variation = ', best_pix_scale_variation

	; generate the wavelength axis given the best-fit pixel scale and pixel scale variation
	; make an array longer than needed to account for shifting
	lambda_axis = flame_generate_lambda_axis(N_pixels=n_elements(sky) + 2*abs(best_delta), $
		lambda_central=lambda_central, delta_lambda=best_pix_scale, delta_delta_lambda=best_pix_scale_variation)

	; apply horizontal shift found via cross correlation
	lambda_axis = shift( lambda_axis, -best_delta )

	; trim the leftover part
	lambda_axis = lambda_axis[ abs(best_delta) : -abs(best_delta)-1 ]

	print, 'lambda axis range = ', lambda_axis[0], lambda_axis[-1]

	; interpolate model flux onto new wavelength axis
	new_model_flux = interpol( model_flux, model_lambda, lambda_axis)
	new_model_flux /= max(new_model_flux, /nan)

	cgplot, lambda_axis, sky-0.05, charsize=1, thick=3, xtit='wavelength (micron)', title=title
	cgplot, lambda_axis, new_model_flux, color='red', /overplot, thick=3
	cgtext, 0.75, 0.80, 'observed sky', charsize=1, /normal
	cgtext, 0.75, 0.75, 'model sky', charsize=1, /normal, color='red'

	; return the wavelength solution
	lambda_axis_output = lambda_axis

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_wavecal_approximate, fuel=fuel, this_slit=this_slit

	;
	; Two steps:
	; First, use a coarse, logarithmic grid with constant pixel scale (uniform wavelength solution)
	; Second, use a fine grid centered on the coarse value of pixel scale, and explore also the
	; pix_scale_variation parameter for non-uniform wavelength solution
	;

  ; take the first frame
  slit_filename = (*this_slit.filenames)[0]
  print, 'Using the central pixel rows of ', slit_filename

	; load the model sky spectrum
	readcol, fuel.util.sky_emission_filename, model_lambda, model_flux	; lambda in micron

	; select the range of interest for this particular slit
	w_slit = where( model_lambda GT this_slit.approx_wavelength_lo $
		and model_lambda LT this_slit.approx_wavelength_hi, /null )
	model_lambda = model_lambda[w_slit]
	model_flux = model_flux[w_slit]

	; normalize sky model
	model_flux /= max(model_flux, /nan)

	; for a tidy extrapolation, first and last elements of model should be zeros
	model_flux[0:5] = 0.
	model_flux[-5:*] = 0.


	; load the slit image
	;---------------------

	; read in slit
	im = readfits(slit_filename, hdr)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(im))[2]

	; sky spectrum: extract along the central 5 pixels
  ; NB: NEED TO MAKE SURE THAT THIS IS, INDEED, SKY, AND NOT THE OBJECT
	sky = median(im[*, N_spatial_pix/2-3 : N_spatial_pix/2+2], dimension=2)

	; get rid of NaNs, which create problems for the cross-correlation
	sky[where(~finite(sky), /null)] = 0

	; create the x-axis in pixel coordinates
	pix_axis = dindgen( n_elements(sky) )

	; trim edge to avoid problems
	sky[0:10] = 0.
	sky[-11:-1] = 0.

	; normalize the observed spectrum
	sky /= max(sky, /nan)

	; start PS file
  ps_filename = file_dirname(slit_filename, /mark_directory) + 'wavelength_solution_estimate.ps'
	cgPS_open, ps_filename, /nomatch


	print, ''
	print, 'FIRST STEP: rough estimate of lambda and delta_lambda'
	print, '-----------------------------------------------------'

	; estimated pixel scale from the approximate wavelength range
	estimated_pix_scale = (this_slit.approx_wavelength_hi-this_slit.approx_wavelength_lo) $
		/ double(n_elements(sky))

	; estimated central wavelength
	estimated_lambda = 0.5*(this_slit.approx_wavelength_hi+this_slit.approx_wavelength_lo)

	; pixel scale for the model
	model_pix_scale = abs( median(model_lambda - shift(model_lambda,-1)) )

	; let's smooth observed and model sky so that we can get an approximate match
	smoothing_sigma = estimated_lambda / (2.36 * fuel.input.rough_wavecal_R)
	sky_sm = gauss_smooth(sky, smoothing_sigma / estimated_pix_scale )
	model_flux_sm = gauss_smooth(model_flux, smoothing_sigma / model_pix_scale)

	; WARNING: THIS CAUSES A PROBLEM WITH LRIS, BECAUSE 100 pixels are not enough for a continuum-level smoothing!
	; ; remove continuum variation from the observed sky (important in the K band)
	; sky_sm -= median(sky, 100)

	; renormalize spectra
	sky_sm /= max(sky_sm)
	model_flux_sm /= max(model_flux_sm)

	; show spectra before cross-correlation
	cgplot, this_slit.approx_wavelength_lo + estimated_pix_scale * dindgen(n_elements(sky_sm)), sky_sm-0.05, $
		charsize=1, thick=3, xtit='wavelength (micron)', title=(strsplit(slit_filename,'/', /extract))[-1] + ' - initial guess'
	cgplot, model_lambda, model_flux_sm, color='red', /overplot, thick=3
	cgtext, 0.75, 0.80, 'observed sky', charsize=1, /normal
	cgtext, 0.75, 0.75, 'model sky', charsize=1, /normal, color='red'

	; there are two parameters that define the wavelength axis: central pixel scale and variation in the pixel scale
	; first, we assume a constant pixel scale, in micron per pixels:
	;coarse_pix_scale_grid = 10^(-4.5 + 2.5*dindgen(1000)/999.)
	coarse_pix_scale_grid = estimated_pix_scale * (10.0+dindgen(500))/500.0*2.0

	flame_wavecal_crosscorr, sky=sky_sm, model_lambda=model_lambda, model_flux=model_flux_sm, $
		 approx_lambda_central=median(model_lambda), pix_scale_grid=coarse_pix_scale_grid, $
		 lambda_axis=lambda_axis, title='after first cross-correlation'

	; extract central lambda and pixel scale
	coarse_lambda_central = median(lambda_axis)
	coarse_pixel_scale = abs( median(lambda_axis - shift(lambda_axis,-1)) )

	; estimate the error on the pixel scale by using the step size in the grid
	w_coarse_pixel_scale = ( where( coarse_pix_scale_grid GE coarse_pixel_scale, /null) )[0]
	if w_coarse_pixel_scale eq 0 or w_coarse_pixel_scale eq n_elements(coarse_pix_scale_grid)-1 then $
		message, 'Cross correlation selected values at the edge of the grid!'
	coarse_pixel_scale_error = coarse_pix_scale_grid[ w_coarse_pixel_scale + 1 ] - coarse_pix_scale_grid[ w_coarse_pixel_scale - 1 ]

	print, 'Results of first coarse fit:'
	print, 'central wavelength: ', coarse_lambda_central, ' micron'
	print, 'pixel scale: ', coarse_pixel_scale, ' +/- ', coarse_pixel_scale_error, ' micron per pixel'

	print, ''
	print, 'SECOND STEP: refine wavelength solution by using non-uniform wavelength axis'
	print, '-----------------------------------------------------------------------------'

	; set the size of the grid
	N_pix_scale = 50
	N_pix_scale_variation = 50

	; make the grid
	; for the pixel scale, bracket the value found in the coarse fit, +/- 25% of its value
	pix_scale_grid = coarse_pixel_scale + 0.25*coarse_pixel_scale * (-0.5 + dindgen(N_pix_scale)/double(N_pix_scale-1) )

	; the other parameter is the relative variation in pixel scale. This is a guess but should include all realistic cases.
	pix_scale_variation_grid = -2d-4 + 4d-4 * dindgen(N_pix_scale_variation)/double(N_pix_scale_variation-1)

	; cross-correlate the observed and model sky spectra and find a good wavelength solution
	flame_wavecal_crosscorr, sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
		 approx_lambda_central=coarse_lambda_central, pix_scale_grid=pix_scale_grid, pix_scale_variation_grid=pix_scale_variation_grid, $
		 lambda_axis=lambda_axis, title='after second cross-correlation'


	print, ''
	print, 'THIRD STEP: further refine wavelength solution'
	print, '-----------------------------------------------------------------------------'


	; set the size of the grid
	N_pix_scale = 20
	N_pix_scale_variation = 20


	; make the grid
	; for the pixel scale, bracket the value found in the second step, +/- 15% of its value
	i_mid = n_elements(lambda_axis)/2
 	approx_pixel_scale = lambda_axis[i_mid+1] - lambda_axis[i_mid]
	pix_scale_grid = approx_pixel_scale + 0.15*approx_pixel_scale * (-0.5 + dindgen(N_pix_scale)/double(N_pix_scale-1) )

	; the other parameter is the relative variation in pixel scale
	approx_pix_scale_variation = (lambda_axis[i_mid+1]-lambda_axis[i_mid]) / (lambda_axis[i_mid] - lambda_axis[i_mid-1]) - 1d
	pix_scale_variation_grid = approx_pix_scale_variation + 0.15*approx_pix_scale_variation * (-0.5 + dindgen(N_pix_scale_variation)/double(N_pix_scale_variation-1) )

	; cross-correlate the observed and model sky spectra and find a good wavelength solution
	flame_wavecal_crosscorr, sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
		 approx_lambda_central=median(lambda_axis), pix_scale_grid=pix_scale_grid, pix_scale_variation_grid=pix_scale_variation_grid, $
		 lambda_axis=lambda_axis, title='after third cross-correlation'

  cgPS_close

	; return the approximate wavelength axis
	return, lambda_axis

END




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_wavecal_rough, fuel=fuel

	start_time = systime(/seconds)

  print, ' '
  print, 'flame_wavecal_rough'
  print, '*******************'
  print, ' '

	; extract the slits structures
	slits = fuel.slits

	; loop through all slits
	for i_slit=0, n_elements(slits)-1 do begin

		this_slit = fuel.slits[i_slit]

	  print, 'Rough wavelength calibration for slit ', strtrim(this_slit.number,2), ' - ', this_slit.name
		print, ' '

		rough_wavecal = flame_wavecal_approximate( fuel=fuel, this_slit=this_slit)
    *(fuel.slits[i_slit].rough_wavecal) = rough_wavecal

  endfor

	print, 'It took ', systime(/seconds) - start_time, ' seconds'

END
