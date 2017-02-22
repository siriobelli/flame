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


PRO flame_wavecal_crosscorr, observed_sky=observed_sky, model_lambda=model_lambda, model_flux=model_flux, $
	 approx_lambda_0=approx_lambda_0, pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, a3_grid=a3_grid, $
   R_smooth = R_smooth, plot_title=plot_title, $
	 wavecal_coefficients=wavecal_coefficients

   ; inputs:
   ; observed_sky (1D spectrum)
   ; model_lambda, model_flux (1D spectrum and corresponding wavelength axis)
   ; approx_lambda_0 (scalar, approximate wavelength of first pixel of observed_sky)
   ; pix_scale_grid (array of pixel scale values to be considered)
   ; a2_grid (optional, array of a2 values to be considered)
   ; a3_grid (optional, array of a3 values to be considered)
   ; R_smooth (optional, scalar, spectral resolution R for smoothing)
   ; plot_title (string)
   ;
   ; outputs:
   ; wavecal_coefficients (array, best-fit poly coefficients to describe the lambda axis)

	;
	; Estimate the wavelength solution by comparing the observed sky spectrum with a model sky spectrum
	; The wavelength solution is parameterized as a polynomial function of the pixel coordinate
  ; Loop through a grid of pixel scale values and higher order coefficients and
  ; at each step cross-correlate the spectra to find the zero-th coefficient.
	; Return the polynomial coefficients.
	;


  ; if some of the coefficient grids are not specified, then fix them to 0
  if ~keyword_set(a2_grid) then a2_grid = [0.0]
  if ~keyword_set(a3_grid) then a3_grid = [0.0]

  ; setup observed spectrum
  ;-----------------------------------------------------------------------------

  ; name change to avoid overriding argument when calling this procedure
  sky = observed_sky

  ; number of pixels
  N_skypix = n_elements(sky)

  ; expected wavelength range
  expected_lambdarange = [approx_lambda_0, approx_lambda_0 + median(pix_scale_grid)*N_skypix]
  expected_lambdacen = 0.5 * (expected_lambdarange[0] + expected_lambdarange[1])

  ; extended wavelength range
  extended_waverange = expected_lambdarange + $
    [-0.5, 0.5] * (expected_lambdarange[1]-expected_lambdarange[0])

  ; if necessary, smooth spectrum
  if keyword_set(R_smooth) then begin

    ; calculate the sigma for smoothing given the spectral resolution R
  	smoothing_sigma = expected_lambdacen / (2.36 * R_smooth)

    ; smooth observed sky
  	sky = gauss_smooth(sky, smoothing_sigma / median(pix_scale_grid) )

  endif

  ; trim edge to avoid problems
  sky[0:10] = 0.
  sky[-11:-1] = 0.

  ; normalize spectrum
  sky -= median(sky)
  sky /= max(sky, /nan)

  ; get rid of NaNs, which create problems for the cross-correlation
  sky[where(~finite(sky), /null)] = 0


  ; setup model spectrum
  ;-----------------------------------------------------------------------------

  ; select the wavelength range of interest
  w_interest = where( model_lambda GT extended_waverange[0] $
  and model_lambda LT extended_waverange[1], /null )
  if w_interest EQ !null then $
    message, 'the sky model does not cover the required wavelength range'
  model_l = model_lambda[w_interest]
  model_f = model_flux[w_interest]

  ; if necessary, smooth spectrum
  if keyword_set(R_smooth) then $
    model_f = gauss_smooth(model_f, smoothing_sigma / (model_l[1]-model_l[0]) )

  ; normalize spectrum
  model_f -= median(model_f)
  model_f /= max(model_f, /nan)

  ; get rid of NaNs, which create problems for the cross-correlation
  model_f[where(~finite(model_f), /null)] = 0

  ; for a tidy extrapolation, first and last elements of model should be zeros
  model_f[0:5] = 0.
  model_f[-5:*] = 0.


	; show spectra before cross-correlation
  ;-----------------------------------------------------------------------------
  erase
  approx_coeff = [approx_lambda_0, median(pix_scale_grid), median(a2_grid), median(a3_grid) ]
	cgplot, poly(indgen(N_skypix), approx_coeff), sky, $
		charsize=1, thick=3, xtit='', title=plot_title, layout=[1,2,1]
	cgplot, model_l, model_f, color='red', /overplot, thick=3
	cgtext, 0.75, 0.88, 'observed sky', charsize=1, /normal
	cgtext, 0.75, 0.84, 'model sky', charsize=1, /normal, color='red'


  ; cross-correlation
  ;-----------------------------------------------------------------------------

  ; dimensionality of the grid
  N1 = n_elements(pix_scale_grid)
  N2 = n_elements(a2_grid)
  N3 = n_elements(a3_grid)

	; wavelength shift values to probe
	lag = -n_elements(sky)/2 + indgen(n_elements(sky))	; in pixels

  ; these tables will store the result of the cross correlation
  ; note: need to store the parameters backwards (N3, N2, N1) so that even
  ; when there are no a2 and a3 we still have a multidimensional array, e.g.[1,1,500]
  cc_table = dblarr(N3, N2, N1)
  delta_table = intarr(N3, N2, N1)

	print, 'Finding wavelength solution via cross-correlation...'

	for i1=0, N1-1 do $
    for i2=0, N2-1 do $
      for i3=0, N3-1 do begin

    ; show progress as percentage of loops completed
		print, strjoin(replicate(string(8B),10)) + $
      strtrim( round( 100. * float(i3+i2*N3+i1*N2*N3)/float(N1*N2*N3-1) ) , 2) + '%' , format='(a,$)'

    this_coeff = [ approx_lambda_0, pix_scale_grid[i1], a2_grid[i2], a3_grid[i3] ]

		; this is the new wavelength axis, with the chosen coefficients
		;lambda_axis = flame_generate_lambda_axis(N_pixels=n_elements(sky), lambda_central=lambda_central, delta_lambda=pix_scale_grid[i], delta_delta_lambda=pix_scale_variation_grid[j])
    lambda_axis = poly(indgen(N_skypix), this_coeff )

		; need to interpolate the model onto the new wavelength axis
		new_model_f = interpol( model_f, model_l, lambda_axis)

		; find the best wavelength shift
		crosscorr = c_correlate(sky, new_model_f, lag)
		max_cc = max(crosscorr, ind, /nan)
		cc_table[i3, i2, i1] = max_cc
		delta_table[i3, i2, i1] = lag[ind]

	endfor

	; find the best point in the grid, the one with the highest value of cross-correlation
	max_cc = max(cc_table, ind)

  ; check that we are within the boundary
	if ind eq 0 or ind eq n_elements(max_cc)-1 then $
		print, 'WARNING: Cross correlation selected value at the edge of the grid!'

	; wavelength shift in pixels
	best_delta = delta_table[ind]

  ; best-fit coefficients
	ind3d = array_indices(cc_table, ind)
  best_coefficients = [ approx_lambda_0, pix_scale_grid[ind3d[2]], a2_grid[ind3d[1]], a3_grid[ind3d[0]] ]

	print, ''
  print, 'pixel shift = ', best_delta
  print, 'pixel scale = ', best_coefficients[1]

  ; applying the shift best_delta to the polynomial we find another polynomial function:
  new_coefficients = [ poly(best_delta, best_coefficients), $
    best_coefficients[1] - 2*best_coefficients[2]*best_delta + 3*best_coefficients[3]*best_delta^2, $
    best_coefficients[2] -3*best_coefficients[3]*best_delta, $
    best_coefficients[3] ]

  ; generate the final best-fit wavelength axis
  lambda_axis = poly(indgen(N_skypix), new_coefficients)
	print, 'lambda axis range = ', lambda_axis[0], lambda_axis[-1]

	; interpolate model flux onto new wavelength axis
	new_model_f = interpol( model_f, model_l, lambda_axis)
	new_model_f /= max(new_model_f, /nan)

  ; show spectra with new wavelength calibration
	cgplot, lambda_axis, sky, charsize=1, thick=3, xtit='wavelength (micron)', title=title, layout=[1,2,2]
	cgplot, lambda_axis, new_model_f, color='red', /overplot, thick=3

	; return the wavelength solution
	wavecal_coefficients = new_coefficients

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

	; load the slit image
	;---------------------

	; read in slit
	im = readfits(slit_filename, hdr)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(im))[2]

	; sky spectrum: extract along the central 5 pixels
  ; NB: NEED TO MAKE SURE THAT THIS IS, INDEED, SKY, AND NOT THE OBJECT
	sky = median(im[*, N_spatial_pix/2-3 : N_spatial_pix/2+2], dimension=2)

	; start PS file
  ps_filename = file_dirname(slit_filename, /mark_directory) + 'wavelength_solution_estimate.ps'
	cgPS_open, ps_filename, /nomatch

	print, ''
	print, 'FIRST STEP: rough estimate of lambda and delta_lambda'
	print, '-----------------------------------------------------'

	; estimated pixel scale from the approximate wavelength range
	estimated_pix_scale = (this_slit.approx_wavelength_hi-this_slit.approx_wavelength_lo) $
		/ double(n_elements(sky))

	; first, we assume a constant pixel scale, in micron per pixels:
	;pix_scale_grid = estimated_pix_scale * (10.0+dindgen(200))/200.0*3.0
  pix_scale_grid = estimated_pix_scale * 10.0^(-0.66 + 1.32*dindgen(100)/99.0)

  flame_wavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
	 approx_lambda_0=this_slit.approx_wavelength_lo, pix_scale_grid=pix_scale_grid, $
   R_smooth = fuel.input.rough_wavecal_R, plot_title='first cross-correlation', $
	 wavecal_coefficients=wavecal_coefficients

	print, ''
	print, 'SECOND STEP: refine wavelength solution by using non-uniform wavelength axis'
	print, '-----------------------------------------------------------------------------'

	; set the size of the grid
	N1 = 11
	N2 = 21

	; make the grid
	; for the pixel scale, bracket the value found in the coarse fit, +/- 15% of its value
	pix_scale_grid = wavecal_coefficients[1] *( 0.85 + 0.30*dindgen(N1)/double(N1-1) )

	; make the grid for the second order polynomial coefficient
	a2_grid = wavecal_coefficients[1] * 10.0^(-7.0 + 6.0 * dindgen(N2/2)/double(N2/2-1) )
  a2_grid = [reverse(-a2_grid), a2_grid]

  flame_wavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
	 approx_lambda_0=wavecal_coefficients[0], pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, $
   R_smooth = fuel.input.rough_wavecal_R, plot_title='second cross-correlation', $
	 wavecal_coefficients=wavecal_coefficients

	print, ''
	print, 'THIRD STEP: further refine wavelength solution using 3rd order polynomial'
	print, '-----------------------------------------------------------------------------'

  	; set the size of the grid
  	N1 = 11
  	N2 = 11
    N3 = 21

  	; make the grid
  	; for the pixel scale, bracket the value found in the previous fit, +/- 5% of its value
  	pix_scale_grid = wavecal_coefficients[1] * ( 0.95 + 0.10 * dindgen(N1)/double(N1-1) )

  	; same thing for the second order coefficient (+/- 50%)
  	a2_grid = wavecal_coefficients[2] * ( 0.5 + 1.0* dindgen(N2)/double(N2-1) )

    ; make the grid for the third-order coefficient
    a3_grid = wavecal_coefficients[2] * 10.0^(-7.0 + 6.0 * dindgen(N3/2)/double(N3/2-1) )
    a3_grid = [reverse(-a3_grid), a3_grid]

    flame_wavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
  	 approx_lambda_0=wavecal_coefficients[0], pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, a3_grid=a3_grid, $
     R_smooth = fuel.input.rough_wavecal_R, plot_title='third cross-correlation', $
  	 wavecal_coefficients=wavecal_coefficients

  cgPS_close

	; return the approximate wavelength axis
	return, poly(indgen(n_elements(sky)), wavecal_coefficients)

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
