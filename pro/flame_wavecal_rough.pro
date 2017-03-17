;
; First pass at a rough wavelength calibration.
; For each slit, load only the first frame, extract the spectrum
; of the central few pixel rows, and cross-correlate that
; with a model spectrum of the sky until you find a reasonable match
;


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_crosscorr, observed_sky=observed_sky, model_lambda=model_lambda, model_flux=model_flux, $
	 lambda0_range=lambda0_range, pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, a3_grid=a3_grid, $
   R_smooth = R_smooth, plot_title=plot_title, $
	 wavecal_coefficients=wavecal_coefficients

   ; inputs:
   ; observed_sky (1D spectrum)
   ; model_lambda, model_flux (1D spectrum and corresponding wavelength axis)
   ; lambda0_range (array, the two values bracketing the wavelength of the first pixel)
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

  ; if necessary, smooth spectrum
  if keyword_set(R_smooth) then begin

    ; calculate the sigma for smoothing given the spectral resolution R
    expected_lambdacen = mean(lambda0_range) + 0.5*median(pix_scale_grid)*N_skypix
  	smoothing_sigma = expected_lambdacen / (2.36 * R_smooth)

    ; smooth observed sky
  	sky = gauss_smooth(sky, smoothing_sigma / median(pix_scale_grid) )

  endif

  ; trim edge to avoid problems
  sky[0:10] = !values.d_nan
  sky[-11:-1] = !values.d_nan
	sky[where(sky eq 0.0, /null)] = !values.d_nan

  ; normalize spectrum
  sky -= median(sky)
  sky /= max(sky, /nan)

  ; let's assume the mid-point lambda as the reference one
  lambda_ref = mean(lambda0_range)


  ; setup model spectrum
  ;-----------------------------------------------------------------------------

  ; name change to avoid overriding argument when calling this procedure
  model_l = model_lambda
  model_f = model_flux

  ; if necessary, smooth spectrum
  if keyword_set(R_smooth) then $
    model_f = gauss_smooth(model_f, smoothing_sigma / (model_l[1]-model_l[0]) )

  ; normalize spectrum
  model_f -= median(model_f)
  model_f /= max(model_f, /nan)


	; show spectra before cross-correlation
  ;-----------------------------------------------------------------------------
  erase
  approx_coeff = [lambda_ref, median(pix_scale_grid), median(a2_grid), median(a3_grid) ]
	cgplot, poly(indgen(N_skypix), approx_coeff), sky, $
		charsize=1, thick=3, xtit='', title=plot_title, layout=[1,2,1]
	cgplot, model_l, model_f, color='red', /overplot, thick=3
  cgtext, 0.5, 0.93, 'before cross-correlation', charsize=0.8, /normal, alignment=0.5
	cgtext, 0.75, 0.93, 'observed sky', charsize=0.8, /normal
	cgtext, 0.75, 0.91, 'model sky', charsize=0.8, /normal, color='red'


  ; cross-correlation
  ;-----------------------------------------------------------------------------

  ; dimensionality of the grid
  N1 = n_elements(pix_scale_grid)
  N2 = n_elements(a2_grid)
  N3 = n_elements(a3_grid)

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

    ; coefficients for this loop, assuming the reference lambda
    this_coeff = [ lambda_ref, pix_scale_grid[i1], a2_grid[i2], a3_grid[i3] ]

    ; calculate the lambda axis corresponding to these coefficients
    lambda_axis = poly(indgen(N_skypix), this_coeff )

    ; make a new model lambda that is linear and uses the current pixel scale
    this_model_l = model_l[0] + this_coeff[1]*dindgen( (model_l[-1]-model_l[0])/this_coeff[1] )

    ; interpolate the model onto the new wavelength axis
	  this_model_f = interpol( model_f, model_l, this_model_l)

    ; we now need to linearize the observed lambda axis (does nothing if a2=a3=0)
    lambda_axis_linear = lambda_axis[0] + $
      this_coeff[1]*dindgen( (lambda_axis[-1]-lambda_axis[0])/this_coeff[1] )

    ; interpolate the observed sky on a linear axis wavelength
    sky_linear = interpol( sky, lambda_axis, lambda_axis_linear)

    ; identify the indices corresponding to the sky spectrum assuming the reference lambda
    w_start = (where(this_model_l GE lambda_ref, /null))[0]
    w_end = w_start + N_skypix-1
    if w_end GE n_elements(this_model_f) then begin
      print, 'sky model does not cover enough wavelength range'
      continue
    endif

    ; calculate the cross-correlation
    ; ----------------------

    ; wavelength shift values to probe
    lag = indgen( (lambda0_range[1]-lambda0_range[0]) / this_coeff[1] ) ; in pixels

    ; but with zero corresponding to the reference lambda, at the center of the range
    lag -= n_elements(lag)/2

    ; array that will contain the cross-correlation values
    cross = dblarr(n_elements(lag))

    ; loop through all the lag values
    for k=0, n_elements(lag)-1 do begin

      ; shift the model spectrum according to the lag
      y = shift(this_model_f, -lag[k])

      ; select the region of interest
      y = y[w_start: w_end]

      ; normalize model
      y -= median(y)
      y /= max(y, /nan)

      ; cgplot, sky_linear, title=number_formatter(this_coeff[1], decimals=4)
      ; cgplot, y, /overplot, color='red'
      ; wait, 0.000001

      ; calculate cross-correlation (total of product divided by variance)
      cross[k] = total( sky_linear * y, /nan ) / sqrt( total(sky_linear^2, /nan) * total(y^2, /nan) )

    endfor

    ; find the lag that maximizes the cross-correlation
		max_cc = max(cross, ind, /nan)

    ; keep only the lag and cc value at the peak
		cc_table[i3, i2, i1] = max_cc
		delta_table[i3, i2, i1] = lag[ind]

	endfor

	; find the best point in the grid, the one with the highest value of cross-correlation
	max_cc = max(cc_table, ind)

	; wavelength shift in pixels
	best_delta = delta_table[ind]

  ; best-fit coefficients
  if (size(cc_table))[0] eq 2 then begin  ; necessary for when pix_scale is fixed
    ind2d = array_indices(cc_table, ind)
    ind3d = [ind2d[0], ind2d[1], 0]
  endif else $
	 ind3d = array_indices(cc_table, ind)
  best_coefficients = [ lambda_ref, pix_scale_grid[ind3d[2]], a2_grid[ind3d[1]], a3_grid[ind3d[0]] ]

	print, ''
  print, 'cross correlation at peak = ', max_cc
  print, 'pixel shift = ', best_delta
  print, 'pixel scale = ', best_coefficients[1]
	print, 'a2 = ', best_coefficients[2]

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
  cgtext, 0.5, 0.43, 'after cross-correlation', charsize=0.8, /normal, alignment=0.5

	; return the wavelength solution
	wavecal_coefficients = new_coefficients

END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_wavecal_rough_oneslit, fuel=fuel, this_slit=this_slit

	;
	; Two steps:
	; First, use a coarse, logarithmic grid with constant pixel scale (uniform wavelength solution)
	; Second, use a fine grid centered on the coarse value of pixel scale, and explore also the
	; pix_scale_variation parameter for non-uniform wavelength solution
	;

	; load the slit image
	;---------------------

  ; take the first frame
  slit_filename = (*this_slit.filenames)[0]
  print, 'Using the central pixel rows of ', slit_filename

	; read in slit
	im = readfits(slit_filename, hdr)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(im))[2]

	; sky spectrum: extract along the central 5 pixels
  ; NB: NEED TO MAKE SURE THAT THIS IS, INDEED, SKY, AND NOT THE OBJECT
	sky = median(im[*, N_spatial_pix/2-3 : N_spatial_pix/2+2], dimension=2)

  range_pixel_scale = this_slit.range_pixel_scale
  range_start_lambda = this_slit.range_lambda0


  ; load the model sky spectrum
	;---------------------
  readcol, fuel.util.sky_emission_filename, model_lambda, model_flux	; lambda in micron

  ; cut out a reasonable range
  wide_range = [ range_start_lambda[0], range_start_lambda[1] + range_pixel_scale[1]*n_elements(sky) ]
  w_reasonable = where(model_lambda GT wide_range[0] $
    and model_lambda LT wide_range[1], /null)
  if w_reasonable EQ !null then $
      message, 'the sky model does not cover the required wavelength range'
  model_lambda = model_lambda[w_reasonable]
  model_flux = model_flux[w_reasonable]


  ; first step: smoothed spectrum, find zero-point and pixel scale
	;---------------------

	; start PS file
  ps_filename = file_dirname(slit_filename, /mark_directory) + 'wavelength_solution_estimate.ps'
	cgPS_open, ps_filename, /nomatch

	print, ''
	print, 'FIRST STEP: rough estimate of lambda and delta_lambda'
	print, '-----------------------------------------------------'

	; first, we assume a constant pixel scale, in micron per pixels:
  pix_scale_grid = range_pixel_scale[0] + $
    (range_pixel_scale[1]-range_pixel_scale[0])*dindgen(51)/50.0

  flame_wavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
	 lambda0_range=range_start_lambda, pix_scale_grid=pix_scale_grid, $
   R_smooth = fuel.input.rough_wavecal_R[0], plot_title='first: find pixel scale and zero-point', $
	 wavecal_coefficients=wavecal_coefficients


  ; second step: smoothed spectrum, add second-order variations
	;---------------------

	print, ''
	print, 'SECOND STEP: rough estimate of second-order variation'
	print, '-----------------------------------------------------'

	; set the size of the grid
	N1 = 20
	N2 = 40

	; make the grid
  ; for the pixel scale, bracket the value found in the coarse fit, +/- 20% of its value
	pix_scale_grid = wavecal_coefficients[1] *( 0.80 + 0.40*dindgen(N1)/double(N1-1) )

  ; we assume that the pixel scale does not vary by more than a factor of 2
  ; across the full spectrum. This gives us the extreme negative values for a2:
  a2_ref = wavecal_coefficients[1] / (4d*n_elements(sky) + 1d)

  ; then we make a generic grid, from 1/1000 to 1 with logarithmic spacing
  log_grid = 10.0^( -3.0 + 3.0*dindgen(N2/2)/double(N2/2-1))

  ; the grid for a2 needs to include 0, and be logarithmic but also positive and negative
  a2_grid = [ -a2_ref * reverse(log_grid), 0.0, a2_ref * log_grid ]

  ; there should not be a large shift in wavelength now - allow up to 30% of the full lambda range
	fullrange = wavecal_coefficients[1]*n_elements(sky)
  lambda0_range = wavecal_coefficients[0] + fullrange*[-0.3,0.3]

	print, 'lambda0: ', lambda0_range
	print, 'a1: ', pix_scale_grid[0], pix_scale_grid[-1]
	print, 'a2: ', a2_grid[0], a2_grid[-1]
	print, a2_grid

  flame_wavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
	 lambda0_range=lambda0_range, pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, $
   R_smooth = fuel.input.rough_wavecal_R[1], plot_title='second: use second-order polynomial', $
	 wavecal_coefficients=wavecal_coefficients


	print, ''
	print, 'THIRD STEP: refine wavelength solution by using non-uniform wavelength axis'
	print, '-----------------------------------------------------------------------------'

	; set the size of the grid
	N1 = 21
	N2 = 41

	; make the grid
  ; for the pixel scale, bracket the value found in the previous fit, +/- 10% of its value
	pix_scale_grid = wavecal_coefficients[1] *( 0.90 + 0.20*dindgen(N1)/double(N1-1) )

	; for the a2 grid, bracket the value found in the previous fit, +/- a factor of 3
  a2_grid =  wavecal_coefficients[2] *( 0.3 + 2.7*dindgen(N1)/double(N1-1) )

  ; there should not be a large shift in wavelength now
	fullrange = wavecal_coefficients[1]*n_elements(sky)
  lambda0_range = wavecal_coefficients[0] + fullrange*[-0.1,0.1]


  flame_wavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
	 lambda0_range=lambda0_range, pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, $
   R_smooth = fuel.input.rough_wavecal_R[2], plot_title='third: use second-order polynomial', $
	 wavecal_coefficients=wavecal_coefficients

	; print, ''
	; print, 'THIRD STEP: further refine wavelength solution using 3rd order polynomial'
	; print, '-----------------------------------------------------------------------------'
  ;
  ; 	; set the size of the grid
  ; 	N1 = 3
  ; 	N2 = 61
  ;   N3 = 61
  ;
  ; 	; make the grid
  ; 	; for the pixel scale, bracket the value found in the previous fit, +/- 5% of its value
  ; 	pix_scale_grid = wavecal_coefficients[1] * ( 0.95 + 0.10 * dindgen(N1)/double(N1-1) )
  ;   pix_scale_grid = [ wavecal_coefficients[1] ]
  ;
  ; 	; same thing for the second order coefficient (+/- 50%)
  ; 	a2_grid = wavecal_coefficients[2] * ( 0.5 + 1.0* dindgen(N2)/double(N2-1) )
  ;
  ;   ; make the grid for the third-order coefficient
  ;   a3_grid = wavecal_coefficients[2] * 10.0^(-7.0 + 6.0 * dindgen(N3/2)/double(N3/2-1) )
  ;   a3_grid = [reverse(-a3_grid), a3_grid]
  ;
  ;   flame_wavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
  ; 	 approx_lambda_0=wavecal_coefficients[0], pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, a3_grid=a3_grid, $
  ;    plot_title='third: use third-order polynomial', $
  ; 	 wavecal_coefficients=wavecal_coefficients

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

		rough_wavecal = flame_wavecal_rough_oneslit( fuel=fuel, this_slit=this_slit)
    *(fuel.slits[i_slit].rough_wavecal) = rough_wavecal

  endfor

	print, 'It took ', systime(/seconds) - start_time, ' seconds'

END
