;
; First pass at a rough wavelength calibration.
; For each slit, load only the first frame, extract the spectrum
; of the central few pixel rows, and cross-correlate that
; with a model spectrum of the sky until you find a reasonable match
;


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_roughwavecal_crosscorr, observed_sky=observed_sky, model_lambda=model_lambda, model_flux=model_flux, $
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

	; workaround so that the median() command always works
	if n_elements(pix_scale_grid) eq 1 then pix_scale_grid = [pix_scale_grid]

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
    expected_lambdacen = mean(lambda0_range)
  	smoothing_sigma = expected_lambdacen / (2.36 * double(R_smooth))

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

		; check that the wavelength solution is not double valued - in that case skip this loop
		if this_coeff[3] NE 0.0 then begin
			x_sol1 = ( -this_coeff[2] + sqrt(thiscoeff[2]^2 - 3*this_coeff[1]*this_coeff[3]) ) / (3.0*this_coeff[3])
			x_sol2 = ( -this_coeff[2] - sqrt(thiscoeff[2]^2 - 3*this_coeff[1]*this_coeff[3]) ) / (3.0*this_coeff[3])
			if x_sol1 GT 0.0 and x_sol1 LT N_skypix then continue
			if x_sol2 GT 0.0 and x_sol2 LT N_skypix then continue
		endif else $
			if this_coeff[2] NE 0.0 then begin
				x_sol = -0.5*this_coeff[1]/this_coeff[2]
				if x_sol GT 0.0 and x_sol LT N_skypix then continue
			endif

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

	print, ' '

	; check if the best-fit values are at the edge of the grids
	if n_elements(pix_scale_grid) GT 1 then if ind3d[2] eq 0 or ind3d[2] eq n_elements(pix_scale_grid) -1 then $
		print, '******* WARNING! best-fit value for the pixel scale is at the edge of the grid! *******'

	if n_elements(a2_grid) GT 1 then if ind3d[1] eq 0 or ind3d[1] eq n_elements(a2_grid) -1 then $
		print, '******* WARNING! best-fit value for a2 is at the edge of the grid! *******'

	if n_elements(a3_grid) GT 1 then if ind3d[0] eq 0 or ind3d[0] eq n_elements(a3_grid) -1 then $
		print, '******* WARNING! best-fit value for a3 is at the edge of the grid! *******'

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

FUNCTION flame_roughwavecal_solution, fuel=fuel, this_slit=this_slit, sky=sky, $
		 model_lambda=model_lambda, model_flux=model_flux, wavecal_coefficients=wavecal_coefficients, $
		 wavecal_rough_R=wavecal_rough_R

	range_delta_lambda = this_slit.range_delta_lambda
  range_start_lambda = this_slit.range_lambda0

	if keyword_set(wavecal_rough_R) then begin
		R_smooth_1 = wavecal_rough_R[0]
		R_smooth_2 = wavecal_rough_R[1]
		R_smooth_3 = wavecal_rough_R[2]
	endif


	  ; first step: smoothed spectrum, find zero-point and pixel scale
		;---------------------

		print, ''
		print, 'STEP 1: rough estimate of lambda and delta_lambda'

		; set the size of the grid
		N1 = 51

		; first, we assume a constant pixel scale, in micron per pixels:
	  pix_scale_grid = range_delta_lambda[0] + $
	    (range_delta_lambda[1]-range_delta_lambda[0])*dindgen(N1)/double(N1-1)

	  flame_roughwavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
		 lambda0_range=range_start_lambda, pix_scale_grid=pix_scale_grid, $
	   R_smooth = R_smooth_1, plot_title='step 1: find pixel scale and zero-point', $
		 wavecal_coefficients=wavecal_coefficients


	  ; second step: smoothed spectrum, add second-order variations
		;---------------------

		print, ''
		print, 'STEP 2: rough estimate of second-order variation'

		; set the size of the grid
		N1 = 40
		N2 = 60

		; make the grid
	  ; for the pixel scale, bracket the value found in the coarse fit, +/- 30% of its value
		pix_scale_grid = wavecal_coefficients[1] *( 0.70 + 0.60*dindgen(N1)/double(N1-1) )

	  ; we assume that the pixel scale does not vary by more than a factor of 2
	  ; across the full spectrum. This gives us the extreme negative values for a2:
	  a2_ref = 1.5 * wavecal_coefficients[1] / (4d*n_elements(sky) + 1d)

	  ; then we make a generic grid, from 1/1000 to 1 with logarithmic spacing
	  log_grid = 10.0^( -3.0 + 3.0*dindgen(N2/2)/double(N2/2-1))

	  ; the grid for a2 needs to include 0, and be logarithmic but also positive and negative
	  a2_grid = [ -a2_ref * reverse(log_grid), 0.0, a2_ref * log_grid ]

	  ; there should not be a large shift in wavelength now - allow up to 15% of the full lambda range
		fullrange = wavecal_coefficients[1]*n_elements(sky)
	  lambda0_range = wavecal_coefficients[0] + fullrange*[-0.15,0.15]

		; print, 'lambda0: ', lambda0_range
		; print, 'a1: ', pix_scale_grid[0], pix_scale_grid[-1]
		; print, 'a2: ', a2_grid[0], a2_grid[-1]

	  flame_roughwavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
		 lambda0_range=lambda0_range, pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, $
	   R_smooth = R_smooth_2, plot_title='step 2: use second-order polynomial', $
		 wavecal_coefficients=wavecal_coefficients


		if ~fuel.settings.wavecal_rough_split then begin


			print, ''
			print, 'STEP 3: refine wavelength solution by using non-uniform wavelength axis'

			; set the size of the grid
			N1 = 41
			N2 = 61

			; make the grid
		  ; for the pixel scale, bracket the value found in the previous fit, +/- 20% of its value
			pix_scale_grid = wavecal_coefficients[1] *( 0.80 + 0.40*dindgen(N1)/double(N1-1) )

			; for the a2 grid, bracket the value found in the previous fit, +/- a factor of 5
		  a2_grid =  wavecal_coefficients[2] *( 0.2 + 4.8*dindgen(N1)/double(N1-1) )

		  ; there should not be a large shift in wavelength now
			fullrange = wavecal_coefficients[1]*n_elements(sky)
		  lambda0_range = wavecal_coefficients[0] + fullrange*[-0.05,0.05]

		  flame_roughwavecal_crosscorr, observed_sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
			 lambda0_range=lambda0_range, pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, $
		   R_smooth = R_smooth_3, plot_title='step 3: use second-order polynomial', $
			 wavecal_coefficients=wavecal_coefficients


			wavecal_solution = poly(dindgen(n_elements(sky)), wavecal_coefficients)
			return, wavecal_solution

		endif


		print, ''
		print, 'STEP 3: find wavelength solution for the two halves separately'


		; set the size of the grid
		N1 = 41
		N2 = 61

		; make the grid
	  ; for the pixel scale, bracket the value found in the previous fit, +/- 20% of its value
		pix_scale_grid = wavecal_coefficients[1] *( 0.80 + 0.40*dindgen(N1)/double(N1-1) )

		; make the a2 grid; assume the sign does not change
		a2_ref = wavecal_coefficients[1] / (4d*n_elements(sky) + 1d)
		; log_grid = 10.0^( -3.0 + 3.0*dindgen(N2/2)/double(N2/2-1))
		; a2_grid = [ -a2_ref * reverse(log_grid), 0.0, a2_ref * log_grid ]
	  a2_grid =  wavecal_coefficients[2] * 3.0*dindgen(N2/2)/double(N2/2-1)

	  ; there should not be a large shift in wavelength now
		fullrange = wavecal_coefficients[1]*n_elements(sky)
	  lambda0_range = wavecal_coefficients[0] + fullrange*[-0.05,0.05]

		; split sky into two halves
		x_c = n_elements(sky)/2
		sky_1 = sky[0 : x_c-1]
		sky_2 = sky[x_c : *]

		; for the second half, need to adjust the parameters
		pix_scale_grid_2 = pix_scale_grid + 2.0*wavecal_coefficients[2]
		lambda0_range_2 = lambda0_range + x_c * $
		 (wavecal_coefficients[1]+2.0*wavecal_coefficients[2]) + wavecal_coefficients[2]*x_c^2

		; cross-correlate the first half
	  flame_roughwavecal_crosscorr, observed_sky=sky_1, model_lambda=model_lambda, model_flux=model_flux, $
		 lambda0_range=lambda0_range, pix_scale_grid=pix_scale_grid, a2_grid=a2_grid, $
	   R_smooth = R_smooth_3, plot_title='step 3: first half', $
		 wavecal_coefficients=wavecal_coefficients_1

		 ; cross-correlate the second half
 	  flame_roughwavecal_crosscorr, observed_sky=sky_2, model_lambda=model_lambda, model_flux=model_flux, $
 		 lambda0_range=lambda0_range_2, pix_scale_grid=pix_scale_grid_2, a2_grid=a2_grid, $
 	   R_smooth = R_smooth_3, plot_title='step 3: second half', $
 		 wavecal_coefficients=wavecal_coefficients_2

		 ; calculate the two wavelength solutions for the two chunks
		 wavecal_solution_1 = poly(dindgen(n_elements(sky_1)), wavecal_coefficients_1)
		 wavecal_solution_2 = poly(dindgen(n_elements(sky_2)), wavecal_coefficients_2)

		 ; merge the two halves
		 wavecal_solution = [ wavecal_solution_1, wavecal_solution_2 ]

		 ; smooth out the potential edge in the middle
		 wavecal_solution = smooth(wavecal_solution, 100)

		return, wavecal_solution

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_roughwavecal_oneslit, fuel=fuel, this_slit=this_slit, rough_skyflux=rough_skyflux
	;
	; Three steps:
	; First, use a coarse, logarithmic grid with constant pixel scale (uniform wavelength solution)
	; Second and third, use a fine grid centered on the coarse value of pixel scale, and explore also the
	; pix_scale_variation parameter for non-uniform wavelength solution
	;

	; load the slit image
	;---------------------

  ; take the first frame
  slit_filename = (this_slit.cutouts.filename)[0]
  print, 'Using the central pixel rows of ', slit_filename

	; read in slit
	im = readfits(slit_filename, hdr)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(im))[2]

	; calculate the distance of each pixel row from the spatial center of the slit
	dist_to_center = abs( indgen(N_spatial_pix) - N_spatial_pix/2 )

	; identify a sky region near the center of the slit
	;---------------------

	; extract profile
	profile = median(im, dimension=1)

	; detect positive flux fluctuations, probably due to bright objects
	w_objects = where( profile - median(profile) GT 3.0* stddev(profile, /nan), /null)

	; exclude these pixels by assigning them an arbitrary high distance
	if w_objects NE !NULL then dist_to_center[w_objects] = 999

	; identify the most-central five pixel rows that are not contaminated by bright objects
	w_skyregion = (sort(dist_to_center))[0:4]

	; extract the sky region
	sky_region = []
	for i=0, n_elements(w_skyregion)-1 do sky_region = [ [sky_region], [im[*,w_skyregion[i]]] ]

	; extract the spectrum from the sky region
	sky = median(sky_region, dimension=2)


	; filter the sky spectrum: find the running minimum pixel value
	;--------------------------------------------------------------
	if fuel.settings.wavecal_rough_smooth_window GT 0 then begin

		; make a 2D matrix where at each value of x-pixel you have a column with all the neighhboring pixel values
		matrix = []
		for i_shift = -fuel.settings.wavecal_rough_smooth_window/2, fuel.settings.wavecal_rough_smooth_window/2 do $
		 	matrix = [ [matrix], [shift(sky, i_shift)]]

		; for each x-pixel take the minimum of all the neighboring pixel values
		sky_minfilter = min(matrix, /nan, dimension=2)

		; subtract the continuum
		sky -= sky_minfilter

		; crop the edges
		sky[0:fuel.settings.wavecal_rough_smooth_window/2] = 0
		sky[-fuel.settings.wavecal_rough_smooth_window/2-1:*] = 0

	endif


  ; load the model sky spectrum
	;---------------------
  readcol, fuel.settings.sky_emission_filename, model_lambda, model_flux	; lambda in micron

  ; cut out a reasonable range
  wide_range = [ this_slit.range_lambda0[0], this_slit.range_lambda0[1] + this_slit.range_delta_lambda[1]*n_elements(sky) ]
  w_reasonable = where(model_lambda GT wide_range[0] $
    and model_lambda LT wide_range[1], /null)
  if w_reasonable EQ !null then $
      message, 'the sky model does not cover the required wavelength range'
  model_lambda = model_lambda[w_reasonable]
  model_flux = model_flux[w_reasonable]

	; start PS file
  ps_filename = file_dirname(slit_filename, /mark_directory) + 'rough_wavelength_calibration.ps'
	cgPS_open, ps_filename, /nomatch

	; find the wavelength solution
	wavecal_solution = flame_roughwavecal_solution(fuel=fuel, this_slit=this_slit, sky=sky, $
		model_lambda=model_lambda, model_flux=model_flux, wavecal_rough_R = fuel.settings.wavecal_rough_R)


  cgPS_close

	; output the observed sky spectrum used for the wavecal
	rough_skyflux = sky

	; return the approximate wavelength axis
	return, wavecal_solution

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_roughwavecal_oneslit_witharcs, fuel=fuel, this_slit=this_slit, rough_arcflux=rough_arcflux
	;
	; Same as above, but use arcs instead of skylines
	;

	; load the arc spectrum
	;---------------------

  ; filename for cutout of arc spectrum
	arc_filename = this_slit.arc_cutout.filename
  print, 'Using the central pixel rows of ', arc_filename

	; read in arc spectrum
	im = readfits(arc_filename, hdr)

	; how many pixels on either directions
	N_spatial_pix = (size(im))[2]
	N_spectral_pix = (size(im))[1]

	; extract the spectrum
	arc_spectrum = median(im[*,N_spatial_pix/2-2:N_spatial_pix/2+2], dimension=2)
	arc_xcoord = dindgen(N_spectral_pix)

	; smooth spectrum
	arc_spectrum = gauss_smooth(arc_spectrum, 3, /nan)

	; define a flux threshold for bright lines: standard deviation of the spectrum without the extreme points
	sorted_values = arc_spectrum[sort(arc_spectrum)]
	threshold = 0.5*stddev(sorted_values[N_spectral_pix/10 : -N_spectral_pix/10], /nan)

	; automatically identify peaks
	w_peaks = where( arc_spectrum GT shift(arc_spectrum, 1) and arc_spectrum GT shift(arc_spectrum, 2) and $
		arc_spectrum GT shift(arc_spectrum, 3) and arc_spectrum GT shift(arc_spectrum, -1) and $
		arc_spectrum GT shift(arc_spectrum, -2) and arc_spectrum GT shift(arc_spectrum, -3) and $
		arc_spectrum - 0.5*(shift(arc_spectrum, 3)+shift(arc_spectrum, -3)) GT threshold, /null )

	print, n_elements(w_peaks), ' bright lines identified in the arc spectrum'

	; make array that will contain the x coordinate of each arc line
	arc_lines = []

	; fit a Gaussian to each line
	for i_line=0,n_elements(w_peaks)-1 do begin

		; select the region to fit - +/- 10 pixels
		w_fit = where( abs(arc_xcoord-w_peaks[i_line]) LE 10, /null )

		; check that the region is within the observed range
		if w_fit eq !NULL then continue

    ; check that there actually is signal and it's not just a bunch of NaNs
    if n_elements( where( finite(arc_spectrum[w_fit]), /null ) ) LE 5 then continue

		; error handling for the gaussian fitting
		catch, error_gaussfit
		if error_gaussfit ne 0 then begin
			print, 'GAUSSFIT ERROR STATUS: ' + strtrim(error_gaussfit,2)
			catch, /cancel
			continue
		endif

		; estimate parameters of the Gaussian
		est_peak = arc_spectrum[w_peaks[i_line]]
		est_center = w_peaks[i_line]
		est_sigma = 2.0
		est_cont = min( median( arc_spectrum[w_fit], 3) , /nan)

		; Gaussian fit
		!NULL = gaussfit( arc_xcoord[w_fit], arc_spectrum[w_fit], gauss_param, nterms=4, $
			estimates=[est_peak, est_center, est_sigma, est_cont], sigma=gauss_err, chisq=chisq )

		; check that chi square makes sense
		if ~finite(chisq) then continue

		; check that the peak of the Gaussian is positive
		if gauss_param[0] LT 0.0 then continue

		; check that the SNR is high
		if gauss_param[0] LT 5.0*gauss_err[0] then continue

		; check that the center of the Guassian is in the observed range
		if gauss_param[1] LT min(arc_xcoord[w_fit]) or gauss_param[1] GT max(arc_xcoord[w_fit]) then continue

		; check that the Gaussian width makes sense
		if gauss_param[2] LT 0.3 or gauss_param[2] GT 20.0 then continue

		; add to the stack
		arc_lines = [arc_lines, gauss_param[1]]

	endfor

	; make a simple spectrum from the identified line list, assuming all lines have equal intensity
	arc_simplespectrum = arc_xcoord*0.0
	simple_sigma = 2.0	; assume a sigma of two pixels
	for i=0, n_elements(arc_lines)-1 do arc_simplespectrum += exp(-0.5*(arc_xcoord-arc_lines[i])^2/simple_sigma^2)


  ; load the model spectrum of the arc lamp
	; ------------------------------------------

	; load the linelist
  readcol, fuel.util.intermediate_dir + 'linelist_arcs.txt', all_lines

	; select a reasonable range
  wide_range = [ this_slit.range_lambda0[0] - this_slit.range_delta_lambda[1] * 0.2*N_spectral_pix, this_slit.range_lambda0[1] + this_slit.range_delta_lambda[1] * 1.2*N_spectral_pix]
	all_lines = all_lines[ where(all_lines GT wide_range[0] and all_lines LT wide_range[1], /null) ]

	; make a simple theoretical spectrum from the line list (assuming all lines have equal intensity)
	model_lambda = wide_range[0] + (wide_range[1]-wide_range[0]) * dindgen(N_spectral_pix)/double(N_spectral_pix)
	model_flux = model_lambda*0.0
	for i=0, n_elements(all_lines)-1 do model_flux[ (where(model_lambda GE all_lines[i]))[0] ] = 1
	model_flux[-1] = 0.0	; for when where() returned -1

	; start PS file
  ps_filename = file_dirname( arc_filename, /mark_directory ) + 'rough_wavelength_calibration_arc.ps'
	cgPS_open, ps_filename, /nomatch

	; find the wavelength solution
	wavecal_solution = flame_roughwavecal_solution(fuel=fuel, this_slit=this_slit, sky=arc_simplespectrum, $
		model_lambda=model_lambda, model_flux=model_flux)

	cgPS_close

	; output the observed sky spectrum used for the wavecal
	rough_arcflux = arc_spectrum

	; return the approximate wavelength axis
	return, wavecal_solution

END




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_roughwavecal, fuel

	flame_util_module_start, fuel, 'flame_roughwavecal'


	; extract the slits structures
	slits = fuel.slits

	; loop through all slits
	for i_slit=0, n_elements(slits)-1 do begin

		if fuel.slits[i_slit].skip then continue

		print, ''
		print, ''
	  print, 'Rough wavelength calibration for slit ', strtrim(fuel.slits[i_slit].number,2), ' - ', fuel.slits[i_slit].name
		print, '--------------------------------------------------------------------'
		print, ''

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

		; deal with arcs
		if fuel.util.arc.n_frames GT 0 then begin

			; find rough calibration using integrated arc spectrum
			rough_arclambda = flame_roughwavecal_oneslit_witharcs( fuel=fuel, this_slit=fuel.slits[i_slit], rough_arcflux=rough_arcflux)
    	*(fuel.slits[i_slit].rough_arclambda) = rough_arclambda
			*(fuel.slits[i_slit].rough_arcflux) = rough_arcflux

		endif else begin

			; find rough calibration using integrated sky spectrum
			rough_skylambda = flame_roughwavecal_oneslit( fuel=fuel, this_slit=fuel.slits[i_slit], rough_skyflux=rough_skyflux)
	    *(fuel.slits[i_slit].rough_skylambda) = rough_skylambda
			*(fuel.slits[i_slit].rough_skyflux) = rough_skyflux

		endelse


  endfor


  flame_util_module_end, fuel

END
