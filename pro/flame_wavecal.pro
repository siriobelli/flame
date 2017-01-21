
;
; NOTE: calibration gives wavelengths in vacuum!!
;
;

FUNCTION flame_poly_surface, X, Y, P
	;
	; Used to fit a 2D wavelength solution to a slit using MPFIT
	;

	if n_elements(p) eq 4 then return, p[0] + p[1]*x + p[2]*y + p[3]*x^2

	if n_elements(p) eq 5 then return, p[0] + p[1]*x + p[2]*y + p[3]*x^2 + p[4]*x*y

	if n_elements(p) eq 6 then return, p[0] + p[1]*x + p[2]*y + p[3]*x^2 + p[4]*x*y + p[5]*x^3

	if n_elements(p) eq 7 then return, p[0] + p[1]*x + p[2]*y + p[3]*x^2 + p[4]*x*y + p[5]*y^2 + p[6]*x^3

	return, !values.d_nan


END

; ---------------------------------------------------------------------------------------------------------------------------


PRO flame_wavecal_writeds9, OH_lines, filename=filename
;
; write a ds9 region file with all the OH line detections
;

  ; extract the wavelength of all identifications
  line_lambdas = OH_lines[*,0]
  uniq_lambdas = line_lambdas[UNIQ(line_lambdas, SORT(line_lambdas))]

  ; open file
  openw, lun, filename, /get_lun

  ; write header
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=1 include=1 source=1'
  printf, lun, 'image'

  for i_line=0, n_elements(uniq_lambdas)-1 do begin

  	; the wavelength of this OH line
  	this_lambda = uniq_lambdas[i_line]

  	; select the x and y coordinates for this line
  	w_thisline = where(OH_lines[*,0] eq this_lambda)
  	this_x = OH_lines[w_thisline, 1]
 	this_y = OH_lines[w_thisline, 2]

 	; sort points by y coordinate 
 	wsort = sort(this_y)
 	this_x = this_x[wsort]
 	this_y = this_y[wsort]
 	
    ; concatenate points
    ;all_x = [top_x, reverse(bottom_x)]
    ;all_y = [top_y, reverse(bottom_y)]
    all_x = [ this_x, reverse(this_x) ]
    all_y = [ this_y, reverse(this_y) ]

    ; make the string with all the points
    all_points = ''
    for i=0,n_elements(all_x)-2 do all_points += strtrim(all_x[i],2) + ',' + cgnumber_formatter(all_y[i], decimals=1) + ','
    ; add the last two points without the final comma
    all_points += strtrim(all_x[-1],2) + ',' + cgnumber_formatter(all_y[-1], decimals=1)

    ; write the line corresponding to this slit
    printf, lun, 'polygon(' + all_points + ') # text={' + cgnumber_formatter(this_lambda, decimals=5) + '}'

  endfor

  ; close file
  free_lun, lun


END



; ---------------------------------------------------------------------------------------------------------------------------



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


; ---------------------------------------------------------------------------------------------------------------------------


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


; ---------------------------------------------------------------------------------------------------------------------------


FUNCTION flame_wavecal_skylines, x=x, y=y, $
	wavelength_axis_guess=wavelength_axis_guess, wavecal_settings=wavecal_settings, $
	skylines_pixel=skylines_pixel, skylines_wavelength=skylines_wavelength, extra_plot=extra_plot

	;
	; Given the x coordinate in pixels and the sky spectrum, it finds an accurate wavelength solution
	; by fitting the sky emission lines and comparing them to a line list. It outputs the wavelength axis.
	; It needs an approximate wavelength solution to identify the sky emission lines.
	; It also outputs the coordinates of the sky lines, both in pixels and in wavelength units
	;

		; these arrays will contain the fitted position of the sky lines
		skylines_wavelength = []
		skylines_pixel = []
		skylines_chisq = []

		; identify the OH lines that are in this wavelength range
		w_lines = where(wavecal_settings.line_list GT min(wavelength_axis_guess, /nan) $
			AND wavecal_settings.line_list LT max(wavelength_axis_guess, /nan), /null )

		; make sure there are OH lines here
		if w_lines EQ !NULL then message, 'Wavelength range does not contain OH lines?!'

		; keep only the OH lines of interest
		line_list = wavecal_settings.line_list[w_lines]

		; estimate width of OH lines
		approximate_linewidth_A = median(wavelength_axis_guess) / wavecal_settings.instrument_resolution
		approximate_linewidth_pix = approximate_linewidth_A / ( wavelength_axis_guess[3] - wavelength_axis_guess[2] )
			
		; fit a Gaussian to every sky line
		for i_line=0,n_elements(line_list)-1 do begin

			; select the region to fit
			w_fit = where( abs(wavelength_axis_guess-line_list[i_line]) LT 4.0*approximate_linewidth_A, /null )

			; check that the region is within the observed range
			if w_fit eq !NULL then continue

			; error handling for the gaussian fitting
			catch, error_gaussfit
			if error_gaussfit ne 0 then begin
				print, 'GAUSSFIT ERROR STATUS: ' + strtrim(error_gaussfit,2)
				catch, /cancel
				continue
			endif 

			; estimate parameters of the Gaussian
			est_peak = max( median( y[w_fit], 3) , /nan)
			est_center = w_fit[ n_elements(w_fit)/2 ]
			est_sigma = 0.5* approximate_linewidth_pix
			est_cont = min( median( y[w_fit], 3) , /nan)

			; Gaussian fit
			junk = gaussfit( x[w_fit], y[w_fit], gauss_param, nterms=4, $
				estimates=[est_peak, est_center, est_sigma, est_cont], sigma=gauss_err, chisq=chisq )

			; check that chi square makes sense
			if ~finite(chisq) then continue

			; check that the peak of the Gaussian is positive
			if gauss_param[0] LT 0.0 then continue

			; check that the SNR is high
			if gauss_param[0] LT 5.0*gauss_err[0] then continue

			; check that the center of the Guassian is in the observed range
			if gauss_param[1] LT min(x[w_fit]) or gauss_param[1] GT max(x[w_fit]) then continue

			; check that the Gaussian width makes sense
			if gauss_param[2] LT 0.5 or gauss_param[2] GT 0.5*n_elements(w_fit) then continue

			skylines_wavelength = [ skylines_wavelength, line_list[i_line] ]
			skylines_pixel = [ skylines_pixel, gauss_param[1] ]
			skylines_chisq = [ skylines_chisq, chisq ]

		endfor
		
		; if too few lines were found, then no reliable wavelength solution exists
		if n_elements(skylines_pixel) LT wavecal_settings.Nmin_lines then return, replicate(!values.d_nan, n_elements(x))

		; fit a polynomial to the skyline positions
		wavesol_coeff = poly_fit( skylines_pixel, skylines_wavelength, wavecal_settings.poly_degree )
		
		; calculate polynomial solution
		poly_wl = poly(x, wavesol_coeff)

		; properly handle regions with no information
		if where( ~finite(y), /null ) NE !NULL then begin	; check if there are NaNs in the input sky spectrum
			nosky_regions = label_region( ~finite(y) )		; this will have zero everywhere and N in the Nth "region" of NaNs
			nosky_regions[0] = nosky_regions[1]			; boundary issue
			nosky_regions[-1] = nosky_regions[-2]		; boundary issue
			if nosky_regions[1] ne 0 then poly_wl[ where(nosky_regions eq nosky_regions[1]) ] = !values.d_NaN	; void the region at the beginning
			if nosky_regions[-2] ne 0 then poly_wl[ where(nosky_regions eq nosky_regions[-2]) ] = !values.d_NaN	; void the region at the end
		endif

		; plot the spectrum
		erase
		cgplot, x, y, charsize=1, xsty=1, xtit='pixel coordinate', ytit='sky flux', layout=[1,2,1], _extra=extra_plot
		for i_line=0, n_elements(skylines_pixel)-1 do cgplot, skylines_pixel[i_line] + [0,0], [-1,1], /overplot, color='red'

		; show the result of Gaussian fitting
		cgplot, skylines_pixel, skylines_wavelength, /ynozero, xra=[x[0], x[-1]], xsty=1, psym=16, color='red', symsize=0.7, $
			xtit='pixel coordinate', ytit='expected wavelength', charsize=1, layout=[1,2,2]

		; show the polynomial fit
		cgplot, x, poly_wl, color='blue', /overplot

		; and here is our wavelength calibration!
		return, poly_wl

END


; ---------------------------------------------------------------------------------------------------------------------------

PRO flame_wavecal_accurate, slit_filename=slit_filename, $
	 wavecal_settings=wavecal_settings, approx_lambda_axis=approx_lambda_axis, $
	 OH_lines=OH_lines, wavelength_solution=wavelength_solution

	print, ' '
	print, 'Wavelength solution for ', slit_filename
	print, '*************************************************************************************************'
	print, ' '

	cgPS_open, flame_util_replace_string(slit_filename, '.fits', '_wavecal.ps'), /nomatch


	; load the model sky spectrum
	;-----------------------------------------------------------------------------------------------------------

	model_lambda = wavecal_settings.model_lambda
	model_flux = wavecal_settings.model_flux

	; load the slit image
	;-----------------------------------------------------------------------------------------------------------

	; read in slit
	im = readfits(slit_filename, hdr)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(im))[2]

	; how many pixels on the wavelength direction
	N_lambda_pix = (size(im))[1]

	; make a fine adjustment to the approximate wavelength solution via crosscorrelation using the central spectrum
	;--------------------------------------------------------------------------------------------------------------

	; sky spectrum: extract along the central 5 pixels
	sky = median(im[*, N_spatial_pix/2-3 : N_spatial_pix/2+2], dimension=2)

	; create the x-axis in pixel coordinates
	pix_axis = dindgen( n_elements(sky) )
	
	; trim edge to avoid problems
	sky[0:10] = 0.
	sky[-11:-1] = 0.

	; normalize the observed spectrum
	sky /= max(sky, /nan)

	; get rid of NaNs, which create problems for the cross-correlation
	sky[where(~finite(sky), /null)] = 0

	; set the size of the grid
	N_pix_scale = 20
	N_pix_scale_variation = 20

	; make the grid
	; for the pixel scale, bracket the value found in the approximate fit, +/- 15% of its value
	i_mid = n_elements(approx_lambda_axis)/2
 	approx_pixel_scale = approx_lambda_axis[i_mid+1] - approx_lambda_axis[i_mid]
	pix_scale_grid = approx_pixel_scale + 0.15*approx_pixel_scale * (-0.5 + dindgen(N_pix_scale)/double(N_pix_scale-1) )

	; the other parameter is the relative variation in pixel scale
	approx_pix_scale_variation = (approx_lambda_axis[i_mid+1]-approx_lambda_axis[i_mid]) / (approx_lambda_axis[i_mid] - approx_lambda_axis[i_mid-1]) - 1d
	pix_scale_variation_grid = approx_pix_scale_variation + 0.15*approx_pix_scale_variation * (-0.5 + dindgen(N_pix_scale_variation)/double(N_pix_scale_variation-1) )

	; cross-correlate the observed and model sky spectra and find a good wavelength solution
	flame_wavecal_crosscorr, sky=sky, model_lambda=model_lambda, model_flux=model_flux, $
		 approx_lambda_central=median(approx_lambda_axis), pix_scale_grid=pix_scale_grid, pix_scale_variation_grid=pix_scale_variation_grid, $
		 lambda_axis=lambda_axis_central_row, title=(strsplit(slit_filename,'/', /extract))[-1]

	; now fit all the pixel rows
	;--------------------------------------------------------------------------------------------------------------

	; create the 2D array that will contain the wavelength value for each pixel
	wavelength_solution = im
	wavelength_solution[*] = 0.

	; also, create 1D arrays that will contain the information for each OH line found
	OH_xpixel = []
	OH_wavelength = []
	OH_ypixel = []

	; start from the central row and go up until the top row, then start from center and go down
	row_number = indgen(N_spatial_pix)
	sorted_rows = [ row_number[N_spatial_pix/2: N_spatial_pix-1] , reverse(row_number[0:N_spatial_pix/2-1]) ]

	; as initial guess for the wavelength axis, use what we found for the central rows
	wavelength_axis_guess = lambda_axis_central_row

	print, 'Fitting individual sky lines for every pixel row...'

	; loop through all the pixel rows
	for counter=0, N_spatial_pix-1 do begin

		; index of the row we are considering now
		i_row = sorted_rows[counter]

		; skip 3 pixels at the edges to avoid the noisy part
		if i_row LT 3 or N_spatial_pix-1-i_row LT 3 then continue

		; print info on the row
		print, 'row ' + strtrim(i_row, 2) + ' ', format='(a,$)'

		; extract this pixel row from the slit
		this_row = im[*, i_row]

		; if needed, to gain signal, extract this row and sum it to the four neighbor rows
		if wavecal_settings.use_more_pixelrows then $
			this_row = total(im[*, i_row-2:i_row+2], 2)

		; normalize
		this_row /= max(this_row, /nan)

		; fit the emission lines and find the wavelength solution
		wavelength_axis_for_this_row = flame_wavecal_skylines( x=pix_axis, y=this_row, $
			wavelength_axis_guess=wavelength_axis_guess, wavecal_settings=wavecal_settings, $
			skylines_pixel=skylines_pixel, skylines_wavelength=skylines_wavelength, $
			extra_plot={title:'row ' + strtrim(i_row,2)} )

		; save the wavelength of all the OH lines found
		if n_elements(skylines_pixel) GT 0 then begin
			OH_xpixel = [OH_xpixel, skylines_pixel]
			OH_wavelength = [OH_wavelength, skylines_wavelength]
			OH_ypixel = [OH_ypixel, replicate(i_row, n_elements(skylines_pixel))]
		endif

		; save the wavelength solution
		wavelength_solution[*, i_row] = wavelength_axis_for_this_row

		; if a solution was found, then use it as the initial guess for the next row
		if where(finite(wavelength_axis_for_this_row), /null) NE !NULL then $
			wavelength_axis_guess = wavelength_axis_for_this_row

		; but check that we are not at the top edge of the slit: in that case use the central row
		if i_row eq max(sorted_rows)-3 then wavelength_axis_guess = lambda_axis_central_row

	endfor

	print, ''

	; take the median wavelength solution as a reference
	wavelength_solution_reference = median(wavelength_solution, dimension=2)

	; calculate the typical shift in wavelength for each row from the reference solution. Take the median wl of the central 5 pixels
	wavelength_shift = dblarr(N_spatial_pix)
	for i_row=0, N_spatial_pix-1 do $
		wavelength_shift[i_row] = $
			median(wavelength_solution[ N_lambda_pix/2-2 : N_lambda_pix/2+2, i_row]) - $
			median(wavelength_solution_reference[ N_lambda_pix/2-2 : N_lambda_pix/2+2 ] )

	; the 3 pixels at each edge are not being fit
	wavelength_shift[0:2] = 0.0
	wavelength_shift[-3:-1] = 0.0

	; plot the shift as a function of vertical position
	cgplot, 1d4*wavelength_shift, psym=-16, thick=3, charsize=1, $
		xtit = 'pixel position along the vertical (spatial) axis', ytit='wavelength shift from the reference pixel row (angstrom)'
	cgplot, [-1d4, 1d4], [0,0], /overplot, thick=2

	; plot the individual detections on a 2D view of the slit
	erase
	cgplot, OH_xpixel, OH_ypixel, psym=16, xtit='x pixel', ytitle='y pixel on this slit', $
		title='OH line detections', charsize=1, layout=[1,2,1], symsize=0.5

	cgPS_close


	; output the coordinates of the OH lines 
	OH_lines = [ [OH_wavelength], [OH_xpixel], [OH_ypixel] ]

	; write a ds9 region file with the identified OH lines 
	flame_wavecal_writeds9, OH_lines, filename =  $
		flame_util_replace_string(slit_filename, '.fits', '_OHlines.reg')


END


; ---------------------------------------------------------------------------------------------------------------------------


PRO flame_wavecal_approximate, slit_filename=slit_filename, this_slit=this_slit, $
	wavecal_settings=wavecal_settings, approx_lambda_axis=approx_lambda_axis
	;
	; Two steps:
	; First, use a coarse, logarithmic grid with constant pixel scale (uniform wavelength solution)
	; Second, use a fine grid centered on the coarse value of pixel scale, and explore also the 
	; pix_scale_variation parameter for non-uniform wavelength solution
	;

	; load the model sky spectrum
	;---------------------

	model_lambda = wavecal_settings.model_lambda
	model_flux = wavecal_settings.model_flux
	
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

	; sky spectrum: extract along the central 5 pixels 			; NB: NEED TO MAKE SURE THAT THIS IS, INDEED, SKY, AND NOT THE OBJECT
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
	filename_pieces = strsplit(slit_filename, '/', /extract) ; necessary for finding the slit directory
	filename_pieces[-1] = 'wavelength_solution_estimate.ps'
	ps_filename =  strjoin(filename_pieces, '/')
	cgPS_open, ps_filename, /nomatch


	print, ''
	print, 'FIRST STEP: rough estimate of lambda and delta_lambda'
	print, '-----------------------------------------------------'

	; let's smooth observed and model sky so that we can get an approximate match
	sky_sm = gauss_smooth(sky, wavecal_settings.smoothing_length)
	model_flux_sm = gauss_smooth(model_flux, 3)

	; remove continuum variation from the observed sky (important in the K band)
	sky_sm -= median(sky, 100)

	; renormalize spectra
	sky_sm /= max(sky_sm)
	model_flux_sm /= max(model_flux_sm)

	; estimated pixel scale from the approximate wavelength range
	estimated_pix_scale = (this_slit.approx_wavelength_hi-this_slit.approx_wavelength_lo) $
		/ double(n_elements(sky_sm))

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

	; return the approximate wavelength axis
	approx_lambda_axis = lambda_axis

	cgPS_close

END




; ---------------------------------------------------------------------------------------------------------------------------


PRO flame_wavecal_output_grid, wavelength_solution=wavelength_solution, $
	OH_lines=OH_lines, slit=slit

	; fit a 3rd degree polynomial to the 2D wavelength solution array in order to smooth it
	; 1. this could be more robust
	; 2. Is there a better way to set the output wavelength grid?


	; how many pixels on the spatial direction
	N_spatial_pix = (size(wavelength_solution))[2]

	; how many pixels on the wavelength direction
	N_lambda_pix = (size(wavelength_solution))[1]

	; unpack the OH line detection array
	OH_wavelength = OH_lines[*,0]
	OH_xpixel = OH_lines[*,1]
	OH_ypixel = OH_lines[*,2]

	; guess starting parameters for smooth solution
	start_params = [ min(wavelength_solution, /nan), $
		(max(wavelength_solution, /nan) - min(wavelength_solution, /nan)) / double(N_lambda_pix) , $
		(median(wavelength_solution[*,-5:-1]) - median(wavelength_solution[*,0:4])) / double(N_spatial_pix) , $
		1d-5, $
		1d-5, $
		1d-5]

	; fit a 3rd degree polynomial to all the OH lines found at each pixel row 
	fit_params = mpfit2dfun('flame_poly_surface', $
		OH_xpixel, OH_ypixel, OH_wavelength, replicate(1.0, n_elements(OH_xpixel)), start_params, /quiet)

	; generate the smooth wavelength solution using the best-fit parameters
	x_coordinate = dindgen(N_lambda_pix) # replicate(1.0, N_spatial_pix)
	y_coordinate = replicate(1.0, N_lambda_pix) # dindgen(N_spatial_pix)
	wavelength_solution_smooth = flame_poly_surface(x_coordinate, y_coordinate, fit_params)

	; find the wavelength range
	lambda_min_smooth = min(wavelength_solution_smooth, /nan)
	lambda_max_smooth = max(wavelength_solution_smooth, /nan)

	; find the median delta lambda
 	diff_lambda_smooth = abs( wavelength_solution_smooth - shift(wavelength_solution_smooth,1) )
 	diff_lambda_smooth = diff_lambda_smooth[5:-5,*]	; avoid edge effects
 	lambda_delta_smooth = median(diff_lambda_smooth)

 	; find the rounded value for the pixel scale, on a logarithimc scale
 	; (the idea here is to try to have the same value for slightly different datasets)
 	lambda_delta_out = 10.0^( round(alog10(lambda_delta_smooth)*10.0)/10.0 )

 	; define lambda range, being a little conservative
 	lambda_min = lambda_min_smooth - 10.0*lambda_delta_out
 	lambda_max = lambda_max_smooth + 10.0*lambda_delta_out

 	; save output grid to the slit structure
 	slit.outlambda_min = 10.0^( floor(alog10(lambda_min)*100.0)/100.0 )
	slit.outlambda_delta = lambda_delta_out
 	slit.outlambda_Npix = round( (lambda_max - lambda_min) / lambda_delta_out + 0.5 )

END 


; ---------------------------------------------------------------------------------------------------------------------------


PRO flame_wavecal_2D_calibration, filename=filename, slit=slit, OH_lines=OH_lines

; This routine calculates the 2D wavelength solution and y-rectification.
; These are two mappings from the observed pixel coordinates to the rectified grid
; lambdax, gamma, where lambdax is a pixel grid linear in lambda and gamma is the 
; vertical distance to the edge of the slit (taking warping into account)
; The result of this routine is a pair of matrices, Klambda and Kgamma, saved in fuel.slits,
; that can be used to rectify the image via poly_2D()
;

	; polynomial degree for image warping
	degree=3

	; read in file to calibrate
	im = mrdfits(filename, 0, header)

	; read dimensions of the image
	N_imx = (size(im))[1]
	N_imy = (size(im))[2]

	; find the minimum y value for the bottom edge of the slit
	bottom_edge = poly(indgen(N_imx), slit.bottom_poly)
	ymin_edge = min(bottom_edge)

	; this is the y-coordinate of the bottom pixel row in the cutout
	first_pixel =  ceil(ymin_edge)	

	; OH line coordinates
	OH_lambda = OH_lines[*,0]
	OH_x = OH_lines[*,1]
	OH_y = OH_lines[*,2]

	; output lambda axis
	lambda_0 = slit.outlambda_min
	delta_lambda = slit.outlambda_delta

	; translate every OH detection into the new coordinate system
	OH_lambdax = (OH_lambda - lambda_0)/delta_lambda
	OH_gamma = OH_y + first_pixel - poly(OH_x, slit.bottom_poly)

	; indices of the pixels we want to use - start with using all of them
	w_goodpix = lindgen(n_elements(OH_x))
	Ngoodpix = n_elements(w_goodpix)+1

	; loops are used to throw away outliers and make polywarp more robust
	WHILE n_elements(w_goodpix) LT Ngoodpix DO BEGIN

		; save old number of good pixels
		Ngoodpix = n_elements(w_goodpix)

		; calculate transformation Klambda,Kgamma from (x,y) to (lambda, gamma)
		polywarp, OH_lambdax[w_goodpix], OH_gamma[w_goodpix], $
			OH_x[w_goodpix], OH_Y[w_goodpix], degree, Klambda, Kgamma, /double, status=status

		; check that polywarp worked
		if status NE 0 then message, 'polywarp did not find a good solution'

		; calculate the model lambda given x,y where x,y are arrays
		lambda_modelx = fltarr(n_elements(OH_x))
		for i=0,degree do for j=0,degree do lambda_modelx +=  Klambda[i,j] * OH_x^j * OH_y^i
		lambda_model = lambda_0 + lambda_modelx*delta_lambda

		discrepancy = OH_lambda - lambda_model
		;w_outliers = where( abs(discrepancy/delta_lambda) GT 1.5, complement=w_goodpix, /null)
		w_outliers = where( abs(discrepancy) GT 2.0*stddev(discrepancy), complement=w_goodpix, /null)
		print, 'Outliers found: ', n_elements(w_outliers)

	ENDWHILE

	; now find the inverse transformation, using only the good pixels
	; calculate transformation Kx,Ky from (lambda, gamma) to (x,y)
	polywarp, OH_x[w_goodpix], OH_y[w_goodpix], $
		OH_lambdax[w_goodpix], OH_gamma[w_goodpix], degree, Kx, Ky, /double, status=status
	
	; save into slit structure
	this_rectification = {Klambda:Klambda, Kgamma:Kgamma, Kx:Kx, Ky:Ky}
	*slit.rectification = [ *slit.rectification, this_rectification]

END

; ---------------------------------------------------------------------------------------------------------------------------


PRO flame_wavecal_init, fuel=fuel, wavecal_settings=wavecal_settings
;
; initialize settings and load things
; outputs everything into the wavecal_settings structure
;

	; load line list
	readcol, fuel.util.linelist_filename, line_list

	; approximate value of R needed to estimate the width of OH lines (assuming one arcsec slit width!)
	instrument_resolution = (*fuel.instrument).resolution_slit1arcsec

	; the degree of the polynomial used to describe the wavelength solution
	poly_degree = 3
	
	; minimum number of OH lines for a reliable wavelength solution
	Nmin_lines = 6 
	
	; normally, extract individual pixel rows and detect OH lines
	use_more_pixelrows = 0

	; read in sky model
	readcol, fuel.util.sky_emission_filename, model_lambda, model_flux	; lambda in micron
	
	; create the wavecal_settings structure
	wavecal_settings = { $
		line_list : line_list, $
		instrument_resolution : instrument_resolution, $
		poly_degree : poly_degree, $
		Nmin_lines : Nmin_lines, $
		model_lambda : model_lambda, $
		model_flux : model_flux, $
		smoothing_length : fuel.input.wavecal_approx_smooth, $
		use_more_pixelrows : use_more_pixelrows $
	}



END


; ---------------------------------------------------------------------------------------------------------------------------


PRO flame_wavecal, fuel=fuel, verbose=verbose

	start_time = systime(/seconds)

	; extract the slits structures
	slits = *fuel.slits

	; avoid printing too much stuff on the terminal
	quiet_state = !QUIET
	if keyword_set(verbose) then !QUIET = 0 else !QUIET = 1

	; initialization: create wavecal_settings structure
	flame_wavecal_init, fuel=fuel, wavecal_settings=wavecal_settings

	; loop through all slits
	for i_slit=0, n_elements(slits)-1 do begin 

		this_slit = (*fuel.slits)[i_slit]

		print, '**********************************************************************'
		print, 'Wavecal for slit ', strtrim(this_slit.number,2), ' - ', this_slit.name
		print, '**********************************************************************'

		; take the first frame
		reference_filename = (*this_slit.filenames)[0]

		print, 'Using the central pixel rows of ', reference_filename

		flame_wavecal_approximate, slit_filename=reference_filename, this_slit=this_slit, $
			wavecal_settings=wavecal_settings, approx_lambda_axis = approx_lambda_axis

		; make sure there are no old rectification coefficients left over
		this_slit.rectification = ptr_new(/allocate_heap)

		for i_frame=0, n_elements(*slits[i_slit].filenames)-1 do begin

			flame_wavecal_accurate, slit_filename=(*slits[i_slit].filenames)[i_frame], $
				wavecal_settings=wavecal_settings, approx_lambda_axis=approx_lambda_axis, $
				OH_lines=OH_lines, wavelength_solution=wavelength_solution

			flame_wavecal_output_grid, wavelength_solution=wavelength_solution, $
				OH_lines=OH_lines, slit=this_slit

			flame_wavecal_2D_calibration, filename=(*slits[i_slit].filenames)[i_frame], $
				slit=this_slit, OH_lines=OH_lines

			; update the slit structure with the output wavelength grid of the last frame
			(*fuel.slits)[i_slit] = this_slit

		endfor

	endfor


	; revert to original !QUIET state
	!QUIET = quiet_state

	print, 'It took ', systime(/seconds) - start_time, ' seconds'

END
