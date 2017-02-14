;
; Wavelength calibration using OH sky lines.
; For each slit, extract the spectrum of one pixel row, starting
; at the center and assuming the rough calibration. Identify and fit
; the sky emission lines for each row, and then find a polynomial
; warping that describes the 2D transformation from observed frame
; to lambda-calibrated and vertically-rectified frame.
;




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


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


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


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

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



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


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


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


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_wavecal_skylines, fuel=fuel, x=x, y=y, $
	wavelength_axis_guess=wavelength_axis_guess, $
	skylines_pixel=skylines_pixel, skylines_wavelength=skylines_wavelength, extra_plot=extra_plot

	;
	; Given the x coordinate in pixels and the sky spectrum, it finds an accurate wavelength solution
	; by fitting the sky emission lines and comparing them to a line list. It outputs the wavelength axis.
	; It needs an approximate wavelength solution to identify the sky emission lines.
	; It also outputs the coordinates of the sky lines, both in pixels and in wavelength units
	;

  ; settings:
	; the degree of the polynomial used to describe the wavelength solution
	poly_degree = 3

	; minimum number of OH lines for a reliable wavelength solution
	Nmin_lines = 6

	; approximate value of R needed to estimate the width of OH lines (assuming one arcsec slit width!)
  instrument_resolution = fuel.instrument.resolution_slit1arcsec

	; load line list
	readcol, fuel.util.linelist_filename, line_list, format='D', /silent

	; these arrays will contain the fitted position of the sky lines
	skylines_wavelength = []
	skylines_pixel = []
	skylines_chisq = []

	; identify the OH lines that are in this wavelength range
	w_lines = where(line_list GT min(wavelength_axis_guess, /nan) $
		AND line_list LT max(wavelength_axis_guess, /nan), /null )

	; make sure there are OH lines here
	if w_lines EQ !NULL then message, 'Wavelength range does not contain OH lines?!'

	; keep only the OH lines of interest
	line_list = line_list[w_lines]

	; estimate width of OH lines
	approximate_linewidth_A = median(wavelength_axis_guess) / instrument_resolution
	approximate_linewidth_pix = approximate_linewidth_A / ( wavelength_axis_guess[3] - wavelength_axis_guess[2] )

	; fit a Gaussian to every sky line
	for i_line=0,n_elements(line_list)-1 do begin

		; select the region to fit
		w_fit = where( abs(wavelength_axis_guess-line_list[i_line]) LT 4.0*approximate_linewidth_A, /null )

		; check that the region is within the observed range
		if w_fit eq !NULL then continue

    ; check that there actually is signal and it's not just a bunch of NaNs
    if n_elements( where( finite(y[w_fit]), /null ) ) LE 5 then continue

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
		!NULL = gaussfit( x[w_fit], y[w_fit], gauss_param, nterms=4, $
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
	if n_elements(skylines_pixel) LT Nmin_lines then return, replicate(!values.d_nan, n_elements(x))

	; fit a polynomial to the skyline positions
	wavesol_coeff = poly_fit( skylines_pixel, skylines_wavelength, poly_degree )

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


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_oneslit, fuel=fuel, slit_filename=slit_filename, $
	 approx_lambda_axis=approx_lambda_axis, $
	 OH_lines=OH_lines, wavelength_solution=wavelength_solution

	print, ' '
	print, 'Wavelength solution for ', slit_filename
	print, '*************************************************************************************************'
	print, ' '

	cgPS_open, flame_util_replace_string(slit_filename, '.fits', '_wavecal.ps'), /nomatch


	; load the slit image
	;-----------------------------------------------------------------------------------------------------------

	; read in slit
	im = readfits(slit_filename, hdr)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(im))[2]

	; how many pixels on the wavelength direction
	N_lambda_pix = (size(im))[1]

	; create the x-axis in pixel coordinates
	pix_axis = dindgen( N_lambda_pix )

	; fit all the pixel rows
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

	; as initial guess for the wavelength axis, use what found during the rough wavecal
	wavelength_axis_guess = approx_lambda_axis

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

		; normalize
		this_row /= max(this_row, /nan)

		; fit the emission lines and find the wavelength solution
		wavelength_axis_for_this_row = flame_wavecal_skylines( fuel=fuel, $
      x=pix_axis, y=this_row, $
			wavelength_axis_guess=wavelength_axis_guess, $
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
		if i_row eq max(sorted_rows)-3 then wavelength_axis_guess = approx_lambda_axis

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

	; get rid of wildly wrong rows
	wavelength_shift[where(abs(wavelength_shift - median(wavelength_shift)) $
		GT 5.0 *stddev(wavelength_shift, /nan), /null)] = !values.d_NaN

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


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_clean, slit=slit, index=i_slit
  ;
  ; clean cutout from cosmic rays
  ; cosmic rays are found by comparing the observed frame
  ; to the median of the frames in its vicinity
  ;

  ; how many frames are we using for the cosmic ray identification?
  Nframes = 3

  ; read in the cutout file names
  filenames = *slit.filenames

  ; select the adjacent ones
  indices = indgen(n_elements(filenames))

  ; distance from the index of interest
  distance = abs( indices - i_slit )

  ; sort indices by distance
  closest_indices = indices[sort(distance)]

  ; pick the closest Nframes (including the current one of course)
  w_touse = closest_indices[0:Nframes-1]

  ; read all the frames
  stack = []
  for i=0, Nframes-1 do stack = [ [[stack]] , [[ readfits(filenames[w_touse[i]]) ]] ]

  ; this is the cleaned frame (no CRs)
  clean = median(stack, dimension=3)

  ; now take the current frame
  current = readfits(filenames[i_slit], header)

  ; and subtract the cleaned version from it
  sub = current-clean

  ; calculate properties of sub frame
  median = median(sub)
  sigma = stddev(sub, /nan)   ; this is a conservative sigma, because it includes all the cosmic rays

  ; identify cosmic rays as positive fluctuations
  w_cr = where( sub-median GT 3.0*sigma, /null)

  ; make mask with flaggeed pixels
  mask = current
  mask[*] = 0.0
  mask[w_cr] = 1.0

  ; since the flagging is very conservative,
  ; let's grow the mask by 1 pixel along the four directions
  mask_grow = mask + shift(mask, [1,0]) + shift(mask, [-1,0]) + shift(mask, [0,1]) + shift(mask, [0,-1])
  w_tomask = where(mask_grow GT 0.0, /null)

  ; make clean cutout
  clean = current
  clean[w_tomask] = !values.d_NaN

  ; rename the original cutout for archival reasons
  file_move, filenames[i_slit], flame_util_replace_string(filenames[i_slit], '.fits', '_raw.fits'), /overwrite

  ; write clean cutout
  writefits, filenames[i_slit], clean, header


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************




PRO flame_wavecal_accurate, fuel=fuel

	start_time = systime(/seconds)

  print, ' '
  print, 'flame_wavecal_accurate'
  print, '**********************'
  print, ' '

  ; avoid printing too much stuff (especially from GAUSSFIT)
  quiet_state = !QUIET
  !QUIET = 1

	; extract the slits structures
	slits = fuel.slits

	; loop through all slits
	for i_slit=0, n_elements(slits)-1 do begin

		this_slit = fuel.slits[i_slit]

	  print, 'Accurate wavelength calibration for slit ', strtrim(this_slit.number,2), ' - ', this_slit.name
		print, ' '

		; make sure there are no old rectification coefficients left over
		this_slit.rectification = ptr_new(/allocate_heap)

		for i_frame=0, n_elements(*slits[i_slit].filenames)-1 do begin

      ; if needed, clean the cutout from cosmic rays
      if fuel.input.clean_individual_frames then $
        flame_wavecal_clean, slit=this_slit, index=i_frame

			flame_wavecal_oneslit, fuel=fuel, slit_filename=(*slits[i_slit].filenames)[i_frame], $
				approx_lambda_axis=*this_slit.rough_wavecal, $
				OH_lines=OH_lines, wavelength_solution=wavelength_solution

			flame_wavecal_output_grid, wavelength_solution=wavelength_solution, $
				OH_lines=OH_lines, slit=this_slit

			flame_wavecal_2D_calibration, filename=(*slits[i_slit].filenames)[i_frame], $
				slit=this_slit, OH_lines=OH_lines

			; update the slit structure with the output wavelength grid of the last frame
			fuel.slits[i_slit] = this_slit

		endfor

	endfor

	; revert to original !QUIET state
	!QUIET = quiet_state
	print, 'It took ', systime(/seconds) - start_time, ' seconds'

END
