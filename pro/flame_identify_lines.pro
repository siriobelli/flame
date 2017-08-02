;
; For each slit & frame, extract the spectrum of one pixel row, starting
; at the center and assuming the rough calibration. Identify and fit
; the sky emission lines for each row.
;





;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_identify_writeds9, speclines, filename=filename
;
; write a ds9 region file with all the OH line detections
;

  ; extract the wavelength of all identifications
  line_lambdas = speclines.lambda
  uniq_lambdas = line_lambdas[UNIQ(line_lambdas, SORT(line_lambdas))]

  ; open file
  openw, lun, filename, /get_lun

  ; write header
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=1 include=1 source=1'
  printf, lun, 'image'

  for i_line=0, n_elements(uniq_lambdas)-1 do begin

  	; the wavelength of this OH line
  	this_lambda = uniq_lambdas[i_line]

  	; select the x and y coordinates for this line
  	w_thisline = where(speclines.lambda eq this_lambda)
  	this_x = speclines[w_thisline].x
 		this_y = speclines[w_thisline].y

   	; sort points by y coordinate
   	wsort = sort(this_y)
   	this_x = this_x[wsort]
   	this_y = this_y[wsort]

    ; concatenate points
    all_x = [ this_x, reverse(this_x) ]
    all_y = [ this_y, reverse(this_y) ]

		; in ds9, the first pixel is (1,1), not (0,0)
		all_x += 1
		all_y += 1

    ; make the string with all the points
    all_points = ''
    for i=0,n_elements(all_x)-2 do all_points += strtrim(all_x[i],2) + ',' + cgnumber_formatter(all_y[i], decimals=1) + ','
    ; add the last two points without the final comma
    all_points += strtrim(all_x[-1],2) + ',' + cgnumber_formatter(all_y[-1], decimals=1)

    ; the color for this line (red if we are using it for wavecal, otherwise black)
    if speclines[w_thisline[0]].trust_lambda eq 1 then color ='red' $
      else color='black'

    ; write the region corresponding to this line
    printf, lun, 'polygon(' + all_points + ') # text={' + cgnumber_formatter(this_lambda, decimals=5) + '} color = ' + color

  endfor

  ; close file
  free_lun, lun


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_identify_fitskylines, fuel=fuel, x=x, y=y, $
	approx_wavecal=approx_wavecal, linewidth=linewidth, $
	line_list=line_list_in, line_trust=line_trust_in, $
	speclines=speclines, wavecal=wavecal, plot_title=plot_title

	;
	; Given a 1D sky spectrum, it finds an accurate wavelength solution
	; by fitting the sky emission lines and comparing them to a line list. It outputs
	; the wavelength axis and the Gaussian parameters of each sky line.
	; It needs an approximate wavelength solution and linewidth to identify the sky emission lines.
	;
	;
	; x: (input) array with the pixel x-coordinate
	; y: (input) array with the observed sky in one pixel row
	; approx_wavecal: (input) array with the approximate wavelength for each x
	; linewidth: (input and output) approximate sigma of unresolved sky lines, which will be updated. NB:in pixels
	; line_list: (input) array of expected wavelengths of usable OH lines
  ; line_trust: (input) array of flags: 1 if the line can be used for the wavelength calibration. If not specified, use all lines
	; speclines: (output) array of structures with parameters for each OH line
	; wavecal: (output) array with the wavelength solution
	; plot_title : (input) string to print as title of the plot
	;

  ; if line_trust is not provided, then trust all lines
  if ~keyword_set(line_trust_in) then $
    line_trust_in = replicate(1, n_elements(line_list_in))

	; convert linewidth to micron
	linewidth_um = linewidth * median( approx_wavecal - shift(approx_wavecal, 1) )

	; identify the OH lines that are in this wavelength range
	w_lines = where(line_list_in GT min(approx_wavecal, /nan) $
		AND line_list_in LT max(approx_wavecal, /nan), /null )

	; make sure there are OH lines here
	if w_lines EQ !NULL then begin
    print, 'Warning: wavelength range does not contain OH lines'
    speclines = !NULL
    wavecal = !NULL
    return
  endif

	; keep only the OH lines of interest
	line_list = line_list_in[w_lines]
  line_trust = line_trust_in[w_lines]

  ; check that there are at least three trustable lines
  if n_elements(where(line_trust eq 1, /null)) LT 3 then $
    message, 'linelist contains less than three lines that can be used for wavelength calibration'

	; make the array that will contain the result of the fitting
	speclines = []

	; fit a Gaussian to every sky line
	for i_line=0,n_elements(line_list)-1 do begin

		; select the region to fit
		w_fit = where( abs(approx_wavecal-line_list[i_line]) LT 0.5*fuel.settings.identify_lines_linefit_window*linewidth_um, /null )

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
		est_sigma = linewidth
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
		if gauss_param[2] LT linewidth/10.0 or gauss_param[2] GT linewidth*10.0 then continue

		; make OHline structure
		this_OHline = { lambda: line_list[i_line], $
			x: gauss_param[1], $
			y: -1.0, $
			sigma: gauss_param[2], $
			peak: gauss_param[0], $
			chisq: chisq, $
      trust_lambda:line_trust[i_line] }

		; add to the stack
		speclines = [speclines, this_OHline]

	endfor

  ; did we find any speclines at all?
  if n_elements(speclines) eq 0 then return

  ; select only the lines we can use for the wavelength calibration
  speclines_trust = speclines[where(speclines.trust_lambda eq 1, /null)]
  speclines_donttrust = speclines[where(speclines.trust_lambda eq 0, /null)]

	; if too few lines were found, then no reliable wavelength solution exists
	if n_elements(speclines_trust) LT fuel.settings.identify_lines_Nmin_lines then begin
		speclines = !NULL
		return
	endif

  ; set the degree for the polynomial fit - if there are few lines, decrease the degree
  poly_degree = fuel.settings.identify_lines_poly_degree
  if poly_degree GT (n_elements(speclines_trust)+1)/3 then poly_degree = (n_elements(speclines_trust)+1)/3

  ; fit a polynomial to the skyline positions using only the lines we can trust
  wavesol_coeff = robust_poly_fit( speclines_trust.x, speclines_trust.lambda, poly_degree, /DOUBLE )

  ; ; calculate residuals
  ; residuals = speclines.lambda-poly(speclines.x, wavesol_coeff)
  ; sigma_res = stddev(residuals, /nan)
  ; w_outliers = where(abs(residuals) GT 3.0*sigma_res, /null)
  ; if n_elements(w_outliers) GT 0 then print, speclines[w_outliers].lambda

	; calculate polynomial solution
	poly_wl = poly(x, wavesol_coeff)

	; charsize
	ch = 0.8

	; panel 1: plot the spectrum
	erase
	cgplot, x, y, charsize=ch, xsty=1, xtit='', ytit='observed flux', title=plot_title, $
		position = [0.15, 0.70, 0.95, 0.96], xtickformat="(A1)", xra=[x[0], x[-1]], /nodata

	; show the OH lines that were identified
	for i_line=0, n_elements(speclines_trust)-1 do cgplot, speclines_trust[i_line].x + [0,0], [-2,2]*max(abs(y)), /overplot, color='red'
  for i_line=0, n_elements(speclines_donttrust)-1 do cgplot, speclines_donttrust[i_line].x + [0,0], [-2,2]*max(abs(y)), /overplot, color='blk4'

	; show the spectrum on top, for clarity
	cgplot, x, y, /overplot

	; panel 2: show the wavelength solution
	cgplot, speclines_trust.x, speclines_trust.lambda, /ynozero, xra=[x[0], x[-1]], xsty=1, psym=16, color='red', symsize=0.7, $
		xtit='', ytit='expected wavelength', charsize=ch, $
		/noerase, position = [0.15, 0.50, 0.95, 0.70], xtickformat="(A1)"

  if n_elements(speclines_donttrust) GT 0 then $
    cgplot, speclines_donttrust.x, speclines_donttrust.lambda, /overplot, psym=16, color='blk4', symsize=0.7

	; show the polynomial fit
	cgplot, x, poly_wl, color='blue', /overplot

	; panel 3: show the residuals
	cgplot, speclines_trust.x, 1d4 * (speclines_trust.lambda-poly(speclines_trust.x, wavesol_coeff)), /ynozero, xra=[x[0], x[-1]], xsty=1, psym=16, color='red', symsize=0.7, $
		ytit='residuals (' + string("305B) + ')', charsize=ch, $
		/noerase, position = [0.15, 0.30, 0.95, 0.50], xtickformat="(A1)"

  if n_elements(speclines_donttrust) GT 0 then $
    cgplot, speclines_donttrust.x, 1d4 * (speclines_donttrust.lambda-poly(speclines_donttrust.x, wavesol_coeff)), /overplot, psym=16, color='blk4', symsize=0.7

  cgplot, [x[0], x[-1]], [0,0], /overplot, thick=3, linestyle=2


	; panel 4: plot the line widths
	cgplot, speclines_trust.x, speclines_trust.sigma, /ynozero, xra=[x[0], x[-1]], xsty=1, psym=16, color='red', symsize=0.7, $
		xtit='pixel coordinate', ytit='line width (pixel)', charsize=ch, $
		/noerase, position = [0.15, 0.10, 0.95, 0.30]

  if n_elements(speclines_donttrust) GT 0 then $
    cgplot, speclines_donttrust.x, speclines_donttrust.sigma, /overplot, psym=16, color='blk4', symsize=0.7

	; show median value of line width
	cgplot, [x[0], x[-1]], [0,0]+median(speclines.sigma), /overplot, thick=3, linestyle=2

	; update value of linewidth
	linewidth = median(speclines.sigma)

	; output wavelength solution
	wavecal = poly_wl

END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_identify_find_speclines, fuel=fuel, filename=filename, $
	 slit=slit, rough_lambda=rough_lambda, rough_flux=rough_flux, linelist_filename=linelist_filename, $
	 speclines=speclines, wavelength_solution=wavelength_solution

	print, ' '
	print, 'Identifying lines in ', filename

	cgPS_open, flame_util_replace_string(filename, '.fits', '_speclines.ps'), /nomatch


	; load the slit image
	;-----------------------------------------------------------------------------------------------------------

	; read in slit
	im = readfits(filename, hdr)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(im))[2]

	; how many pixels on the wavelength direction
	N_lambda_pix = (size(im))[1]

	; create the x-axis in pixel coordinates
	pix_axis = dindgen( N_lambda_pix )


  ; find a zeroth-order solution starting with the rough_wavecal
	;--------------------------------------------------------------------------------------------------------------

  ; setup the order to use for the pixel rows:
  ; vertical coordinate of each row
	row_number = indgen(N_spatial_pix)

  ; check if we are stacking more than one pixel row
  if fuel.settings.identify_lines_stack_rows LT 0 then message, 'fuel.settings.identify_lines_stack_rows can only be positive!'
  N_stack = fuel.settings.identify_lines_stack_rows

  ; if stacking more than one pixel row, then cut the edges
  if N_stack ne 0 then $
    row_number = row_number[N_stack:-1-N_stack]

  ; number of rows that will be extracted
  N_rows = n_elements(row_number)

  ; start from the central row and go up until the top, then start from center and go down
	sorted_rows = [ row_number[N_rows/2: N_rows-1] , reverse(row_number[0:N_rows/2-1]) ]

	; identify first row of bottom half
	i0_bottom = row_number[N_rows/2-1]

  ; extract the sky spectrum from a bunch of central pixel rows
  central_skyspec = median( im[*, sorted_rows[0]:sorted_rows[2]], dimension=2 )

  ; median smooth it to make sure the lines look like Gaussians (especially for very tilted slits)
  central_skyspec = median(central_skyspec, 5)

  ; clean from NaNs because they don't work with the cross correlation
  central_skyspec_clean = central_skyspec
  central_skyspec_clean[ where(~finite(central_skyspec), /null) ] = 0.0

  ; take the reference spectrum from the rough wavecal
  ref_skyspec = median(rough_flux, 3)
  ref_skyspec[ where(~finite(ref_skyspec), /null) ] = 0.0

  ; measure the overall shift since the rough wavecal may have
  ; been obtained on different spatial positions or different frames
  lag = indgen(100)-50 ; up to 50 pixels each direction
  crosscorr = c_correlate( central_skyspec_clean, ref_skyspec, lag)
  max_crosscorr = max( crosscorr, max_ind, /nan)
  delta = -lag[max_ind]
  print, 'shifting the central pixel row by ' + strtrim(delta, 2) + ' pixels'

  ; apply the shift to the rough wavecal to obtain a good starting guess
  approx_lambda_axis = shift(rough_lambda, delta)
  if delta GT 0 then approx_lambda_axis[0:delta] = !values.d_nan
  if delta LT 0 then approx_lambda_axis[-1-abs(delta):-1] = !values.d_nan


  ; identify the spectral lines
	;--------------------------------------------------------------------------------------------------------------

  ; load line list
	readcol, linelist_filename, line_list, line_trust, format='D,I', /silent

  ; calcolate typical wavelength step of one pixel
  lambda_step = median( approx_lambda_axis - shift(approx_lambda_axis, 1) )

  ; approximate sky line width (assuming one arcsec slit width)
	approximate_linewidth_um = median(approx_lambda_axis) / (2.36 * fuel.instrument.resolution_slit1arcsec)
	linewidth = approximate_linewidth_um / lambda_step ; in pixel

  ; start with a larger linewidth, for a generous range where the line could be
  linewidth *= 2.0

  ; use the shifted rough wavelength calibration as starting solution
  lambda_axis = approx_lambda_axis

  ; number of identified lines
  Nlines = 0
  delta_Nlines = 1

  ; loop number
  i_loop = 0

  print, 'Identifying the spectral lines in the central rows...'

  ; re-iteratively identify lines and improve wavelength solution
  while delta_Nlines GT 0 and i_loop LT 10 do begin

    ; for this initial fits, use only the lines that we can trust
    line_list_trust = line_list[where(line_trust eq 1, /null)]

    ; fit the emission lines and find the wavelength solution
  	flame_identify_fitskylines, fuel=fuel, x=pix_axis, y=central_skyspec, $
  		approx_wavecal=lambda_axis, linewidth=linewidth, $
  		line_list=line_list_trust, $
  		speclines=speclines_thisloop, wavecal=lambda_axis_output, plot_title='central rows / ' + strtrim(i_loop,2)

    ; update the wavelength solution
    lambda_axis = lambda_axis_output

    ; compare the number of identified lines to previous loop
    delta_Nlines = n_elements(speclines_thisloop) - Nlines
    Nlines = n_elements(speclines_thisloop)
    i_loop++

  endwhile

  w_l = where(line_list_trust GT lambda_axis[4] and line_list_trust LT lambda_axis[-5], /null)
  print, 'Identified ', strtrim(Nlines, 2), ' out of ', strtrim(n_elements(w_l), 2), ' lines present in the line list.'


  ; fit all the pixel rows
	;--------------------------------------------------------------------------------------------------------------

	; create the 2D array that will contain the wavelength value for each pixel
	wavelength_solution = im
	wavelength_solution[*] = !values.d_nan

	; as initial guess for the wavelength axis, use what we just found
	wavelength_axis_guess = lambda_axis

	; create the empty array of speclines structures
	speclines = []

	print, 'Fitting individual sky lines for every pixel row...'

	; loop through all the pixel rows
	for counter=0, N_rows-1 do begin

		; index of the row we are considering now
		i_row = sorted_rows[counter]

		; print info on the row
		print, 'row ' + strtrim(i_row, 2) + ' ', format='(a,$)'

		; extract this pixel row from the slit (if needed, stack more than one row to increase SNR)
		if N_stack eq 0 then $
      this_row = im[*, i_row] else $
      this_row = median(im[*, i_row-N_stack:i_row+N_stack], dimension=2)

		; if this is the first row of the bottom half, then use the wavelength solution found for the first row
		if i_row eq i0_bottom then wavelength_axis_guess = wavelength_solution_0

		; fit the emission lines and find the wavelength solution
		flame_identify_fitskylines, fuel=fuel, x=pix_axis, y=this_row, $
			approx_wavecal=wavelength_axis_guess, linewidth=linewidth, $
			line_list=line_list, line_trust=line_trust, $
			speclines=speclines_thisrow, wavecal=wavelength_axis_for_this_row, plot_title='row '+strtrim(i_row,2)

		; if sky lines were not found, then skip to next row
		if n_elements(speclines_thisrow) EQ 0 then begin
      print, ' (lines not found) ', format='(a,$)'
      continue
    endif

		; set the y coordinate for the OH lines
		speclines_thisrow.y = i_row

		; save the OH lines from this row
		speclines = [ speclines, speclines_thisrow ]

		; save the wavelength solution
		wavelength_solution[*, i_row] = wavelength_axis_for_this_row

		; if a solution was found
		if where(finite(wavelength_axis_for_this_row), /null) NE !NULL then begin

      ; use it as the initial guess for the next row
			wavelength_axis_guess = wavelength_axis_for_this_row

      ; if not already done, save the "zero" solution for when doing the bottom half of the slit
      if n_elements(wavelength_solution_0) eq 0 then $
        wavelength_solution_0 = wavelength_axis_for_this_row

    endif

	endfor

	print, ''

  ; check that at least some lines were found
  if n_elements(speclines) eq 0 then message, 'no lines were identified in this cutout!'

	cgPS_close

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_identify_output_grid, wavelength_solution=wavelength_solution, slit=slit

;
; Set the wavelength grid that will be used for the rectified frame,
; in such a way that it is a good fit to the observed wavelength range
;

  ; smoothing length, in pixels
  beta = 5
	if (size(wavelength_solution))[2] LE 5 then beta = 3

  ; apply median filtering to the 2D wavelength solution
  wavelength_solution_smooth = median(wavelength_solution, beta)

  ; get rid of the edge values that cannot be properly smoothed
  wavelength_solution_smooth = wavelength_solution_smooth[beta:-1-beta, beta:-1-beta]

	; find the median lambda step
 	diff_lambda = abs( wavelength_solution_smooth - shift(wavelength_solution_smooth, 1, 0) )
 	lambda_delta = double(median(diff_lambda))

  ; find the median wavelength of the first and last pixel
  lambda_first_pixel = median(wavelength_solution_smooth[0,*])
  lambda_last_pixel = median(wavelength_solution_smooth[-1,*])

  ; number of pixels along the vertical direction
  N_pix_y = (size(wavelength_solution))[2]

	; determine the min and max values
  ; (go a few pixel beyond the values found above, to account for tilted slits)
	lambda_min = lambda_first_pixel - 0.5*N_pix_y*lambda_delta
	lambda_max = lambda_last_pixel + 0.5*N_pix_y*lambda_delta

 	; find the rounded values, on a logarithmic scale
 	; (the idea here is to try to have the same value for slightly different datasets)
  ; and save output grid to the slit structure
	slit.outlambda_min = 10.0^( floor(alog10(lambda_min)*100.0)/100.0 )
	slit.outlambda_delta = 10.0^( round(alog10(lambda_delta)*20.0d)/20.0d )
  slit.outlambda_Npix = round( (lambda_max - slit.outlambda_min) / slit.outlambda_delta )

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_identify_lines, fuel

	flame_util_module_start, fuel, 'flame_identify_lines'


  ; avoid printing too much stuff (especially from GAUSSFIT)
  quiet_state = !QUIET
  !QUIET = 1

	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

    if fuel.slits[i_slit].skip then continue

	  print, 'Finding emission lines for calibration for slit ', strtrim(fuel.slits[i_slit].number,2), ' - ', fuel.slits[i_slit].name
		print, ' '

		; handle errors by ignoring that slit
		if fuel.settings.debugging eq 0 then begin
  		catch, error_status
  		if error_status ne 0 then begin
  			print, ''
        print, ''
  	    print, '**************************'
  	    print, '***       WARNING      ***'
  	    print, '**************************'
        print, 'Error found. Skipping slit ' + strtrim(fuel.slits[i_slit].number,2), ' - ', fuel.slits[i_slit].name
        print, ''
        fuel.slits[i_slit].skip = 1
  			catch, /cancel
  			continue
  		endif
    endif


		; this slit
		this_slit = fuel.slits[i_slit]


    ; arcs ---------------------------------------------------------------------

		; filename of the cutout
    arc_filename = this_slit.arc_cutout.filename_step1

		; identify and measure the speclines
		flame_identify_find_speclines, fuel=fuel, filename=arc_filename, slit=this_slit, $
    rough_lambda=*this_slit.rough_arclambda, rough_flux=*this_slit.rough_arcflux, $
    linelist_filename=fuel.util.intermediate_dir + 'linelist_arcs.txt', $
			speclines=speclines, wavelength_solution=wavelength_solution

		; save the speclines in the slit structure
		*this_slit.arc_cutout.speclines = speclines

				; write a ds9 region file with the identified speclines
		flame_identify_writeds9, speclines, filename=flame_util_replace_string(arc_filename, '.fits', '_speclines.reg')


    ; skylines -----------------------------------------------------------------

		for i_frame=0, n_elements(fuel.slits[i_slit].cutouts)-1 do begin

				; filename of the cutout
				slit_filename = fuel.slits[i_slit].cutouts[i_frame].filename_step1

				; identify and measure the speclines
				flame_identify_find_speclines, fuel=fuel, filename=slit_filename, slit=this_slit, $
        rough_lambda=*this_slit.rough_skylambda, rough_flux=*this_slit.rough_skyflux, $
        linelist_filename=fuel.settings.linelist_filename, $
					speclines=speclines, wavelength_solution=wavelength_solution

				; save the speclines in the slit structure
				*this_slit.cutouts[i_frame].speclines = speclines

				; write a ds9 region file with the identified speclines
				flame_identify_writeds9, speclines, filename=flame_util_replace_string(slit_filename, '.fits', '_speclines.reg')

				; use the pixel-by-pixel wavelength solution OF THE FIRST FRAME to set the output grid in wavelength
				if i_frame eq 0 then begin
					flame_identify_output_grid, wavelength_solution=wavelength_solution, slit=this_slit
					fuel.slits[i_slit] = this_slit
				endif

		endfor

	endfor

	; revert to original !QUIET state
	!QUIET = quiet_state


  flame_util_module_end, fuel

END
