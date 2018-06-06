;
; For each slit & frame, extract the spectrum of one pixel row, starting
; at the center and assuming the rough calibration. Identify and fit
; the sky emission lines for each row.
;





;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_findlines_writeds9, speclines, filename=filename
;
; write a ds9 region file with all the emission line detections
;

  ; extract the wavelength of all identifications
  line_lambdas = speclines.lambda
  uniq_lambdas = line_lambdas[UNIQ(line_lambdas, SORT(line_lambdas))]

  ; open file
  openw, lun, filename, /get_lun

  ; write header
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=1 include=1 source=1'
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

    ; the color for this line (red if we are using it for wavecal, otherwise blue)
    if speclines[w_thisline[0]].trust_lambda eq 1 then color ='red' $
      else color='blue'

		; in ds9, the first pixel is (1,1), not (0,0)
		this_x += 1
		this_y += 1

    ; for each detection make a point region
    for i=0, n_elements(this_x)-1 do $
      printf, lun, 'point(' + strtrim(this_x[i], 2) + ',' + strtrim(this_y[i], 2) + ') # point=cross color = ' + color

    ; for the top detection, make an extra region with the line wavelength in the text
    !NULL = max(this_y, w_top, /nan)
    printf, lun, 'point(' + strtrim(this_x[w_top], 2) + ',' + strtrim(this_y[w_top], 2) + $
      ') # point=cross color = ' + color + ' text={' + cgnumber_formatter(this_lambda, decimals=5) + '}'

  endfor

  ; close file
  free_lun, lun


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_findlines_fitskylines, x=x, y=y, $
	approx_wavecal=approx_wavecal, linewidth=linewidth, $
  reflines=reflines, check_shift=check_shift, Nmin_lines=Nmin_lines, fit_window=fit_window, $
  poly_degree=poly_degree, verbose=verbose, $
	speclines=speclines, wavecal=wavecal, plot_title=plot_title

	;
  ; Given a 1D sky spectrum, identify and measure a set of emission lines.
  ; Output the wavelength axis and the Gaussian parameters of each line.
	; It needs an approximate wavelength solution and linewidth and a list
  ; of emission line wavelengths.
	;
	; x: (input) array with the pixel x-coordinate
	; y: (input) array with the observed sky in one pixel row
	; approx_wavecal: (input) array with the approximate wavelength for each x
	; linewidth: (input and output) approximate sigma of unresolved sky lines, which will be updated. NB:in pixels
  ; reflines: (input) array of structures to be used as reference for the line identification:
  ;           reflines.lambda: wavelength of the lines to be identified
  ;           reflines.trust_lambda: if zero, consider the wavelength value for this line as approximate
  ;           reflines.x, measured x-coordinates for the emission lines in the reference spectrum
  ; check_shift: (input, optional) if set, the x-positions of the identified lines are compared to the x positions of the reference lines.
  ; Nmin_lines: (input, optional) minimum number of lines for the fit to be valid (default: 3)
  ; fit_window: (input, optional) window to be used for fitting each lines (default: 8 times the linewidth)
  ; poly_degree: (input, optional) polynomial degree to be used for the wavelength solution (default: 5)
  ; verbose: (input, optional) if set, and the fit fails, info on what went wrong is printed
	; speclines: (output) array of structures with parameters for each OH line
	; wavecal: (output) array with the wavelength solution
	; plot_title : (input) string to print as title of the plot
	;

  ; set default values
  if ~keyword_set(Nmin_lines) then Nmin_lines = 3
  if ~keyword_set(fit_window) then fit_window = 8.0*linewidth
  if ~keyword_set(poly_degree) then poly_degree = 5

	; convert linewidth to micron
  conversion_to_um = median( approx_wavecal - shift(approx_wavecal, 1) )
	linewidth_um = linewidth * conversion_to_um

	; identify the emission lines that are in this wavelength range
	w_lines = where(reflines.lambda GT min(approx_wavecal, /nan) $
		AND reflines.lambda LT max(approx_wavecal, /nan), /null )

	; make sure there are some lines here
	if w_lines EQ !NULL then begin
    print, 'Warning: wavelength range does not contain emission lines'
    speclines = !NULL
    wavecal = !NULL
    return
  endif

	; keep only the emission lines of interest
	line_list = reflines[w_lines].lambda
  line_trust = reflines[w_lines].trust_lambda

  ; check that there are at least three trustable lines
  if n_elements(where(line_trust eq 1, /null)) LT 3 then $
    message, 'linelist contains less than three lines that can be used for wavelength calibration'

	; make the array that will contain the result of the fitting
	speclines = []

	; fit a Gaussian to every sky line
	for i_line=0,n_elements(line_list)-1 do begin

		; select the region to fit
		w_fit = where( abs(approx_wavecal-line_list[i_line]) LT 0.5*fit_window*conversion_to_um and $
      finite(y), /null )

    ; check that there actually is signal and it's not just a bunch of NaNs or it's outside the range
    if n_elements(w_fit) LE 5 then continue

		; error handling for the gaussian fitting
		catch, error_gaussfit
		if error_gaussfit ne 0 then begin
			print, 'GAUSSFIT ERROR STATUS: ' + strtrim(error_gaussfit,2)
			catch, /cancel
			continue
		endif

		; estimate parameters of the Gaussian
		est_peak = max( median( y[w_fit], 3), ind_max , /nan)
		;est_center = w_fit[ n_elements(w_fit)/2 ]
    est_center = x[w_fit[ind_max]]
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

    ; cgplot, x[w_fit], y[w_fit], charsize=1, psym=10
    ; color='red'
    ; if line_trust[i_line] eq 0 then color='blk6'
    ; cgplot, gauss_param[1]+[0,0], [0, max(y)], /overplot, color=color
    ; print, 0.5*fit_window*conversion_to_um
    ; print, fit_window

	endfor

  ; did we find any speclines at all?
  if n_elements(speclines) eq 0 then begin
    if keyword_set(verbose) then print, 'No spectral lines were found!'
    return
  endif

  ; compare the detections to the reference positions,
  ; to make sure that the fit did not jump to an adjacent emission line
  if keyword_set(check_shift) then begin

    ; find the reference x position for all the detected lines
    epsilon = 0.01* min( abs(line_list-shift(line_list,1)) , /nan)  ; maximum distance for matching float numbers
    match, reflines.lambda, speclines.lambda, w_reflines, w_speclines, epsilon=epsilon

    ; check that lines are not too far from the (shifted) reference position
    xshift = speclines[w_speclines].x - reflines[w_reflines].x
    w_ok = where( abs( xshift - median([xshift], /even) ) LT 5.0*linewidth, /null)

    ; select only the good lines
    if w_ok EQ !NULL then return
    speclines = speclines[w_ok]

  endif

  ; select only the lines we can use for the wavelength calibration
  speclines_trust = speclines[where(speclines.trust_lambda eq 1, /null)]
  speclines_donttrust = speclines[where(speclines.trust_lambda eq 0, /null)]

  ; charsize for the plots
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

	; if too few lines were found, then no reliable wavelength solution exists
	if n_elements(speclines_trust) LT Nmin_lines then begin
    if keyword_set(verbose) then begin
      print, 'Only ', n_elements(speclines_trust), ' lines were found,'
      print, 'minimum allowed: ', Nmin_lines
    endif
		speclines = !NULL
		return
	endif

  ; set the degree for the polynomial fit - if there are few lines, decrease the degree
  if poly_degree GT (n_elements(speclines_trust)+1)/3 then poly_degree = (n_elements(speclines_trust)+1)/3

  ; fit a polynomial to the skyline positions using only the lines we can trust
  wavesol_coeff = robust_poly_fit( speclines_trust.x, speclines_trust.lambda, poly_degree, /DOUBLE )

	; calculate polynomial solution
	poly_wl = poly(x, wavesol_coeff)

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
		xtit='pixel coordinate', ytit='line width (sigma, in pixel)', charsize=ch, $
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


PRO flame_findlines_find_speclines, fuel=fuel, filename=filename, $
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
  if fuel.settings.findlines_stack_rows LT 0 then message, 'fuel.settings.findlines_stack_rows can only be positive!'
  N_stack = fuel.settings.findlines_stack_rows

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

  ; make a dummy speclines structure with the line list
  reflines_initial = replicate({lambda:0.0, trust_lambda:0, x:0.0}, n_elements(line_list))
  reflines_initial.lambda = line_list
  reflines_initial.trust_lambda = line_trust

  ; calculate typical wavelength step of one pixel, in um
  lambda_step = median( approx_lambda_axis - shift(approx_lambda_axis, 1) )

  ; approximate sky line width
	approximate_linewidth_um = median(approx_lambda_axis) / (2.36 * slit.approx_R)
	linewidth = approximate_linewidth_um / lambda_step ; in pixel

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

    ; fit the emission lines and find the wavelength solution
  	flame_findlines_fitskylines, x=pix_axis, y=central_skyspec, $
  		approx_wavecal=lambda_axis, linewidth=linewidth, reflines=reflines_initial, $
      Nmin_lines=3, fit_window=10.0*linewidth, poly_degree=fuel.settings.findlines_poly_degree, verbose=1, $
  		speclines=speclines_thisloop, wavecal=lambda_axis_output, plot_title='central rows / ' + strtrim(i_loop,2)

    ; check that lines were identified
    if n_elements(speclines_thisloop) eq 0 then begin
      cgps_close
      message, 'Not enough lines were identified in the spectrum extracted from the central rows!'
    endif

    ; update the wavelength solution
    lambda_axis = lambda_axis_output

    ; compare the number of identified lines to previous loop
    delta_Nlines = n_elements(speclines_thisloop) - Nlines
    Nlines = n_elements(speclines_thisloop)
    i_loop++

  endwhile

  w_l = where(reflines_initial.lambda GT lambda_axis[4] and reflines_initial.lambda LT lambda_axis[-5], /null)
  print, 'Identified ', strtrim(Nlines, 2), ' out of ', strtrim(n_elements(w_l), 2), ' lines present in the line list.'

  ; consider only the lines that have been identified here and use them as reference
  reflines = speclines_thisloop


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
		flame_findlines_fitskylines, x=pix_axis, y=this_row, $
			approx_wavecal=wavelength_axis_guess, linewidth=linewidth, $
			reflines=reflines, check_shift=1, poly_degree=fuel.settings.findlines_poly_degree, $
      Nmin_lines=fuel.settings.findlines_Nmin_lines, fit_window=fuel.settings.findlines_linefit_window*linewidth, $
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



PRO flame_findlines_diagnostics, fuel, slit

  ;
  ; Produce diagnostic plots to check whether the lines were properly identified
  ; in all frames. Assumes that sky lines were used
  ;

  ; number of frames
  Nfr = n_elements(slit.cutouts)

  ; start plot
	cgPS_open, file_dirname(slit.cutouts[0].filename, /mark_directory) + 'summary_line_identification.ps', /nomatch


	; take care of the frame numbers for the x axis titles (see also flame_diagnostics)
	;----------------------------------------------------------------------------------

  ; these are the values used to label the x axis
  xtickname = strtrim(fix(fuel.diagnostics.frame_num), 2)
  xtickv = indgen(Nfr)
  xminor = 1

  ; if there are too many frames, only label some of them
  while n_elements(xtickv) GT 18 do begin

    ; make the array of subindices used to select every other element
    if n_elements(xtickname) mod 2 eq 0 then $
      subset = 2*indgen(n_elements(xtickname)/2) else $   ; if odd
      subset = 2*indgen((1+n_elements(xtickname))/2)      ; if even, select also the last element

    ; keep only every other element for the labeling of the axis
    xtickname = xtickname[ subset ]
    xtickv = xtickv[ subset ]

    ; therefore, need to double the number of minor tick marks between two major marks
    xminor *= 2

  endwhile


	; read in and stack together all the speclines from all cutouts
	;----------------------------------------------------------------------------------

  ; read in the first specline
  template = (*slit.cutouts[0].speclines)[0]

  ; expand the structure by adding the frame index
  new_template = create_struct('frame', 0, template)

  ; set all fields to zero, just to be sure
  for i_tag=0, n_tags(new_template)-1 do new_template.(i_tag) = 0.0

  ; make the new, extended array of speclines
  speclines_plus = []

  ; for each cutout, extract the speclines and add them to the array
  for i_frame=0, Nfr-1 do begin

    ; speclines from this frame
    speclines_thisframe = *slit.cutouts[i_frame].speclines

    ; copy into new structure
    speclines_plus_thisframe = replicate(new_template, n_elements(speclines_thisframe))
    struct_assign, speclines_thisframe, speclines_plus_thisframe

    ; set the frame index
    speclines_plus_thisframe.frame = i_frame

    ; save in final array
    speclines_plus = [speclines_plus, speclines_plus_thisframe]

  endfor



  ; PLOT 1: number of identified lines for each frame
	;----------------------------------------------------------------------------------

  ; empty arrays with number of identified lines
  count_frame = indgen(Nfr)
  count_pixels = intarr(Nfr)
  count_lines = intarr(Nfr)
  count_lines_50p = intarr(Nfr)
  count_lines_80p = intarr(Nfr)
  count_lines_90p = intarr(Nfr)

  ; loop through each frame and count identified lines
  for i_frame=0, Nfr-1 do begin

    ; all speclines for this frame
    speclines = *slit.cutouts[i_frame].speclines

    ; identify unique wavelengths
    line_lambda = speclines.lambda
    uniq_lambda = line_lambda[UNIQ(line_lambda, SORT(line_lambda))]

    ; of these, which ones have at least 50% of pixels identified?
    line_pixels = intarr(n_elements(uniq_lambda))
    for i_line=0, n_elements(uniq_lambda)-1 do $
      line_pixels[i_line] = n_elements( where(line_lambda eq uniq_lambda[i_line], /null) )

    count_pixels[i_frame] = n_elements(speclines)
    count_lines[i_frame] = n_elements(uniq_lambda)
    count_lines_50p[i_frame] = n_elements(where(line_pixels GE 0.50*slit.height))
    count_lines_80p[i_frame] = n_elements(where(line_pixels GE 0.80*slit.height))
    count_lines_90p[i_frame] = n_elements(where(line_pixels GE 0.90*slit.height))

  endfor

  ; plot number of identified lines
	cgplot, count_frame, count_lines, psym=-16, symsize=1, charsize=1, color='black', $
		ytit='number of identified lines', $
		xra=[-1, Nfr], /xsty, xtit='frame number', $
    xtickv=xtickv, xticks=n_elements(xtickv)-1, xtickname=xtickname, xminor=xminor

  ; overplot number of lines with at least xx% of pixels identified
  cgplot, count_frame, count_lines_50p, psym=-16, symsize=1, charsize=1, color='red8', /overplot
  cgplot, count_frame, count_lines_80p, psym=-16, symsize=1, charsize=1, color='red6', /overplot
  cgplot, count_frame, count_lines_90p, psym=-16, symsize=1, charsize=1, color='red4', /overplot

  ; legend
  cgtext, 0.2, 0.30, 'number of lines detected', charsize=1, /normal, color='black'
  cgtext, 0.2, 0.26, 'number of lines covering more than 50% of the slit', charsize=1, /normal, color='red8'
  cgtext, 0.2, 0.22, 'number of lines covering more than 80% of the slit', charsize=1, /normal, color='red6'
  cgtext, 0.2, 0.18, 'number of lines covering more than 90% of the slit', charsize=1, /normal, color='red4'

  ; check if some of the frames are suspicious
  meanclip, count_lines_50p, mean_50p, sigma_50p, clipsig=2.0
  w_out_50p = where( mean_50p - count_lines_50p GT (3.0*sigma_50p > 2.0), /null)

  meanclip, count_lines_90p, mean_90p, sigma_90p, clipsig=2.0
  w_out_90p = where( mean_90p - count_lines_90p GT (3.0*sigma_90p > 2.0), /null)

  if w_out_50p NE !NULL or w_out_90p NE !NULL then begin
    print, '*********************'
    print, ''
    print, 'WARNING'
    print, 'Some frames may have bad line identification: '
    if w_out_50p NE !NULL then print, 'low number of lines identified: ', strtrim(fix(fuel.diagnostics[w_out_50p].frame_num), 2)
    if w_out_90p NE !NULL then print, 'low number of lines fully traced:', strtrim(fix(fuel.diagnostics[w_out_90p].frame_num), 2)
    print, ''
    print, '*********************'
  endif


  ; unique lines and their properties
	;----------------------------------------------------------------------------------

  ; identify unique lines
  line_lambda = speclines_plus.lambda
  uniq_lambda = line_lambda[UNIQ(line_lambda, SORT(line_lambda))]

  ; for each line, count the number of frames and detections
  uniq_framecount = intarr( n_elements(uniq_lambda) )
  uniq_detections = intarr( n_elements(uniq_lambda) )
  for i_line=0, n_elements(uniq_lambda)-1 do begin
    frames = speclines_plus[where(speclines_plus.lambda eq uniq_lambda[i_line])].frame
    uniq_framecount[i_line] = n_elements(frames[UNIQ(frames, SORT(frames))])
    uniq_detections[i_line] = n_elements(speclines_plus[where(speclines_plus.lambda eq uniq_lambda[i_line])])
  endfor

  ; write out detailed list of line detections
  forprint, uniq_lambda, uniq_detections, uniq_framecount, $
    textout=file_dirname(slit.cutouts[0].filename, /mark_directory) + 'summary_line_identification.txt', $
    comment = '# wavelength    number of detections    frames in which it is detected'

  ; identify possibly problematic lines
  w_probl = where(uniq_detections LT 0.33*median(uniq_detections), /null)
  if w_probl NE !NULL then begin
    print, '***'
    print, 'These lines have fewer detections; you may want to check them:'
    forprint, uniq_lambda[w_probl]
    print, '***'
  endif


  ; PLOT 2: properties for three representative lines
	;----------------------------------------------------------------------------------

  ; get the approximate lambda axis
  if n_elements(*slit.rough_skylambda) ne 0 then $
    lambda_range = *slit.rough_skylambda else $
    lambda_range = *slit.rough_arclambda

  ; identify "best" line in each third of the slit
  bestlines = [-1d, -1d, -1d]
  bins = lambda_range[0] + [0.0, 0.33, 0.66, 1.0] * (lambda_range[-1]-lambda_range[0])

  for j=0,2 do begin

    ; select all lines in this third of the slit
    w_third = where(uniq_lambda GE bins[j] AND uniq_lambda LT bins[j+1], /null)
    if w_third eq !NULL then continue

    ; choose the line with the highest number of identifications
    w_line = w_third[where(uniq_detections[w_third] eq max(uniq_detections[w_third]), /null)]

    ; save the best line in this third
    bestlines[j] = uniq_lambda[w_line[0]]

  endfor

  ; check that lines exist; otherwise duplicate them
  w_nogood = where(bestlines LT 0.0, /null, complement=w_good)
  if w_nogood NE !NULL then bestlines[w_nogood] = bestlines[w_good[0]]

  ; select speclines for each of the best lines (labeled A, B, and C)
  speclines_A = speclines_plus[where(speclines_plus.lambda eq bestlines[0], /null)]
  speclines_B = speclines_plus[where(speclines_plus.lambda eq bestlines[1], /null)]
  speclines_C = speclines_plus[where(speclines_plus.lambda eq bestlines[2], /null)]

  ; make a 2D array with all the identification in all frames
  lineshapes_A = dblarr(Nfr, max(speclines_plus.y)) + !values.d_NaN
  lineshapes_A[speclines_A.frame, speclines_A.y] = speclines_A.x
  lineshapes_B = dblarr(Nfr, max(speclines_plus.y)) + !values.d_NaN
  lineshapes_B[speclines_B.frame, speclines_B.y] = speclines_B.x
  lineshapes_C = dblarr(Nfr, max(speclines_plus.y)) + !values.d_NaN
  lineshapes_C[speclines_C.frame, speclines_C.y] = speclines_C.x

  ; construct the median profiles
  median_A = median(lineshapes_A, dimension=1)
  median_B = median(lineshapes_B, dimension=1)
  median_C = median(lineshapes_C, dimension=1)

  ; calculate distance to median profile
  distance_to_median_A = lineshapes_A - median_A ## replicate(1, Nfr)
  distance_to_median_B = lineshapes_B - median_B ## replicate(1, Nfr)
  distance_to_median_C = lineshapes_C - median_C ## replicate(1, Nfr)


  ; show all line identifications and line A, B, C
  cgplot, speclines_plus.x, speclines_plus.y, psym=16, symsize=0.5, color='blk5', $
    title='all line identifications', charsize=1, xtit='x', ytit='y'
  cgplot, median_A, indgen(max(speclines_plus.y)), /overplot, color='forest green', thick=5
  cgplot, median_B, indgen(max(speclines_plus.y)), /overplot, color='red', thick=5
  cgplot, median_C, indgen(max(speclines_plus.y)), /overplot, color='blue', thick=5

  erase
  ; plot the shape of the A line for each frame
  cgplot, speclines_A.x, speclines_A.y, /ynozero, /nodata, charsize=1.5, layout=[3,1,1], $
    xtitle='x', ytitle='y', title='line A: ' + strtrim(bestlines[0], 2) + ' um'
  for i=0,Nfr-1 do cgplot, lineshapes_A[i,*], indgen(max(speclines_plus.y)), /overplot, thick=2
  cgplot, median_A, indgen(max(speclines_plus.y)), /overplot, color='forest green', thick=5

  ; plot the shape of the B line for each frame
  cgplot, speclines_B.x, speclines_B.y, /ynozero, /nodata, charsize=1.5, layout=[3,1,2], $
    xtitle='x', ytitle='y', title='line B: ' + strtrim(bestlines[1], 2) + ' um'
  for i=0,Nfr-1 do cgplot, lineshapes_B[i,*], indgen(max(speclines_plus.y)), /overplot, thick=2
  cgplot, median_B, indgen(max(speclines_plus.y)), /overplot, color='red', thick=5

  ; plot the shape of the C line for each frame
  cgplot, speclines_C.x, speclines_C.y, /ynozero, /nodata, charsize=1.5, layout=[3,1,3], $
    xtitle='x', ytitle='y', title='line C: ' + strtrim(bestlines[2], 2) + ' um'
  for i=0,Nfr-1 do cgplot, lineshapes_C[i,*], indgen(max(speclines_plus.y)), /overplot, thick=2
  cgplot, median_C, indgen(max(speclines_plus.y)), /overplot, color='blue', thick=5

  erase
  ; show distance from median profile
  tot_x = [indgen(Nfr), indgen(Nfr), indgen(Nfr)] ; these are just to set the area to plot
  tot_y = [ mean(distance_to_median_A, dimension=2, /nan), $
    mean(distance_to_median_B, dimension=2, /nan), mean(distance_to_median_C, dimension=2, /nan)]
  cgplot, tot_x, tot_y, /nodata, charsize=0.8, $
    ytit='average of distance to median profile (pix)', layout=[1,2,1], $
		xra=[-1, Nfr], /xsty, xtit='frame number', $
    xtickv=xtickv, xticks=n_elements(xtickv)-1, xtickname=xtickname, xminor=xminor
  cgplot, indgen(Nfr), mean(distance_to_median_A, dimension=2, /nan), psym=-16, /overplot, color='forest green'
  cgplot, indgen(Nfr), mean(distance_to_median_B, dimension=2, /nan), psym=-16, /overplot, color='red'
  cgplot, indgen(Nfr), mean(distance_to_median_C, dimension=2, /nan), psym=-16, /overplot, color='blue'
  cgplot, [-1,Nfr], [0,0], /overplot, linestyle=2

  ; show stddev
  tot_y = [ stddev(distance_to_median_A, dimension=2, /nan), $
    stddev(distance_to_median_B, dimension=2, /nan), stddev(distance_to_median_C, dimension=2, /nan)]
  cgplot, tot_x, tot_y, /nodata, charsize=0.8, $
    ytit='std. dev. of distance to median profile (pix)', layout=[1,2,2], $
    xra=[-1, Nfr], /xsty, xtit='frame number', $
    xtickv=xtickv, xticks=n_elements(xtickv)-1, xtickname=xtickname, xminor=xminor
  cgplot, indgen(Nfr), stddev(distance_to_median_A, dimension=2, /nan), psym=-16, /overplot, color='forest green'
  cgplot, indgen(Nfr), stddev(distance_to_median_B, dimension=2, /nan), psym=-16, /overplot, color='red'
  cgplot, indgen(Nfr), stddev(distance_to_median_C, dimension=2, /nan), psym=-16, /overplot, color='blue'


  cgps_close

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_findlines_output_grid, fuel=fuel, wavelength_solution=wavelength_solution, slit=slit

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

 	; round this value on a logarithmic scale - the idea here is to try
 	; to have the same value for slightly different datasets
  outlambda_delta = 10.0^( round(alog10(lambda_delta)*20.0d)/20.0d )

  ; however, this value can also be directly specified by the user
  if fuel.settings.lambda_step GT 0.0 then begin
    if fuel.settings.lambda_step GT 0.01 then message, 'settings.lambda_step should be specified in units of micron per pixel'
    outlambda_delta = fuel.settings.lambda_step
  endif

  ; find the median wavelength of the first and last pixel
  lambda_first_pixel = median(wavelength_solution_smooth[0,*])
  lambda_last_pixel = median(wavelength_solution_smooth[-1,*])

  ; number of pixels along the vertical direction
  N_pix_y = (size(wavelength_solution))[2]

	; determine the min and max values
  ; (go a few pixel beyond the values found above, to account for tilted slits)
	lambda_min = lambda_first_pixel - 0.5*N_pix_y*lambda_delta
	lambda_max = lambda_last_pixel + 0.5*N_pix_y*lambda_delta

  ; calculate how many pixels (at the given wavelength scale)
  ; are between 1 micron and lambda_min
  Npix_1um = (lambda_min-1.0) / outlambda_delta

  ; round this down to a number that is multiple of 100
  Npix_1um_rounded = floor(Npix_1um/100.0)*100

  ; take that distance from 1um as the initial wavelength for the output grid
  outlambda_min = 1.0 + Npix_1um_rounded * outlambda_delta

  ; this is the number of pixels needed to get to lambda_max
  pixel_needed = (lambda_max-outlambda_min) / outlambda_delta

  ; round it up to a multiple of 100
  outlambda_Npix = ceil(pixel_needed/100.0)*100

  ; save output grid to the slit structure
	slit.outlambda_min = outlambda_min
	slit.outlambda_delta = outlambda_delta
  slit.outlambda_Npix = outlambda_Npix

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_findlines, fuel

	flame_util_module_start, fuel, 'flame_findlines'


  ; avoid printing too much stuff (especially from GAUSSFIT)
  quiet_state = !QUIET
  !QUIET = 1

	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

    if fuel.slits[i_slit].skip then continue

	  print, 'Finding emission lines for calibration for slit ', strtrim(fuel.slits[i_slit].number,2), ' - ', fuel.slits[i_slit].name
		print, ' '

		; handle errors by ignoring that slit
		if fuel.settings.stop_on_error eq 0 then begin
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

    if fuel.util.arc.n_frames GT 0 then begin

  		; filename of the cutout
      arc_filename = this_slit.arc_cutout.filename

  		; identify and measure the speclines
  		flame_findlines_find_speclines, fuel=fuel, filename=arc_filename, slit=this_slit, $
      rough_lambda=*this_slit.rough_arclambda, rough_flux=*this_slit.rough_arcflux, $
      linelist_filename=fuel.settings.linelist_arcs_filename, $
  			speclines=speclines, wavelength_solution=wavelength_solution

  		; save the speclines in the slit structure
  		*this_slit.arc_cutout.speclines = speclines

  		; write a ds9 region file with the identified speclines
  		flame_findlines_writeds9, speclines, filename=flame_util_replace_string(arc_filename, '.fits', '_speclines.reg')

			; use the pixel-by-pixel wavelength solution of the arc frame to set the output grid in wavelength
			flame_findlines_output_grid, fuel=fuel, wavelength_solution=wavelength_solution, slit=this_slit
			fuel.slits[i_slit] = this_slit


    endif else begin

    ; skylines -----------------------------------------------------------------

  		for i_frame=0, n_elements(fuel.slits[i_slit].cutouts)-1 do begin

  				; filename of the cutout
  				slit_filename = fuel.slits[i_slit].cutouts[i_frame].filename

  				; identify and measure the speclines
  				flame_findlines_find_speclines, fuel=fuel, filename=slit_filename, slit=this_slit, $
          rough_lambda=*this_slit.rough_skylambda, rough_flux=*this_slit.rough_skyflux, $
          linelist_filename=fuel.settings.linelist_sky_filename, $
  					speclines=speclines, wavelength_solution=wavelength_solution

  				; save the speclines in the slit structure
  				*this_slit.cutouts[i_frame].speclines = speclines

  				; write a ds9 region file with the identified speclines
  				flame_findlines_writeds9, speclines, filename=flame_util_replace_string(slit_filename, '.fits', '_speclines.reg')

  				; use the pixel-by-pixel wavelength solution OF THE FIRST FRAME to set the output grid in wavelength
  				if i_frame eq 0 then begin
  					flame_findlines_output_grid, fuel=fuel, wavelength_solution=wavelength_solution, slit=this_slit
  					fuel.slits[i_slit] = this_slit
  				endif

  		endfor

      ; diagnostic plot to check the line identification
      flame_findlines_diagnostics, fuel, this_slit

    endelse

	endfor

	; revert to original !QUIET state
	!QUIET = quiet_state


  flame_util_module_end, fuel

END
