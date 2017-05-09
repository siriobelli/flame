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


PRO flame_wavecal_2D_calibration, filename=filename, slit=slit, OHlines=OHlines, $
		wavecal_accurate=wavecal_accurate, diagnostics=diagnostics, this_diagnostics=this_diagnostics

; This routine calculates the 2D wavelength solution and y-rectification.
; These are two mappings from the observed pixel coordinates to the rectified grid
; lambdax, gamma, where lambdax is a pixel grid linear in lambda and gamma is the
; vertical distance to the edge of the slit (taking into account warping and also vertical drift)
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
	first_pixel = ceil(ymin_edge)

	; OH line coordinates
	OH_lambda = OHlines.lambda
	OH_x = OHlines.x
	OH_y = OHlines.y

	; output lambda axis
	lambda_0 = slit.outlambda_min
	delta_lambda = slit.outlambda_delta

	; calculate vertical offset
	w_this_offset = where(diagnostics.offset_pos eq this_diagnostics.offset_pos)
	ref_diagnostics = diagnostics[w_this_offset[0]]
	vertical_offset = this_diagnostics.position - floor(ref_diagnostics.position)

	; translate every OH detection into the new coordinate system
	OH_lambdax = (OH_lambda - lambda_0)/delta_lambda
	OH_gamma = OH_y + first_pixel - poly(OH_x, slit.bottom_poly) - vertical_offset

	; indices of the pixels we want to use - start with using all of them
	w_goodpix = lindgen(n_elements(OH_x))
	Ngoodpix = n_elements(w_goodpix)+1

	; check that we have enough points to calculate warping polynomial
	if n_elements(w_goodpix) LT (degree+1.0)^2 then message, 'not enough data points for polywarp'

	; loops are used to throw away outliers and make polywarp more robust
	WHILE n_elements(w_goodpix) LT Ngoodpix AND n_elements(w_goodpix) GE (degree+1.0)^2  DO BEGIN

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

	; finally, output the actual wavelength calibration as a 2D array
	wavecal_accurate = im * 0.0

	; order of polynomial
	Nord = (size(Klambda))[1]
	xexp  = findgen(Nord)
	yexp  = findgen(Nord)

	for ix=0.0,N_imx-1 do $
		for iy=0.0,N_imy-1 do $
			wavecal_accurate[ix,iy] = lambda_0 + delta_lambda * total(((iy)^xexp # (ix)^yexp ) * Klambda)

	; write the accurate solution to a FITS file
	writefits, flame_util_replace_string(filename, '.fits', '_wavecal_smooth.fits'), wavecal_accurate, hdr


END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_wavecal_output_grid, wavelength_solution=wavelength_solution, $
	OHlines=OHlines, slit=slit

	; fit a 3rd degree polynomial to the 2D wavelength solution array in order to smooth it
	; 1. this could be more robust
	; 2. Is there a better way to set the output wavelength grid?


	; how many pixels on the spatial direction
	N_spatial_pix = (size(wavelength_solution))[2]

	; how many pixels on the wavelength direction
	N_lambda_pix = (size(wavelength_solution))[1]

	; unpack the OH line detection array
	OH_wavelength = OHlines.lambda
	OH_xpixel = OHlines.x
	OH_ypixel = OHlines.y

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


PRO flame_wavecal_illum_correction, OHlines=OHlines, filename=filename, rectification=rectification, slit=slit
	;
	; use the OH lines to dervie and apply an illumination correction
	; along the spatial slit axis
	;

	;
	; WARNING: THE ILLUMINATION-CORRECTED FRAMES ARE NOT ACTUALLY USED FOR NOW
	;

	; read in slit
	im = readfits(filename, hdr)
	N_pixel_x = (size(im))[1]
	N_pixel_y = (size(im))[2]

	; sort by wavelength
	sorted_lambdas = OHlines[ sort(OHlines.lambda) ].lambda

	; find unique wavelengths
	lambdas = sorted_lambdas[uniq(sorted_lambdas)]

	; calculate Gussian fluxes
	OHflux = sqrt(2.0*3.14) * OHlines.peak * OHlines.sigma

	; for each line, normalize flux to the median
	OHnorm = OHflux * 0.0
	for i_line=0, n_elements(lambdas)-1 do begin
		w_thisline = where(OHlines.lambda eq lambdas[i_line], /null)
		OHnorm[w_thisline] = OHflux[w_thisline] / median(OHflux[w_thisline])
	endfor

	; extract the rectification matrix for gamma
	Kgamma = rectification.Kgamma

	; order of polynomial
	Nord = (size(Kgamma))[1]
	xexp  = findgen(Nord)
	yexp  = findgen(Nord)

	; calculate the gamma coordinate of each OHline detection
	OHgamma = dblarr(n_elements(OHlines))
	for i_line=0, n_elements(OHlines)-1 do $
		OHgamma[i_line] = total(((OHlines[i_line].y)^xexp # (OHlines[i_line].x)^yexp ) * Kgamma)

	; sort by gamma
	sorted_gamma = OHgamma[sort(OHgamma)]
	sorted_illum = OHnorm[sort(OHgamma)]

	; do not consider measurements that are more than a factor of three off
	w_tofit = where(sorted_illum GT 0.33 and sorted_illum LT 3.0, /null)

	; fit polynomial to the illumination correction as a function of gamma
	poly_coeff = poly_fit(sorted_gamma[w_tofit], sorted_illum[w_tofit], 10)

	; set the boundaries for a meaningful correction
	gamma_min = sorted_gamma[3]
	gamma_max = sorted_gamma[-4]

	; scatter plot of the illumination (show all OH lines)
	cgplot, sorted_gamma[w_tofit], sorted_illum[w_tofit], psym=3, /ynozero, charsize=1.2, $
		xtitle='gamma coordinate', ytitle='Illumination', title='Illumination correction'

	; overplot the smooth illumination
	x_axis = gamma_min + (gamma_max-gamma_min)*dindgen(200)/199.
	cgplot, x_axis, poly(x_axis, poly_coeff), $
		/overplot, color='red', thick=3

	; overplot flat illumination
	cgplot, [gamma_min - 0.5*gamma_max , gamma_max*1.5], [1,1], /overplot, linestyle=2, thick=3

	; calculate the gamma coordinate for each observed pixel
	gamma_coordinate = im * 0.0
	for ix=0.0,N_pixel_x-1 do $
		for iy=0.0,N_pixel_y-1 do $
			gamma_coordinate[ix,iy] = total(((iy)^xexp # (ix)^yexp ) * Kgamma)

	; calculate the illumination correction at each pixel
	illumination_correction = poly(gamma_coordinate, poly_coeff)

	; set the correction to NaN when outside the boundary
	illumination_correction[where(gamma_coordinate LT gamma_min OR $
		gamma_coordinate GT gamma_max, /null)] = !values.d_NaN

	; apply illumination correction
	im /= illumination_correction

	; write out the illumination-corrected cutout
  writefits, flame_util_replace_string(filename, '.fits', '_illumcorr.fits'), im, hdr


END




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_writeds9, OHlines, filename=filename
;
; write a ds9 region file with all the OH line detections
;

  ; extract the wavelength of all identifications
  line_lambdas = OHlines.lambda
  uniq_lambdas = line_lambdas[UNIQ(line_lambdas, SORT(line_lambdas))]

  ; open file
  openw, lun, filename, /get_lun

  ; write header
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global color=red dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=1 include=1 source=1'
  printf, lun, 'image'

  for i_line=0, n_elements(uniq_lambdas)-1 do begin

  	; the wavelength of this OH line
  	this_lambda = uniq_lambdas[i_line]

  	; select the x and y coordinates for this line
  	w_thisline = where(OHlines.lambda eq this_lambda)
  	this_x = OHlines[w_thisline].x
 		this_y = OHlines[w_thisline].y

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

    ; write the line corresponding to this slit
    printf, lun, 'polygon(' + all_points + ') # text={' + cgnumber_formatter(this_lambda, decimals=5) + '}'

  endfor

  ; close file
  free_lun, lun


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_wavecal_fitskylines, x=x, y=y, $
	approx_wavecal=approx_wavecal, linewidth=linewidth, $
	line_list=line_list_in, $
	OHlines=OHlines, wavecal=wavecal, plot_title=plot_title

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
	; OHlines: (output) array of structures with parameters for each OH line
	; wavecal: (output) array with the accurate wavelength solution
	; plot_title : (input) string to print as title of the plot
	;


  ; settings:
	; the degree of the polynomial used to describe the wavelength solution
	poly_degree = 3

	; minimum number of OH lines for a reliable wavelength solution
	Nmin_lines = 6

	; convert linewidth to micron (assuming linear wavelength solution)
	linewidth_um = linewidth * (approx_wavecal[3]-approx_wavecal[2])

	; identify the OH lines that are in this wavelength range
	w_lines = where(line_list_in GT min(approx_wavecal, /nan) $
		AND line_list_in LT max(approx_wavecal, /nan), /null )

	; make sure there are OH lines here
	if w_lines EQ !NULL then message, 'Wavelength range does not contain OH lines?!'

	; keep only the OH lines of interest
	line_list = line_list_in[w_lines]

	; make the array that will contain the result of the fitting
	OHlines = []

	; fit a Gaussian to every sky line
	for i_line=0,n_elements(line_list)-1 do begin

		; select the region to fit
		w_fit = where( abs(approx_wavecal-line_list[i_line]) LT 6.0*linewidth_um, /null )

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
			chisq: chisq }

		; add to the stack
		OHlines = [OHlines, this_OHline]

	endfor

	; if too few lines were found, then no reliable wavelength solution exists
	if n_elements(OHlines) LT Nmin_lines then begin
		OHlines = !NULL
		return
	endif

	; fit a polynomial to the skyline positions
	wavesol_coeff = poly_fit( OHlines.x, OHlines.lambda, poly_degree )

	; calculate polynomial solution
	poly_wl = poly(x, wavesol_coeff)

	; properly handle regions with no information
	; if where( ~finite(y), /null ) NE !NULL then begin	; check if there are NaNs in the input sky spectrum
	; 	nosky_regions = label_region( ~finite(y) )		; this will have zero everywhere and N in the Nth "region" of NaNs
	; 	nosky_regions[0] = nosky_regions[1]			; boundary issue
	; 	nosky_regions[-1] = nosky_regions[-2]		; boundary issue
	; 	if nosky_regions[1] ne 0 then poly_wl[ where(nosky_regions eq nosky_regions[1]) ] = !values.d_NaN	; void the region at the beginning
	; 	if nosky_regions[-2] ne 0 then poly_wl[ where(nosky_regions eq nosky_regions[-2]) ] = !values.d_NaN	; void the region at the end
	; endif

	; panel 1: plot the spectrum
	erase
	cgplot, x, y, charsize=1, xsty=1, xtit='', ytit='sky flux', title=plot_title, $
		position = [0.15, 0.69, 0.95, 0.96], xtickformat="(A1)", xra=[x[0], x[-1]]
	for i_line=0, n_elements(OHlines)-1 do cgplot, OHlines[i_line].x + [0,0], [-2,2]*max(abs(y)), /overplot, color='red'

	; panel 2: show the result of Gaussian fitting
	cgplot, OHlines.x, OHlines.lambda, /ynozero, xra=[x[0], x[-1]], xsty=1, psym=16, color='red', symsize=0.7, $
		xtit='', ytit='expected wavelength', charsize=1, $
		/noerase, position = [0.15, 0.42, 0.95, 0.69], xtickformat="(A1)"

	; show the polynomial fit
	cgplot, x, poly_wl, color='blue', /overplot

	; panel 3: plot the line widths
	cgplot, OHlines.x, OHlines.sigma, /ynozero, xra=[x[0], x[-1]], xsty=1, psym=16, color='red', symsize=0.7, $
		xtit='pixel coordinate', ytit='line width (pixel)', charsize=1, $
		/noerase, position = [0.15, 0.15, 0.95, 0.42]

	; show median value of line width
	cgplot, [-1, 2*max(OHlines.x)], [0,0]+median(OHlines.sigma), /overplot, thick=3, linestyle=2

	; update value of linewidth
	linewidth = median(OHlines.sigma)

	; output wavelength solution
	wavecal = poly_wl

END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_find_OHlines, fuel=fuel, slit_filename=slit_filename, $
	 approx_lambda_axis=approx_lambda_axis, $
	 OHlines=OHlines, wavelength_solution=wavelength_solution

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

	; load line list
	readcol, fuel.util.linelist_filename, line_list, format='D', /silent

	; approximate sky line width (assuming one arcsec slit width)
	approximate_linewidth_um = median(approx_lambda_axis) / (2.36 * fuel.instrument.resolution_slit1arcsec)
	linewidth = approximate_linewidth_um / ( approx_lambda_axis[3] - approx_lambda_axis[2] )

	; start from the central row and go up until the top row, then start from center and go down
	row_number = indgen(N_spatial_pix)
	sorted_rows = [ row_number[N_spatial_pix/2: N_spatial_pix-1] , reverse(row_number[0:N_spatial_pix/2-1]) ]

	; identify first row of bottom half
	i0_bottom = row_number[N_spatial_pix/2-1]

	; as initial guess for the wavelength axis, use what found during the rough wavecal
	wavelength_axis_guess = approx_lambda_axis

	; create the empty array of OHlines structures
	OHlines = []

	print, 'Fitting individual sky lines for every pixel row...'

	; loop through all the pixel rows
	for counter=0, N_spatial_pix-1 do begin

		; index of the row we are considering now
		i_row = sorted_rows[counter]

		; print info on the row
		print, 'row ' + strtrim(i_row, 2) + ' ', format='(a,$)'

		; extract this pixel row from the slit
		this_row = im[*, i_row]

		; if this is the first row of the bottom half, then use the initial wavelength guess
		if i_row eq i0_bottom then wavelength_axis_guess = approx_lambda_axis

		; fit the emission lines and find the wavelength solution
		flame_wavecal_fitskylines, x=pix_axis, y=this_row, $
			approx_wavecal=wavelength_axis_guess, linewidth=linewidth, $
			line_list=line_list, $
			OHlines=OHlines_thisrow, wavecal=wavelength_axis_for_this_row, plot_title='row '+strtrim(i_row,2)

		; if sky lines were not found, then skip to next row
		if n_elements(OHlines_thisrow) EQ 0 then continue

		; set the y coordinate for the OH lines
		OHlines_thisrow.y = i_row

		; save the OH lines from this row
		OHlines = [ OHlines, OHlines_thisrow ]

		; save the wavelength solution
		wavelength_solution[*, i_row] = wavelength_axis_for_this_row

		; if a solution was found, then use it as the initial guess for the next row
		if where(finite(wavelength_axis_for_this_row), /null) NE !NULL then $
			wavelength_axis_guess = wavelength_axis_for_this_row

	endfor

	print, ''

	cgPS_close

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_plots, wavelength_solution=wavelength_solution, OHlines=OHlines

	; get the dimensions
	N_lambda_pix = (size(wavelength_solution))[1]
	N_spatial_pix = (size(wavelength_solution))[2]

	; -------------------------------------------------------
	; plot the individual detections on a 2D view of the slit
	cgplot, OHlines.x, OHlines.y, psym=16, xtit='x pixel', ytitle='y pixel on this slit', $
	title='OH line detections', charsize=1, layout=[1,2,1], symsize=0.5

	; show the line widths
	cgplot, OHlines.x, OHlines.sigma, psym=16, xtit='x pixel', ytitle='line width (pixel)', $
	charsize=1, layout=[1,2,2], symsize=0.5, yra = median(OHlines.sigma)*[0.5, 1.5]
	cgplot, [0, 2*max(OHlines.x)], median(OHlines.sigma) + [0,0], /overplot, thick=3, linestyle=2

	; -------------------------------------------------------
	; plot the shift in wavelength as a function of spatial position
	erase

	; take the median wavelength solution as a reference
	wavelength_solution_reference = median(wavelength_solution, dimension=2)

	; calculate the typical shift in wavelength for each row from the reference solution. Take the median wl of the central 5 pixels
	wavelength_shift = dblarr(N_spatial_pix)
	for i_row=0, N_spatial_pix-1 do $
	wavelength_shift[i_row] = $
	median(wavelength_solution[ N_lambda_pix/2-2 : N_lambda_pix/2+2, i_row]) - $
	median(wavelength_solution_reference[ N_lambda_pix/2-2 : N_lambda_pix/2+2 ] )

	; cut the pixels near the edge
	wavelength_shift[0:2] = !values.d_nan
	wavelength_shift[-3:-1] = !values.d_nan

	; ; get rid of wildly wrong rows
	; wavelength_shift[where(abs(wavelength_shift - median(wavelength_shift)) $
	; 	GT 5.0 *stddev(wavelength_shift, /nan), /null)] = !values.d_NaN

	; plot the shift as a function of vertical position
	cgplot, 1d4*wavelength_shift, psym=-16, thick=3, charsize=1, $
	xtit = 'pixel position along the vertical (spatial) axis', ytit='wavelength shift from the reference pixel row (angstrom)'
	cgplot, [-1d4, 1d4], [0,0], /overplot, thick=2


END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_onecutout, fuel=fuel, i_slit=i_slit, i_frame=i_frame, $
	guess_lambda_axis=guess_lambda_axis

	; filename of the cutout
	slit_filename = (*fuel.slits[i_slit].filenames)[i_frame]

	; this slit
	this_slit = fuel.slits[i_slit]

	; identify and measure the OH lines
	flame_wavecal_find_OHlines, fuel=fuel, slit_filename=slit_filename, $
		approx_lambda_axis=guess_lambda_axis, $
		OHlines=OHlines, wavelength_solution=wavelength_solution

	; write a ds9 region file with the identified OH lines
	flame_wavecal_writeds9, OHlines, filename=flame_util_replace_string(slit_filename, '.fits', '_OHlines.reg')

	; write a FITS file with the pixel-by-pixel wavelength solution
	writefits, flame_util_replace_string(slit_filename, '.fits', '_wavecal.fits'), $
		wavelength_solution, headfits(slit_filename)

	; show plots of the wavelength calibration and OH line identification
	cgPS_open, flame_util_replace_string(slit_filename, '.fits', '_plots.ps'), /nomatch
	flame_wavecal_plots, wavelength_solution=wavelength_solution, OHlines=OHlines

	flame_wavecal_output_grid, wavelength_solution=wavelength_solution, $
		OHlines=OHlines, slit=this_slit

	flame_wavecal_2D_calibration, filename=slit_filename, $
		slit=this_slit, OHlines=OHlines, wavecal_accurate=wavecal_accurate, $
		diagnostics=fuel.diagnostics, this_diagnostics=(fuel.diagnostics)[i_frame]

	; calculate and apply the illumination correction
		flame_wavecal_illum_correction, OHlines=OHlines, filename=slit_filename, $
		rectification = (*this_slit.rectification)[i_frame], slit=this_slit

	cgPS_close

	; update the slit structure with the output wavelength grid of the last frame
	fuel.slits[i_slit] = this_slit

	; use the accurate wavecal of the central pixel row as guess for the next frame
	guess_lambda_axis = wavecal_accurate[ * , (size(wavecal_accurate))[2]/2 ]


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

		; the initial guess for the first frame is from the rough wavelength calibration
		guess_lambda_axis = *this_slit.rough_wavecal

		for i_frame=0, n_elements(*slits[i_slit].filenames)-1 do begin

			flame_wavecal_onecutout, fuel=fuel, i_slit=i_slit, i_frame=i_frame, $
		 	 guess_lambda_axis=guess_lambda_axis

		endfor

	endfor

	; revert to original !QUIET state
	!QUIET = quiet_state
	print, 'It took ', systime(/seconds) - start_time, ' seconds'

END
