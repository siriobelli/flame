;
; Wavelength calibration.
; For each slit & frame, find a polynomial
; warping that describes the 2D transformation from observed frame
; to lambda-calibrated and vertically-rectified frame, using
; the emission line measurements from flame_identify_lines
;

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************

FUNCTION lambda_calibration, coefficients, speclines=speclines, lambdax=lambdax, lambda_polyorder=lambda_polyorder
	;
	; This function is used for finding the best-fit coefficients
	; describing the wavelength of each pixek in the observed frame
	;
	; The coefficients are: [P00, P01, P02, P03, P10, P11, P12, P13, P20, ..... PNN]
	; and the lambdax calibration is of the form SUM(Pij*x^i*y^j)
	;

	; check that the total dimension of the coefficients matches the polynomial order
	if n_elements(coefficients) NE (lambda_polyorder[0]+1)*(lambda_polyorder[1]+1) then message, 'dimension of coefficients is wrong!'

	; transform the coefficients from 1D to 2D
	Klambda = reform(coefficients, lambda_polyorder[1]+1, lambda_polyorder[0]+1)

	; given the coefficients, calculate the predicted normalized lambda for each specline
	predicted_lambdax = dblarr(n_elements(speclines))*0.0
	for i=0,lambda_polyorder[0] do for j=0,lambda_polyorder[1] do predicted_lambdax += double(Klambda[j,i]) * double(speclines.x)^i * double(speclines.y)^j

	; construct array with deviation in lambda for each line
	dev = dblarr(n_elements(speclines))

	; different way to calculate dev for the two types of lines
	w_trust = where(speclines.trust_lambda eq 1, /null, compl=w_donttrust)

	; for the lines that can be used for the wavelength calibration,
	; calculate the deviation from the theoretical lambda
	dev[w_trust] = lambdax[w_trust] - predicted_lambdax[w_trust]

	; for the lines with unreliable lambda, calculate the deviation from the median
	; value of all the other detections of the same line (this helps with the rectification)
	if n_elements(w_donttrust) GT 0 then $
		for i_line=0, n_elements(w_donttrust)-1 do $
			dev[i_line] = predicted_lambdax[i_line] - median(predicted_lambdax[where(speclines.lambda eq speclines[i_line].lambda, /null)])

	; return deviation of true lambdax from predicted value
	return, lambdax - predicted_lambdax

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_2D_calibration, fuel=fuel, slit=slit, cutout=cutout, $
		diagnostics=diagnostics, this_diagnostics=this_diagnostics

;
; This routine calculates the 2D wavelength solution and y-rectification.
; These are two mappings from the observed pixel coordinates to the rectified grid
; lambdax, gamma, where lambdax is a pixel grid linear in lambda and gamma is the
; vertical distance to the edge of the slit (taking into account warping and also vertical drift)
; The result of this routine is a pair of coefficient sets, lambda_coeff and gamma_coeff, saved in fuel.slits,
; that can be used to rectify the image
;

	print, ''
	print, 'Accurate 2D wavelength solution for ', cutout.filename_step1

	; read in file to calibrate
	im = mrdfits(cutout.filename_step1, 0, header, /silent)

	; read dimensions of the image
	N_imx = (size(im))[1]
	N_imy = (size(im))[2]

	; find the minimum y value for the bottom edge of the slit
	bottom_edge = poly(indgen(N_imx), slit.bottom_poly)
	ymin_edge = min(bottom_edge)

	; this is the y-coordinate of the bottom pixel row in the cutout
	first_pixel = ceil(ymin_edge)

	; specline coordinates
	speclines = *cutout.speclines
	OH_lambda = speclines.lambda
	OH_x = speclines.x
	OH_y = speclines.y

	; output lambda axis
	lambda_0 = slit.outlambda_min
	delta_lambda = slit.outlambda_delta

	; calculate vertical offset
	w_this_offset = where(diagnostics.offset_pos eq this_diagnostics.offset_pos)
	ref_diagnostics = diagnostics[w_this_offset[0]]
	vertical_offset = this_diagnostics.position - floor(ref_diagnostics.position)

	; translate every OH detection into the new coordinate system
	OH_lambdax = (OH_lambda - lambda_0)/delta_lambda
	;OH_gamma = OH_y + first_pixel - poly(OH_x, slit.bottom_poly) - vertical_offset


	; find the gamma coefficients analytically -------------------------------------------------------

	; make the 4x4 matrix
	Kgamma = dblarr(4, 4)

	; by definition, dgamma/dy = 1
	Kgamma[1,0] = 1.0

	; use the definition of the bottom edge of the slit to build the gamma matrix
	Kgamma[0,*] = -slit.bottom_poly

	; but the zero is now at the edge of the cutout, and is corrected for the drift
	Kgamma[0,0] += first_pixel - vertical_offset


	; guess the lambda coefficients -------------------------------------------------------

	; the degrees of the 2D polynomial that describes lambda as a function of (x,y) (plus 1)
	polyorder_x = fuel.util.wavesolution_order_x
	polyorder_y = fuel.util.wavesolution_order_y
	if polyorder_x LT 1 or polyorder_x GT 10 then message, 'lambda_polyorder_x is out of range'
	if polyorder_y LT 1 or polyorder_y GT 10 then message, 'lambda_polyorder_y is out of range'

	; guess starting coefficients by fitting a polynomial to the lines at the bottom of the slit
	sorted_y = speclines[sort(speclines.y)].y
	y_threshold = sorted_y[0.1*n_elements(sorted_y)]
	w_tofit = where(speclines.y LE y_threshold, /null)
	guess_coeff1d = robust_poly_fit(speclines[w_tofit].x, (speclines[w_tofit].lambda-slit.outlambda_min)/slit.outlambda_delta, polyorder_x)

	; scale by orders of magnitude for the higher order terms
	starting_coefficients_l = dblarr( polyorder_y+1, polyorder_x+1 )
	for i=0,polyorder_y do starting_coefficients_l[i,*] = guess_coeff1d * (1d-3)^i
	starting_coefficients_l = starting_coefficients_l[*]


	; fit the observed speclines and find the best-fit lambda coefficients --------------------------------------------

	; number of speclines we are using
	Ngoodpix = n_elements(speclines)+1

	; check that we have enough points to calculate warping polynomial
	if n_elements(speclines) LT 5 then message, 'not enough data points for a wavelength solution'

	; loops are used to throw away outliers and make the fit more robust
	WHILE n_elements(speclines) LT Ngoodpix AND n_elements(speclines) GE 5 DO BEGIN

		; save old number of good speclines
		Ngoodpix = n_elements(speclines)

		args = {speclines:speclines, lambdax:OH_lambdax, lambda_polyorder:[polyorder_x, polyorder_y] }

		; fit the data and find the coefficients for the lambda calibration
		lambda_coeff = mpfit('lambda_calibration', starting_coefficients_l, functargs=args, $
			bestnorm=bestnorm_l, best_resid=best_resid_l, /quiet, status=status_l)

		; check that mpfit worked
		if status_l LE 0 then message, 'mpfit did not find a good solution'

		w_outliers = where( abs(best_resid_l) GT 3.0*stddev(best_resid_l), complement=w_goodpix, /null)
		print, strtrim( n_elements(w_outliers), 2) + ' outliers rejected. ', format='(a,$)'

		; keep only the non-outliers
		speclines = speclines[w_goodpix]
		OH_lambdax = OH_lambdax[w_goodpix]

	ENDWHILE

	print, ''

	; convert the coefficients to a 2D matrix
	Klambda = reform(lambda_coeff, polyorder_y+1, polyorder_x+1)

; ---------------------------------------------------------------------------------------

	; save into slit structure - copy also the lambda_min and lambda_step parameters
	*cutout.rectification = {Klambda:Klambda, Kgamma:Kgamma, $
		lambda_min:slit.outlambda_min, lambda_delta:slit.outlambda_delta}

	; finally, output the actual wavelength calibration as a 2D array

	; create empty frame that will contain the wavelength solution
	wavelength_solution = im * 0.0

	; apply the polynomial transformation to calculate (lambda, gamma) at each point of the 2D grid
	for ix=0, N_imx-1 do $
		for iy=0, N_imy-1 do begin
			flame_util_transform_direct, *cutout.rectification, x=ix, y=iy, lambda=lambda, gamma=gamma
			wavelength_solution[ix, iy] = lambda
		endfor

	; write the accurate solution to a FITS file
	writefits, flame_util_replace_string(cutout.filename_step1, '.fits', '_wavecal_2D.fits'), wavelength_solution, hdr


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_illum_correction, fuel=fuel, i_slit=i_slit, i_frame=i_frame
	;
	; use the speclines to derive and apply an illumination correction
	; along the spatial slit axis
	;

	cutout = fuel.slits[i_slit].cutouts[i_frame]

	speclines = *cutout.speclines

	; read in slit
	im = mrdfits(cutout.filename_step1, 0, hdr, /silent)
	im_sigma = mrdfits(cutout.filename_step1, 1, /silent)

	; cutout dimensions
	N_pixel_x = (size(im))[1]
	N_pixel_y = (size(im))[2]

	; sort by wavelength
	sorted_lambdas = speclines[ sort(speclines.lambda) ].lambda

	; find unique wavelengths
	lambdas = sorted_lambdas[uniq(sorted_lambdas)]

	; calculate Gussian fluxes
	OHflux = sqrt(2.0*3.14) * speclines.peak * speclines.sigma

	; for each line, normalize flux to the median
	OHnorm = OHflux * 0.0
	for i_line=0, n_elements(lambdas)-1 do begin
		w_thisline = where(speclines.lambda eq lambdas[i_line], /null)
		OHnorm[w_thisline] = OHflux[w_thisline] / median(OHflux[w_thisline])
	endfor

	; calculate the gamma coordinate of each OHline detection
	flame_util_transform_direct, *cutout.rectification, x=speclines.x, y=speclines.y, $
		lambda=OHlambda, gamma=OHgamma

	; sort by gamma
	sorted_gamma = OHgamma[sort(OHgamma)]
	sorted_illum = OHnorm[sort(OHgamma)]

	; do not consider measurements that are more than a factor of three off
	w_tofit = where(sorted_illum GT 0.33 and sorted_illum LT 3.0, /null)

	; fit polynomial to the illumination correction as a function of gamma
	poly_coeff = poly_fit(sorted_gamma[w_tofit], sorted_illum[w_tofit], 8)

	; set the boundaries for a meaningful correction
	gamma_min = sorted_gamma[3]
	gamma_max = sorted_gamma[-4]

	; scatter plot of the illumination (show all OH lines)
	cgplot, sorted_gamma[w_tofit], sorted_illum[w_tofit], psym=3, /ynozero, charsize=1.2, $
		xtitle='gamma coordinate', ytitle='Illumination'

	; overplot the smooth illumination
	x_axis = gamma_min + (gamma_max-gamma_min)*dindgen(200)/199.
	cgplot, x_axis, poly(x_axis, poly_coeff), $
		/overplot, color='red', thick=3

	; overplot flat illumination
	cgplot, [gamma_min - 0.5*gamma_max , gamma_max*1.5], [1,1], /overplot, linestyle=2, thick=3

	; overplot the limit to the correction (25%)
	cgplot, [gamma_min - 0.5*gamma_max , gamma_max*1.5], 0.75+[0,0], /overplot, linestyle=2, thick=1
	cgplot, [gamma_min - 0.5*gamma_max , gamma_max*1.5], 1.25+[0,0], /overplot, linestyle=2, thick=1

	; if we do not have to apply the illumination correction, then we are done
	if ~fuel.util.illumination_correction then begin
		fuel.slits[i_slit].cutouts[i_frame].filename_step2 = fuel.slits[i_slit].cutouts[i_frame].filename_step1
		return
	endif

	; calculate the gamma coordinate for each observed pixel, row by row
	gamma_coordinate = im * 0.0
	for i_row=0, N_pixel_y-1 do begin
		flame_util_transform_direct, *cutout.rectification, x=dindgen(N_pixel_x), y=replicate(i_row, N_pixel_x), $
			lambda=lambda_row, gamma=gamma_row
		gamma_coordinate[*,i_row] = gamma_row
	endfor

	; calculate the illumination correction at each pixel
	illumination_correction = poly(gamma_coordinate, poly_coeff)

	; set the correction to NaN when outside the boundary
	illumination_correction[where(gamma_coordinate LT gamma_min OR $
		gamma_coordinate GT gamma_max, /null)] = !values.d_NaN

	; set the correction to NaN if it's more than 25%
	illumination_correction[where(illumination_correction GT 1.25 or $
		illumination_correction LT 0.75, /null)] = !values.d_NaN

	; apply illumination correction
	im /= illumination_correction
	im_sigma /= illumination_correction

	; filename for the illumination-corrected cutout
	illcorr_filename = flame_util_replace_string(cutout.filename_step1, '_corr', '_illcorr')

	; write out the illumination-corrected cutout
  writefits, illcorr_filename, im, hdr
	writefits, illcorr_filename, im_sigma, /append

	; save the filename of the illumination-corrected frame in the cutout structure as the step2 frame
	fuel.slits[i_slit].cutouts[i_frame].filename_step2 = illcorr_filename


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_plots, slit=slit, cutout=cutout

	speclines = *cutout.speclines

	; -------------------  select two example lines  -------------------------------

		; sort by wavelength
		sorted_lambdas = speclines[ sort(speclines.lambda) ].lambda

		; find unique wavelengths
		lambdas = sorted_lambdas[uniq(sorted_lambdas)]

		; for each wavelength find number of detections, median x position and flux
		lambdas_det = intarr(n_elements(lambdas))
		lambdas_x = fltarr(n_elements(lambdas))
		lambdas_peak = fltarr(n_elements(lambdas))
		for i=0, n_elements(lambdas)-1 do begin
			w0 = where(speclines.lambda eq lambdas[i], /null)
			lambdas_det[i] = n_elements(w0)
			if lambdas_det[i] LT 2 then continue
			lambdas_x[i] = median(speclines[w0].x)
			lambdas_peak[i] = median(speclines[w0].peak)
		endfor

		; find the center of the x-distribution of the speclines
		x_center = 0.5*(max(lambdas_x) - min(lambdas_x))

		; split into left and right sides
		w_left = where(lambdas_x LT x_center, complement=w_right, /null)

		; select lines that have a good number of detections
		w_left_num = cgsetintersection(w_left, where(lambdas_det GE median(lambdas_det[w_left]), /null) )
		w_right_num = cgsetintersection(w_right, where(lambdas_det GE median(lambdas_det[w_right]), /null) )

		; pick the brightest ones
		!NULL = max(lambdas_peak[w_left_num], w_chosen, /nan )
		ind_line1 = w_left_num[w_chosen]
		!NULL = max(lambdas_peak[w_right_num], w_chosen, /nan )
		ind_line2 = w_right_num[w_chosen]

		; select all the detections for these two lines
		lambda1 = lambdas[ind_line1]
		lambda2 = lambdas[ind_line2]
		w_line1 = where(speclines.lambda eq lambda1, /null)
		w_line2 = where(speclines.lambda eq lambda2, /null)

		color1 = 'blu4'
		color2 = 'red4'

	; -------------------------------------------------------
	; plot the individual detections on a 2D view of the slit
	cgplot, speclines.x, speclines.y, psym=16, xtit='x pixel', ytitle='y pixel on this slit', $
	title='OH line detections', charsize=1, layout=[1,2,1], symsize=0.5

	cgplot, speclines[w_line1].x, speclines[w_line1].y, /overplot, psym=16, symsize=0.6, color=color1
	cgplot, speclines[w_line2].x, speclines[w_line2].y, /overplot, psym=16, symsize=0.6, color=color2


	; show the line widths
	cgplot, speclines.x, speclines.sigma, psym=16, xtit='x pixel', ytitle='line width (pixel)', $
	charsize=1, layout=[1,2,2], symsize=0.5, yra = median(speclines.sigma)*[0.5, 1.5]

	cgplot, speclines[w_line1].x, speclines[w_line1].sigma, /overplot, psym=16, symsize=0.6, color=color1
	cgplot, speclines[w_line2].x, speclines[w_line2].sigma, /overplot, psym=16, symsize=0.6, color=color2

	cgplot, [0, 2*max(speclines.x)], median(speclines.sigma) + [0,0], /overplot, thick=3, linestyle=2


	; -------------------------------------------------------
	; plot the residuals of the wavelength solution

	; calculate the wavelength solution at the location of the speclines
	flame_util_transform_direct, *cutout.rectification, x=speclines.x, y=speclines.y, lambda=lambda_model, gamma=gamma_model

	; show the residuals
	cgplot, speclines.x, 1d4 * (speclines.lambda-lambda_model), psym=16, $
		xtit='x pixel', ytitle='delta lambda (angstrom)', $
		title='Residuals of the wavelength solution (stddev: ' + $
			cgnumber_formatter( stddev( 1d4 * (speclines.lambda-lambda_model), /nan), decimals=3) + $
			' angstrom)', charsize=1, thick=3

	cgplot, speclines[w_line1].x, 1d4 * (speclines[w_line1].lambda-lambda_model[w_line1]), $
		psym=16, /overplot, color=color1
	cgplot, speclines[w_line2].x, 1d4 * (speclines[w_line2].lambda-lambda_model[w_line2]), $
			psym=16, /overplot, color=color2

	cgplot, [0, 2*max(speclines.x)], [0,0], /overplot, thick=3, linestyle=2



	; -------------------------------------------------------
	; plot the shift in wavelength as a function of spatial position

	; measured coordinates for the selected lines
	x1 = speclines[w_line1].x
	y1 = speclines[w_line1].y
	x2 = speclines[w_line2].x
	y2 = speclines[w_line2].y

	; find the reference pixel row
	y_ref = median(y1)
	x1_ref = x1[ (sort(abs(y1-y_ref)))[0] ]
	x2_ref = x2[ (sort(abs(y2-y_ref)))[0] ]

	; calculate relative offset to improve plot clarity
	offset = stddev([x1-x1_ref, x2-x2_ref], /nan)

	; plots
	cgplot, [x1-x1_ref, x2-x2_ref+offset], [y1, y2], /nodata, charsize=1, $
		xtit = 'horizontal shift from central row (pixels)', ytit='vertical pixel coordinate'
	cgplot, x1-x1_ref, y1, psym=9, /overplot, color=color1
	cgplot, x2-x2_ref+offset, y2, psym=9, /overplot, color=color2
	;cgplot, [-1d4, 1d4], y_ref+[0,0], /overplot, thick=2

	; legend
	cgtext, 0.85, 0.25, 'line at ' + strtrim(lambda1, 2) + ' um', /normal, alignment=1.0, color=color1, charsize=1
	cgtext, 0.85, 0.20, 'line at ' + strtrim(lambda2, 2) + ' um', /normal, alignment=1.0, color=color2, charsize=1
	cgtext, 0.85, 0.15, '(offset by ' + cgnumber_formatter(offset, decimals=2) + ' pixels)', /normal, alignment=1.0, color=color2, charsize=1


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_wavecal_accurate, fuel

	flame_util_module_start, fuel, 'flame_wavecal_accurate'


	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		if fuel.slits[i_slit].skip then continue

	  print, 'Accurate wavelength calibration for slit ', strtrim(fuel.slits[i_slit].number,2), ' - ', fuel.slits[i_slit].name
		print, ' '

		; handle errors by ignoring that slit
		if fuel.util.debugging eq 0 then begin
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


		for i_frame=0, n_elements(fuel.slits[i_slit].cutouts)-1 do begin

				; this slit
				this_slit = fuel.slits[i_slit]

				; the speclines measured for this slit
				speclines = *this_slit.cutouts[i_frame].speclines

				; calculate the polynomial transformation between observed and rectified frame
				flame_wavecal_2D_calibration, fuel=fuel, slit=this_slit, cutout=this_slit.cutouts[i_frame], $
					diagnostics=fuel.diagnostics, this_diagnostics=(fuel.diagnostics)[i_frame]

				; show plots of the wavelength calibration and specline identification
				cgPS_open, flame_util_replace_string(fuel.slits[i_slit].cutouts[i_frame].filename_step1, '.fits', '_plots.ps'), /nomatch
				flame_wavecal_plots, slit=this_slit, cutout=this_slit.cutouts[i_frame]

				; calculate and apply the illumination correction
				flame_wavecal_illum_correction, fuel=fuel, i_slit=i_slit, i_frame=i_frame

				cgPS_close


		endfor

	endfor


  flame_util_module_end, fuel

END
