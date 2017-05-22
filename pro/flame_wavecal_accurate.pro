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


PRO flame_wavecal_2D_calibration, slit=slit, cutout=cutout, $
		diagnostics=diagnostics, this_diagnostics=this_diagnostics

; This routine calculates the 2D wavelength solution and y-rectification.
; These are two mappings from the observed pixel coordinates to the rectified grid
; lambdax, gamma, where lambdax is a pixel grid linear in lambda and gamma is the
; vertical distance to the edge of the slit (taking into account warping and also vertical drift)
; The result of this routine is a pair of matrices, Klambda and Kgamma, saved in fuel.slits,
; that can be used to rectify the image via poly_2D()
;

	print, ''
	print, 'Accurate 2D wavelength solution for ', cutout.filename

	; polynomial degree for image warping
	degree=3

	; read in file to calibrate
	im = mrdfits(cutout.filename, 0, header, /silent)

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
	*cutout.rectification = {Klambda:Klambda, Kgamma:Kgamma, Kx:Kx, Ky:Ky}

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
	writefits, flame_util_replace_string(cutout.filename, '.fits', '_wavecal_2D.fits'), wavecal_accurate, hdr


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_wavecal_illum_correction, slit=slit, cutout=cutout
	;
	; use the speclines to dervie and apply an illumination correction
	; along the spatial slit axis
	;

	;
	; WARNING: THE ILLUMINATION-CORRECTED FRAMES ARE NOT ACTUALLY USED FOR NOW
	;

	speclines = *cutout.speclines

	; read in slit
	im = readfits(cutout.filename, hdr)
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

	; extract the rectification matrix for gamma
	Kgamma = (*cutout.rectification).Kgamma

	; order of polynomial
	Nord = (size(Kgamma))[1]
	xexp  = findgen(Nord)
	yexp  = findgen(Nord)

	; calculate the gamma coordinate of each OHline detection
	OHgamma = dblarr(n_elements(speclines))
	for i_line=0, n_elements(speclines)-1 do $
		OHgamma[i_line] = total(((speclines[i_line].y)^xexp # (speclines[i_line].x)^yexp ) * Kgamma)

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
  writefits, flame_util_replace_string(cutout.filename, '.fits', '_illumcorr.fits'), im, hdr


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
		w_left_num = intersect(w_left, where(lambdas_det GE median(lambdas_det[w_left]), /null) )
		w_right_num = intersect(w_right, where(lambdas_det GE median(lambdas_det[w_right]), /null) )

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

		color1 = 'blue'
		color2 = 'red'

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

	; extract the rectification matrix for lambda for this frame
	Klambda = (*cutout.rectification).Klambda

	; order of polynomial
	Nord = (size(Klambda))[1]

	xexp  = findgen(Nord)
	yexp  = findgen(Nord)
	lambda_model = dblarr(n_elements(speclines))
	for i_OH=0, n_elements(speclines)-1 do lambda_model[i_OH] = $
		slit.outlambda_min + slit.outlambda_delta * total(((speclines[i_OH].y)^xexp # (speclines[i_OH].x)^yexp ) * Klambda)

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
	cgplot, x1-x1_ref, y1, psym=16, /overplot, color=color1
	cgplot, x2-x2_ref+offset, y2, psym=16, /overplot, color=color2
	;cgplot, [-1d4, 1d4], y_ref+[0,0], /overplot, thick=2

	; legend
	cgtext, 0.85, 0.25, 'line at ' + strtrim(lambda1, 2) + ' um', /normal, alignment=1.0, color=color1, charsize=1
	cgtext, 0.85, 0.20, 'line at ' + strtrim(lambda2, 2) + ' um', /normal, alignment=1.0, color=color2, charsize=1
	cgtext, 0.85, 0.15, '(offset by ' + strtrim(offset, 2) + ')', /normal, alignment=1.0, color=color2, charsize=1

	; -------------------------------------------------------
	; show the predicted position from the accurate 2D wavelength solution

	; extract the rectification matrix for going from lambda,gamma to x,y to for this frame
	Kx = (*cutout.rectification).Kx
	Ky = (*cutout.rectification).Ky
	Nord = (size(Kx))[1]
	xexp  = findgen(Nord)
	yexp  = findgen(Nord)

	; generate the locus of points with fixed lambda and variable gamma
	gamma_array = dindgen(max(y1))

	; line1
	x1_model = dblarr(n_elements(gamma_array))
	for i=0, n_elements(x1_model)-1 do x1_model[i] = $
		total(((gamma_array[i])^xexp # ( (lambda1-slit.outlambda_min)/slit.outlambda_delta )^yexp ) * Kx)

	y1_model = dblarr(n_elements(gamma_array))
	for i=0, n_elements(y1_model)-1 do y1_model[i] = $
		total(((gamma_array[i])^xexp # ( (lambda1-slit.outlambda_min)/slit.outlambda_delta )^yexp ) * Ky)

	cgplot, x1_model-x1_model[ (sort(abs(y1_model-y_ref)))[0] ], y1_model, /overplot, color=color1, thick=3


	; line2
	x2_model = dblarr(n_elements(gamma_array))
	for i=0, n_elements(x2_model)-1 do x2_model[i] = $
		total(((gamma_array[i])^xexp # ( (lambda2-slit.outlambda_min)/slit.outlambda_delta )^yexp ) * Kx)

	y2_model = dblarr(n_elements(gamma_array))
	for i=0, n_elements(y2_model)-1 do y2_model[i] = $
		total(((gamma_array[i])^xexp # ( (lambda2-slit.outlambda_min)/slit.outlambda_delta )^yexp ) * Ky)

	cgplot, x2_model-x2_model[ (sort(abs(y2_model-y_ref)))[0] ] + offset, y2_model, /overplot, color=color2, thick=3


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_wavecal_accurate, fuel

		start_time = systime(/seconds)

	  print, ''
	  print, '-------------------------------------'
	  print, '---    flame_wavecal_accurate     ---'
	  print, '-------------------------------------'
	  print, ''


  ; avoid printing too much stuff (especially from GAUSSFIT)
  quiet_state = !QUIET
  !QUIET = 1

	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

	  print, 'Accurate wavelength calibration for slit ', strtrim(fuel.slits[i_slit].number,2), ' - ', fuel.slits[i_slit].name
		print, ' '

		for i_frame=0, n_elements(fuel.slits[i_slit].cutouts)-1 do begin

				; filename of the cutout
				slit_filename = fuel.slits[i_slit].cutouts[i_frame].filename

				; this slit
				this_slit = fuel.slits[i_slit]

				; the speclines measured for this slit
				speclines = *this_slit.cutouts[i_frame].speclines

				; calculate the polynomial transformation between observed and rectified frame
				flame_wavecal_2D_calibration, slit=this_slit, cutout=this_slit.cutouts[i_frame], $
					diagnostics=fuel.diagnostics, this_diagnostics=(fuel.diagnostics)[i_frame]

				; show plots of the wavelength calibration and specline identification
				cgPS_open, flame_util_replace_string(slit_filename, '.fits', '_plots.ps'), /nomatch
				flame_wavecal_plots, slit=this_slit, cutout=this_slit.cutouts[i_frame]

				; calculate and apply the illumination correction
					flame_wavecal_illum_correction, slit=this_slit, cutout=this_slit.cutouts[i_frame]

				cgPS_close


		endfor

		; make summary plots, including all frames for one slit

	endfor

	; revert to original !QUIET state
	!QUIET = quiet_state


	print, ''
  print, 'flame_wavecal_accurate took ', $
    cgnumber_formatter( systime(/seconds) - start_time, decimals=2), ' seconds'
  print, ''

END
