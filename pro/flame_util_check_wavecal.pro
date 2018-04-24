PRO flame_util_check_wavecal, slit=slit, diagnostics=diagnostics, ascii_filename=ascii_filename
;
; Utility to plot information about the speclines and the quality of the wavelength
; solution for all cutouts relative to one slit. Used both in flame_wavecal and in
; flame_checkdata. Does not plot anything if arcs were used.
;
; If ascii_filename is specified, then a list of the lines used and the wavelength
; residuals is written to a text file
;

	; read in and stack together all the relevant info for each single speclines, for all frames
	;----------------------------------------------------------------------------------

	line_frame = []
	line_lambda = []
	line_width = []
	line_resid = []

	; number of frames
	Nfr = n_elements(slit.cutouts)

	; loop through each frame
	for i_frame=0, Nfr-1 do begin
		speclines = *slit.cutouts[i_frame].speclines

		; check that there are speclines measured
		if speclines EQ !NULL then continue

		; select only lines we can trust
		speclines = speclines[where(speclines.trust_lambda, /null)]

		; check that there are speclines measured
		if speclines EQ !NULL then continue

		; calculate residuals
		lambda = flame_util_transform_coord(speclines.x, speclines.y, *slit.cutouts[i_frame].lambda_coeff )

		line_frame = [line_frame, replicate(i_frame, n_elements(speclines))]
		line_lambda = [line_lambda, speclines.lambda]
		line_width = [line_width, speclines.sigma]
		line_resid = [line_resid, 1d4*(lambda-speclines.lambda)]	; in angstrom

	endfor

	; total number of speclines for all frames
	Nlines = n_elements(line_frame)

	; if there are no speclines (because arcs were used), then we are done
	if Nlines eq 0 then return


	; take care of the frame numbers for the x axis titles (see also flame_diagnostics)
	;----------------------------------------------------------------------------------

	; slightly perturb the frame number to improve clarity in the plot
	line_frame_pert = line_frame + 0.1*randomn(seed, Nlines)

  ; these are the values used to label the x axis
  xtickname = strtrim(fix(diagnostics.frame_num), 2)
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


	; plot 1: residuals vs wavelength
	;----------------------------------------------------------------------------------

	; calculate reasonable range for plot
	resid_sorted = line_resid[sort(line_resid)]
	resid_range = [ resid_sorted[0.05*Nlines], resid_sorted[0.95*Nlines] ]
	resid_range += [-1.0, 1.0] * 0.5*(resid_range[1]-resid_range[0])

	; plot distribution of residuals
	cgplot, line_lambda, line_resid, psym=16, symsize=0.4, charsize=1, color='blk4', $
		yra=resid_range, /ysty, ytit='residuals (angstrom)', xtit='wavelength (micron)'

	; identify different lines
  uniq_lambda = line_lambda[UNIQ(line_lambda, SORT(line_lambda))]
	uniq_res_median = fltarr(n_elements(uniq_lambda))
	uniq_res_mad = fltarr(n_elements(uniq_lambda))
	uniq_detections = lonarr(n_elements(uniq_lambda))

	; overplot median for each emission line
	for i_line=0, n_elements(uniq_lambda)-1 do begin

		; select speclines that belong to this frame
		line_resid_0 = line_resid[ where(line_lambda eq uniq_lambda[i_line], /null) ]
		uniq_detections[i_line] = n_elements(line_resid_0)

		; find percentiles to plot
		resid_sorted = line_resid_0[sort(line_resid_0)]
		top68 = resid_sorted[0.84*n_elements(resid_sorted)]
		bottom_68 = resid_sorted[0.16*n_elements(resid_sorted)]
		uniq_res_median[i_line] = median(line_resid_0) ; i.e., 50th percentile

		; find median absolute deviation
		; NOTE: these are the deviations from the median of the residuals, not from zero!
		uniq_res_mad[i_line] = median(abs( line_resid_0-uniq_res_median[i_line] ))

		; overplot MAD
		cgplot, [ uniq_lambda[i_line] ], [ uniq_res_median[i_line] ], $
			/overplot, psym=3, color='blue', $
			/err_clip, err_yhigh=uniq_res_mad[i_line], err_ylow=uniq_res_mad[i_line]

		; overplot median and 68 percentile
		cgplot, [ uniq_lambda[i_line] ], [ uniq_res_median[i_line] ], $
			/overplot, psym=16, color='red', $
			/err_clip, err_yhigh=top68-uniq_res_median[i_line], err_ylow=uniq_res_median[i_line]-bottom_68

	endfor

	cgplot, [0.5*min(line_lambda), 2.0*max(line_lambda)], [0, 0], /overplot, thick=3, linestyle=2



	; plot 2: line widths vs frame number
	;----------------------------------------------------------------------------------

	; calculate reasonable range for plot
	width_sorted = line_width[sort(line_width)]
	width_range = [ width_sorted[0.05*Nlines], width_sorted[0.95*Nlines] ]
	width_range += [-1.0, 1.0] * 0.5*(width_range[1]-width_range[0])

	; plot distribution of linewidths
	cgplot, line_frame_pert, line_width, psym=16, symsize=0.4, charsize=1, color='blk4', $
		yra=width_range, /ysty, ytit='line width (raw pixels)', $
		xra=[-1, Nfr], /xsty, xtit='frame number', $
    xtickv=xtickv, xticks=n_elements(xtickv)-1, xtickname=xtickname, xminor=xminor

	; overplot median for each frame
	for i_frame=0, Nfr-1 do begin

		; select speclines that belong to this frame
		line_width_0 = line_width[ where(line_frame eq i_frame, /null) ]

		; find percentiles to plot
		width_sorted = line_width_0[sort(line_width_0)]
		top68 = width_sorted[0.84*n_elements(width_sorted)]
		bottom_68 = width_sorted[0.16*n_elements(width_sorted)]

		; median absolute deviation
		mad = median( abs(line_width_0-median(line_width_0)) )

		; overplot MAD
		cgplot, [ i_frame ], [ median(line_width_0) ], $
		/overplot, psym=3, color='blue', $
		/err_clip, err_yhigh=mad, err_ylow=mad

		; overplot median and 68 percentile
		cgplot, [ i_frame ], [ median(line_width_0) ], $
		/overplot, psym=16, color='red', $
		/err_clip, err_yhigh=top68-median(line_width_0), err_ylow=median(line_width_0)-bottom_68

	endfor


	; plot 3: residuals vs frame number
	;----------------------------------------------------------------------------------

	; calculate reasonable range for plot
	resid_sorted = line_resid[sort(line_resid)]
	resid_range = [ resid_sorted[0.05*Nlines], resid_sorted[0.95*Nlines] ]
	resid_range += [-1.0, 1.0] * 0.5*(resid_range[1]-resid_range[0])

	; plot distribution of residuals
	cgplot, line_frame_pert, line_resid, psym=16, symsize=0.4, charsize=1.2, color='blk4', $
		xthick=4, ythick=4, $
		yra=resid_range, /ysty, ytit='residuals (angstrom)', $
		xra=[-1, Nfr], /xsty, xtit='frame number', $
    xtickv=xtickv, xticks=n_elements(xtickv)-1, xtickname=xtickname, xminor=xminor

	; median absolute deviation for each frame
	mad = dblarr(Nfr)

	; overplot median for each frame
	for i_frame=0, Nfr-1 do begin

		; select speclines that belong to this frame
		line_resid_0 = line_resid[ where(line_frame eq i_frame, /null) ]

		; median absolute deviation
		mad = median ( abs(line_resid_0) )

		; find percentiles to plot
		resid_sorted = line_resid_0[sort(line_resid_0)]
		top68 = resid_sorted[0.84*n_elements(resid_sorted)]
		bottom_68 = resid_sorted[0.16*n_elements(resid_sorted)]

		; overplot MAD
		cgplot, [ i_frame ], [ 0 ], /overplot, psym=3, color='blue', $
		/err_clip, err_yhigh=mad, err_ylow=mad, err_thick=4

		; overplot median and 68 percentile
		cgplot, [ i_frame ], [ median(line_resid_0) ], /overplot, psym=16, color='red', $
		/err_clip, err_yhigh=top68-median(line_resid_0), err_ylow=median(line_resid_0)-bottom_68, err_thick=6

	endfor

	cgplot, [-1, Nfr], [0, 0], /overplot, thick=3, linestyle=2


	; write text file with line list
	;----------------------------------------------------------------------------------

	if n_elements(ascii_filename) EQ 1 then $
		forprint, uniq_lambda, uniq_detections, uniq_res_median, uniq_res_mad, format='F12,I12,F12.4,F12.4', $
     	textout=ascii_filename, $
     	comment = '# wavelength    number of detections    median residual (angstrom)    scatter (MAD, angstrom)'


END
