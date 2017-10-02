

PRO flame_skysub_oneframe, fuel=fuel, cutout=cutout

	slit_filename = cutout.filename
	if cutout.illcorr_applied then slit_filename = flame_util_replace_string(cutout.filename, '_corr', '_illcorr')
	print, 'Sky subtraction for ', slit_filename

	; read in the slit image
	slit_image = mrdfits(slit_filename, 0, header, /silent)
	slit_image_sigma = mrdfits(slit_filename, 1, /silent)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(slit_image))[2]

	; how many pixels on the wavelength direction
	N_lambda_pix = (size(slit_image))[1]

	; create empty frame that will contain the wavelength solution
	wavelength_solution = dblarr(N_lambda_pix, N_spatial_pix)

	; create 2D arrays containing the observed coordinates of each pixel
	x_2d = indgen(N_lambda_pix) # replicate(1, N_spatial_pix)
	y_2d = replicate(1, N_lambda_pix) # indgen(N_spatial_pix)

	; create 2D array containing wavelength at each pixel
	wavelength_solution = flame_util_transform_coord(x_2d, y_2d, *cutout.lambda_coeff )


	;**************
	;*  B-spline  *
	;**************

	; set the order of the B-spline
	bspline_nord = 4

	; set the fraction of pixels to reject
	reject_fraction = fuel.settings.skysub_reject_fraction
	if reject_fraction LT 0.0 or reject_fraction GE 1.0 then $
		message, 'fuel.settings.skysub_reject_fraction must be between 0 and 1'

	; make image with the pixel y coordinate
	pixel_ycoord_2d = replicate(1, N_lambda_pix ) # indgen( N_spatial_pix )

	; exclude bad pixels and the top and bottom 2 pixel rows
	w_good = where(finite(slit_image) and finite(wavelength_solution) AND $
		pixel_ycoord_2d GT 2 and pixel_ycoord_2d LT N_spatial_pix-3, /null, complement=w_bad )

	; work in 1D
	pixel_wavelength = wavelength_solution[w_good]
	pixel_flux = slit_image[w_good]
	pixel_sigma = slit_image_sigma[w_good]
	pixel_ycoord = pixel_ycoord_2d[w_good]

	; sort by wavelength
	w_sort = sort(pixel_wavelength)
	pixel_wavelength = pixel_wavelength[w_sort]
	pixel_flux = pixel_flux[w_sort]
	pixel_sigma = pixel_sigma[w_sort]
	pixel_ycoord = pixel_ycoord[w_sort]

	; for each object calculate the distance from the running median
	pixel_delta = pixel_flux - median(pixel_flux, N_spatial_pix)

  ; calculate local average
	running_average = smooth(pixel_delta, N_spatial_pix*2, edge_truncate=1)

  ; calculate local average of squares
	running_average_square = smooth(pixel_delta^2, N_spatial_pix*2, edge_truncate=1)

	; calculate running variance of the pixel_delta values
	pixel_variance =  running_average_square - running_average^2

	; get the wavelengths of all the pixels in the central row as reference
	wavepoints = wavelength_solution[*,N_spatial_pix/2]

	; set the breakpoints (or nodes) for the B-spline
	if fuel.settings.skysub_bspline_oversample EQ 1.0 then breakpoints = wavepoints else begin
		if fuel.settings.skysub_bspline_oversample GT 1.0 then begin
			oversample = fix(fuel.settings.skysub_bspline_oversample)
			breakpoints = rebin(wavepoints, oversample*n_elements(wavepoints))
		endif else begin
			undersample = fix(1.0/fuel.settings.skysub_bspline_oversample)
			breakpoints = rebin(wavepoints, n_elements(wavepoints)/undersample)
		endelse
	endelse

	; set the size of the bin used for rejection, as a multiple of the breakpoint size
	bin_size = fuel.settings.skysub_reject_window * median( abs(breakpoints-shift(breakpoints, 1)) )

	; make a mask array, 1 if the pixel is not to be used
	pixel_mask = bytarr( n_elements(pixel_flux) )

	; check the number of loops for rejection
	if fuel.settings.skysub_reject_loops LT 0 or fuel.settings.skysub_reject_loops GT 15 then $
		message, 'skysub_reject_loops must be between 0 and 15'

	i_loop = 0
  while i_loop LT fuel.settings.skysub_reject_loops do begin

		; select good pixels
		pixel_wavelength_ok = pixel_wavelength[where(pixel_mask eq 0, /null)]
		pixel_flux_ok = pixel_flux[where(pixel_mask eq 0, /null)]

		; calculate B-spline model of the sky, using only good pixels
		sset = bspline_iterfit(pixel_wavelength_ok, pixel_flux_ok, nord=bspline_nord, $
			fullbkpt=breakpoints, $
			outmask=outmask, upper=5.0, lower=5.0)

		; calculate the absolute deviations of all pixels
		pixel_deviations = abs(pixel_flux - bspline_valu(pixel_wavelength, sset))

		; remove the pixels with the largest absolute deviations, in each bin
		bin_start = min(pixel_wavelength, /nan)
		while bin_start LT max(pixel_wavelength, /nan) do begin

			; select all points within this bin
			w_bin = where(pixel_wavelength GE bin_start and pixel_wavelength LT bin_start+bin_size, /null)
			if w_bin EQ !NULL then begin
				bin_start += bin_size
				continue
			endif

			; select the deviations of the good pixels in this bin
			pixel_deviations_bin = pixel_deviations[cgsetintersection(w_bin, where(pixel_mask eq 0, /null))]

			; sort the absolute value of the deviations for the good pixels within this bin
			sorted_deviations = pixel_deviations_bin[ sort(pixel_deviations_bin) ]

			; calculate the threshold corresponding to the set percentile level
			alpha = fuel.settings.skysub_reject_fraction
			max_deviation = sorted_deviations[ (1.0-alpha) * (n_elements(sorted_deviations)-1) ]

			; mask the pixels in this bin that are outliers
			pixel_mask[w_bin[where(pixel_deviations[w_bin] GT max_deviation, /null)]] = 1

			; advance to next bin
			bin_start += bin_size

		endwhile

		i_loop += 1

	endwhile

	; select points that are not masked out
	w_good = where(pixel_mask eq 0, /null)

	; print fraction of outliers rejected
	print, cgnumber_formatter( (1.0 - n_elements(w_good)/double(n_elements(pixel_mask)))*100.0, decimals=2) + '% of pixels were discarded'

	; set the range for the plot
	rel_range = fuel.settings.skysub_plot_range
	xrange=breakpoints[ [n_elements(breakpoints)*rel_range[0], n_elements(breakpoints)*rel_range[1]] ]

	; plot all pixels in a small wavelength range
	cgplot, pixel_wavelength[w_good], pixel_flux[w_good], psym=3, color='blk3', /nodata, $
		xtit='wavelength (micron)', xra=xrange, title=(strsplit(slit_filename,'/', /extract))[-1], /ynozero

	; plot all pixels in gray
	cgplot, pixel_wavelength, pixel_flux, psym=16, color='blk3', /overplot

	; plot good points in black
	cgplot, pixel_wavelength[w_good], pixel_flux[w_good], psym=16, /overplot, color='black'

	; calculate B-spline model of the sky - using the final selection of pixels
	sset = bspline_iterfit(pixel_wavelength[w_good], pixel_flux[w_good], nord=bspline_nord, $
		fullbkpt=breakpoints, $
		outmask=outmask, upper=2.0, lower=2.0)

	; wavelength axis finely sampled
	wl_axis = min(pixel_wavelength) + (max(pixel_wavelength) - min(pixel_wavelength)) * dindgen( N_lambda_pix * 10 ) / double(N_lambda_pix * 10)

	; overplot the B-spline model
	cgplot, wl_axis, bspline_valu(wl_axis, sset), /overplot, color='red'

	; ; show pixels that were masked out
	; if where(~outmask, /null) NE !NULL then $
	; 	cgplot, pixel_wavelength[where(~outmask, /null)], pixel_flux[where(~outmask, /null)], /overplot, psym=16, color='blue'

	; generate sky model for the whole slit
	sky_model = bspline_valu(wavelength_solution, sset)

	; avoid extrapolating the sky model outside the observed wavelength range
	w_extrap = where(wavelength_solution LT min(pixel_wavelength) or wavelength_solution GT max(pixel_wavelength), /null)
	if w_extrap NE !NULL then sky_model[w_extrap] = !values.d_nan

	; add a frame of NaNs to avoid extrapolation when rectifying
	sky_model[0,*] = !values.d_nan
	sky_model[-1,*] = !values.d_nan
	sky_model[*,0] = !values.d_nan
	sky_model[*,-1] = !values.d_nan

	; save sky model
	writefits, flame_util_replace_string(slit_filename, '.fits', '_skymodel.fits'), sky_model

	; save sky-subtracted image
	writefits, flame_util_replace_string(slit_filename, '.fits', '_skysub.fits'), slit_image - sky_model, header
	writefits, flame_util_replace_string(slit_filename, '.fits', '_skysub.fits'), slit_image_sigma, /append

END



; ---------------------------------------------------------------------------------------------------------------------------





PRO flame_skysub, fuel

	flame_util_module_start, fuel, 'flame_skysub'


 	; loop through all the slits & frames
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		if fuel.slits[i_slit].skip then continue

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

		for i_frame=0, fuel.util.science.n_frames-1 do $
			flame_skysub_oneframe, fuel=fuel, cutout=fuel.slits[i_slit].cutouts[i_frame]

	endfor


  flame_util_module_end, fuel

END
