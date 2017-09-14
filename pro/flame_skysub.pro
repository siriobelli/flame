

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
	; 2D B-spline *
	;**************

	; set the order of the B-spline
	bspline_nord = fuel.settings.skysub_bspline_order

	; set the fraction of bright pixels to reject
	reject_fraction = fuel.settings.skysub_reject_fraction
	if reject_fraction LT 0.0 or reject_fraction GT 1.0 then $
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

	; use the wavelengths of all the pixels in the central row as breakpoints (or nodes) for the B-spline
	breakpoints = wavelength_solution[*,N_spatial_pix/2]

	; set the range for the plot
	rel_range = fuel.settings.skysub_plot_range
	xrange=breakpoints[ [n_elements(breakpoints)*rel_range[0], n_elements(breakpoints)*rel_range[1]] ]

	; plot all pixels in a small wavelength range, in gray
	cgplot, pixel_wavelength, pixel_flux, psym=3, xtit='wavelength (micron)', $
		xra=xrange, color='blk3', $
		title=(strsplit(slit_filename,'/', /extract))[-1]

	; calculate B-spline model of the sky
	sset = bspline_iterfit(pixel_wavelength, pixel_flux, nord=bspline_nord, $
		fullbkpt=breakpoints, invvar=1.0/pixel_variance, $
		outmask=outmask, upper=3.0, lower=3.0)

	; calculate the deviations from the model
	pixel_deviations = pixel_flux -  bspline_valu(pixel_wavelength, sset)

	; remove half of the pixels, selecting the ones with the largest positive deviations

	; make a mask array, 1 if the pixel is not to be used
	pixel_mask = bytarr( n_elements(pixel_flux) )

	for i_bin=0, n_elements(breakpoints)-2 do begin

		; select the points within this bin
		w_bin = where(pixel_wavelength GE breakpoints[i_bin] and pixel_wavelength LT breakpoints[i_bin+1], /null)
		if w_bin EQ !NULL then continue

		; sort the deviation values for the points within this bin
		sorted_deviations = pixel_deviations[w_bin[ sort(pixel_deviations[w_bin]) ]]

		; calculate the threshold corresponding to the set percentile
		threshold_deviation = sorted_deviations[ (1.0-reject_fraction) * (n_elements(sorted_deviations)-1) ]

		; mask the points above the threshold
		pixel_mask[w_bin[ where(pixel_deviations[w_bin] GT threshold_deviation) ]] = 1

	endfor

	; select points that are not masked out
	w_good = where(pixel_mask eq 0, /null)

	; plot them in black
	cgplot, pixel_wavelength[w_good], pixel_flux[w_good], psym=3, /overplot, color='black'

	; calculate B-spline model of the sky - now removing the masked points
	sset = bspline_iterfit(pixel_wavelength[w_good], pixel_flux[w_good], nord=bspline_nord, $
		fullbkpt=breakpoints, $
		outmask=outmask, upper=2.0, lower=2.0)

	; wavelength axis finely sampled
	wl_axis = min(pixel_wavelength) + (max(pixel_wavelength) - min(pixel_wavelength)) * dindgen( N_lambda_pix * 10 ) / double(N_lambda_pix * 10)

	; overplot the B-spline model at each row
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
