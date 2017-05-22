

PRO flame_skysub_oneframe, slit_filename=slit_filename, rectification=rectification, $
	lambda_0=lambda_0, delta_lambda=delta_lambda

	print, 'Sky subtraction for ', slit_filename

	; read in the slit image
	slit_image = mrdfits(slit_filename, 0, header, /silent)
	slit_image_sigma = mrdfits(slit_filename, 1, /silent)

	; how many pixels on the spatial direction
	N_spatial_pix = (size(slit_image))[2]

	; how many pixels on the wavelength direction
	N_lambda_pix = (size(slit_image))[1]

	; calculate the wavelength solution using the rectification matrices

	; make the arrays that will have the new coordinate of each pixel
	lambda = slit_image * 0.0
	gamma = lambda

	; order of polynomial
	Nord = (size(rectification.Klambda))[1]
	xexp  = findgen(Nord)
	yexp  = findgen(Nord)

	for ix=0.0,N_lambda_pix-1 do $
		for iy=0.0,N_spatial_pix-1 do $
			lambda[ix,iy] = lambda_0 + delta_lambda * total(((iy)^xexp # (ix)^yexp ) * rectification.Klambda)

	for ix=0.0,N_lambda_pix-1 do $
		for iy=0.0,N_spatial_pix-1 do $
			gamma[ix,iy] = total(((iy)^xexp # (ix)^yexp ) * rectification.Kgamma)

	wavelength_solution = lambda

	;**************
	; 2D B-spline *
	;**************

	pixel_ycoord_2d = replicate(1, N_lambda_pix ) # indgen( N_spatial_pix )

	w_good = where(finite(slit_image) and finite(wavelength_solution) AND $
		pixel_ycoord_2d GT 2 and pixel_ycoord_2d LT N_spatial_pix-3, /null, complement=w_bad )

	pixel_wavelength = wavelength_solution[w_good]
	pixel_flux = slit_image[w_good]
	pixel_ycoord = pixel_ycoord_2d[w_good]

	; simple Poisson error - NOT USING THIS ANYMORE DURING ITERATIVE BSPLINE FITTING
;	inverse_variance = 1.0/pixel_flux
;	sigma = sqrt(pixel_flux)

	; use the wavelengths of all the pixels in the central row as breakpoints (or nodes) for the B-spline
	breakpoints = wavelength_solution[*,N_spatial_pix/2]

	; plot all pixels in a small wavelength range (from 1/4 to 2/4 of entire observed range)
	cgplot, pixel_wavelength, pixel_flux, psym=3, xtit='wavelength (micron)', $
		xra=breakpoints[ [n_elements(breakpoints)*1/4, n_elements(breakpoints)*2/4] ], $
		title=(strsplit(slit_filename,'/', /extract))[-1]

	; sset = bspline_iterfit(pixel_wavelength, pixel_flux, nord=4, $
	; 	fullbkpt=breakpoints, x2=pixel_ycoord, npoly=1, $
	; 	outmask=outmask, upper=4.0, lower=4.0 )

	sset = bspline_iterfit(pixel_wavelength, pixel_flux, nord=4, $
		fullbkpt=breakpoints, $
		outmask=outmask, upper=4.0, lower=4.0 )

	; wavelength axis finely sampled
	wl_axis = min(pixel_wavelength) + (max(pixel_wavelength) - min(pixel_wavelength)) * dindgen( N_lambda_pix * 10 ) / double(N_lambda_pix * 10)

	; overplot the B-spline model at each row
	; cgplot, wl_axis, bspline_valu(wl_axis, sset, x2=0.0001), /overplot, color='red' ; bspline_valu bug: setting i_row=0 is interpreted as keyword not set
	; for i_row=1, N_spatial_pix-1 do cgplot, wl_axis, bspline_valu(wl_axis, sset, x2=i_row), /overplot, color='red'
	cgplot, wl_axis, bspline_valu(wl_axis, sset), /overplot, color='red'

	; show pixels that were masked out
	if where(~outmask, /null) NE !NULL then $
		cgplot, pixel_wavelength[where(~outmask, /null)], pixel_flux[where(~outmask, /null)], /overplot, psym=16, color='blue'

	; generate sky model for the whole slit
	;sky_model = bspline_valu(wavelength_solution, sset, x2=pixel_ycoord_2d)
	sky_model = bspline_valu(wavelength_solution, sset)

	; save sky model
	writefits, flame_util_replace_string(slit_filename, '.fits', '_skymodel.fits'), sky_model

	; save sky-subtracted image
	writefits, flame_util_replace_string(slit_filename, '.fits', '_skysub.fits'), slit_image - sky_model, header
	writefits, flame_util_replace_string(slit_filename, '.fits', '_skysub.fits'), slit_image_sigma, /append

END



; ---------------------------------------------------------------------------------------------------------------------------





PRO flame_skysub, fuel

		start_time = systime(/seconds)

	  print, ''
	  print, '-------------------------------------'
	  print, '---         flame_skysub          ---'
	  print, '-------------------------------------'
	  print, ''


	; extract the slits structures
	slits = fuel.slits

 	; loop through all the slits
	for i_slit=0, n_elements(slits)-1 do begin

		this_slit = fuel.slits[i_slit]

		for i_frame=0, n_elements(this_slit.cutouts)-1 do begin

			this_cutout = this_slit.cutouts[i_frame]

			slit_filename = this_cutout.filename
			rectification = *(this_cutout.rectification)
			lambda_0 = this_slit.outlambda_min
			delta_lambda = this_slit.outlambda_delta

			flame_skysub_oneframe, slit_filename=slit_filename, $
				rectification=rectification, lambda_0=lambda_0, delta_lambda=delta_lambda

		endfor

	endfor


	print, ''
  print, 'flame_skysub took ', $
    cgnumber_formatter( systime(/seconds) - start_time, decimals=2), ' seconds'
  print, ''

END
