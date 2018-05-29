

PRO flame_util_measure_seeing, spec_filename, psfilename=psfilename, pixel_scale=pixel_scale
;
; Measure and plot the effective seeing from the final 2D spectrum.
; Also track seeing and vertical position as a function of wavelength.
;
; spec_filename (input): name of a FITS file output by flame (e.g. slit01-12345_A-B.fits)
; psfilename (input): name of the output ps file with plots of seeing and position
; pixel_scale (input): spatial pixel scale in units of arcsec/pixel
;


  cgPS_open, psfilename, /nomatch

	; load 2D spectrum
	; --------------------------------------

	; load file
	ref_spec = mrdfits(spec_filename, 0, header, /silent)

	; get the wavelength calibration from the header
 	lambda_unit = strlowcase( strtrim(sxpar(header, 'CUNIT1'), 2) )
	lambda_axis = sxpar(header, 'CRVAL1') + sxpar(header,'CDELT1') * $
		( findgen(sxpar(header,'NAXIS1')) )


 	; median profile and fit
 	; -----------------------------

	; starting guess: seeing is 1 arcsec
	approx_seeing = 1.0							; FWHM, in arcsec
	approx_seeing /= pixel_scale		; FWHM, in pixels

	; make the y-axis in units of gamma
	yaxis_cutout = dindgen( (size(ref_spec))[2] ) + sxpar(header, 'CRVAL2')

	; ; select only the region around the reference star (which is by definition at gamma=0)
	; w_region = where(abs(yaxis_cutout) LT 7.0*approx_seeing, /NULL)
	; if n_elements(w_region) LT 4 then message, 'Could not detect the trace of the reference star'

	; extract the spatial profile
	ref_profile_cutout = median(ref_spec, dimension=1)

	; use the full profile
	ref_profile = ref_profile_cutout
	yaxis = yaxis_cutout

	; select the positive points (to exclude the negative A-B traces)
	w_fit = where(ref_profile GE 0.0 and finite(ref_profile), /null)
	if w_fit EQ !NULL then w_fit = where( finite(ref_profile), /null)

	; fit a gaussian to the integrated profile
	est_peak = max(ref_profile, ind_peak)
	est_coeff = [ est_peak, yaxis[ind_peak], approx_seeing/2.36, 0.0]
  fit_result = gaussfit(yaxis[w_fit], ref_profile[w_fit], ref_coeff, nterms=4, $
    estimates = est_coeff, $
    chisq=chisq, sigma=coeff_err)

	; calculate seeing in arcsec
	median_seeing = 2.355 * ref_coeff[2] * pixel_scale
	print, 'The final effective seeing calculated from the reference star is ' + $
	 	cgnumber_formatter(median_seeing, decimals=2) + ' arcsec.'
	print, ''

	; plot the median profile
	cgplot, yaxis, ref_profile, charsize=1, xtit='y pixel coordinate', color='blk3', $
		ytit='median flux', title='reference star: median profile (seeing = ' + $
		cgnumber_formatter(median_seeing, decimals=2) + ' arcsec)'

	; overplot the points that have been used
	cgplot, yaxis[w_fit], ref_profile[w_fit], psym=16, /overplot

	; overplot Gaussian fit
	xaxis = yaxis[0] + dindgen(300)/299.0*(yaxis[-1]-yaxis[0])
	cgplot, xaxis, ref_coeff[0] * exp( -0.5* ( (xaxis-ref_coeff[1])/ref_coeff[2] )^2 ) + ref_coeff[3], $
		/overplot, color='red', thick=3


	; fit profile as a function of wavelength
	; ----------------------------------------

	; set the bin size, in pixels, along the wavelength direction
	binsize = 100
	starting_pixel = 0

	; dimensions
	N_pixel_x = (size(ref_spec))[1]
	N_pixel_y = (size(ref_spec))[2]

	; empty arrays for seeing and position measurement
	coord_x = []
	seeing = []
	center = []

	while starting_pixel LT N_pixel_x do begin

    ; extract the bin
    end_pixel = min([starting_pixel + binsize - 1, N_pixel_x-1])
    cutout_bin = ref_spec[starting_pixel : end_pixel, *]

    ; spatial profile
    profile = median(cutout_bin, dimension=1)

		; range to fit
		wfit = where( abs(yaxis_cutout-ref_coeff[1]) LT 3.0*ref_coeff[2] and finite(profile), /null )

		; if there's nothing to fit, skip this bin
		if n_elements(wfit) GT 4 then begin

			; fit a Gaussian
			fit_result = gaussfit(yaxis_cutout[wfit], profile[wfit], coeff, nterms=4, $
	    	estimates=ref_coeff, $
	    	chisq=chisq, sigma=coeff_err)

			; save the result of the fit
			coord_x = [coord_x, 0.5*(starting_pixel+end_pixel)]
			seeing = [seeing, 2.355 * coeff[2] * pixel_scale]
			center = [center, coeff[1]]

		endif

    ; advance to next bin
    starting_pixel += binsize

  endwhile

	; plot seeing as a function of wavelength
	cgplot, coord_x, seeing, thick=4, charsize=1, $
		ytit='seeing (arcsec)', $
		yra=[min(seeing, /nan)-0.1, max(seeing, /nan)+0.1], $
		xra=[0, N_pixel_x], xsty=1+8, xthick=4, ythick=4, $
		position = [0.15, 0.55, 0.9, 0.9], xtickformat='(A1)'

	; top x-axis with pixel coordinate
	cgaxis, xaxis=1, xra = [ 0 , N_pixel_x ], xsty=1, charsize=1, $
		xtit='x-coordinate (pixel)'

	; show the median seeing
	cgplot, [0, N_pixel_x], [0,0]+median_seeing, /overplot, thick=3, linestyle=2

	; make center relative to initial position
	center -= center[0]

	; show centering as a function of wavelength
	cgplot, coord_x, center, thick=4, charsize=1, $
		xtit='x-coordinate (pixel)', ytit='center position (pixel)', $
		yra=[min(center, /nan)-1.0, max(center, /nan)+1.0], $
		xra=[0, N_pixel_x], xsty=1+4, xthick=4, ythick=4, $
		position = [0.15, 0.15, 0.9, 0.55], /noerase

	; bottom x-axis with wavelength
	cgaxis, xaxis=0, xra = [ lambda_axis[0] , lambda_axis[-1] ], xsty=1, charsize=1, $
		xtit='wavelength (' + lambda_unit + ')'

	; overplot zero line
	cgplot, [0, N_pixel_x], [0,0], /overplot, thick=3, linestyle=2


	cgPS_close

END
