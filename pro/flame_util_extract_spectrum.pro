;
; flame_util_extract_spectrum, filename, $
;     boxcar_aperture=boxcar_aperture, boxcar_auto=boxcar_auto, $
;     optimal_profile=optimal_profile, optimal_gaussian=optimal_gaussian, $
;     ycrop=ycrop, output_filename=output_filename
;
; This function extracts a 1D spectrum from a 2D spectrum. The input must be a file
; that was output by Flame (i.e., 2D FITS file with six extensions: flux, noise,
;  sigma, sky, exptime, weight)
; It is used by Flame at the end of the data reduction, but can also be used
; as a standalone routine, for example to extract the spectra after having combined
; observations from different nights. The brightest trace in the frame is detected and
; used to define the aperture.
; The output is a FITS file with a structure containing the fields lambda (in angstrom), flux,
; ivar (inverse variance), ivar_empirical, and sky, and is compatible with SpecPro.
; A ps file showing the spatial profile of the trace, extracted spectrum, and SNR is also produced.
;
; filename (input) - a string with the name of the input FITS file
;
; NOTE: one and only one of the following four keywords must be set:
; boxcar_aperture - an array of the form [a, b] specifying the vertical pixel coordinates
;          to be considered for the boxcar extraction (starting from zero, including extrema).
; /boxcar_auto - if set, the brightest trace is automatically detected, and the boxcar aperture
;          is set to the central +/- 2 sigma of the Gaussian fit.
; /optimal_profile - if set, optimal extraction is performed, using as a weight the profile of the brightest object.
; /optimal_gaussian - if set, optimal extraction is performed, using a Gaussian fit
;          to the profile as the weight, instead of the profile itself
;
; ycrop (input, optional) - an array of the form [a, b] specifying the vertical pixel range
;          to be considered (useful if many traces are present in the 2D spectrum)
; output_filename (input, optional) - a string array with the name of the output FITS file
;


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_util_extract_spectrum_boxcar, structure_2d, aperture

	; make nice structure for 1D spectrum
	structure_1d = { $
		lambda: structure_2d.lambda, $
		flux: total(structure_2d.flux[*,min(aperture):max(aperture)], 2, /nan), $
		ivar: 1. / total(1./structure_2d.ivar[*,min(aperture):max(aperture)], 2, /nan), $
		ivar_empirical: 1. / total(1./structure_2d.ivar_empirical[*,min(aperture):max(aperture)], 2, /nan), $
		sky: total(structure_2d.sky[*,min(aperture):max(aperture)], 2, /nan) }

    return, structure_1d

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_util_extract_spectrum_optimal, structure_2d, weight1d

  ; normalize weight by its integral
  weight1d /= total(weight1d, /nan)

  ; extend weight to 2d frame
  weight2d = replicate(1, (size(structure_2d.flux))[1] ) # weight1d

	; make nice structure for 1D spectrum
	structure_1d = { $
		lambda: structure_2d.lambda, $
		flux: total(weight2d*structure_2d.ivar*structure_2d.flux, 2, /nan) / total(weight2d^2*structure_2d.ivar, 2, /nan), $
		ivar: total(weight2d^2*structure_2d.ivar, 2, /nan), $
		ivar_empirical: total(weight2d^2*structure_2d.ivar_empirical, 2, /nan), $
		sky: total(weight2d*structure_2d.sky, 2, /nan) / total(weight2d^2, 2, /nan) }

    return, structure_1d

END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_util_extract_spectrum_gaussfit, profile, y_1d
  ;
  ; given a profile and corresponding pixel coordinates, perform
  ; Gaussian fit and return an array with the best-fit parameters
  ;

  ; cut away non positive pixels
	profile[where(profile LT 0.0, /null)] = 0.0

	; cut away NaNs
	profile[where(~finite(profile), /null)] = 0.0

	; median filter the profile to get rid of noise peaks
	if n_elements(profile) GT 20 then profile_sm = median(profile, 9) else profile_sm = median(profile, 5)

	; estimate parameters of the Gaussian
	est_peak = max(profile_sm, est_center_index)
  ;est_center = 0.5*n_elements(profile)
	;est_center = total(y_1d*profile, /nan) / total(profile, /nan)
	est_sigma = 2.0
	est_cont = median(profile_sm)
	est_param = [est_peak, y_1d[est_center_index], est_sigma, est_cont]

	; Gaussian fit
	gaussian_model = gaussfit( y_1d, profile, gauss_param, nterms=4, $
		estimates=est_param, sigma=gauss_err, chisq=chisq )

	; calculate the rms of the profile, once the peak has been subtracted
	profile_residuals = profile - gaussian_model
	w_trace = where( abs(y_1d-gauss_param[1]) LT 5.0*gauss_param[2], /null)
	if w_trace EQ !NULL then $
	 	profile_rms = stddev(profile_residuals, /nan) else $
		profile_rms = stddev(profile_residuals[w_trace], /nan)

  ; in case of non-detection
	if ~finite(chisq) or chisq LE 0.0 or  $			; check that chi square makes sense
		gauss_param[0] LE 0.0 or $	; check that the peak of the Gaussian is positive
		gauss_param[0] LT 5.0*gauss_err[0] or $ 	; check that the SNR is high
		gauss_param[0] LT 5.0*profile_rms or $ 	; check that the peak is high compared to the noise in the profile
	 	gauss_param[1] LT min(y_1d) or gauss_param[1] GT max(y_1d) or $ 			; check that the center of the Guassian is in the observed range
		gauss_param[2] LT 0.1 or gauss_param[2] GT n_elements(profile) $		 	; check that the Gaussian width makes sense
	then return, !NULL

  return, gauss_param

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************

PRO flame_util_extract_spectrum, filename, $
    boxcar_aperture=boxcar_aperture, boxcar_auto=boxcar_auto, $
    optimal_profile=optimal_profile, optimal_gaussian=optimal_gaussian, $
    ycrop=ycrop, output_filename=output_filename


  ; preliminary operations
  ; ----------------------------------------------------------------------------

  ; number of input files
  if n_elements(filename) NE 1 then message, 'Input must be one and only one filename'

  ; check consistency of keywords
  extraction_type = [keyword_set(boxcar_aperture), keyword_set(boxcar_auto), $
    keyword_set(optimal_profile), keyword_set(optimal_gaussian)]
  if total(abs(extraction_type)) EQ 0 then message, 'Extraction type not specified!'
  if total(abs(extraction_type)) GT 1 then message, 'Only one type of extraction can be specified!'
  if keyword_set(boxcar_aperture) and keyword_set(ycrop) then message, 'Only one among boxcar_aperture and ycrop can be specified!'

  ; set the output filename
  if ~keyword_set(output_filename) then begin
    file_extension = (strsplit(filename, '.', /extract))[-1]
    output_filename = flame_util_replace_string(filename, '.' + file_extension, '_spec1d.fits')
  endif

  ; set filename for the ps file with the profile
  file_extension = (strsplit(output_filename, '.', /extract))[-1]
  ps_filename = flame_util_replace_string(output_filename, '.' + file_extension, '.ps')


	; read in data
	; ----------------------------------------------------------------------------

  ; read in 2d spectrum
  spec2d = mrdfits(filename, 0, header, /silent)

  ; read in 2d theoretical error spectrum
  err2d = mrdfits(filename, 1, /silent)

	; read in 2d empirical error spectrum
	sig2d = mrdfits(filename, 2, /silent)

	; read in 2d sky
  sky2d = mrdfits(filename, 3, /silent)

  ; create ivar images
  ivar2d_th = 1d/err2d^2
	ivar2d_emp = 1d/sig2d^2

  ; make spatial axis
  y_1d = findgen( (size(spec2d))[2] )

  ; if required, crop input data to specified window
  if keyword_set(ycrop) then begin

    ; chek range of specified window
    if min(ycrop) LT 0 or max(ycrop) GE n_elements(spec2d[0,*]) then $
      message, 'the specified ycrop is outside the observed pixel range!'

    ; crop all frames
    spec2d = spec2d[*, min(ycrop) : max(ycrop)]
    ivar2d_th = ivar2d_th[*, min(ycrop) : max(ycrop)]
    ivar2d_emp = ivar2d_emp[*, min(ycrop) : max(ycrop)]
    sky2d = sky2d[*, min(ycrop) : max(ycrop)]
    y_1d = y_1d[min(ycrop) : max(ycrop)]

  endif

  ; read in wavelength axis
	lambda_unit = strlowcase( strtrim(sxpar(header, 'CUNIT1'), 2) )
	lambda_1d = sxpar(header,'CRVAL1') + (findgen(sxpar(header,'NAXIS1')) - sxpar(header,'CRPIX1') + 1d) * sxpar(header,'CDELT1')

	; if angstroms, convert lambda_axis to microns
	if lambda_unit eq 'angstrom' then $
		lambda_1d /= 1e4 $
	else if lambda_unit ne 'micron' then message, lambda_unit + ' not supported!'

	; make nice structure 2D spectra
	structure_2d = { $
		lambda: lambda_1d, $
		flux: spec2d, $
		ivar: ivar2d_th, $
		ivar_empirical: ivar2d_emp, $
		sky: sky2d }


	; derive and plot the spatial profile
	; ----------------------------------------------------------------------------

  cgPS_open, ps_filename, /nomatch

	; obtain the median profile
  profile = median(spec2d, dimension=1)

  ; show spatial profile
	cgplot, y_1d, profile, charsize=1, $
		ytitle='Total flux', $
    position = [0.15, 0.55, 0.95, 0.95]


  ; 1) BOXCAR APERTURE
	; ----------------------------------------------------------------------------

  if keyword_set(boxcar_aperture) then begin

    ; check that the aperture is well defined
    if n_elements(boxcar_aperture) ne 2 then message, 'boxcar_aperture must be an array with two elements!'
    if min(boxcar_aperture) LT 0 or max(boxcar_aperture) GE n_elements(profile) then $
      message, 'the specified boxcar_aperture is outside the observed pixel range!'

    ; plot aperture on top of profile
  	cgplot, [0,0] + min(boxcar_aperture), [-1d5, 1d5], /overplot, color='red', thick=3
  	cgplot, [0,0] + max(boxcar_aperture), [-1d5, 1d5], /overplot, color='red', thick=3

    ; boxcar extraction
    structure_1d = flame_util_extract_spectrum_boxcar(structure_2d, boxcar_aperture)


  ; OTHERWISE, FIND AND FIT THE TRACE OF THE OBJECT
	; ----------------------------------------------------------------------------
  endif else begin

    gauss_param = flame_util_extract_spectrum_gaussfit(profile, y_1d)

    ; if there is no trace, then we are done
    if gauss_param EQ !NULL then begin
      print, 'no object was detected'
      print, ''
      cgps_close
      return
    endif

  	; overplot the Gaussian fit
  	x_axis = min(y_1d) + n_elements(y_1d) * dindgen(300)/299.0
  	cgplot, x_axis, gauss_param[0] * exp( -0.5*( (x_axis-gauss_param[1])/gauss_param[2] )^2 ) + gauss_param[3], $
  		/overplot, color='red'

  endelse


  ; 2) BOXCAR AUTO
	; ----------------------------------------------------------------------------

  if keyword_set(boxcar_auto) then begin

  	; define boxcar aperture (round up to include all relevant pixels)
    auto_aperture = [ floor(gauss_param[1] - 2.0*gauss_param[2]), $
        1 + floor(gauss_param[1] + 2.0*gauss_param[2])]

    ; check that the aperture is well defined
    if min(auto_aperture) LT min(y_1d) or max(auto_aperture) GE max(y_1d) then $
        message, 'the automatic boxcar aperture is outside the observed pixel range!'

  	; overplot aperture
  	cgplot, [0,0] + auto_aperture[0], [-1d5, 1d5], /overplot, color='red3', linestyle=2
  	cgplot, [0,0] + auto_aperture[1], [-1d5, 1d5], /overplot, color='red3', linestyle=2

    ; boxcar extraction
    structure_1d = flame_util_extract_spectrum_boxcar(structure_2d, auto_aperture-y_1d[0])

  endif


  ; 3) OPTIMAL USING PROFILE
	; ----------------------------------------------------------------------------

  if keyword_set(optimal_profile) then begin

    ; make weight from profile
    weight = profile

    ; cut the profile beyond three sigma
	  weight[ where( abs(y_1d-gauss_param[1]) GT 3.0*gauss_param[2], /null ) ] = 0.0

    ; optimal extraction
    structure_1d = flame_util_extract_spectrum_optimal(structure_2d, profile)

  endif


  ; 4) OPTIMAL USING GAUSSIAN FIT
	; ----------------------------------------------------------------------------

  if keyword_set(optimal_gaussian) then begin

    ; make weight from Gaussian fit (without adding the zero point)
    weight = gauss_param[0] * exp( -0.5*( (y_1d-gauss_param[1])/gauss_param[2] )^2 )

    ; cut the profile beyond three sigma
	  weight[ where( abs(y_1d-gauss_param[1]) GT 3.0*gauss_param[2], /null ) ] = 0.0

    ; optimal extraction
    structure_1d = flame_util_extract_spectrum_optimal(structure_2d, weight)

  endif


  ; MAKE PLOTS
	; ----------------------------------------------------------------------------

  ; make fake weight for boxcar extraction, just for plotting purposes
  if keyword_set(boxcar_aperture) then w_ap = where(y_1d GE min(boxcar_aperture) and y_1d LE max(boxcar_aperture), /null)
  if keyword_set(boxcar_auto) then w_ap = where(y_1d GE min(auto_aperture) and y_1d LE max(auto_aperture), /null)
  if keyword_set(boxcar_aperture) or keyword_set(boxcar_auto) then begin
    weight = fltarr(n_elements(y_1d))
    weight[w_ap] = 0.9
  endif

  ; plot weight used for extraction
  cgplot, y_1d, weight, psym=0, thick=3, $
    charsize=1, xtitle='Position along the slit (pixel)', ytitle='Relative weight', $
    position = [0.15, 0.15, 0.95, 0.50], /noerase

	; plot extracted spectrum
	cgplot, structure_1d.lambda, median(structure_1d.flux, 7), $
		charsize=1, ytitle='Extracted flux', $
		title = 'black: observed flux, blue; observed uncertainty', xrange=xrange, /xstyle, $
      position = [0.15, 0.55, 0.95, 0.95]
 	cgplot, structure_1d.lambda, 1.0/sqrt(structure_1d.ivar), /overplot, color='blue'

	; plot SNR
	cgplot, structure_1d.lambda, median(structure_1d.flux*sqrt(structure_1d.ivar), 7), $
		charsize=1, xtitle='Wavelength (um)', ytitle='Signal-to-noise ratio', $
		xrange=xrange, /xstyle, $
      position = [0.15, 0.15, 0.95, 0.50], /noerase

  cgPS_close



	; write FITS file
	; ----------------------------------------------------------------------------

	; convert wavelength to angstrom in order to be compatible with SpecPro
	structure_1d.lambda *= 1d4

	; make new FITS header, with units
	header_output = header
 	sxdelpar, header_output, 'naxis'
	sxdelpar, header_output, 'naxis1'
	sxdelpar, header_output, 'naxis2'
	sxdelpar, header_output, 'bitpix'
	sxaddpar, header_output, 'TUNIT1', 'Angstrom'
	sxaddpar, header_output, 'TUNIT2', 'electron / (pix s)'

	; write structure to FITS file
	mwrfits, structure_1d, output_filename, header_output, /create

	print, ''
  print, '1D spectrum extracted: ' + output_filename
	print, ''



END
