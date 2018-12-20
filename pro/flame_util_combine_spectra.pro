;
; flame_util_combine_spectra, filenames, output_filename=output_filename, difference=difference, $
;    observed_frame=observed_frame, rectified_frame=rectified_frame, nan=nan, $
;    useweights=useweights, usenoise=usenoise, usesigma=usesigma, $
;    alignment_box=alignment_box
;
; This function can be used to combine - without any resampling - two or more 2D spectra, which must be in
; the flame "output" format: FITS files with six extensions each (DATA, NOISE, SIGMA, SKY, EXPTIME, WEIGHT).
; It is used by Flame to combine the A and B frames of the same slit into an A-B;
; to further combine the A-B with the shifted B-A for the same slit (on-slit nodding);
; to combine a slit with another slit (paired-slit nodding).
; It can also be used outside of Flame to combine different observations
; of the same slit taken in different nights, etc.
;
; filenames (input) - a string array with the names of the input FITS files
; output_filename (input, optional) - a string with the name for the output FITS file
; /difference (input, optional) - if set, then the second spectrum is subtracted from the first one
;   (works only if exactly two spectra are input)
; /observed_frame (input, optional) - if set, the spatial alignment is done on the observed frame,
;    i.e. using the YCUTOUT keyword from the FITS header
; /rectified_frame (input, optional) - if set, the spatial alignment is done on the rectified frame,
;    i.e. using the gamma coordinate (from the CRVAL2 keyword in the FITS header)
; /nan (input, optional) - if set, NaNs are ignored when stacking the spectra
; /useweights (input, optional) - if set, the weight map is used when stacking the spectra
; /usenoise (input, optional) - if set, the (theoretical) noise map is used for the weights
; /usesigma (input, optional) - if set, the sigma map is used for the weights
; /alignment_box (input, optional) - the array [x0, y0, x1, y1] (in pixels) defining
;    the region to be used for the vertical alignment (cross-correlation) before combining
;    the spectra. This can be useful when only an emission line is detected. The box
;    should be large enough to contain the emission feature in all frames. Note that
;    if the alignment box is provided, then the /observed_frame and /rectified_frame
;    keywords cannot be set.
;


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************

FUNCTION flame_util_combine_spectra_align, lambda_axis=lambda_axis, $
    cube_data=cube_data, cube_sigma=cube_sigma, alignment_box=alignment_box, output_filename=output_filename

    ; open plot file
    cgPS_open, output_filename, /nomatch

    ; get the dimensions of the cube
    Nx = (size(cube_data))[1]
    Ny = (size(cube_data))[2]
    N  = (size(cube_data))[3]

    ; determine edges of alignment box
    x0 = alignment_box[0] > 0
    x1 = alignment_box[2] < (Nx-1)
    y0 = alignment_box[1] > 0
    y1 = alignment_box[3] < (Ny-1)

    ; extended range for boxcar spectrum
    delta_x = alignment_box[2] - alignment_box[0]
    x0_ext = (alignment_box[0] - 2*delta_x) > 0
    x1_ext = (alignment_box[2] + 2*delta_x) < (Nx-1)

    ; make a long cyclic array of colors
    colors = ['red6', 'blu6', 'grn6', 'red2', 'blu2', 'grn6']
    colors = ['black', colors, colors, colors]

    ; make the SNR array
    cube_snr = cube_data / cube_sigma

    ; show the boxcar extraction of the SNR around the alignment box
    s = stddev( median( total( cube_snr[x0_ext:x1_ext, y0:y1, 0], 2), 3), /nan )
    cgplot, lambda_axis[x0_ext:x1_ext], median( total( cube_snr[x0_ext:x1_ext, y0:y1, 0], 2), 3), $
      title='SNR spectra', charsize=1, xtit='wavelength'
    for i=1, N-1 do cgplot, lambda_axis[x0_ext:x1_ext], median( total( cube_snr[x0_ext:x1_ext, y0:y1, i], 2), 3), $
      /overplot, color=colors[i]

    ; show the error spectra around the alignment box
    cgplot, lambda_axis[x0_ext:x1_ext], total( cube_sigma[x0_ext:x1_ext, y0:y1, 0], 2), $
      title='error spectra', charsize=1, xtit='wavelength'
    for i=1, N-1 do cgplot, lambda_axis[x0_ext:x1_ext], total( cube_sigma[x0_ext:x1_ext, y0:y1, i], 2), $
      /overplot, color=colors[i]

    ; extract spatial profile of the galaxy
    profile = dblarr(y1-y0+1, N)
    for i=0, N-1 do profile[*,i] = total(cube_data[x0:x1,y0:y1,i], 1, /nan)

    ; show galaxy spatial profile
    cgplot, y0 + indgen(y1-y0+1), profile[*,0], title='spatial profile before alignment', $
      xtitle='pixel coordinate along the y axis', charsize=1
    for i=1, N-1 do cgplot, y0 + indgen(y1-y0+1), profile[*,i], /overplot, color=colors[i]

    ; lag for cross-correlation
    d = (y1-y0)/3
    lag = -d + indgen(2*d+1)

    ; cross-correlate the profiles with the first frame
    cc = dblarr(n_elements(lag), N-1)
    for i=1, N-1 do cc[*, i-1] = c_correlate(profile[*,0], profile[*,i], lag)

    ; show the cross-correlation
    cgplot, lag, cc[*,0], psym=-16, color=colors[1], $
      title='cross-correlation', charsize=1, xtit='vertical shift in pixels'
    for i=2, N-1 do cgplot, lag, cc[*,i-1], psym=-16, /overplot, color=colors[i-1]

    ; determine the location of the max cross-correlation for each frame
    shift = intarr(N)
    for i=1, N-1 do begin
      !NULL = max(cc[*,i-1], ind)
      shift[i] = lag[ind]
    endfor
    print, 'Results of the cross-correlation between each frame and the first one:'
    forprint, ' shift: ' + strtrim(shift,2) + ' pixels'

    ; apply the vertical shift to all frames
    new_cube_data = cube_data
    for i=0, N-1 do new_cube_data[*,*,i] = shift( cube_data[*,*,i], 0, -shift[i])

    ; re-extract the spatial profile
    newprofile = dblarr(y1-y0+1, N)
    for i=0, N-1 do newprofile[*,i] = total(new_cube_data[x0:x1,y0:y1,i], 1, /nan)

    ; show galaxy spatial profile
    cgplot, y0 + indgen(y1-y0+1), newprofile[*,0], title='spatial profile after alignment', $
      xtitle='pixel coordinate along the y axis', charsize=1
    for i=1, N-1 do cgplot, y0 + indgen(y1-y0+1), newprofile[*,i], /overplot, color=colors[i]

    ; close plot file
    cgPS_close

    ; return the pixel shifts for each frame
    return, shift


END




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_util_combine_spectra, filenames, output_filename=output_filename, difference=difference, $
    observed_frame=observed_frame, rectified_frame=rectified_frame, nan=nan, $
    useweights=useweights, usenoise=usenoise, usesigma=usesigma, $
    alignment_box=alignment_box


  ; *************************** test input ***************************

  ; number of files to combine
  N = n_elements(filenames)
  if N LT 2 then message, 'At least two file names must be input'

  ; set the output filename
  if ~keyword_set(output_filename) then output_filename='stack.fits'

  ; check the difference case
  if keyword_set(difference) and N ne 2 then $
    message, '/difference option can only be used when combining two files'

  ; check the coordinate system
  if keyword_set(observed_frame) + keyword_set(rectified_frame) + keyword_set(alignment_box) NE 1 then $
    message, 'One and only one among /observed_frame, /rectified_frame, or alignment_box must be specified'

  ; check that the alignment box has the right format
  if n_elements(alignment_box) NE 0 then begin
    if n_elements(alignment_box) NE 4 then $
      message, 'alignment box must be of the form [x0, y0, x1, y1]'
  endif


  ; *************************** load coordinates from headers ***************************

  ; make the arrays for the header info
  naxis1 = intarr(N)
  naxis2 = intarr(N)
  crval1 = fltarr(N)
  crval2 = fltarr(N)
  crpix1 = fltarr(N)
  crpix2 = fltarr(N)
  cdelt1 = fltarr(N)
  cdelt2 = fltarr(N)
  ycutout = fltarr(N)

  ; read all the headers and get the relevant information
  ; (go in reverse order so that at the end the header is that of the first frame)
  for i=N-1, 0, -1 do begin
    header = headfits(filenames[i])
    naxis1[i] = sxpar(header, 'NAXIS1')
    naxis2[i] = sxpar(header, 'NAXIS2')
		lambda_unit = strlowcase( strtrim(sxpar(header, 'CUNIT1'), 2) )
		case lambda_unit of
			'angstrom': begin ; if angstroms, convert lambda_axis to microns
				crval1[i] = sxpar(header, 'CRVAL1')/1e4
				cdelt1[i] = sxpar(header, 'CDELT1')/1e4
			end
			'micron': begin
				crval1[i] = sxpar(header, 'CRVAL1')
				cdelt1[i] = sxpar(header, 'CDELT1')
			end
			else: message, lambda_unit + ' not supported!'
		endcase
    crpix1[i] = sxpar(header, 'CRPIX1')
    crval2[i] = sxpar(header, 'CRVAL2')
    cdelt2[i] = sxpar(header, 'CDELT2')
    crpix2[i] = sxpar(header, 'CRPIX2')
    ycutout[i] = sxpar(header, 'YCUTOUT')
    print, 'reading header for ', filenames[i], $
      '  -  size: ', strtrim( naxis1[i], 2), ' x ', strtrim(naxis2[i], 2)
  endfor


  ; *************************** check headers ***************************

  ; check that the wavelength scale is the same for all files
  if n_elements(uniq(crpix1)) NE 1 then message, 'Wavelength axes are different: CRPIX1 = ' + strjoin(string(crpix1))
  if n_elements(uniq(cdelt1)) NE 1 then message, 'Wavelength axes are different: CDELT1 = ' + strjoin(string(cdelt1))
  if cdelt1[0] eq 0.0 then message, 'Wavelength axes are not defined!'


  ; *************************** set up the common wavelength grid ***************************

  ; find extrema of wavelengths
  lambda_min = min(crval1)
  lambda_max = max(crval1 + naxis1*cdelt1)

  ; for each frame, how many pixels do we leave blank at the beginning of the grid?
  pixel_padding = (crval1-lambda_min)/cdelt1[0]
  print, 'Frames need to be shifted by the following pixel amount:', reverse(pixel_padding)

  ; check that the wavelength grids can be shifted and matched
  if total( abs(pixel_padding-round(pixel_padding)) GT 0.05 ) NE 0.0 then $
    print, 'WARNING! Wavelength grids are misaligned!'

  ; number of pixels in the output wavelength grid
  Nx = round( (lambda_max-lambda_min)/cdelt1[0] )

  ; create the lambda axis
  lambda_axis = lambda_min + (dindgen(Nx) - lambda_min + 1d) * cdelt1[0]

  ; update the header
  sxaddpar, header, 'CRVAL1', lambda_min

  ; define the x-pixel range that will be occupied by each frame in the common wavelength grid
  x0 = round(pixel_padding)
  x1 = round(pixel_padding + naxis1 - 1)


  ; *************************** set up the common spatial grid ***************************

  ;; NB: this assumes that all frames are vertically misplaced by integer number of pixels!!!

  ; spatial alignment in the observed frame, i.e. on the detector
  if keyword_set(observed_frame) then begin

    ; find the extrema of the ycutout values
    ycutout_min = min(ycutout)
    ycutout_max = max(ycutout+naxis2)-1

    ; update the YCUTOUT keyword in the header
    sxaddpar, header, 'YCUTOUT', ycutout_min

    ; number of pixels in the output spatial grid
    Ny = round(ycutout_max - ycutout_min + 1)

    ; define the y-pixel range that will be occupied by each frame in the common wavelength grid
    y0 = round(ycutout-ycutout_min)
    y1 = y0 + naxis2 - 1

  endif

  ; spatial alignment in the rectified frame, i.e. on the sky
  if keyword_set(rectified_frame) then begin

  ; find the extrema of the CRVAL2 (i.e. gamma) values
    gamma_min = min(crval2)
    gamma_max = max(crval2+naxis2)-1

    ; update the CRVAL2 keyword in the header
    sxaddpar, header, 'CRVAL2', gamma_min

    ; number of pixels in the output spatial grid
    Ny = round( gamma_max - gamma_min + 1 )

    ; define the y-pixel range that will be occupied by each frame in the common wavelength grid
    y0 = round(crval2-gamma_min)
    y1 = y0 + naxis2 - 1

  endif

  ; spatial alignment using alignment box - start by simply aligning the detectors
  if keyword_set(alignment_box) then begin

    ; number of pixels in the output spatial grid is simply the number of pixels in the first
    Ny = max(naxis2)

    ; for each frame, simply take all the pixels
    y0 = intarr(N)
    y1 = y0 + naxis2 - 1

  endif

  ; make sure that the type of alignment has been specified
  if n_elements(y0) eq 0 then message, 'Either /observed_frame of /rectified_frame or alignment_box must be set'


  ; *************************** load spectra into cube ***************************

  ; create the 3D arrays
  cube_data = fltarr(Nx, Ny, N) + !values.d_nan
  cube_error = fltarr(Nx, Ny, N) + !values.d_nan
  cube_sigma = fltarr(Nx, Ny, N) + !values.d_nan
  cube_sky = fltarr(Nx, Ny, N) + !values.d_nan
  cube_exptime = fltarr(Nx, Ny, N) + !values.d_nan
  cube_weight = fltarr(Nx, Ny, N) + !values.d_nan
  print, 'The output will be ', strtrim( Nx, 2), ' x ', strtrim(Ny, 2)

  ; loop through all files and load them in
  ; (go in reverse order so that at the end the extension headers are those of the first frame)
  for i_file=N-1, 0, -1 do begin

    print, 'Loading ', filenames[i_file]

    ; read in all the extensions
    im_data =    mrdfits(filenames[i_file], 0, hdr_data, /silent)
  	im_error =   mrdfits(filenames[i_file], 1, hdr_error, /silent)
  	im_sigma =   mrdfits(filenames[i_file], 2, hdr_sigma, /silent)
  	im_sky =     mrdfits(filenames[i_file], 3, hdr_sky, /silent)
  	im_exptime = mrdfits(filenames[i_file], 4, hdr_exptime, /silent)
  	im_weight =  mrdfits(filenames[i_file], 5, hdr_weight, /silent)

    ; load the 2D spectra into the right position of the cube
    cube_data[    x0[i_file]:x1[i_file] , y0[i_file]:y1[i_file] , i_file ] = im_data
    cube_error[   x0[i_file]:x1[i_file] , y0[i_file]:y1[i_file] , i_file ] = im_error
    cube_sigma[   x0[i_file]:x1[i_file] , y0[i_file]:y1[i_file] , i_file ] = im_sigma
    cube_sky[     x0[i_file]:x1[i_file] , y0[i_file]:y1[i_file] , i_file ] = im_sky
    cube_exptime[ x0[i_file]:x1[i_file] , y0[i_file]:y1[i_file] , i_file ] = im_exptime
    cube_weight[  x0[i_file]:x1[i_file] , y0[i_file]:y1[i_file] , i_file ] = im_weight

  endfor


  ; *************************** align spectra with cross-correlation ***************************

  if keyword_set(alignment_box) then begin

    ; name of ps file
    ps_filename = flame_util_replace_string(output_filename, '.fits', '.ps')

    ; determine pixel shifts for each frame
    shift = flame_util_combine_spectra_align(lambda_axis=lambda_axis, $
      cube_data=cube_data, cube_sigma=cube_sigma, alignment_box=alignment_box, output_filename=ps_filename)

    ; apply pixel shift to each frame
    for i=0, N-1 do begin
      cube_data[*,*,i]    = shift( cube_data[*,*,i], 0, -shift[i])
      cube_error[*,*,i]   = shift( cube_error[*,*,i], 0, -shift[i])
      cube_sigma[*,*,i]   = shift( cube_sigma[*,*,i], 0, -shift[i])
      cube_sky[*,*,i]     = shift( cube_sky[*,*,i], 0, -shift[i])
      cube_exptime[*,*,i] = shift( cube_exptime[*,*,i], 0, -shift[i])
      cube_weight[*,*,i]  = shift( cube_weight[*,*,i], 0, -shift[i])
    endfor

  endif


  ; *************************** combine spectra ***************************

  ; for the A-B case the stack is not an average, so there are different rules
  if keyword_set(difference) then begin

    ; flip the sign of the B frame
    cube_data[*,*,1] *= -1d

    ; subtract B from A
    stack_data = total(cube_data, 3, nan=nan)

    ; make the error spectra
    stack_error = sqrt( total( cube_error^2, 3, nan=nan) )
    stack_sigma = sqrt( total( cube_sigma^2, 3, nan=nan) )

    ; use the sky of the A frame
  	stack_sky = cube_sky[*,*,0]

    ; use the exposure time and the weight map of the A frame
	  stack_exptime = cube_exptime[*,*,0]
	  stack_weight = cube_weight[*,*,0]

  endif else begin

    ; weights for the stacking
    weights = cube_data
    weights[*] = 1.0
    if keyword_set(useweights) then weights = cube_weight
    if keyword_set(usenoise) then   weights = 1.0/cube_error^2
    if keyword_set(usesigma) then   weights = 1.0/cube_sigma^2

  	; mean-stack the frames
  	stack_data = total(weights*cube_data, 3, nan=nan) / total(weights, 3, nan=nan)

  	; make the error spectra
  	stack_error = sqrt( total( weights^2 * cube_error^2, 3, nan=nan) ) / total(weights, 3, nan=nan)
  	stack_sigma = sqrt( total( weights^2 * cube_sigma^2, 3, nan=nan) ) / total(weights, 3, nan=nan)

  	; mean-stack the sky model
  	stack_sky = total(weights*cube_sky, 3, nan=nan) / total(weights, 3, nan=nan)

    ; add the exposure time and the weights
  	stack_exptime = total(cube_exptime, 3, nan=nan)
  	stack_weight = total(cube_weight, 3, nan=nan)

  endelse


  ; *************************** write output ***************************


	; write out FITS file
	writefits, output_filename, stack_data, header
  writefits, output_filename, stack_error, hdr_error, /append
  writefits, output_filename, stack_sigma, hdr_sigma, /append
  writefits, output_filename, stack_sky, hdr_sky, /append
  writefits, output_filename, stack_exptime, hdr_exptime, /append
  writefits, output_filename, stack_weight, hdr_weight, /append
	print, output_filename, ' written'


END
