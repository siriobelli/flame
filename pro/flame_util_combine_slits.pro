
PRO flame_util_combine_slits, filenames, output=output, alignment_box=alignment_box, $
    sky_filenames=sky_filenames, noweights=noweights, signs=signs
;
; combine the flame output from different observations of the same object
; if specified, alignment_box sets the region [x0, y0, x1, y1] (in pixels)
; to be used for the vertical alignment before combining the spectra
;
; signs is an arrays of +1 and -1
;

  ; *************************** load headers ***************************

  ; number of files to combine
  N = n_elements(filenames)

  ; check the signs
  if ~keyword_set(signs) then signs = 1 + intarr(N)
  if n_elements(signs) NE N then $
    message, 'signs must have the same number of elements as the array of file names'

  if ~keyword_set(output) then output='combined_frames.fits'

  ; make the arrays for the header info
  naxis1 = intarr(N)
  naxis2 = intarr(N)
  cunit1 = strarr(N)
  crval1 = fltarr(N)
  crval2 = fltarr(N)
  crpix1 = fltarr(N)
  crpix2 = fltarr(N)
  cdelt1 = fltarr(N)

  ; read all the headers and get the relevant information
  for i=0, N-1 do begin
    header = headfits(filenames[i])
    naxis1[i] = sxpar(header, 'NAXIS1')
    naxis2[i] = sxpar(header, 'NAXIS2')
    cunit1[i] = strtrim(sxpar(header, 'CUNIT1'), 2)
    crval1[i] = sxpar(header,'CRVAL1')
    crval2[i] = sxpar(header,'CRVAL2')
    crpix1[i] = sxpar(header,'CRPIX1')
    crpix2[i] = sxpar(header,'CRPIX2')
    cdelt1[i] = sxpar(header,'CDELT1')
    print, 'reading header for ', filenames[i], $
      '  -  size: ', strtrim( naxis1[i], 2), ' x ', strtrim(naxis2[i], 2)
  endfor

  ; *************************** check headers ***************************

  ; check that the units are the same
  if n_elements(uniq(cunit1)) NE 1 then message, 'Wavelength units are different: ', cunit1

  ; check that the wavelength axis is the same for all files
  if n_elements(uniq(crval1)) NE 1 then message, 'Wavelength axes are different: CRVAL1 = ', crval1
  if n_elements(uniq(crpix1)) NE 1 then message, 'Wavelength axes are different: CRPIX1 = ', crpix1
  if n_elements(uniq(cdelt1)) NE 1 then message, 'Wavelength axes are different: CDELT1 = ', cdelt1
  print, 'Wavelength axis are identical for all files.'

  ; find the smallest common dimensions
  Nx = min(naxis1)
  Ny = min(naxis2)
  print, 'The output will be ', strtrim( Nx, 2), ' x ', strtrim(Ny, 2)

  ; *************************** read in spectra ***************************

  ; create the 3D arrays
  cube = fltarr(Nx, Ny, N)
  cube_sigma = fltarr(Nx, Ny, N)
  cube_sky = fltarr(Nx, Ny, N)

  ; create the lambda axis
  lambda_axis = crval1[0] + (findgen(Nx) - crpix1[0] + 1d) * cdelt1[0]

  ; loop through all files and read them in
  for i=0, N-1 do begin

    print, 'Loading ', filenames[i]

    ; read in the spectrum and error spectrum
    spec = signs[i] * mrdfits(filenames[i], 0, /silent)
    spec_sigma = mrdfits(filenames[i], 1, /silent)
    if keyword_set(sky_filenames) then sky = mrdfits(sky_filenames[i], 0, /silent)

    ; trim them and save them in the array
    cube[*,*,i] = spec[ 0:Nx-1 , 0:Ny-1 ]
    cube_sigma[*,*,i] = spec_sigma[ 0:Nx-1 , 0:Ny-1 ]
    if keyword_set(sky_filenames) then cube_sky[*,*,i] = sky[ 0:Nx-1 , 0:Ny-1 ]

  endfor

  ; make the SNR array
  cube_snr = cube / cube_sigma

  ; *************************** align spectra ***************************

  if keyword_set(alignment_box) then begin

    ; open plot file
    cgPS_open, 'flame_combine_observations.ps', /nomatch

    ; determine edges of alignment box
    x0 = alignment_box[0] > 0
    x1 = alignment_box[2] < Nx-1
    y0 = alignment_box[1] > 0
    y1 = alignment_box[3] < Ny-1

    ; extended range for boxcar spectrum
    delta_x = alignment_box[2] - alignment_box[0]
    x0_ext = alignment_box[0] - 2*delta_x > 0
    x1_ext = alignment_box[2] + 2*delta_x < Nx-1

    ; make a long cyclic array of colors
    colors = 'red' + strtrim(indgen(8)+1,2)
    colors = [colors, colors, colors]

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

    ; show the sky spectra around the alignment box (if provided)
    if keyword_set(sky_filenames) then begin
      cgplot, lambda_axis[x0_ext:x1_ext], total( cube_sky[x0_ext:x1_ext, y0:y1, 0], 2), $
        title='sky spectra', charsize=1, xtit='wavelength'
      for i=1, N-1 do cgplot, lambda_axis[x0_ext:x1_ext], total( cube_sky[x0_ext:x1_ext, y0:y1, i], 2), $
        /overplot, color=colors[i]
    endif

    ; extract spatial profile of the galaxy
    profile = dblarr(y1-y0+1, N)
    for i=0, N-1 do profile[*,i] = total(cube[x0:x1,y0:y1,i], 1, /nan)

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
    forprint, filenames, ' shift: ' + strtrim(shift,2) + ' pixels'

    ; apply the vertical shift to all frames
    for i=0, N-1 do cube[*,*,i] = shift( cube[*,*,i], 0, -shift[i])
    for i=0, N-1 do cube_sigma[*,*,i] = shift( cube_sigma[*,*,i], 0, -shift[i])
    if keyword_set(sky_filenames) then $
      for i=0, N-1 do cube_sky[*,*,i] = shift( cube_sky[*,*,i], 0, -shift[i])

    ; re-extract the spatial profile
    newprofile = dblarr(y1-y0+1, N)
    for i=0, N-1 do newprofile[*,i] = total(cube[x0:x1,y0:y1,i], 1, /nan)

    ; show galaxy spatial profile
    cgplot, y0 + indgen(y1-y0+1), newprofile[*,0], title='spatial profile after alignment', $
      xtitle='pixel coordinate along the y axis', charsize=1
    for i=1, N-1 do cgplot, y0 + indgen(y1-y0+1), newprofile[*,i], /overplot, color=colors[i]

    ; close plot file
    cgPS_close

  endif

  ; *************************** combine spectra ***************************

  ; option 1: weighted mean of each pixel
  if ~keyword_set(noweights) then begin
    cube_weight = 1.0/cube_sigma^2
    spec_tot = total( cube*cube_weight , 3, /nan) / total( cube_weight , 3, /nan)
    spec_sigma_tot = 1.0 / sqrt( total( cube_weight , 3, /nan) )
    if keyword_set(sky_filenames) then $
      sky_tot = total( cube_sky*cube_weight , 3, /nan) / total( cube_weight , 3, /nan)
  endif

  ; ; option 2: weighted mean, all pixels in one frame have the same weight
  ; frame_typical_sigma = median( median(cube_sigma, dimension=1), dimension=1)
  ; cube_weight = cube * 0.0
  ; for i=0, N-1 do cube_weight[*,*,i] = 1.0 / frame_typical_sigma[i]^2
  ; spec_tot = total( cube*cube_weight , 3, /nan) / total( cube_weight , 3, /nan)
  ; spec_sigma_tot = sqrt( total( cube_sigma^2 * cube_weight^2 , 3, /nan) / (total( cube_weight , 3, /nan))^2 )

  ; option 3: arithmetic mean
  if keyword_set(noweights) then begin
    spec_tot = mean(cube, dimension=3, /nan)
    spec_sigma_tot = sqrt( total( cube_sigma^2, 3, /nan) ) / float(N)
    if keyword_set(sky_filenames) then $
      sky_tot = mean(cube_sky, dimension=3, /nan)
  endif

  ; read the header of the first input frame
  header0 = headfits(filenames[0])

  ; write spectra
	writefits, output, spec_tot, header0
  writefits, output, spec_sigma_tot, /append
    if keyword_set(sky_filenames) then $
      writefits, flame_util_replace_string(output, '.fits', '_sky.fits'), sky_tot, header0


END
