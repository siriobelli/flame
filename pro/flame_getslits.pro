
; ****************************************************************************************


FUNCTION flame_getslits_findshift, image, approx_top, approx_bottom
  ;
  ; find the shift between the expected slit positions and the real ones
  ;

  ; integrate the 2D spectrum along the wavelength direction, and smooth
  x = median(total(image, 1, /nan), 15)

  ; find the edges of the slits (1d profile where positive peaks mark the beginning and negative peaks mark the end of the slit)
  edges = x - shift(x, 4)
  edges[ where( abs(edges) LT stddev(edges, /nan), /null ) ] = 0
  edges[ where(edges GT 0.0, /null) ] =  1.0
  edges[ where(edges LT 0.0, /null) ] = -1.0

  ; these are the expected edges
  expected_edges = dblarr( n_elements(edges) )
  expected_edges[approx_top] = -1.
  expected_edges[approx_bottom] = 1.

  ; cross-correlate to find the shift between expected and measured slit edges
  lag = indgen(400)-200
  crosscorr = c_correlate(edges, expected_edges, lag)
  max_crosscorr = max( crosscorr, max_ind, /nan)
  delta = lag[max_ind]

  ; return the shift
  return, delta

END


; ****************************************************************************************


FUNCTION flame_getslits_trace_smear, array, filter_length, up=up, down=down
  ;
  ; little pretty function that "smears" a 2D array along the vertical direction.
  ; the direction is set with /up and /down keywords (default is up)
  ;

  ; create the kernel
  if keyword_set(down) then $
    kernel = [ replicate(0, filter_length-1), replicate(1, filter_length) ] $
  else $
    kernel = [ replicate(1, filter_length), replicate(0, filter_length-1) ]

  ; make it vertical
  kernel = transpose(kernel)

  ; smear
  array_smeared = convol(array, kernel, /center)

  return, array_smeared

END


; ****************************************************************************************


FUNCTION flame_getslits_trace_continuum, image_in, approx_edge, top=top, bottom=bottom

  ; check keywords
  if ~keyword_set(top) AND ~keyword_set(bottom) then message, 'Please select either /top or /bottom'
  if keyword_set(top) AND keyword_set(bottom) then message, 'Please select either /top or /bottom'

  ; first, mask out bright sources and sky lines
  image = image_in
  image[where(image GT median(image) + 0.5*stddev(image, /nan), /null)] = !values.d_nan

  mwrfits, image, 'test.fits', /create

  ; from now on, the code assumes that we are interested in finding the top edge of a slit
  ; if instead we want the bottom edge, simply take the negative of the image
  if keyword_set(bottom) then image *= -1.0

  ; read in the x-size of the image
  sz = size(image)
  N_pixel_x = sz[1]

  ; how big, in the y direction, is the cutout?
  cutout_size = 14      ; even number please

  ; for each bin along the spectral direction, detect the edge using the sky background
  binsize = 100
  starting_pixel = N_pixel_x/2
  x_edge = []
  y_edge = []

  ; starting guess for the ycoord
  previous_ycoord = approx_edge

  ; starting from the center (because we assume that's where the approximate values refer to)
  ; first do the right half
  while starting_pixel LT N_pixel_x do begin

    ; extract the bin
    end_pixel = min([starting_pixel + binsize - 1, N_pixel_x-1])
    cutout_bin = image[starting_pixel : end_pixel, previous_ycoord - cutout_size/2: previous_ycoord + cutout_size/2]

    ; spatial profile
    profile = median(cutout_bin, dimension=1)

    ; detect the edge
    derivative = abs(profile - shift(profile,1))
    derivative[0] = 0
    derivative[-1] = 0
    peak = max(derivative, peak_location)

    ; add the amount of pixels left out by the cutout
    peak_location += previous_ycoord - cutout_size/2

    ; save x and y coordinate of the detection
    x_edge = [ x_edge, 0.5*(starting_pixel + end_pixel) ]
    y_edge = [ y_edge, peak_location ]

    ; save the y coordinate for the next bin
    previous_ycoord = peak_location

    ; advance to next bin (to the right)
    starting_pixel += binsize

  endwhile

  ; and then, starting from the center, to the left
  starting_pixel = N_pixel_x/2
  previous_ycoord = approx_edge

  while starting_pixel GT 0 do begin

    ; extract the bin
    end_pixel = min([starting_pixel + binsize - 1, N_pixel_x-1])
    cutout_bin = image[starting_pixel : end_pixel, previous_ycoord - cutout_size/2: previous_ycoord + cutout_size/2]

    ; spatial profile
    profile = median(cutout_bin, dimension=1)

    ; detect the edge
    derivative = abs(profile - shift(profile,1))
    derivative[0] = 0
    derivative[-1] = 0
    peak = max(derivative, peak_location)

    ; add the amount of pixels left out by the cutout
    peak_location += previous_ycoord - cutout_size/2

    ; save x and y coordinate of the detection
    x_edge = [ x_edge, 0.5*(starting_pixel + end_pixel) ]
    y_edge = [ y_edge, peak_location ]

    ; save the y coordinate for the next bin
    previous_ycoord = peak_location

    ; advance to next bin (to the right)
    starting_pixel -= binsize

  endwhile

  ; fill in the missing pixels with NaNs
  y_edge_full = fltarr(N_pixel_x) + !values.d_nan
  y_edge_full[x_edge] = y_edge

  ; return array with the y-coordinates of the edges
  return, y_edge_full


END


; ****************************************************************************************



FUNCTION flame_getslits_trace_skylines, image, approx_edge, top=top, bottom=bottom

  ;
  ; Given a frame containing bright sky emission lines, it traces the edge of a slit.
  ; One of the two keywords /top and /bottom must be specified
  ; The output is a 1D array with the y-coordinate of the edge at each x position.
  ; Pixels with no detected edges are set to NaN
  ;

  ; check keywords
  if ~keyword_set(top) AND ~keyword_set(bottom) then message, 'Please select either /top or /bottom'
  if keyword_set(top) AND keyword_set(bottom) then message, 'Please select either /top or /bottom'

  ; read in the x-size of the image
  sz = size(image)
  N_pixel_x = sz[1]

  ; how big, in the y direction, is the cutout?
  cutout_size = 34      ; even number please

  ; let's extract a cutout centered on the edge
  cutout_bottom_ycoord = approx_edge-cutout_size/2
  cutout = image[ * , cutout_bottom_ycoord : cutout_bottom_ycoord + cutout_size - 1 ]

  ; from now on, the code assumes that we are interested in finding the top edge of a slit
  ; if instead we want the bottom edge, simply flip vertically the cutout
  if keyword_set(bottom) then cutout = reverse(cutout, 2)

  ; get rid of negative values and NaNs
  cutout[ where(cutout LT 0.0 or ~finite(cutout), /null) ] = 0.0

  ; extract a spectrum from the fiducial slit regions
  sky_spectrum = median(cutout[*,0:cutout_size/4], dimension=2)

  ; subtract the continuum from the sky spectrum (particularly important in the K band)
  sky_spectrum_padded = [replicate(0d,200), sky_spectrum, replicate(0d,200)]
  sky_spectrum_continuum = (median(sky_spectrum_padded, 200))[200:-201]
  sky_spectrum -= sky_spectrum_continuum

  ; measure the sky flux in between the OH lines
  mmm, sky_spectrum, level_betweenlines, sigma_betweenlines

  ; consider as a sky line all those pixels that are more than 3 sigma above the normal level of between-lines sky
  w_OH = where(sky_spectrum GT level_betweenlines + 3.0*sigma_betweenlines, /null)
  if w_OH eq !NULL then w_OH = -1 ; it won't be used anyway

  ; for each pixel column, take the ratio of each pixel to the fiducial sky value in the slit
  sky_spectrum_2d = sky_spectrum # replicate(1, cutout_size)
  flux_ratio = cutout / sky_spectrum_2d
  flux_ratio[where(sky_spectrum_2d eq 0.0, /null)] = 1.0

  ; now delete all pixel with flux less than 50% of what is found in the OH line of the corresponding pixel column
  flux_ratio[ where(flux_ratio LT 0.50, /null) ] = 0.0

  ; because of bad pixels, there might be some zeroes on the OH lines. Let's median smooth everything along the x axis
  flux_ratio_sm = median(flux_ratio, 3, dimension=2)

  ; now let's set to 1 all the "bright" pixels, and we end up with a binary array
  binary_mask = flux_ratio_sm
  binary_mask[ where(binary_mask GT 0.0, /null) ] = 1.0

  ; FIND EDGES
  ; here is my definition of an edge: you need at least X consecutive "bright" pixels just below the edge,
  ; and at least X consecutive "dark" pixels above the edge.
  ; set X:
  filter_length = 5

  ; the first condition is: there are only ones among the last X pixels
  ; easy way to do this is to "smear" the binary mask vertically
  binary_mask_smeared_up = flame_getslits_trace_smear( binary_mask, filter_length, /up)

  ; if the sum is less than 1+1+1+... X times, then there was a zero
  condition_A = binary_mask_smeared_up*0.0
  condition_A[where(binary_mask_smeared_up eq 5.0, /null)] = 1.0

  ; the second condition is: there are only zeroes among the next X pixels
  binary_mask_smeared_down = flame_getslits_trace_smear( binary_mask, filter_length, /down)

  ; if the sum is not zero, then there was a one
  condition_B = binary_mask_smeared_down*0.0
  condition_B[where(binary_mask_smeared_down eq 0.0, /null)] = 1.0

  ; finally, at the edge you have one pixel satisying condition A, and the next one satisfying condition_B
  both_conditions = condition_A * shift(condition_B, [0,-1])

  ; select all the candidate "edge" points
  candidate_edge_ind =  where(both_conditions, /null)

  ; split into x and y coordinates
  candidate_edge_xy = array_indices(binary_mask, candidate_edge_ind)
  candidate_edge_x = reform( candidate_edge_xy[0,*] )
  candidate_edge_y = reform( candidate_edge_xy[1,*] )

  ; make axis of x-coordinate
  x_axis = indgen(N_pixel_x)

  ; make axis of y coordinates for the edges
  y_edge = fltarr(N_pixel_x)

  ; at each x coordinate, see if there are edge candidates and if more than one,
  ; take the one with the lowest y coordinate (i.e. the one more toward the slit)
  for i_x=0L,N_pixel_x-1 do begin
    w = where( candidate_edge_x eq i_x, /null )
    if w NE !NULL then y_edge[i_x] = min( candidate_edge_y[w] )
  endfor

  ; exclude those artificial edges at the edge of the cutout
  w_borders = where( y_edge LT filter_length+1 or y_edge GT cutout_size-filter_length-1, /null)
  y_edge[w_borders] = 0.0

  ; count how many unique OH lines are there
  w_uniq = where( w_OH-1 - shift(w_OH,1) NE 0, /null) ; get rid of consecutive pixels
  Nlines = n_elements(w_uniq)

  if Nlines LT 9 then begin

    message, 'I found only ' + strtrim(Nlines,2) + $
      ' OH lines! Edge tracing should be based on sky background.', /informational

  endif else begin

    print, strtrim(Nlines,2) + ' OH lines found. Tracing edge...'
    if keyword_set(top) then print, 'top edge at about ' + strtrim(approx_edge, 2)
    if keyword_set(bottom) then print, 'bottom edge at about ' + strtrim(approx_edge, 2)

    ; select the good edge measurements
    w_ok = cgsetintersection( w_OH, where( y_edge ne 0.0, /null ) )

    ; set all the non-good edge measurements to NaNs
    y_edge[ cgsetdifference( indgen(N_pixel_x), w_ok ) ] = !values.d_NaN

  endelse

  ; transform into real y coordinates (not just within the cutout anymore)
  if keyword_set(bottom) then $
    real_y_edge = cutout_bottom_ycoord + cutout_size - float(y_edge) $
  else $
    real_y_edge = cutout_bottom_ycoord + float(y_edge)

  ; return array with the y-coordinates of the edges
  return, real_y_edge

END


; ****************************************************************************************

FUNCTION flame_getslits_crosscorr, image_in, expected_bottom, expected_top, rectified_image=rectified_image
;
; cross-correlate the pixel rows of an image to find the approximate edges of a slit
; return [bottom, top] pixel coordinates
;

  ; vertical size of the image
  N_pixel_y = (size(image_in))[2]

  ; work with integer numbers
  expected_top = round(expected_top)
  expected_bottom = round(expected_bottom)

  ; cross-correlation does not work well with NaNs
  image = image_in
  image[ where( ~finite(image), /null ) ] = 0.0

  ; fiducial center
  fiducial_center = (expected_top + expected_bottom) / 2

  ; use this as reference (should definitely be within the slit)
  reference_spectrum = median(image[*,fiducial_center-1 : fiducial_center+1], dimension=2)

  ; these are the shifts considered (in pixels, relative to the previous row)
  lag = -15 + indgen(31)

  ; these array will contain the cross-correlation results for each row
  cc_shift = intarr( N_pixel_y )
  cc_peak = fltarr( N_pixel_y )

  ; cross-correlate all the rows from the center going up
  for i_row=fiducial_center, N_pixel_y-1 do begin

    ; reference shift from previous rows
    ref_shift = median(cc_shift[i_row-3:i_row-1])

    ; cross-correlate this row with the reference spectrum
    cc = c_correlate( reference_spectrum, image[*,i_row], lag + ref_shift)

    ; find the peak of the cross-correlation
    cc_peak[i_row] = max(cc, indmax, /nan)
    if ~finite(cc_peak[i_row]) then continue

    ; store this shift and move on to the next row
    cc_shift[i_row] = lag[indmax] + ref_shift

  endfor

    ; again but now from the center going down
    for i_row=fiducial_center, 0, -1 do begin

      ; reference shift from previous rows
      ref_shift = median(cc_shift[i_row+1:i_row+3])

      ; cross-correlate this row with the reference spectrum
      cc = c_correlate( reference_spectrum, image[*,i_row], lag + ref_shift)

      ; find the peak of the cross-correlation
      cc_peak[i_row] = max(cc, indmax)

      ; store this shift and move on to the next row
      cc_shift[i_row] = lag[indmax] + ref_shift

    endfor

  ; smooth out the cross-correlation results in case one row failed
  cc_peak = median(cc_peak, 2)
  cc_shift = median(cc_shift, 2)

  ; using the expected slit edges, select a conservative region that should be inside the slit
  i_row = indgen(N_pixel_y)
  margin = 0.15*(expected_top-expected_bottom)
  w_inslit = where( i_row GT expected_bottom+margin and i_row LT expected_top-margin, /null )

  ; find the top edge
  for i_top = fiducial_center, N_pixel_y-1 do $
    if cc_peak[i_top] LT 0.5 * median(cc_peak[w_inslit]) then break
  i_top -= 1

  ; find the bottom edge
  for i_bottom = fiducial_center, 0, -1 do $
    if cc_peak[i_bottom] LT 0.5 * median(cc_peak[w_inslit]) then break
  i_bottom += 1

  ; rectify the slit
  rectified_image = image_in
  for i=i_bottom, i_top do rectified_image[*,i] = shift(rectified_image[*,i], -cc_shift[i])

  ; return the slit edges
  return, [i_bottom, i_top]


END



; ****************************************************************************************



FUNCTION flame_getslits_update_slit, fuel, old_slit, yshift, slitid_top, slitid_bottom, slit_height, poly_coeff

    ; make the cutout structure
    cutout = { $
      filename: '', $
      rectification: ptr_new(), $
      speclines: ptr_new() }

    ; replicate for each of the frames
    cutouts = replicate(cutout, fuel.util.N_frames)

    ; initialize the pointers
    cutouts.rectification = ptrarr(fuel.util.N_frames, /allocate_heap)
    cutouts.speclines = ptrarr(fuel.util.N_frames, /allocate_heap)

    ; add new fields to slit structure
    new_slit = create_struct( $
      'yshift', yshift, $
      'slitid_top', slitid_top, $
      'slitid_bottom', slitid_bottom, $
      'height', slit_height, $
      'bottom_poly', poly_coeff, $
      'rough_wavecal', ptr_new(/allocate_heap), $
      'cutouts', cutouts, $
      'outlambda_min', 0d, $
      'outlambda_delta', 0d, $
      'outlambda_Npix', 0L )

    ; now we need to merge old and new slit structures

    ; check if new fields are already present in the slits structure
    if tag_exist(old_slit, 'yshift') then begin

      ; in that case, assign new values
      struct_assign, new_slit, old_slit, /nozero

      ; otherwise, append new fields
    endif else $
      new_slit = create_struct( old_slit, new_slit )

    return, new_slit

END

;******************************************************************


PRO flame_getslits_multislit, fuel=fuel

  ; frame to use to find the edges
  frame_filename = fuel.util.master_getslit

  ; read in the frame
  image=readfits(frame_filename, hdr)

  ; if the user specified the slit positions manually, then do not apply a shift
  if fuel.input.slit_position_file ne 'none' then $
    yshift = 0.0 else $
    ; get the overall vertical shift between the expected and the actual slit positions
    yshift = flame_getslits_findshift( image, fuel.slits.approx_top, fuel.slits.approx_bottom )

  ; create array of new slit structures
  slits = []

  ; trace the edges of the slits using the sky emission lines or sky continuum
  for i_slit=0, n_elements(fuel.slits)-1 do begin

    ; read the old slits structure - containing the info from the header
    old_slits_struc = fuel.slits[i_slit]

    print, 'Finding the edges of slit ' + strtrim(old_slits_struc.number,2) + ' - ' + old_slits_struc.name

    ; cross-correlation: find approximate edges of the slit and make rectified image
    approx_edges = flame_getslits_crosscorr( image, old_slits_struc.approx_bottom - yshift, old_slits_struc.approx_top - yshift, rectified_image=rectified_image )

    ; if the slit position were specified manually, then do not use the edges detected automatically
    if fuel.input.slit_position_file ne 'none' then approx_edges = [old_slits_struc.approx_bottom, old_slits_struc.approx_top]

    ; split slit into three chunks and use cross-correlation to find slit edges
    N_pixel_x = (size(image))[1]
    edges_left = flame_getslits_crosscorr( image[0:N_pixel_x/3-1, *], approx_edges[0], approx_edges[1])
    edges_center = flame_getslits_crosscorr( image[N_pixel_x/3 : N_pixel_x*2/3-1, *], approx_edges[0], approx_edges[1])
    edges_right = flame_getslits_crosscorr( image[N_pixel_x*2/3 : -1, *], approx_edges[0], approx_edges[1])

    if keyword_set(fuel.input.use_sky_edge) then begin
      ; identify top and bottom edge using sky background or flat lamp
      slitid_top = flame_getslits_trace_continuum(image, approx_edges[1], /top )
      slitid_bottom = flame_getslits_trace_continuum(image, approx_edges[0], /bottom )
    endif else begin
      ; identify top and bottom edge using OH lines (and in this case use the rectified image)
      slitid_top = flame_getslits_trace_skylines(rectified_image, approx_edges[1], /top )
      slitid_bottom = flame_getslits_trace_skylines(rectified_image, approx_edges[0], /bottom )
    endelse

    ; calculate the slit height
    slit_height = median( slitid_top - slitid_bottom )
    if ~finite(slit_height) then slit_height = median(slitid_top) - median(slitid_bottom)

    ; combine the top and bottom edge measurements together (i.e. subtract slit_height from the top edge)
    x_to_fit = [ where(finite(slitid_top), /null ), where(finite(slitid_bottom), /null) ]
    y_to_fit = [ slitid_top[ where(finite(slitid_top), /null) ] - slit_height, slitid_bottom[where(finite(slitid_bottom), /null)] ]

    ; fit a 3-rd order polynomial to the combined edges
    poly_coeff = robust_poly_fit( x_to_fit, y_to_fit, 3 )

    ; expand the slit structure with new fields and update them
    this_slit = flame_getslits_update_slit( fuel, old_slits_struc, yshift, slitid_top, slitid_bottom, slit_height, poly_coeff)

    ; add to array with the other slits
    slits = [slits, this_slit]

  endfor

  ; save the slit structures in fuel
  new_fuel = { input:fuel.input, util:fuel.util, instrument:fuel.instrument, diagnostics:fuel.diagnostics, slits:slits }
  fuel=new_fuel

END


;******************************************************************


PRO flame_getslits_longslit, fuel=fuel

  ; frame to use to find the edges
  frame_filename = (fuel.util.corrscience_filenames)[0]

  ; read in the frame
  im=readfits(frame_filename, hdr)

  ; read the old slits structure - containing the info from the header
  old_slits_struc = fuel.slits[0]

  ; calculate slit heigh
  slit_height = old_slits_struc.approx_top - old_slits_struc.approx_bottom

  ; expand the slit structure with new fields and update them
  this_slit = flame_getslits_update_slit( fuel, old_slits_struc, !values.d_nan, !values.d_nan, !values.d_nan, slit_height, old_slits_struc.approx_bottom)

  ; save the slit structures in fuel
  new_fuel = { input:fuel.input, util:fuel.util, instrument:fuel.instrument, diagnostics:fuel.diagnostics, slits:this_slit }
  fuel=new_fuel


END


;******************************************************************


PRO flame_getslits_writeds9, fuel=fuel, raw=raw
  ;
  ; write a ds9 region files that shows the slit edges
  ; if /raw is set then the individual "raw" measurements of the slit identifications
  ; are shown, as opposed to the polynomial fit

  ; name of the region file
  if keyword_set(raw) then region_filename = 'slits_raw.reg' else $
    region_filename = 'slits.reg'

  ; read the slits structures
  slits = fuel.slits

  ; number of horizontal pixel in one frame
  N_pix_x = (size( readfits((fuel.util.corrscience_filenames)[0]) ) )[1]

  ; open file
  openw, lun, fuel.input.intermediate_dir + region_filename, /get_lun

  ; write header
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=1 include=1 source=1'
  printf, lun, 'image'

  for i_slit=0, n_elements(slits)-1 do begin


    if keyword_set(raw) then begin    ; show the actual measured points along the slit edges -------------

      ; find the points among the many NaNs
      top_x = where( finite(slits[i_slit].slitid_top), /null )
      top_y = slits[i_slit].slitid_top[top_x]
      bottom_x = where( finite(slits[i_slit].slitid_bottom), /null )
      bottom_y = slits[i_slit].slitid_bottom[bottom_x]

      ; alternate colors for clarity
      color_string = (['green', 'red'])[i_slit mod 2]

      ; radius of each point, in pixels
      radius = '2'

      ; write the line corresponding to each point for the top edge
      if top_x ne !NULL then for i=0, n_elements(top_x)-1 do $
        printf, lun, 'circle(' + strtrim(top_x[i],2) + ',' + strtrim(top_y[i],2) + ',' + radius + $
        ') # color=' + color_string + ' text={SLIT ' + strtrim(slits[i_slit].number,2) + '}'

      ; write the line corresponding to each point for the bottom edge
      if bottom_x ne !NULL then for i=0, n_elements(bottom_x)-1 do $
        printf, lun, 'circle(' + strtrim(bottom_x[i],2) + ',' + strtrim(bottom_y[i],2) + ',' + radius + $
        ') # color=' + color_string + ' text={SLIT ' + strtrim(slits[i_slit].number,2) + '}'


    endif else begin    ; show the polynomial fit to the slit edges --------------------------------

      ; generate points from the polynomial fit
      top_x = [ 8*indgen(N_pix_x/8), N_pix_x-1] ; one point every 8 pixels plus the last pixel
      top_y = poly(top_x, slits[i_slit].bottom_poly) + slits[i_slit].height
      bottom_x = top_x
      bottom_y = poly(bottom_x, slits[i_slit].bottom_poly)

      ; concatenate top and bottom points
      all_x = [top_x, reverse(bottom_x)]
      all_y = [top_y, reverse(bottom_y)]

      ; make the string with all the points
      all_points = ''
      for i=0,n_elements(all_x)-2 do all_points += strtrim(all_x[i],2) + ',' + cgnumber_formatter(all_y[i], decimals=1) + ','
      ; add the last two points without the final comma
      all_points += strtrim(all_x[-1],2) + ',' + cgnumber_formatter(all_y[-1], decimals=1)

      ; alternate colors for clarity
      color_string = (['green', 'red'])[i_slit mod 2]

      ; write the line corresponding to this slit
      printf, lun, 'polygon(' + all_points + ') # color=' + color_string + ' text={SLIT ' + strtrim(slits[i_slit].number,2) + ' - ' + slits[i_slit].name + '}'

    endelse

  endfor

  ; close file
  free_lun, lun

END


;******************************************************************


PRO flame_getslits_write_slitim, fuel=fuel
  ; make an image where the pixels belonging to a slit have the slit number as a value, otherwise zero

  ; read slits structures
  slits = fuel.slits

  ; read in the first science frame to get the right dimensions
  slitim = fix(0 * readfits((fuel.util.corrscience_filenames)[0], hdr))

  ; construct the coordinates for the pixels in the image
  N_pix_x = (size(slitim))[1]
  N_pix_y = (size(slitim))[2]
  x_axis = indgen(N_pix_x)
  y_axis = indgen(N_pix_y)
  pixel_x = x_axis # replicate(1, N_pix_y)
  pixel_y = transpose(y_axis # replicate(1, N_pix_x))

  for i_slit=0,n_elements(slits)-1 do begin

    top_y = (poly(x_axis, slits[i_slit].bottom_poly) + slits[i_slit].height) # replicate(1, N_pix_x)
    bottom_y = poly(x_axis, slits[i_slit].bottom_poly) # replicate(1, N_pix_x)

    w_slit = where( pixel_y LT top_y AND pixel_y GT bottom_y, /null)
    slitim[w_slit] = slits[i_slit].number

  endfor

  writefits, fuel.input.intermediate_dir + fuel.util.slitim_filename, slitim


END


;******************************************************************


PRO flame_getslits, fuel

	start_time = systime(/seconds)

  print, ''
  print, '-------------------------------------'
  print, '---        flame_getslits         ---'
  print, '-------------------------------------'
  print, ''


  ; identify all slits from the data, and write the fuel.slits structures
  if fuel.input.longslit then $
  flame_getslits_longslit, fuel=fuel $
    else $
  flame_getslits_multislit, fuel=fuel

  ; write ds9 region file with the slit traces
  flame_getslits_writeds9, fuel=fuel
  flame_getslits_writeds9, fuel=fuel, /raw

  ; write slitim (FITS file with slit image - just for diagnostics purpose)
  flame_getslits_write_slitim, fuel=fuel

  ; if we are reducing only one slit, then delete all the others
  if fuel.input.reduce_only_oneslit ne 0 and n_elements(fuel.slits) GT 1 then begin
    new_fuel = { input:fuel.input, util:fuel.util, instrument:fuel.instrument, $
      diagnostics:fuel.diagnostics, slits:fuel.slits[fuel.input.reduce_only_oneslit-1] }
    fuel=new_fuel
  endif


	print, ''
  print, 'flame_getslits took ', $
    cgnumber_formatter( systime(/seconds) - start_time, decimals=2), ' seconds'
  print, ''

END
