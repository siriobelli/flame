
; ****************************************************************************************


FUNCTION flame_getslits_findshift, fuel=fuel
  ;
  ; find the shift between the expected slit positions and the real ones
  ; never use the CR-cleaned files
  ;

  ; for longslits this shift makes no sense
  if fuel.input.longslit then return, 0.0

  ; if the user specified the slit positions manually, then do not apply a shift
  if fuel.input.slit_position_file ne 'none' then return, 0.0

  ; read one FITS file
  spec2d = readfits((fuel.util.corrscience_filenames)[0])

  ; integrate the 2D spectrum along the wavelength direction, and smooth
  x = median(total(spec2d, 1, /nan), 15)

  ; find the edges of the slits (1d profile where positive peaks mark the beginning and negative peaks mark the end of the slit)
  edges = x - shift(x, 4)
  edges[ where( abs(edges) LT stddev(edges, /nan), /null ) ] = 0
  edges[ where(edges GT 0.0, /null) ] =  1.0
  edges[ where(edges LT 0.0, /null) ] = -1.0

  ; these are the expected edges
  expected_edges = dblarr( n_elements(edges) )
  expected_edges[fuel.slits.approx_top] = -1.
  expected_edges[fuel.slits.approx_bottom] = 1.

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


FUNCTION flame_getslits_trace_skyedge, image_in, approx_edge, top=top, bottom=bottom

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



FUNCTION flame_getslits_trace_edge, image, approx_edge, top=top, bottom=bottom

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

  ; rectify the OH lines
  ;----------------------

  ; use this as reference (assume is within the slit)
  bottom_row = cutout[*,0]

  ; these are the shifts considered (in pixels)
  lag = -30 + indgen(61)

  ; this array will contain the shift value for each row
  shift = intarr(cutout_size)

  ; cross-correlate all the other rows and find the shift
  for i_row=1, cutout_size-1 do begin

    ; cross-correlate this row with the bottom row
    cc = c_correlate( bottom_row, cutout[*,i_row], lag)

    ; find the peak of the cross-correlation
    !NULL = max(cc, indmax)

    ; store this shift and move on to the next row
    shift[i_row] = lag[indmax]

  endfor

  ; find the differential shift at each row
  differential_shift = shift - shift(shift, 1)
  differential_shift = differential_shift[1:*]

  ; assume the first three rows are definitely within in the slit
  ; and find the row at which the shift is discontinuous
  for i_row=3, n_elements(differential_shift)-1 do begin

    previous_shifts = differential_shift[0:i_row-1]
    tolerance = stddev(previous_shifts) > 1
    if abs( differential_shift[i_row] - median(previous_shifts) ) GT 3.0*tolerance then break

  endfor

  ; was a discontinuity found?
  if i_row LT n_elements(differential_shift) then begin

    ; take the row closest to the discontinuity as the reference for the shift
    shift = shift[i_row] - shift

    ; do not apply shifts above that point
    shift[i_row+1:*] = 0

    ; determine the 'fiducial' number of pixels in the slit
    Npix_slit_fiducial = i_row-1

  endif else Npix_slit_fiducial = 10

  ; rectify the cutout
  for i_row=0, cutout_size-1 do cutout[*,i_row] = shift(cutout[*,i_row], shift[i_row])

  ;----------------------

  ; extract a spectrum from the fiducial slit regions
  sky_spectrum = median(cutout[*,0:Npix_slit_fiducial-1], dimension=2)

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

    print, approx_edge, keyword_set(top), keyword_set(bottom)

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


PRO flame_getslits_trace, image=image, slits=slits, yshift=yshift, poly_coeff=poly_coeff, slit_height=slit_height, use_sky_edge=use_sky_edge

;
; traces the top and bottom edges of a slit in the image and return the slit height and the
; coefficients of a polynomial fit describing the *bottom* edge
;

  ; apply y shift
  expected_top = slits.approx_top - yshift
  expected_bottom = slits.approx_bottom - yshift

  if keyword_set(use_sky_edge) then begin
    ; identify top and bottom edge using sky background
    top_edge = flame_getslits_trace_skyedge(image, expected_top, /top )
    bottom_edge = flame_getslits_trace_skyedge(image, expected_bottom, /bottom )
  endif else begin
    ; identify top and bottom edge using OH lines
    top_edge = flame_getslits_trace_edge(image, expected_top, /top )
    bottom_edge = flame_getslits_trace_edge(image, expected_bottom, /bottom )
  endelse

  ; calculate the slit height
  slit_height = median( top_edge - bottom_edge )
  if ~finite(slit_height) then slit_height = median(top_edge) - median(bottom_edge)

  ; combine the top and bottom edge measurements together (i.e. subtract slit_heigh from the top edge)
  x_to_fit = [ where(finite(top_edge), /null ), where(finite(bottom_edge), /null) ]
  y_to_fit = [ top_edge[ where(finite(top_edge), /null) ] - slit_height, bottom_edge[where(finite(bottom_edge), /null)] ]

  ; fit a 3-rd order polynomial to the combined edges
  poly_coeff = robust_poly_fit( x_to_fit, y_to_fit, 3 )


END



;******************************************************************

PRO flame_getslits_findedges, fuel=fuel, yshift=yshift

  ; frame to use to find the edges
  frame_filename = (fuel.util.corrscience_filenames)[0]

  ; read in the frame
  im=readfits(frame_filename, hdr)

  ; create array of new slit structures
  slits = []

  ; LONGSLIT ---------------------------------------------------------------------
  if fuel.input.longslit then begin

    ; read the old slits structure - containing the info from the header
    old_slits_struc = fuel.slits[0]

    ; make new fields for slit structure
    new_slits_struc = create_struct( $
        'yshift', !values.d_nan, $
        'height', old_slits_struc.approx_top - old_slits_struc.approx_bottom, $
        'bottom_poly', old_slits_struc.approx_bottom, $
        'filenames', ptr_new(/allocate_heap), $
        'rough_wavecal', ptr_new(/allocate_heap), $
        'rectification', ptr_new(/allocate_heap), $
        'outlambda_min', 0d, $
        'outlambda_delta', 0d, $
        'outlambda_Npix', 0L )

    ; merge old and new slits structures

    ; check if new fields are already present in the slits structure
    if tag_exist(old_slits_struc, 'yshift') then begin

      ; in that case, assign new values
      struct_assign, new_slits_struc, old_slits_struc, /nozero
      this_slit = old_slits_struc

    ; otherwise, append new fields
    endif else $
      this_slit = create_struct( old_slits_struc, new_slits_struc )

    ; add to array with the other slits
    slits = [slits, this_slit]

 ; MOS      ---------------------------------------------------------------------
  endif else begin

    ; trace the edges of the slits using the sky emission lines or sky continuum
    for i_slit=0, n_elements(fuel.slits)-1 do begin

      ; read the old slits structure - containing the info from the header
      old_slits_struc = fuel.slits[i_slit]

      ; trace slit
      flame_getslits_trace, image=im, slits=old_slits_struc, yshift=yshift, $
            poly_coeff=poly_coeff, slit_height=slit_height, use_sky_edge=fuel.input.use_sky_edge

        ; add new fields to slit structure
        new_slits_struc = create_struct( $
          'yshift', yshift, $
          'height', slit_height, $
          'bottom_poly', poly_coeff, $
          'filenames', ptr_new(/allocate_heap), $
          'rough_wavecal', ptr_new(/allocate_heap), $
          'rectification', ptr_new(/allocate_heap), $
          'outlambda_min', 0d, $
          'outlambda_delta', 0d, $
          'outlambda_Npix', 0L )

      ; merge old and new slits structures

      ; check if new fields are already present in the slits structure
      if tag_exist(old_slits_struc, 'yshift') then begin

        ; in that case, assign new values
        struct_assign, new_slits_struc, old_slits_struc, /nozero
        this_slit = old_slits_struc

      ; otherwise, append new fields
      endif else $
        this_slit = create_struct( old_slits_struc, new_slits_struc )

      ; add to array with the other slits
      slits = [slits, this_slit]

    endfor

  endelse

  ; save the slit structures in fuel
  new_fuel = { input:fuel.input, util:fuel.util, instrument:fuel.instrument, diagnostics:fuel.diagnostics, slits:slits }
  fuel=new_fuel

END


;******************************************************************


PRO flame_getslits_writeds9, fuel=fuel
  ;
  ; write a ds9 region files that shows the slit edges
  ;

  ; name of the region file
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

    ; generate points from the polynomial fit
    top_x = 8*indgen(N_pix_x/8) ; one point every 8 pixels
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


PRO flame_getslits, fuel=fuel

  ; get the overall vertical shift between the expected and the actual slit positions
  yshift = flame_getslits_findshift( fuel=fuel )

  ; identify all slits from the data, and write the fuel.slits structures
  flame_getslits_findedges, fuel=fuel, yshift=yshift

  ; write ds9 region file with the slit traces
  flame_getslits_writeds9, fuel=fuel

  ; write slitim (FITS file with slit image - just for diagnostics purpose)
  flame_getslits_write_slitim, fuel=fuel

  ; if we are reducing only one slit, then delete all the others
  if fuel.input.reduce_only_oneslit ne 0 and n_elements(fuel.slits) GT 1 then begin
    new_fuel = { input:fuel.input, util:fuel.util, instrument:fuel.instrument, $
      diagnostics:fuel.diagnostics, slits:fuel.slits[fuel.input.reduce_only_oneslit-1] }
    fuel=new_fuel
  endif


END
