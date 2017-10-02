
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


FUNCTION flame_getslits_trace_continuum, image_in, approx_edge, top=top, bottom=bottom

  ; check keywords
  if ~keyword_set(top) AND ~keyword_set(bottom) then message, 'Please select either /top or /bottom'
  if keyword_set(top) AND keyword_set(bottom) then message, 'Please select either /top or /bottom'

  image = image_in

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
    derivative = shift(profile,1) - profile
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
    derivative = shift(profile,1) - profile
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

PRO flame_getslits_trace_emlines, fuel, image, approx_edges, slitid_top=slitid_top, slitid_bottom=slitid_bottom
  ;
  ; Given a frame containing bright emission lines from sky or arcs, it traces the edge of a slit.
  ; approx_edges is an array of two elements indicating the approximate y-coordinates
  ; of the top and bottom edges of a slit
  ;


  ; extract cutout and reference spectrum
  ; ---------------------------------------------------

  ; set the number of pixels above and below the slit to include in the cutout
  ymargin = fuel.settings.trace_slit_ymargin

  ; extract a cutout of the slit
  cutout = image[ * , min(approx_edges) - ymargin : max(approx_edges) + ymargin ]
  Nx = (size(cutout))[1]
  Ny = (size(cutout))[2]

	; calculate the distance of each pixel row from the spatial center of the slit
	dist_to_center = abs( indgen(Ny) - Ny/2 )

	; extract profile
	profile = median(cutout, dimension=1)

	; detect positive flux fluctuations, probably due to bright objects
	w_objects = where( profile - median(profile) GT 3.0* stddev(profile, /nan), /null)

	; exclude these pixels by assigning them an arbitrary high distance
	if w_objects NE !NULL then dist_to_center[w_objects] = 999

	; identify the most-central five pixel rows that are not contaminated by bright objects
	w_skyregion = (sort(dist_to_center))[0:4]

	; extract the sky region
	sky_region = []
	for i=0, n_elements(w_skyregion)-1 do sky_region = [ [sky_region], [cutout[*,w_skyregion[i]]] ]

	; extract the spectrum from the sky region
	central_spectrum = median(sky_region, dimension=2)

  ; identify the continuum (particularly important in the K band)
  central_spectrum_padded = [replicate(0d,200), central_spectrum, replicate(0d,200)]
  central_spectrum_continuum = (median(central_spectrum_padded, 200))[200:-201]

  ; make the reference spectrum: subtract continuum, smooth, trim edges
  ref_spectrum = gauss_smooth(central_spectrum - central_spectrum_continuum, 3, /nan)
  ref_spectrum[0:40]=0.0
  ref_spectrum[-40:-1]=0.0

  ; subtract the continuum from the full slit
  cutout_zeroed = cutout - central_spectrum_continuum # transpose(replicate(1, (size(cutout))[2]))


  ; for each emission line find the top and bottom edges
  ; ---------------------------------------------------

	; define a flux threshold for emission lines: standard deviation of the spectrum without the extreme points
	sorted_values = ref_spectrum[sort(ref_spectrum)]
	threshold = 0.5*stddev(sorted_values[0.1*Nx : -0.3*Nx], /nan)

	; automatically identify peaks
	w_peaks = where( ref_spectrum GT shift(ref_spectrum, 1) and ref_spectrum GT shift(ref_spectrum, 2) and $
		ref_spectrum GT shift(ref_spectrum, 3) and ref_spectrum GT shift(ref_spectrum, -1) and $
		ref_spectrum GT shift(ref_spectrum, -2) and ref_spectrum GT shift(ref_spectrum, -3) and $
		ref_spectrum - 0.5*(shift(ref_spectrum, 3)+shift(ref_spectrum, -3)) GT threshold, /null )

  if n_elements(w_peaks) LT 5 then message, 'Did not find enough emission lines to trace the slit'
  print, n_elements(w_peaks), ' emission lines identified'

  ; make arrays with the final results: x and y coordinates of edge identifications
  ; for both the top and the bottom of the slit
  slit_top_x = []
  slit_top_y = []
  slit_top_height = []
  slit_bottom_x = []
  slit_bottom_y = []
  slit_bottom_height = []

  ; make a binary mask of the cutout to store the pixels that belong to the identified lines
  mask = bytarr(Nx, Ny)

  ; make an array with the height of each line
  line_height = []

  for i_line=0, n_elements(w_peaks)-1 do begin

    ; this is the x coordinate of this emission line
    x_line = w_peaks[i_line]

    ; check that this line has not already been identified with the previous line
    if mask[x_line, Ny/2] eq 1 then continue

    ; minimum flux to be considered as part of the emission line
    min_flux = threshold + 0.5*(ref_spectrum[x_line]-threshold)

    ; check that the peak of the initial pixel is above the min flux (otherwise search2d does not work)
    if cutout_zeroed[x_line, Ny/2] LT min_flux then continue

    ; make a copy of the cutout, using only the pixels around the line
    cutout_thisline = cutout_zeroed
    cutout_thisline[0:x_line-15,*]=min(cutout_zeroed)
    cutout_thisline[x_line+15:-1,*]=min(cutout_zeroed)

    ; find the 2D region corresponding to the emission line
    w_line = search2d( cutout_thisline, x_line, Ny/2, min_flux, max(cutout_zeroed) )

    ; store the region in the mask
    mask[w_line] = 1

    ; transform the 2D indices into 1D indices
    index_y = w_line / Nx
    index_x = w_line - (index_y * Nx)

    ; sort pixels by x index using the histogram command
    h = Histogram(index_x, min=min(index_x), binsize=1, reverse_indices=ri)

    ; check that there are at least five pixel columns in this line
    if n_elements(ri) LE n_elements(w_line)+5 then continue

    ; get the running top and bottom edges of this line
    edge_x = []
    edge_top = []
    edge_bottom = []

    ; for each pixel column identify the top and bottom edges (skip first and last column)
    for x=0,n_elements(ri)-n_elements(w_line)-3 do begin
      w_x = ri[ ri[x]:ri[x+1]-1 ] ; these are the pixels in the regions in this column
      edge_x = [ edge_x, index_x[w_x[0]] ] ; store the x coordinate of this column
      edge_top = [ edge_top, max(index_y[w_x]) ]
      edge_bottom = [ edge_bottom, min(index_y[w_x]) ]
    endfor

    ; only take the outermost pixels and save their x and y coordinates
    top_ref = (edge_top[sort(edge_top)])[-3]      ; use the third pixel from the top as a reference
    w_top = where(edge_top GE top_ref, /null)     ; select all pixels at or above the reference pixel
    slit_top_x = [ slit_top_x, edge_x[w_top] ]    ; store the x coordinate of the selected pixels
    slit_top_y = [ slit_top_y, edge_top[w_top] ]  ; store the y coordinate of the selected pixels

    ; same thing for the bottom edge
    bottom_ref = (edge_bottom[sort(edge_bottom)])[2]
    w_bottom = where(edge_bottom LE bottom_ref, /null)
    slit_bottom_x = [ slit_bottom_x, edge_x[w_bottom] ]
    slit_bottom_y = [ slit_bottom_y, edge_bottom[w_bottom] ]

    ; calculate the slit height along the y axis (i.e. does not account for tilt)
    slit_height = top_ref - bottom_ref

    ; store slit height for this line
    line_height = [line_height, slit_height]

    ; also store slit height in each top and bottom measurement
    slit_top_height = [slit_top_height, replicate(slit_height, n_elements(w_top))]
    slit_bottom_height = [slit_bottom_height, replicate(slit_height, n_elements(w_bottom))]

  endfor

  ; sigma-clip the values of slit heights to find the true heights
  meanclip, line_height, mean_height, sigma_height, clipsig=2.0
  if sigma_height LT 0.5 then sigma_height=0.5

  ; remove detections at the edge of the cutout or with slit heights too far from the mean
  w_ok = where( slit_top_y LT (size(cutout))[2]-1 and abs(slit_top_height-mean_height) LT 3.0*sigma_height, /null )
  slit_top_x = slit_top_x[w_ok]
  slit_top_y = slit_top_y[w_ok]

  ; again but for the bottom edge
  w_ok = where( slit_bottom_y GT 0 and abs(slit_bottom_height-mean_height) LT 3.0*sigma_height, /null )
  slit_bottom_x = slit_bottom_x[w_ok]
  slit_bottom_y = slit_bottom_y[w_ok]

  ; enlarge by one pixel the edges, otherwise they are very strict
  slit_top_y += 1
  slit_bottom_y -= 1


  ; output the identifications for the top and the bottom
  ; ---------------------------------------------------

  slitid_top = dblarr(Nx) + !values.d_nan
  slitid_top[slit_top_x] =  min(approx_edges) - ymargin + slit_top_y  ; add back the y coordinate of the edge of the cutout

  slitid_bottom = dblarr(Nx) + !values.d_nan
  slitid_bottom[slit_bottom_x] =  min(approx_edges) - ymargin + slit_bottom_y


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
      filename: '', $ ; this one has the corr frame
      illcorr_applied: 0, $ ; flag: 1 if illumination correction is applied to cutout.filename
      lambda_coeff: ptr_new(), $
      gamma_coeff: ptr_new(), $
      speclines: ptr_new() }

    ; replicate for each of the frames
    cutouts = replicate(cutout, fuel.util.science.n_frames)

    ; initialize the pointers
    cutouts.lambda_coeff = ptrarr(fuel.util.science.n_frames, /allocate_heap)
    cutouts.gamma_coeff = ptrarr(fuel.util.science.n_frames, /allocate_heap)
    cutouts.speclines = ptrarr(fuel.util.science.n_frames, /allocate_heap)

    ; again for the arc cutout
    arc_cutout = cutout
    arc_cutout.lambda_coeff = ptr_new(/allocate_heap)
    arc_cutout.gamma_coeff = ptr_new(/allocate_heap)
    arc_cutout.speclines = ptr_new(/allocate_heap)

    ; add new fields to slit structure
    new_slit = create_struct( $
      'yshift', yshift, $
      'slitid_top', slitid_top, $
      'slitid_bottom', slitid_bottom, $
      'height', slit_height, $
      'bottom_poly', poly_coeff, $
      'yrange_cutout', [0L,0L], $
      'rough_arclambda', ptr_new(/allocate_heap), $
      'rough_arcflux', ptr_new(/allocate_heap), $
      'rough_skylambda', ptr_new(/allocate_heap), $
      'rough_skyflux', ptr_new(/allocate_heap), $
      'cutouts', cutouts, $
      'arc_cutout', arc_cutout, $
      'outlambda_min', 0d, $
      'outlambda_delta', 0d, $
      'outlambda_Npix', 0L, $
      'output_file', '' )

    ; now we need to merge old and new slit structures

    ; check if new fields are already present in the slits structure
    if tag_exist(old_slit, 'yshift') then begin

      ; in that case, assign new values
      struct_assign, new_slit, old_slit, /nozero
      new_slit = old_slit

      ; otherwise, append new fields
    endif else $
      new_slit = create_struct( old_slit, new_slit )

    return, new_slit

END

;******************************************************************


PRO flame_getslits_multislit, fuel=fuel

  ; frame to use to find the edges
  frame_filename = fuel.util.slitflat.master_file

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

    if fuel.settings.trace_slit_with_emlines eq 0 then begin
      ; identify top and bottom edge using sky background or flat lamp
      slitid_top = flame_getslits_trace_continuum(image, approx_edges[1], /top )
      slitid_bottom = flame_getslits_trace_continuum(image, approx_edges[0], /bottom )
    endif else begin
      ; identify top and bottom edge using emission lines
      flame_getslits_trace_emlines, fuel, image, approx_edges, slitid_top=slitid_top, slitid_bottom=slitid_bottom
    endelse

    ; calculate the slit height
    slit_height = median( slitid_top - slitid_bottom )
    if ~finite(slit_height) then slit_height = median(slitid_top) - median(slitid_bottom)

    ; combine the top and bottom edge measurements together (i.e. subtract slit_height from the top edge)
    x_to_fit = [ where(finite(slitid_top), /null ), where(finite(slitid_bottom), /null) ]
    y_to_fit = [ slitid_top[ where(finite(slitid_top), /null) ] - slit_height, slitid_bottom[where(finite(slitid_bottom), /null)] ]

    ; fit a polynomial to the combined edges
    ; NB: apply offset from the slit flat (should be zero if not using slit flat)
    poly_coeff = robust_poly_fit( x_to_fit, y_to_fit + fuel.input.slitflat_offset, fuel.settings.trace_slit_polydegree )

    ; expand the slit structure with new fields and update them
    this_slit = flame_getslits_update_slit( fuel, old_slits_struc, yshift, slitid_top, slitid_bottom, slit_height, poly_coeff)

    ; add to array with the other slits
    slits = [slits, this_slit]

  endfor

  ; save the slit structures in fuel
  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:fuel.instrument, diagnostics:fuel.diagnostics, slits:slits }
  fuel=new_fuel

END


;******************************************************************


PRO flame_getslits_longslit, fuel=fuel

  ; frame to use to find the edges
  frame_filename = (fuel.util.science.corr_files)[0]

  ; read in the frame
  im=readfits(frame_filename, hdr)

  ; read the old slits structure - containing the info from the header
  old_slits_struc = fuel.slits[0]

  ; calculate slit heigh
  slit_height = old_slits_struc.approx_top - old_slits_struc.approx_bottom

  ; expand the slit structure with new fields and update them
  this_slit = flame_getslits_update_slit( fuel, old_slits_struc, !values.d_nan, !values.d_nan, !values.d_nan, slit_height, old_slits_struc.approx_bottom)

  ; save the slit structures in fuel
  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:fuel.instrument, diagnostics:fuel.diagnostics, slits:this_slit }
  fuel=new_fuel


END


;******************************************************************


PRO flame_getslits_writeds9, fuel=fuel, raw=raw
  ;
  ; write a ds9 region files that shows the slit edges
  ; if /raw is set then the individual "raw" measurements of the slit identifications
  ; are shown, as opposed to the polynomial fit
  ; note: adding +1 to the coordinate to account for different notation of IDL vs ds9

  ; name of the region file
  if keyword_set(raw) then region_filename = 'slits_raw.reg' else $
    region_filename = 'slits.reg'

  ; read the slits structures
  slits = fuel.slits

  ; number of horizontal pixel in one frame
  N_pix_x = (size( readfits((fuel.util.science.corr_files)[0]) ) )[1]

  ; y coordinate of the center of the top edge
  top_edges = dblarr(n_elements(slits))
  for i_slit=0, n_elements(slits)-1 do $
    top_edges[i_slit] = poly(N_pix_x/2, slits[i_slit].bottom_poly)

  ; get the rank order of each slit
  rank = intarr(n_elements(top_edges))
  rank[sort(top_edges)] = indgen(n_elements(top_edges))

  ; assign alternating colors to slits (alternate by position, not slit number)
  slit_color = (['green', 'red'])[rank mod 2]

  ; open file
  openw, lun, fuel.util.intermediate_dir + region_filename, /get_lun

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

      ; write the line corresponding to each point for the top edge
      if top_x ne !NULL then for i=0, n_elements(top_x)-1 do $
        printf, lun, 'point(' + strtrim(top_x[i]+1,2) + ',' + strtrim(top_y[i]+1,2) + $
        ') # point=cross color=' + slit_color[i_slit]

      ; write the line corresponding to each point for the bottom edge
      if bottom_x ne !NULL then for i=0, n_elements(bottom_x)-1 do $
        printf, lun, 'point(' + strtrim(bottom_x[i]+1,2) + ',' + strtrim(bottom_y[i]+1,2) + $
        ') # point=cross color=' + slit_color[i_slit]


    endif else begin    ; show the polynomial fit to the slit edges --------------------------------

      ; generate points from the polynomial fit
      top_x = [ 8*indgen(N_pix_x/8), N_pix_x-1] ; one point every 8 pixels plus the last pixel
      top_y = poly(top_x, slits[i_slit].bottom_poly) + slits[i_slit].height
      bottom_x = top_x
      bottom_y = poly(bottom_x, slits[i_slit].bottom_poly)

      ; concatenate top and bottom points
      all_x = [top_x, reverse(bottom_x)] +1
      all_y = [top_y, reverse(bottom_y)] +1

      ; make the string with all the points
      all_points = ''
      for i=0,n_elements(all_x)-2 do all_points += strtrim(all_x[i],2) + ',' + cgnumber_formatter(all_y[i], decimals=1) + ','
      ; add the last two points without the final comma
      all_points += strtrim(all_x[-1],2) + ',' + cgnumber_formatter(all_y[-1], decimals=1)

      ; write the line corresponding to this slit
      printf, lun, 'polygon(' + all_points + ') # color=' + slit_color[i_slit] + ' text={SLIT ' + strtrim(slits[i_slit].number,2) + ' - ' + slits[i_slit].name + '}'

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
  slitim = fix(0 * readfits((fuel.util.science.corr_files)[0], hdr))

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

  writefits, fuel.util.intermediate_dir + 'slitim.fits', slitim


END


;******************************************************************


PRO flame_getslits, fuel

	flame_util_module_start, fuel, 'flame_getslit'


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

  ; if we are reducing only one slit, then set the skip flag in all the others
  if fuel.input.reduce_only_oneslit ne 0 then begin
    fuel.slits.skip = 1
    fuel.slits[fuel.input.reduce_only_oneslit-1].skip = 0
  endif


  flame_util_module_end, fuel

END
