
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
  cutout_size = 60      ; even number please

  ; of these pixels, we assume that the first ... are certainly inside the slit
  Npix_slit_fiducial = 10

  ; let's extract a cutout centered on the edge
  cutout_bottom_ycoord = approx_edge-cutout_size/2
  cutout = image[ * , cutout_bottom_ycoord : cutout_bottom_ycoord + cutout_size - 1 ]

  ; from now on, the code assumes that we are interested in finding the top edge of a slit
  ; if instead we want the bottom edge, simply flip vertically the cutout 
  if keyword_set(bottom) then cutout = reverse(cutout, 2)

  ; get rid of negative values and NaNs
  cutout[ where(cutout LT 0.0 or ~finite(cutout), /null) ] = 0.0

  ; extract a spectrum from the fiducial slit regions
  sky_spectrum = median(cutout[*,0:Npix_slit_fiducial-1], dimension=2)

  ; measure the sky flux in between the OH lines
  mmm, sky_spectrum, level_betweenlines, sigma_betweenlines

  ; consider as a sky line all those pixels that are more than 4 sigma above the normal level of between-lines sky
  w_OH = where(sky_spectrum GT level_betweenlines + 4.0*sigma_betweenlines, /null)
  if w_OH eq !NULL then w_OH = -1 ; it won't be used anyway

  ; for each pixel column, take the ratio of ech pixel to the fiducial sky value in the slit 
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

  if Nlines LT 15 then begin

    message, 'I found only ' + strtrim(Nlines,2) + $
      ' OH lines! Edge tracing will be based on sky background.', /informational

  endif else begin

    print, strtrim(Nlines,2) + ' OH lines found. Tracing edge...'

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


PRO flame_getslits_trace, image=image, approx_top=approx_top, approx_bottom=approx_bottom, $
  poly_coeff=poly_coeff, slit_height=slit_height
;
; traces the top and bottom edges of a slit in the image and return the slit height and the 
; coefficients of a polynomial fit describing the *bottom* edge
;

  ; identify top edge
  top_edge = flame_getslits_trace_edge(image, approx_top, /top )

  ; identify bottom edge
  bottom_edge = flame_getslits_trace_edge(image, approx_bottom, /bottom )

  ; calculate the slit height
  slit_height = median( top_edge - bottom_edge )

  ; combine the top and bottom edge measurements together (i.e. subtract slit_heigh from the top edge)
  x_to_fit = [ where(finite(top_edge), /null ), where(finite(bottom_edge), /null) ]
  y_to_fit = [ top_edge[ where(finite(top_edge), /null) ] - slit_height, bottom_edge[where(finite(bottom_edge), /null)] ]

  ; fit a 3-rd order polynomial to the combined edges
  poly_coeff = robust_poly_fit( x_to_fit, y_to_fit, 3 )


END


; ****************************************************************************************


FUNCTION flame_getslits_findshift, frame, top, bottom
  ;
  ; find the shift between the expected slit positions and the real ones
  ;

  ; read one FITS file
  spec2d = readfits(frame)
  
  ; integrate the 2D spectrum along the wavelength direction, and smooth
  x = median(total(spec2d, 1, /nan), 15)
  
  ; find the edges of the slits (1d profile where positive peaks mark the beginning and negative peaks mark the end of the slit)
  edges = x - shift(x, 4)
  edges[ where( abs(edges) LT stddev(edges, /nan), /null ) ] = 0
  edges[ where(edges GT 0.0, /null) ] =  1.0
  edges[ where(edges LT 0.0, /null) ] = -1.0
  
  ; these are the expected edges
  expected_edges = dblarr( n_elements(edges) )
  expected_edges[top] = -1.
  expected_edges[bottom] = 1.
  
  ; cross-correlate to find the shift between expected and measured slit edges
  lag = indgen(400)-200
  crosscorr = c_correlate(edges, expected_edges, lag)
  max_crosscorr = max( crosscorr, max_ind, /nan)
  delta = lag[max_ind]
  
  ; return the shift
  return, delta
  
END


;******************************************************************


PRO flame_getslits_writeds9, fuel=fuel
  ;
  ; write a ds9 region files that shows the slit edges
  ;

  ; name of the region file
  region_filename = 'slits.reg'

  ; read the slits structures
  slits = *fuel.slits

  ; number of horizontal pixel in one frame
  N_pix_x = (size( readfits((*fuel.corrscience_files)[0]) ) )[1]

  ; open file
  openw, lun, fuel.intermediate_dir + region_filename, /get_lun

  ; write header
  printf, lun, '# Region file format: DS9 version 4.1'
  printf, lun, 'global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=1 include=1 source=1'
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

    ; write the line corresponding to this slit
    printf, lun, 'polygon(' + all_points + ') # text={SLIT ' + strtrim(slits[i_slit].number,2) + ' - ' + slits[i_slit].name + '}'

  endfor

  ; close file
  free_lun, lun

END


;******************************************************************
      
PRO flame_getslits_findedges, fuel=fuel

  ; read in the frame
  im=readfits((*fuel.corrscience_files)[0], hdr)
  
  ; create array of new slit structures
  slits = []

  ; LONGSLIT ---------------------------------------------------------------------
  if fuel.longslit then begin

    ; read the info about this slit from the header
    info_fromheader = (*fuel.slits_fromheader)[0]

    slits = { $
        number: info_fromheader.number, $
        name: info_fromheader.name, $
        PA: info_fromheader.PA, $
        approx_wavelength_lo: info_fromheader.approx_wavelength_lo, $
        approx_wavelength_hi: info_fromheader.approx_wavelength_hi, $
        yshift: !values.d_nan, $
        height: !values.d_nan, $
        bottom_poly: !values.d_nan, $
        filenames: ptr_new(/allocate_heap), $
        rectification: ptr_new(/allocate_heap) }

 ; MOS      ---------------------------------------------------------------------
  endif else begin

    ; compare the expected position with the measured ones and obtain rough shift
    yshift = flame_getslits_findshift( (*fuel.corrscience_files)[0], $
      (*fuel.slits_fromheader).approx_top, (*fuel.slits_fromheader).approx_bottom )
 
    ; trace the edges of the slits using the sky emission lines
    for i_slit=0, n_elements((*fuel.slits_fromheader))-1 do begin

    ; read the info about this slit from the header
    info_fromheader = (*fuel.slits_fromheader)[i_slit]

    ; trace slit
    flame_getslits_trace, image=im, approx_top=info_fromheader.approx_top-yshift, $
      approx_bottom=info_fromheader.approx_bottom-yshift, $
          poly_coeff=poly_coeff, slit_height=slit_height

      ; add new fields to slit structure
      this_slit = { $
        number: info_fromheader.number, $
        name: info_fromheader.name, $
        PA: info_fromheader.PA, $
        approx_wavelength_lo: info_fromheader.approx_wavelength_lo, $
        approx_wavelength_hi: info_fromheader.approx_wavelength_hi, $
        yshift: yshift, $
        height: slit_height, $
        bottom_poly: poly_coeff, $
        filenames: ptr_new(/allocate_heap), $
        rectification: ptr_new(/allocate_heap) }

      slits = [slits, this_slit]

      endfor

  endelse

  ; save the slit structures in fuel
  *fuel.slits = slits

END


;******************************************************************


PRO flame_getslits_write_slitim, fuel=fuel
  ; make an image where the pixels belonging to a slit have the slit number as a value, otherwise zero

  ; read slits structures
  slits = *fuel.slits

  ; read in the first science frame to get the right dimensions
  slitim = fix(0 * readfits((*fuel.corrscience_files)[0], hdr))

  ; construct the coordinates for the pixels in the image
  N_pix_x = (size(slitim))[2]
  N_pix_y = (size(slitim))[1]
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
  
  writefits, fuel.intermediate_dir + fuel.slitim_filename, slitim
  

END 


;******************************************************************

PRO flame_getslits_cutout_extract, slitim, slit_structure, science_filenames, output_filenames

  ; loop through science frames
  for i_frame=0, n_elements(science_filenames)-1 do begin

    ; read in science frame 
    im = readfits(science_filenames[i_frame], header)

    ; select pixels belonging to this slit
    w_slit = where(slitim eq slit_structure.number, /null)
    if w_slit eq !NULL then message, slit_structure.name + ': slit not found in slitim!'

    ; create a mask that selects only pixels belonging to the slit 
    slit_mask = im
    slit_mask[*] = !values.d_nan
    slit_mask[w_slit] = 1.0

    ; create a new frame where everything outside the slit is a Nan 
    im_masked = im * slit_mask

    ; convert indices to 2D
    w_slit2d = array_indices(slitim, w_slit)

    ; calculate upper and lower limits
    max_y = max(w_slit2d[1,*])
    min_y = min(w_slit2d[1,*])

    ; extract the slit as a rectangle
    this_slit = im_masked[ * , min_y:max_y]

    ; write this out
    writefits, output_filenames[i_frame], this_slit, header

  endfor

END


;******************************************************************


PRO flame_getslits_cutout, fuel=fuel

  ; extract slits structure
  slits = *fuel.slits

  ; read in the slitim image
  slitim = readfits(fuel.intermediate_dir + fuel.slitim_filename)

  print,'Slits: ', n_elements(slits)
  for i_slit=0, n_elements(slits)-1 do begin
  
    print,'Working on slit ', slits[i_slit].name, ' - ', slits[i_slit].number
    
    ; create directory
    slitdir = fuel.intermediate_dir + 'slit' + string(slits[i_slit].number,format='(I02)') + '/'
    spawn,'rm -rf ' + slitdir
    if file_test(slitdir) eq 0 then spawn, 'mkdir ' + slitdir

    ; file names for the cutouts
    output_filenames = strarr(fuel.n_frames)
    for i_frame=0, fuel.n_frames-1 do begin
      naked_filename = ( strsplit((*fuel.corrscience_files)[i_frame], '/', /extract) )[-1]
      output_filenames[i_frame] = flame_util_replace_string( slitdir + naked_filename, '.fits', '_slit' + string(slits[i_slit].number,format='(I02)') + '.fits'  )
    endfor

    ; extract slit
    print,'*** Cutting out slit ', slits[i_slit].name
    flame_getslits_cutout_extract, slitim, slits[i_slit], (*fuel.corrscience_files), output_filenames

    ; add filenames to the slit structure
    *slits[i_slit].filenames = output_filenames

  endfor

END
  

;******************************************************************

     
PRO flame_getslits, fuel=fuel

  ; identify all slits from the data, and write the fuel.slits structures
  flame_getslits_findedges, fuel=fuel

  ; write ds9 region file with the slit traces
  flame_getslits_writeds9, fuel=fuel

  ; write slitim (FITS file with slit image)
  flame_getslits_write_slitim, fuel=fuel

  ; if we are reducing only one slit, then delete all the others 
  if fuel.reduce_only_oneslit ne 0 then $
    *fuel.slits = (*fuel.slits)[fuel.reduce_only_oneslit-1]

  ; cutout slits 
  flame_getslits_cutout, fuel=fuel

END