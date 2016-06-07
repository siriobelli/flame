
FUNCTION flame_trace_slits_smear, array, filter_length, up=up, down=down
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



FUNCTION flame_trace_edge, image, approx_edge, top=top, bottom=bottom 

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
  binary_mask_smeared_up = flame_trace_slits_smear( binary_mask, filter_length, /up)

  ; if the sum is less than 1+1+1+... X times, then there was a zero 
  condition_A = binary_mask_smeared_up*0.0
  condition_A[where(binary_mask_smeared_up eq 5.0, /null)] = 1.0

  ; the second condition is: there are only zeroes among the next X pixels
  binary_mask_smeared_down = flame_trace_slits_smear( binary_mask, filter_length, /down)

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





PRO flame_trace_slit, image=image, approx_top=approx_top, approx_bottom=approx_bottom, $
  poly_coeff=poly_coeff, slit_height=slit_height
;
; traces the top and bottom edges of a slit in the image and return the slit height and the 
; coefficients of a polynomial fit describing the *bottom* edge
;

  ; identify top edge
  top_edge = flame_trace_edge(image, approx_top, /top )

  ; identify bottom edge
  bottom_edge = flame_trace_edge(image, approx_bottom, /bottom )

  ; calculate the slit height
  slit_height = median( top_edge - bottom_edge )

  ; combine the top and bottom edge measurements together (i.e. subtract slit_heigh from the top edge)
  x_to_fit = [ where(finite(top_edge), /null ), where(finite(bottom_edge), /null) ]
  y_to_fit = [ top_edge[ where(finite(top_edge), /null) ] - slit_height, bottom_edge[where(finite(bottom_edge), /null)] ]

  ; fit a 3-rd order polynomial to the combined edges
  poly_coeff = robust_poly_fit( x_to_fit, y_to_fit, 3 )


END

