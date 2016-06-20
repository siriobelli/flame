

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


PRO flame_getslits_writeds9, slits, filename=filename
  ;
  ; write a ds9 region files that shows the slit edges
  ;

  ; number of horizontal pixel in one frame
  N_pix_x = 2048

  ; open file
  openw, lun, filename, /get_lun

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
      

PRO flame_getslits, fuel=fuel

  ; read in the frame
  im=readfits((*fuel.corrscience_files)[0], hdr)
  
  ; create array of new slit structures
  new_slits = []

  ; LONGSLIT ---------------------------------------------------------------------
  if fuel.longslit then begin
    	
    new_slits = create_struct( (*fuel.slits)[0], $
    	'yshift', !values.d_nan, $
    	'bottom_poly', !values.d_nan, $
    	'filenames', ptr_new(/allocate_heap), $
    	'rectification', ptr_new(/allocate_heap) )

 ; MOS      ---------------------------------------------------------------------
  endif else begin

    ; compare the expected position with the measured ones and obtain rough shift
    yshift = flame_getslits_findshift( (*fuel.corrscience_files)[0], $
    	(*fuel.slits).approx_top, (*fuel.slits).approx_bottom )
 
    ; trace the edges of the slits using the sky emission lines
    for i_slit=0, n_elements((*fuel.slits))-1 do begin

		flame_trace_slit, image=im, approx_top=(*fuel.slits)[i_slit].approx_top-yshift, $
			approx_bottom=(*fuel.slits)[i_slit].approx_bottom-yshift, $
    	    poly_coeff=poly_coeff, slit_height=slit_height

    	; add new fields to slit structure
    	this_slit = create_struct( (*fuel.slits)[i_slit], $
    		'yshift', yshift, $
    		'height', slit_height, $
    		'bottom_poly', poly_coeff, $
    		'filenames', ptr_new(/allocate_heap), $
    		'rectification', ptr_new(/allocate_heap) )

	    new_slits = [new_slits, this_slit]

	endfor


  endelse

  ; replace the slit structures in fuel
  *fuel.slits = new_slits

  ; write ds9 region file with the slit traces
  flame_getslits_writeds9, new_slits, filename=fuel.intermediate_dir+'slits.reg'

  ; make an image where the pixels belonging to a slit have the slit number as a value, otherwise zero
  slitim=im
  slitim[*]=0

  ; construct the coordinates for the pixels in the image
  N_pix_x = (size(im))[2]
  N_pix_y = (size(im))[1]
  x_axis = indgen(N_pix_x)
  y_axis = indgen(N_pix_y)
  pixel_x = x_axis # replicate(1, N_pix_y)
  pixel_y = transpose(y_axis # replicate(1, N_pix_x))

  for i_slit=0,n_elements(new_slits)-1 do begin

    top_y = (poly(x_axis, new_slits[i_slit].bottom_poly) + new_slits[i_slit].height) # replicate(1, N_pix_x)
    bottom_y = poly(x_axis, new_slits[i_slit].bottom_poly) # replicate(1, N_pix_x)

    w_slit = where( pixel_y LT top_y AND pixel_y GT bottom_y, /null)
    slitim[w_slit] = new_slits[i_slit].number
   
  endfor
  
  writefits, fuel.intermediate_dir + fuel.slitim_filename, slitim
  


END