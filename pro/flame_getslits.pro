

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
    flame_trace_slit, image=im, approx_top=info_fromheader.approx_top-yshift, $
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