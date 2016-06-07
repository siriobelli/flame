PRO flame_do_slitim_readheader, header, pixel_scale=pixel_scale, bottom=bottom, top=top, slit_num=slit_num, slit_name=slit_name, maxwidth_arcsec=maxwidth_arcsec
;
; return the pixel coordinates of slits as from the header
;

  ; array of structures that will contain the info for each slit
  slit_hdr = []

  i_slit=1
  while fxpar( header, 'TGT' + string(i_slit, format='(I02)')  + 'NAM', missing='NONE' ) NE 'NONE' do begin

    slitnum = string(i_slit, format='(I02)')

    this_slit_hdr = {number : i_slit, $
    name : strtrim( fxpar( header, 'TGT' + slitnum + 'NAM' ), 2), $
    shape : fxpar( header, 'MOS' + slitnum + 'SHA'), $
    width_arcsec : float(fxpar( header, 'MOS' + slitnum + 'WAS')), $
    length_arcsec : float(fxpar( header, 'MOS' + slitnum + 'LAS')), $
    angle : float(fxpar( header, 'MOS' + slitnum + 'PA')), $
    width_mm : float(fxpar(header, 'MOS' + slitnum + 'WMM')), $
    length_mm : float(fxpar(header, 'MOS' + slitnum + 'LMM')), $
    x_mm : float(fxpar(header, 'MOS' + slitnum + 'XPO')), $
    y_mm : float(fxpar(header, 'MOS' + slitnum + 'YPO')) }

    slit_hdr = [ slit_hdr, this_slit_hdr ]

    i_slit++

  endwhile

  ; exclude reference slits
  slit_hdr = slit_hdr( where( strtrim(slit_hdr.name,2) ne 'refslit', /null) )

  ; exclude alignment boxes
  slit_hdr = slit_hdr[where( slit_hdr.width_arcsec LT maxwidth_arcsec, /null) ]

  ; calculate the conversion between arcsec and mm
  mmtoarcsec = median(slit_hdr.length_arcsec/slit_hdr.length_mm)
  

  ; convert from coordinates in mm to coordinates in pixels
  ;---------------------------------------------------------
  N_pixel_y = 2048

  deltay_arcsec = slit_hdr.y_mm * mmtoarcsec  ; how many arcsec from the mask center; positive going down
  deltay_pixels = deltay_arcsec / pixel_scale ; how many pixels from the mask center; positive going down
  y_pixels = N_pixel_y / 2 - deltay_pixels    ; y-pixel position, usual convention

  ; find the height in pixels of each slit
  slitheight_pixels = slit_hdr.length_arcsec / pixel_scale

  ; output interesting values
  bottom = y_pixels - 0.5*slitheight_pixels
  top = y_pixels + 0.5*slitheight_pixels
  slit_num = slit_hdr.number
  slit_name = slit_hdr.name

END


;******************************************************************


FUNCTION flame_do_slitim_findshift, frame, tpx, bpx
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
  expected_edges[tpx] = 1.
  expected_edges[bpx] = -1.
  
  ; cross-correlate to find the shift between expected and measured slit edges
  lag = indgen(400)-200
  crosscorr = c_correlate(edges, expected_edges, lag)
  max_crosscorr = max( crosscorr, max_ind, /nan)
  delta = lag[max_ind]
  
  ; return the shift
  return, delta
  
END


;******************************************************************

PRO flame_do_slitim_writeds9, slits, filename=filename
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

pro flame_do_slitim, fuel=fuel

  ; read in the frame
  im=readfits((*fuel.corrscience_files)[0], hdr)

  ; LONGSLIT ---------------------------------------------------------------------
  if fuel.longslit then begin

    ; get the slit edges from the user
    slit_top = fuel.longslit_edge[1]
    slit_bottom = fuel.longslit_edge[0]

    ; create slit structure 
    slits = {number:1, name:'longslit', height:slit_top-slit_bottom, $
      bottom_poly:[slit_bottom], filenames:ptr_new(/allocate_heap), rectification:ptr_new(/allocate_heap)}
   
  ; MOS      ---------------------------------------------------------------------
  endif else begin
  
    ; get the slit positions from the header
    flame_do_slitim_readheader, hdr, pixel_scale=fuel.pixel_scale, $
      bottom=bpx, top=tpx, slit_num=slit_num, slit_name=slit_name, $
      maxwidth_arcsec=fuel.maxwidth_arcsec
    
    ; compare the expected position with the measured ones and apply rough shift
    shift = flame_do_slitim_findshift((*fuel.corrscience_files)[0], bpx, tpx)
    bpx -= shift
    tpx -= shift
    
    ; create array of slit structures
    slits = []
    ; trace the edges of the slits using the sky emission lines
    for i_slit=0, n_elements(bpx)-1 do begin

      flame_trace_slit, image=im, approx_top=tpx[i_slit], approx_bottom=bpx[i_slit], $
        poly_coeff=poly_coeff, slit_height=slit_height

      this_slit = {number:slit_num[i_slit], name:slit_name[i_slit], height:slit_height, $
        bottom_poly:poly_coeff, filenames:ptr_new(/allocate_heap), rectification:ptr_new(/allocate_heap)}
      slits = [slits, this_slit]

    endfor

    ; if we are reducing only one slit, then exclude all the others 
    if fuel.reduce_only_oneslit ne 0 then slits = slits[fuel.reduce_only_oneslit-1]

  endelse

  ; save the slits in the fuel structure 
  *fuel.slits = slits

  ; write ds9 region file with the slit traces
  flame_do_slitim_writeds9, slits, filename=fuel.intermediate_dir+'slits.reg'

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

  for i_slit=0,n_elements(slits)-1 do begin

    top_y = (poly(x_axis, slits[i_slit].bottom_poly) + slits[i_slit].height) # replicate(1, N_pix_x)
    bottom_y = poly(x_axis, slits[i_slit].bottom_poly) # replicate(1, N_pix_x)

    w_slit = where( pixel_y LT top_y AND pixel_y GT bottom_y, /null)
    slitim[w_slit] = slits[i_slit].number
   
  endfor
  
  writefits, fuel.intermediate_dir + fuel.slitim_filename, slitim
  
end
