
PRO flame_initialize_luci_slits, header, pixel_scale=pixel_scale, $
  slit_num=slit_num, slit_name=slit_name, slit_PA=slit_PA, $
  bottom=bottom, top=top, target=target, wavelength_lo=wavelength_lo, wavelength_hi=wavelength_hi
;
; read the header of a LUCI science frame
; and for each slit finds or calculate the slit number, name, slit PA, 
; bottom pixel position, top pixel position, target pixel position, 
; and wavelength at the center of the slit
;

  ; only tweakable settings: maximum width in arcsec allowed for science slits. 
  ; slits wider than this value are considered alignment boxes
  maxwidth_arcsec = 2.0

  ; determine appropriate value of conversion factor delta_wavel in micron/mm
  band = strtrim(fxpar(header, 'FILTER2'),2)
  case band of
    'J': delta_wavel = 0.00049
    'H': delta_wavel = 0.00066
    'K': delta_wavel = 0.001075
    else: message, 'I do not have the wavelength scale for this band yet'
  endcase

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
  N_pixel_y = fxpar(header, 'NAXIS2')

  deltay_arcsec = slit_hdr.y_mm * mmtoarcsec  ; how many arcsec from the mask center; positive going down
  deltay_pixels = deltay_arcsec / pixel_scale ; how many pixels from the mask center; positive going down
  y_pixels = N_pixel_y / 2 - deltay_pixels    ; y-pixel position, usual convention

  ; find the height in pixels of each slit
  slitheight_pixels = slit_hdr.length_arcsec / pixel_scale

  ; output interesting values
  slit_num = slit_hdr.number
  slit_name = slit_hdr.name
  slit_PA = SLIT_hdr.angle
  bottom = y_pixels - 0.5*slitheight_pixels
  top = y_pixels + 0.5*slitheight_pixels
  target = y_pixels
  ; rough wavelength range; mask is roughly 300mm across
  wavelength_lo = fxpar(header, 'GRATWLEN') - slit_hdr.x_mm * delta_wavel - 150.0*delta_wavel
  wavelength_hi = fxpar(header, 'GRATWLEN') - slit_hdr.x_mm * delta_wavel + 150.0*delta_wavel


END



;******************************************************************



PRO flame_initialize_luci, fuel=fuel
  ;
  ; LUCI-specific routine that initializes the fuel structure after
  ; the user has changed the default values.
  ;

  ; check and setup directory structure
  if file_test(fuel.intermediate_dir) eq 0 then spawn, 'mkdir ' + fuel.intermediate_dir
  if file_test(fuel.output_dir) eq 0 then spawn, 'mkdir ' + fuel.output_dir

  ; set the number of science frames
  fuel.N_frames = file_lines(fuel.science_filelist)

  ; create file name for corrected science frames
  readcol, fuel.science_filelist, science_filenames, format='A'
  corrscience_files = strarr(fuel.N_frames)
  for i_frame=0, fuel.N_frames-1 do begin
    components = strsplit( science_filenames[i_frame], '/', /extract)
    out_filename = components[-1] ; keep just the filename
    out_filename = flame_util_replace_string( out_filename, '.fits', '_corr.fits' )
    print, i_frame
    corrscience_files[i_frame] = fuel.intermediate_dir + out_filename
  endfor

  ; save it to fuel
  *fuel.corrscience_files = corrscience_files

  ; determine what band we are in
  science_header = headfits(science_filenames[0])
  fuel.band =  strtrim(fxpar(science_header, 'FILTER2'),2)

  ; read in the pixel scale 
  fuel.pixel_scale = fxpar(science_header, 'PIXSCALE')  ; arcsec/pixel

  ; read in read-out noise
  fuel.readnoise = fxpar(science_header, 'RDNOISE')   ; e-/read 

  ; read in gain
  fuel.gain = fxpar(science_header, 'GAIN')   ; e-/adu

  ; set the wavelength grid for the final, resampled spectra (unless it's already been specified)
  if fuel.output_lambda_0 eq 0d then $
    case fuel.band of

      'J': begin
        fuel.output_lambda_0 = 1.13d           ; wavelength for the first pixel, in micron
        fuel.output_lambda_delta = 0.75d-4     ; micron per pixel
        fuel.output_lambda_Npix =  3200        ; number of pixels (determines wavelength range)
      end

      'K': begin
        fuel.output_lambda_0 = 1.90d           ; wavelength for the first pixel, in micron
        fuel.output_lambda_delta = 1.5d-4      ; micron per pixel
        fuel.output_lambda_Npix =  4200        ; number of pixels (determines wavelength range)
      end

      else: message, fuel.band + ' band not supported yet'

    endcase

  ; read central wavelength 
  central_wavelength = fxpar(science_header, 'GRATWLEN')

  ; read in the slits parameters from the FITS header
  ; LONGSLIT ---------------------------------------------------------------------
  if fuel.longslit then begin
    ; get the slit edges from the user

    ; create slit structure 
    slits = { $
      number:1, $
      name:'longslit', $
      PA:!values.d_nan, $
      approx_bottom:fuel.longslit_edge[0], $
      approx_top:fuel.longslit_edge[1], $
      approx_target:mean(fuel.longslit_edge), $
      approx_wavelength_lo:central_wavelength, $
      approx_wavelength_hi:central_wavelength }

  ; MOS      ---------------------------------------------------------------------
  endif else begin

    flame_initialize_luci_slits, science_header, pixel_scale=fuel.pixel_scale, $
      slit_num=slit_num, slit_name=slit_name, slit_PA=slit_PA, $
      bottom=bottom, top=top, target=target, wavelength_lo=wavelength_lo, wavelength_hi=wavelength_hi

    ; create array of slit structures
    slits = []
    ; trace the edges of the slits using the sky emission lines
    for i_slit=0, n_elements(bottom)-1 do begin

      this_slit = { $
        number:slit_num[i_slit], $
        name:slit_name[i_slit], $
        PA:slit_PA[i_slit], $
        approx_bottom:bottom[i_slit], $
        approx_top:top[i_slit], $
        approx_target:target[i_slit], $
        approx_wavelength_lo:wavelength_lo[i_slit], $
        approx_wavelength_hi:wavelength_hi[i_slit] }

      slits = [slits, this_slit]

    endfor

  endelse

  ; save the slits in the fuel structure - this is the "slits_fromheader" structure, not "slits"!
  *fuel.slits_fromheader = slits


END