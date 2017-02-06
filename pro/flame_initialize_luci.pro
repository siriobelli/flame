


FUNCTION flame_initialize_luci_resolution, instrument
;
; return the estimated spectral resolution for a slit width of 0.75 arcsec
;
  
  ; tabulated values for R - note they assume a slit width of 0.5 arcsec

  ; grating G210
  if (strsplit(instrument.grating, /extract))[0] eq 'G210' then case instrument.grating_order of
    '2': R = 5000. ; K 
    '3': R = 5900. ; H 
    '4': R = 5800. ; J 
    '5': R = 5400. ; z
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase
       
  ; grating G200
  if (strsplit(instrument.grating, /extract))[0] eq 'G200' then case instrument.grating_order of
    '1': R = 2250. ; HK 
    '2': R = 2250. ; zJ 
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase
       
  ; grating G150
  if (strsplit(instrument.grating, /extract))[0] eq 'G150' then case instrument.grating_order of
    '2': R = 4150. ; Ks 
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase
  
  ; check whether no grating has been found 
  if R EQ !NULL then message, 'grating ' + instrument.grating + ' not supported'

  ;
  ;  ARE WE SURE THAT R IS INDIPENDENT ON THE CAMERA??
  ;

  ; output the resolution for a 1"-wide slit:
  R_1arcsec = R / 2.0 

  return, R_1arcsec

END


;******************************************************************



FUNCTION flame_initialize_luci_longslit, header, instrument=instrument, input=input

  ; rough wavelength range (slit should be central, x ~ 162 mm)
    lambda_range = $
     flame_initialize_luci_waverange(instrument, 162.0)

    ; create slit structure 
    slits = { $
      number:1, $
      name:'longslit', $
      PA:!values.d_nan, $
      approx_bottom:input.longslit_edge[0], $
      approx_top:input.longslit_edge[1], $
      approx_target:mean(input.longslit_edge), $
      approx_wavelength_lo:lambda_range[0], $
      approx_wavelength_hi:lambda_range[1] }

  return, slits

END 



;******************************************************************



FUNCTION flame_initialize_luci_waverange, instrument, slit_xmm
;
; return the estimated wavelength range for a slit
;


  ; look up the wavelength range (in um) covered by the detector: -------------------

  ; grating G210
  if (strsplit(instrument.grating, /extract))[0] eq 'G210' then case instrument.grating_order of
    '2': wavelength_range = 0.328 ; K 
    '3': wavelength_range = 0.202 ; H 
    '4': wavelength_range = 0.150 ; J 
    '5': wavelength_range = 0.124 ; z
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase
       
  ; grating G200
  if (strsplit(instrument.grating, /extract))[0] eq 'G200' then case instrument.grating_order of
    '1': wavelength_range = 0.880 ; HK 
    '2': wavelength_range = 0.440 ; zJ 
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase
       
  ; grating G150
  if (strsplit(instrument.grating, /extract))[0] eq 'G150' then case instrument.grating_order of
    '2': wavelength_range = 0.533 ; Ks 
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase
  
  ; check whether no grating has been found 
  if wavelength_range EQ !NULL then message, 'grating ' + instrument.grating + ' not supported'

  ; finally, if the camera is not the default N1.8, then scale by the ratio of pixel scales
  if (strsplit(instrument.camera, /extract))[0] NE 'N1.8' then $
    if (strsplit(instrument.camera, /extract))[0] eq 'N3.75' then wavelength_range *= 0.47 else $
      message, 'camera ' + instrument.camera + ' not supported'

; -------------------------------------------------------------------------------------

  ; rough wavelength range; mask is roughly 2*162mm across; x=0mm means centered slit.
  lambda_min = instrument.central_wavelength - (162.0 + slit_xmm) / 162.0 * wavelength_range / 2.0
  lambda_max = lambda_min + wavelength_range


  return, [lambda_min, lambda_max]

END


;******************************************************************

FUNCTION flame_initialize_luci_settings, science_header
;
; read the LUCI settings from the FITS header of a science frame
; and create the instrument structure
;

  ; read things from header - - - - - - - - - - - - - - - - - - - 

  ; read instrument name - useful to discriminate between LUCI1 and LUCI2
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)

  ; read grating name
  grating = strtrim(fxpar(science_header, 'GRATNAME'), 2)

  ; read in grating order
  grating_order = strtrim(fxpar(science_header, 'GRATORDE'), 2)

  ; read central wavelength 
  central_wavelength = fxpar(science_header, 'GRATWLEN')

  ; read camera
  camera = fxpar(science_header, 'CAMERA')

  ; read pixel scale 
  pixel_scale = fxpar(science_header, 'PIXSCALE')  ; arcsec/pixel

  ; read filters
  filter1 = strtrim(fxpar(science_header, 'FILTER1'), 2)
  filter2 = strtrim(fxpar(science_header, 'FILTER2'), 2)

  ; read in read-out noise
  readnoise = fxpar(science_header, 'RDNOISE')   ; e-/read 

  ; read in gain
  gain = fxpar(science_header, 'GAIN')   ; e-/adu



  ; calculate things from hard-coded numbers - - - - - - - - - - - - - - - - - - - 

  ; linearity correction: the polynomial coefficients describing the transformation
  linearity_correction = [0.0d, 1.0d, 4.155d-6]

  ; calibration files for when the user doesn't have them - - - - - - - - - - - - - - - - - - - 
  default_badpixel_mask = 'default_badpixel_mask_' + instrument_name + '.fits'
  default_flat_field = 'default_flat_field_' + instrument_name + '.fits'
  

  ; create the instrument structure - - - - - - - - - - - - - - - - - - - 
  instrument = { $
    instrument_name: instrument_name, $
    grating: grating, $
    grating_order: grating_order, $  
    central_wavelength: central_wavelength, $
    camera: camera, $
    pixel_scale: pixel_scale, $
    filter1: filter1, $
    filter2: filter2, $
    readnoise: readnoise, $
    gain: gain, $
    resolution_slit1arcsec: 0.0, $
    linearity_correction: linearity_correction, $
    default_badpixel_mask: default_badpixel_mask, $
    default_flatfield: default_flat_field $
    }

  ; now use the instrument structure to calculate the spectral resolution 
  instrument.resolution_slit1arcsec = flame_initialize_luci_resolution(instrument)

  return, instrument


END



;******************************************************************



FUNCTION flame_initialize_luci_slits, header, instrument=instrument

;
; read the header of a LUCI science frame
; and for each slit finds or calculate the slit number, name, slit PA, 
; bottom pixel position, top pixel position, target pixel position, 
; and wavelength at the center of the slit. It return the slits structure.
;

  ; only tweakable settings: maximum width in arcsec allowed for science slits. 
  ; slits wider than this value are considered alignment boxes
  maxwidth_arcsec = 2.0

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
  deltay_pixels = deltay_arcsec / instrument.pixel_scale ; how many pixels from the mask center; positive going down
  y_pixels = N_pixel_y / 2 - deltay_pixels    ; y-pixel position, usual convention

  ; find the height in pixels of each slit
  slitheight_pixels = slit_hdr.length_arcsec / instrument.pixel_scale

  ; output interesting values
  slit_num = slit_hdr.number
  slit_name = slit_hdr.name
  slit_PA = slit_hdr.angle
  bottom = y_pixels - 0.5*slitheight_pixels
  top = y_pixels + 0.5*slitheight_pixels
  target = y_pixels


  ; create array of slit structures
  slits = []
  ; trace the edges of the slits using the sky emission lines
  for i_slit=0, n_elements(bottom)-1 do begin

    ; calculate approximate wavelength range
    lambda_range = flame_initialize_luci_waverange(instrument, slit_hdr[i_slit].x_mm)

    this_slit = { $
      number:slit_num[i_slit], $
      name:slit_name[i_slit], $
      PA:slit_PA[i_slit], $
      approx_bottom:bottom[i_slit], $
      approx_top:top[i_slit], $
      approx_target:target[i_slit], $
      approx_wavelength_lo:lambda_range[0], $
      approx_wavelength_hi:lambda_range[1] }

    slits = [slits, this_slit]

  endfor

return, slits

END




;******************************************************************



PRO flame_initialize_luci, fuel=fuel
  ;
  ; LUCI-specific routine that initializes the fuel.instrument structure
  ;

  ; read FITS header of first science frame
  science_header = headfits(fuel.util.science_filenames[0])

  ; read the instrument settings from the header 
  instrument = flame_initialize_luci_settings(science_header)

  ; now read in the slits parameters from the FITS header

  ; LONGSLIT ---------------------------------------------------------------------
  if fuel.input.longslit then begin

    ; get the slit edges from the user
    slits = flame_initialize_luci_longslit( header, instrument=instrument, input=fuel.input)


  ; MOS      ---------------------------------------------------------------------
  endif else begin

    ; get all the info from the header
    slits = flame_initialize_luci_slits( science_header, instrument=instrument)

  endelse

  ; save both instrument and slits in the fuel structure
  new_fuel = { input:fuel.input, util:fuel.util, instrument:instrument, diagnostics:fuel.diagnostics, slits:slits }
  fuel=new_fuel    
  

END