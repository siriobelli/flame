
;******************************************************************

FUNCTION flame_initialize_lris_settings, science_header
;
; read the LRIS settings from the FITS header of a science frame
; and create the instrument structure
;

  ; read things from header - - - - - - - - - - - - - - - - - - - 

  ; read instrument name - useful to discriminate between blue and red channel
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)

  ; read grating name
  grating = strtrim(fxpar(science_header, 'GRANAME'), 2)

  ; read in grating angle
  grating_angle = fxpar(science_header, 'GRANGLE')

  ; read grism name
  grism = strtrim(fxpar(science_header, 'GRISNAME'), 2)

  ; read central wavelength for longslit
  central_wavelength = fxpar(science_header, 'WAVELEN')

  ; read central wavelength for slitmasks
  central_wavelength_MOS = fxpar(science_header, 'MSWAVE')

  ; read filters
  blue_filter = strtrim(fxpar(science_header, 'BLUFILT'), 2)
  red_filter = strtrim(fxpar(science_header, 'REDFILT'), 2)
  
  ; ; read in read-out noise
  ; readnoise = fxpar(science_header, 'RDNOISE')   ; e-/read 


  ; ; calculate things from hard-coded numbers - - - - - - - - - - - - - - - - - - - 

  ; gain - for now let's set it one
  gain = 1.0   ; e-/adu

  ; pixel scale 
  pixel_scale = 0.135  ; arcsec/pixel

  ; linearity correction: the polynomial coefficients describing the transformation
  linearity_correction = [0.0d, 1.0d]   ; CCDs are almost exactly linear

  ; calibration files for when the user doesn't have them - - - - - - - - - - - - - - - - - - - 
  default_badpixel_mask = 'none'
  default_flat_field = 'none'
  

  ; create the instrument structure - - - - - - - - - - - - - - - - - - - 
  instrument = { $
    instrument_name: instrument_name, $
    grating: grating, $
    grating_angle: grating_angle, $  
    grism: grism, $
    central_wavelength: central_wavelength, $
    central_wavelength_MOS: central_wavelength_MOS, $
    blue_filter: blue_filter, $
    red_filter: red_filter, $
    gain: gain, $
    pixel_scale: pixel_scale, $
    linearity_correction: linearity_correction, $
    default_badpixel_mask: default_badpixel_mask, $
    default_flatfield: default_flat_field $
    }

  ; ; now use the instrument structure to calculate the spectral resolution 
  ; instrument.resolution_slit1arcsec = flame_initialize_luci_resolution(instrument)

  return, instrument


END



;******************************************************************




FUNCTION flame_initialize_lris_slits, header, instrument=instrument, slit_y=slit_y
;
; For LRIS, need to start with the list of y-pixel coordinate of the approximate slit edges (slit_y).
; Then for each slit finds the slit number,
; bottom pixel position, top pixel position, target pixel position, 
; and wavelength at the center of the slit. It returns the slits structure.
;
  
  ; sort slit_y 
  slit_y = slit_y[ sort(slit_y) ]

  ; create array of slit structures
  slits = []
  ; trace the edges of the slits using the sky emission lines
  for i_slit=0, n_elements(slit_y)/2-1 do begin

  ;   ; calculate approximate wavelength range
  ;   lambda_range = flame_initialize_luci_waverange(instrument, slit_hdr[i_slit].x_mm)

    this_slit = { $
      number:i_slit+1, $
      name:'', $
      PA:!values.d_NaN, $
      approx_bottom:slit_y[2*i_slit], $
      approx_top:slit_y[2*i_slit+1], $
      approx_target: 0.5 * (slit_y[2*i_slit] + slit_y[2*i_slit+1]), $
      approx_wavelength_lo:!values.d_NaN, $
      approx_wavelength_hi:!values.d_NaN }

    slits = [slits, this_slit]

  endfor

  ; ; calculate the height of each slit
  ; slit_height = slits.approx_top - slits.approx_bottom

  ; ; rank them, starting from the smallest one
  ; s = sort(slit_height)

  ; ; assume the Nboxes smaller slits are the alignment boxes
  ; box_height = slit_height[s[Nboxes-1]]

  ; ; remove alignment boxes   
  ; slits = slits[where(slit_height GT box_height, /null)]

  ; ; now fill in the slit numbers
  ; slits.number = indgen(n_elements(slits)) + 1

return, slits

END




;******************************************************************



PRO flame_initialize_lris, fuel=fuel
  ;
  ; LRIS-specific routine that initializes the fuel.instrument structure
  ;

  ; read FITS header of first science frame
  science_header = headfits(fuel.util.science_filenames[0])

  ; read the instrument settings from the header 
  instrument = flame_initialize_lris_settings(science_header)

  ; read the approximate slit positions inserted manually
  if ~file_test(fuel.input.slit_position_file) then $
    message, 'input.slit_position_file is required for LRIS data'
  readcol, fuel.input.slit_position_file, slit_x, slit_y

  ; now read in the slits parameters from the FITS header

  ; LONGSLIT ---------------------------------------------------------------------
  if fuel.input.longslit then begin

    message, 'longslit not supported yet'

    ; get the slit edges from the user
   ; slits = flame_initialize_luci_longslit( header, instrument=instrument, input=fuel.input)


  ; MOS      ---------------------------------------------------------------------
  endif else begin

    ; get all the info from the header
    slits = flame_initialize_lris_slits( science_header, instrument=instrument, slit_y=slit_y )

  endelse

  ; save both instrument and slits in the fuel structure
  new_fuel = { input:fuel.input, util:fuel.util, instrument:instrument, diagnostics:fuel.diagnostics, slits:slits }
  fuel=new_fuel    
  

END