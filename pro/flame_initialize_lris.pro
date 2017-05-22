
;******************************************************************

FUNCTION flame_initialize_lris_settings, science_header
;
; read the LRIS settings from the FITS header of a science frame
; and create the instrument structure
;

  ; read things from header - - - - - - - - - - - - - - - - - - -

  ; read instrument name - useful to discriminate between blue and red channel
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)

  ; read detector name - useful to discriminate between blue and red channel
  detector = strtrim(fxpar(science_header, 'DETECTOR'), 2)

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
  readnoise = 5.0 ; actually LRIS blue has 4.0

  ; ; calculate things from hard-coded numbers - - - - - - - - - - - - - - - - - - -

  ; gain: since it is different for the four amplifiers,
  ; it needs special handling with from readmhdufits.pro
  gain = 1.0   ; e-/adu

  ; pixel scale
  pixel_scale = 0.135  ; arcsec/pixel

  ; linearity correction: the polynomial coefficients describing the transformation
  linearity_correction = [0.0d, 1.0d]   ; CCDs are almost exactly linear

  ; calibration files for when the user doesn't have them - - - - - - - - - - - - - - - - - - -
  default_badpixel_mask = 'none'
  default_dark = 'none'
  default_pixelflat = 'none'
  default_illumflat = 'none'
  default_arc = 'none'

  ; create the instrument structure - - - - - - - - - - - - - - - - - - -
  instrument = { $
    instrument_name: instrument_name, $
    detector: detector, $
    grating: grating, $
    grating_angle: grating_angle, $
    grism: grism, $
    central_wavelength: central_wavelength, $
    central_wavelength_MOS: central_wavelength_MOS, $
    blue_filter: blue_filter, $
    red_filter: red_filter, $
    gain: gain, $
    readnoise: readnoise, $
    pixel_scale: pixel_scale, $
    resolution_slit1arcsec: 2000.0, $
    linearity_correction: linearity_correction, $
    default_badpixel_mask: default_badpixel_mask, $
    default_dark: default_dark, $
    default_pixelflat: default_pixelflat, $
    default_illumflat: default_illumflat, $
    default_arc: default_arc $
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

  ; NEED TO calculate these values for different settings
  range_lambda0 = [0.4, 0.7]
  range_pixel_scale = [5d-5, 10d-5]


    this_slit = { $
      number:i_slit+1, $
      name:'', $
      PA:!values.d_NaN, $
      approx_bottom:slit_y[2*i_slit], $
      approx_top:slit_y[2*i_slit+1], $
      approx_target: 0.5 * (slit_y[2*i_slit] + slit_y[2*i_slit+1]), $
      width_arcsec:!values.d_NaN, $
      approx_R:instrument.resolution_slit1arcsec, $
      range_lambda0:range_lambda0, $
      range_pixel_scale:range_pixel_scale }

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



FUNCTION flame_initialize_lris_gain, instrument
  ;
  ; calculates tabulated gain (for each of the 4 amplifiers)
  ; according to the settings saved in the instrument structure
  ;

  ; LRIS - RED
  if (strsplit(instrument.detector, /extract))[0] eq 'LRIS-R' then begin

    gain = [1.255, 1.180, 1.191, 1.162]

  ; LRIS - BLUE
  endif else begin

    gain = [1.55,	1.56,	1.63, 1.70]

  endelse

  ; make a structure that will be read by readmhdufits.pro
  vidinp = ['VidInp1','VidInp2','VidInp3','VidInp4']
  gaindata = replicate({GAINDATA, vidinp:'', gain:0.}, 4)
  for i=0,3 do gaindata[i] = {GAINDATA, vidinp:vidinp[i], gain:gain[i]}

  return, gaindata

  ; ALSO WHAT ABOUT CCDGAIN AND CCDSPEED???

END


;******************************************************************



PRO flame_initialize_lris, fuel
  ;
  ; LRIS-specific routine that initializes the fuel.instrument structure
  ;

    print, ''
    print, 'Initializing LRIS data reduction'
    print, ''

  ; -----------------------------------------------------------------
  ; need to use the optical sky spectrum
  ; HANDLE THIS IN A BETTER WAY PLEASE
  fuel.util.sky_emission_filename = $
    fuel.util.flame_data_dir + 'sky_emission_model_optical.dat'

  ; -----------------------------------------------------------------

    ; read FITS header of first science frame
    science_header = headfits(fuel.util.science_filenames[0])

    ; read the instrument settings from the header
    instrument = flame_initialize_lris_settings(science_header)

  ; -----------------------------------------------------------------
  ; we need to convert the LRIS frames into the "normal" format

  ; the conversion also applies the gain correction to each amplifier
  gaindata = flame_initialize_lris_gain(instrument)

  ; create directory where we will store the new FITS files
  lris_dir = fuel.input.intermediate_dir + 'data_LRIS/'
  if file_test(lris_dir) then file_delete, lris_dir, /recursive
  file_mkdir, lris_dir

  ; create new file names for the science frames
  new_filenames = lris_dir + file_basename(fuel.util.science_filenames)

  ; use the readmhdufits.pro routine provided by Keck
  for i=0, fuel.util.N_frames-1 do begin
    image = readmhdufits(fuel.util.science_filenames[i], header=header, gaindata=gaindata)
    image = transpose(image)
    mwrfits, image, new_filenames[i], header
  endfor

  ; update the file names in the fuel structure
  fuel.util.science_filenames = new_filenames

  ; -----------------------------------------------------------------

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
