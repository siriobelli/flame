;
; minimal initialization for MOSFIRE - its only purpose is to
; run the first two steps of flame for on-the-fly diagnostics
;
;******************************************************************

FUNCTION flame_initialize_mosfire_settings, science_header
;
; read the MOSFIRE settings from the FITS header of a science frame
; and create the instrument structure
;

  ; read things from header - - - - - - - - - - - - - - - - - - -

  ; read instrument name
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)


  ; calculate things from hard-coded numbers - - - - - - - - - - - - - - - - - - -

  pixel_scale = 0.18  ; arcsec / pixel

  ; linearity correction: the polynomial coefficients describing the transformation
  linearity_correction = [0.0d, 1.0d]   ; do not have linearity correction


  ; calibration files for when the user doesn't have them - - - - - - - - - - - - - - - - - - -
  default_badpixel_mask = 'none'
  default_dark = 'none'
  default_pixelflat = 'none'
  default_illumflat = 'none'
  default_arc = 'none'


  ; create the instrument structure - - - - - - - - - - - - - - - - - - -
  instrument = { $
    instrument_name: instrument_name, $
    pixel_scale: pixel_scale, $
    resolution_slit1arcsec: 3000.0, $
    linearity_correction: linearity_correction, $
    trim_edges: 4, $
    default_badpixel_mask: default_badpixel_mask, $
    default_dark: default_dark, $
    default_pixelflat: default_pixelflat, $
    default_illumflat: default_illumflat, $
    default_arc: default_arc $
    }

  return, instrument


END



;******************************************************************



FUNCTION flame_initialize_mosfire_slits, header, instrument=instrument, input=input

  return, {number:0, name:'none'}

END




;******************************************************************



FUNCTION flame_initialize_mosfire, input
  ;
  ; MOSFIRE-specific routine that initializes the fuel structure
  ;

  print, ''
  print, 'Initializing MOSFIRE data reduction'
  print, ''


  ; first, create the fuel structure
  fuel = flame_util_create_fuel(input)


  ; ---------------------   INSTRUMENT structure   --------------------------------------

  ; read FITS header of first science frame
  science_header = headfits(fuel.util.science.raw_files[0])

  ; read the instrument settings from the header
  instrument = flame_initialize_mosfire_settings(science_header)


  ; now read in the slits parameters from the FITS header
  slits = flame_initialize_mosfire_slits( science_header, instrument=instrument, input=fuel.input)


  ; ---------------------   SETTINGS   --------------------------------------

  ; do not split the spectrum into two when doing the rough wavecal
  fuel.settings.wavecal_rough_split = 0

  ; Use the R3000 line list
  fuel.settings.linelist_filename = fuel.util.flame_data_dir + 'line_list_R3000.dat'


  ; -------------------------------------------------------------------------------


  ; save both instrument and slits in the fuel structure
  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:instrument, diagnostics:fuel.diagnostics, slits:slits }

  ; return the fuel structure
  return, new_fuel

END
