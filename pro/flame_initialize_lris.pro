
;******************************************************************

FUNCTION flame_initialize_lris_red, science_header
;
; initialize LRIS - RED SIDE
;
; read the LRIS settings from the FITS header of a science frame
; and create the instrument structure
;

  ; read things from header - - - - - - - - - - - - - - - - - - -

  ; read instrument name - useful to discriminate between blue and red channel
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)

  ; check that these are LRIS RED data
  if instrument_name NE 'LRIS' then $
    message, 'instrument name keyword ' + instrument_name + ': expected LRIS instead'

  ; the header will show LRISBLUE for the blue side but only LRIS for the red side. Fix this:
  instrument_name = 'LRISRED'

  ; read detector name - useful to discriminate between blue and red channel
  detector = strtrim(fxpar(science_header, 'DETECTOR'), 2)

  ; read grating name
  grating = strtrim(fxpar(science_header, 'GRANAME'), 2)

  ; read in grating angle
  grating_angle = fxpar(science_header, 'GRANGLE')

  ; read central wavelength for longslit
  central_wavelength = fxpar(science_header, 'WAVELEN') * 1d-4

  ; read central wavelength for slitmasks
  central_wavelength_MOS = fxpar(science_header, 'MSWAVE') * 1d-4

  ; read filter
  filter = strtrim(fxpar(science_header, 'REDFILT'), 2)

  ; read in read-out noise
  readnoise = 5.0


  ; calculate things from hard-coded numbers - - - - - - - - - - - - - - - - -

  ; pixel scale
  pixel_scale = 0.135  ; arcsec/pixel

  ; linearity correction: the polynomial coefficients describing the transformation
  linearity_correction = [0.0d, 1.0d]   ; CCDs are almost exactly linear

  ; gain: since it is different for the four amplifiers,
  ; it needs special handling with readmhdufits.pro

  ; set the gain keyword to one
  gain = 1.0   ; e-/adu

  ; then make gaindata structure
  vidinp = ['VidInp1','VidInp2','VidInp3','VidInp4']
  gain_values = [1.255, 1.180, 1.191, 1.162]
  gaindata = replicate({GAINDATA, vidinp:'', gain:0.}, 4)
  for i=0,3 do gaindata[i] = {GAINDATA, vidinp:vidinp[i], gain:gain_values[i]}


  ; look up tables for dispersion elements - - - - - - - - - - - - - - - - - - -

  case grating of
    '150/7500': dispersion = 3.0d-4
    '300/5000': dispersion = 1.59d-4
    '400/8500': dispersion = 1.16d-4
    '600/5000': dispersion = 0.80d-4
    '600/7500': dispersion = 0.80d-4
    '600/10000': dispersion = 0.80d-4
    '831/8200': dispersion = 0.58d-4
    '900/5500': dispersion = 0.53d-4
    '1200/7500': dispersion = 0.40d-4
    '1200/9000': dispersion = 0.40d-4
    else: message, 'grating ' + grating + 'not supported.'
  endcase

  ; calculate the total observed range in micron
  observed_range = 4096*dispersion

  ; estimate initial wavelength for a slit at the center of the mask
  lambda0 = central_wavelength_MOS - 0.5*observed_range

  ; the FWHM is approximately six pixels, and assume lambda~7000A
  resolution_slit1arcsec = 0.7 / (6.0*dispersion)


  ; effect of binning - - - - - - - - - - - - - - - - - - -

  ; read binning
  binning = fix( strsplit( strtrim(fxpar(science_header, 'BINNING'), 2), ',', /extract) )
  spatial_binning = binning[0]
  spectral_binning = binning[1]

  ; spatial binning: change pixel scale
  if spatial_binning ne 1 then pixel_scale *= float(spatial_binning)

  ; spectral binning: change the dispersion value
  if spectral_binning ne 1 then dispersion *= float(spectral_binning)


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
    central_wavelength: central_wavelength, $
    central_wavelength_MOS: central_wavelength_MOS, $
    filter: filter, $
    gain: gain, $
    gaindata: gaindata, $
    readnoise: readnoise, $
    pixel_scale: pixel_scale, $
    resolution_slit1arcsec: resolution_slit1arcsec, $
    dispersion:dispersion, $
    lambda0:lambda0, $
    linearity_correction: linearity_correction, $
    spatial_binning: spatial_binning, $
    spectral_binning: spectral_binning, $
    default_badpixel_mask: default_badpixel_mask, $
    default_dark: default_dark, $
    default_pixelflat: default_pixelflat, $
    default_illumflat: default_illumflat, $
    default_arc: default_arc $
    }


  return, instrument


END



;******************************************************************




;******************************************************************

FUNCTION flame_initialize_lris_blue, science_header
;
; initialize LRIS - BLUE SIDE
;
; read the LRIS settings from the FITS header of a science frame
; and create the instrument structure
;

  ; read things from header - - - - - - - - - - - - - - - - - - -

  ; read instrument name
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)

  ; check that these are LRIS BLUE data
  if instrument_name NE 'LRISBLUE' then $
    message, 'instrument name keyword ' + instrument_name + ': expected LRISBLUE instead'

  ; read detector name - useful to discriminate between blue and red channel
  detector = strtrim(fxpar(science_header, 'DETECTOR'), 2)

  ; read grism name
  grism = strtrim(fxpar(science_header, 'GRISNAME'), 2)

  ; read filter
  filter = strtrim(fxpar(science_header, 'BLUFILT'), 2)

  ; read in read-out noise
  readnoise = 4.0 ; CHECK THIS


  ; calculate things from hard-coded numbers - - - - - - - - - - - - - - - - -

  ; pixel scale
  pixel_scale = 0.135  ; arcsec/pixel

  ; linearity correction: the polynomial coefficients describing the transformation
  linearity_correction = [0.0d, 1.0d]   ; CCDs are almost exactly linear

  ; gain: since it is different for the four amplifiers,
  ; it needs special handling with readmhdufits.pro

  ; set the gain keyword to one
  gain = 1.0   ; e-/adu

  ; then make gaindata structure
  vidinp = ['VidInp1','VidInp2','VidInp3','VidInp4']
  gain_values = [1.55,	1.56,	1.63, 1.70]
  gaindata = replicate({GAINDATA, vidinp:'', gain:0.}, 4)
  for i=0,3 do gaindata[i] = {GAINDATA, vidinp:vidinp[i], gain:gain_values[i]}


  ; look up tables for dispersion elements - - - - - - - - - - - - - - - - - - -

  case grism of
    '300/5000': begin
      dispersion = 1.43d-4
      lambda0 = 0.2210
    end
    '400/3400': begin
      dispersion = 1.09d-4
      lambda0 = 0.1760
    end
    '600/4000': begin
      dispersion = 0.63d-4
      lambda0 = 0.3300
    end
    '1200/3400': begin
      dispersion = 0.24d-4
      lambda0 = 0.3010
    end
    else: message, 'grism ' + grism + 'not supported.'
  endcase

  ; the FWHM is approximately six unbinned pixels, and assume lambda~4000A
  resolution_slit1arcsec = 0.4 / (6.0*dispersion)


  ; effect of binning - - - - - - - - - - - - - - - - - - -

  ; read binning
  binning = fix( strsplit( strtrim(fxpar(science_header, 'BINNING'), 2), ',', /extract) )
  spatial_binning = binning[0]
  spectral_binning = binning[1]

  ; spatial binning: change pixel scale
  if spatial_binning ne 1 then pixel_scale *= float(spatial_binning)

  ; spectral binning: change the dispersion value
  if spectral_binning ne 1 then dispersion *= float(spectral_binning)


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
    grism: grism, $
    filter: filter, $
    gain: gain, $
    gaindata: gaindata, $
    readnoise: readnoise, $
    pixel_scale: pixel_scale, $
    resolution_slit1arcsec: resolution_slit1arcsec, $
    dispersion:dispersion, $
    lambda0:lambda0, $
    linearity_correction: linearity_correction, $
    spatial_binning: spatial_binning, $
    spectral_binning: spectral_binning, $
    default_badpixel_mask: default_badpixel_mask, $
    default_dark: default_dark, $
    default_pixelflat: default_pixelflat, $
    default_illumflat: default_illumflat, $
    default_arc: default_arc $
    }


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

  ; use the approximate wavelength parameters, and add a bit of uncertainty
  range_lambda0 = instrument.lambda0 * [0.8, 1.20]  ; this depends on the horizontal position of the slit
  range_delta_lambda = instrument.dispersion * [0.8, 1.2]

    this_slit = { $
      number:i_slit+1, $
      name:'NN', $
      skip:0, $
      PA:!values.d_NaN, $
      approx_bottom:slit_y[2*i_slit], $
      approx_top:slit_y[2*i_slit+1], $
      approx_target: 0.5 * (slit_y[2*i_slit] + slit_y[2*i_slit+1]), $
      width_arcsec:!values.d_NaN, $
      approx_R:instrument.resolution_slit1arcsec, $
      range_lambda0:range_lambda0, $
      range_delta_lambda:range_delta_lambda }

    slits = [slits, this_slit]

  endfor

return, slits

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_initialize_lris_arcs, fuel

  ; read first frame
  arc_hdr = headfits(fuel.util.arc.raw_files[0])

  ; old files did not have an explicit list of lamps
  !NULL = sxpar(arc_hdr, 'MERCURY', count=count)

  if count eq 0 then begin  ; old files

    lamp_name = ['MERCURY', 'NEON', 'ARGON', 'CADMIUM', 'ZINC']
    lamp_file = ['arc_Hg_air.txt', 'arc_Ne_air.txt', 'arc_Ar_air.txt', 'arc_Cd_air.txt', $
      'arc_Zn_air.txt' ]

    ; read lamp status from header
    lamp_status = sxpar(arc_hdr, 'LAMPS', count=count)

    ; check that the lamps keyword is in fact existing
    if count eq 0 then message, 'Could not find arc lamps in the FITS header'

    ; transform from string of 1,0 to array
    lamp_status = strsplit(lamp_status, ',', /extract)

    ; delete last entry which is the halogen lamp
    lamp_status = lamp_status[0:4]

  endif else begin  ; new files

    lamp_name = ['MERCURY', 'NEON', 'ARGON', 'CADMIUM', 'ZINC', 'KRYPTON', 'XENON', 'FEARGON']
    lamp_file = ['arc_Hg_air.txt', 'arc_Ne_air.txt', 'arc_Ar_air.txt', 'arc_Cd_air.txt', $
      'arc_Zn_air.txt', 'arc_Kr_air.txt', 'arc_Xe_air.txt', 'arc_FeNe_air.txt' ]

    ; read lamp status from header
    lamp_status_string = strarr(n_elements(lamp_name))
    for i=0, n_elements(lamp_name)-1 do lamp_status_string[i] = strtrim(sxpar(arc_hdr, lamp_name[i]), 2)

    ; from string of 'on', 'off' to integer
    lamp_status = intarr(n_elements(lamp_status_string))
    lamp_status[where(lamp_status_string eq 'on', /null)] = 1

  endelse

  print, ''
  print, 'Lamps used for the arc calibrations:'
  forprint, lamp_name, lamp_status, format='(A10,A7)'

  ; check that some are on
  w_on = where(lamp_status eq 1, /null)
  if w_on eq !NULL then message, 'arc frame was taken with no arc lamp on!'

  ; load line lists for the lamps that were used
  all_lines_air = []
  for i=0, n_elements(w_on)-1 do begin
    print, 'Loading line list ' + fuel.util.flame_data_dir + lamp_file[w_on[i]]
    readcol, fuel.util.flame_data_dir + lamp_file[w_on[i]], arc_linelist
    all_lines_air = [all_lines_air, arc_linelist]
  endfor

	; sort them by wavelength
	all_lines_air = all_lines_air[sort(all_lines_air)]

	; convert to vacuum
	airtovac, all_lines_air, all_lines_vac

	; convert to micron
	all_lines = all_lines_vac*1d-4

	; write out the linelist
  forprint, all_lines, replicate(1, n_elements(all_lines)), $
	 	textout=fuel.util.intermediate_dir + 'linelist_arcs.txt', comment='#  arc_lines  trust'

  print, ''
  print, 'Arc line list written to ', fuel.util.intermediate_dir + 'linelist_arcs.txt'


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_lris, input
  ;
  ; LRIS-specific routine that initializes the fuel structure
  ;

  print, ''
  print, 'Initializing LRIS data reduction'
  print, ''

  ; ---------------------   create fuel structure   --------------------------------------


  fuel = flame_util_create_fuel(input)


  ; ---------------------   INSTRUMENT structure   --------------------------------------

  ; read FITS header of first science frame
  science_header = headfits(fuel.util.science.raw_files[0])

  ; read instrument name - useful to discriminate between blue and red channel
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)

  ; initialize the instrument structure
  case instrument_name of

    'LRISBLUE': instrument = flame_initialize_lris_blue(science_header)

    'LRIS': instrument = flame_initialize_lris_red(science_header)

    else: message, 'instrument keyword: ' + instrument_name + $
      ' ; cannot determine if it is LRIS RED or LRIS BLUE'

  endcase


  ; ---------------------   SETTINGS   --------------------------------------

  ; need to use the optical sky spectrum
  fuel.settings.sky_emission_filename = $
    fuel.util.flame_data_dir + 'sky_emission_model_optical.dat'

  ; use the sky background to trace the slit edges
  fuel.settings.trace_slit_with_emlines = 0

  ; how much margin to allow for detecting slit vertical distortion
  fuel.settings.trace_slit_ymargin = 30

  ; identify cosmic rays using L.A.Cosmic in each science frame
  fuel.settings.clean_individual_frames = 0

  ; split the spectrum into two when doing the rough wavecal
  fuel.settings.roughwavecal_split = 0


  ; ---------------------   ARC LINE LIST    -----------------------------------

  if fuel.util.arc.n_frames gt 0 then flame_initialize_lris_arcs, fuel


  ; -----------------------------------------------------------------
  ;        convert the LRIS frames into the "normal" format
  ; -----------------------------------------------------------------

  ; the conversion also applies the gain correction to each amplifier
  gaindata = instrument.gaindata

  ; create directory where we will store the new FITS files
  lris_dir = fuel.util.intermediate_dir + 'data_LRIS' + path_sep()
  if file_test(lris_dir) then file_delete, lris_dir, /recursive
  file_mkdir, lris_dir

  ; store and update file names for the science frames
  old_filenames = fuel.util.science.raw_files
  new_filenames = lris_dir + file_basename(fuel.util.science.raw_files)
  fuel.util.science.raw_files = lris_dir + file_basename(fuel.util.science.raw_files)

  ; same thing for dark frames, if provided
  if fuel.util.dark.n_frames gt 0 then begin
    old_filenames = [old_filenames, fuel.util.dark.raw_files]
    new_filenames = [new_filenames, lris_dir + file_basename(fuel.util.dark.raw_files)]
    fuel.util.dark.raw_files =  lris_dir + file_basename(fuel.util.dark.raw_files)
  endif

  ; same thing for arc frames, if provided
  if fuel.util.arc.n_frames gt 0 then begin
    old_filenames = [old_filenames, fuel.util.arc.raw_files]
    new_filenames = [new_filenames, lris_dir + file_basename(fuel.util.arc.raw_files)]
    fuel.util.arc.raw_files =  lris_dir + file_basename(fuel.util.arc.raw_files)
  endif

  ; same thing for pixelflat frames, if provided
  if fuel.util.pixelflat.n_frames gt 0 then begin
    old_filenames = [old_filenames, fuel.util.pixelflat.raw_files]
    new_filenames = [new_filenames, lris_dir + file_basename(fuel.util.pixelflat.raw_files)]
    fuel.util.pixelflat.raw_files =  lris_dir + file_basename(fuel.util.pixelflat.raw_files)
  endif

  ; same thing for illumflat frames, if provided
  if fuel.util.illumflat.n_frames gt 0 then begin
    old_filenames = [old_filenames, fuel.util.illumflat.raw_files]
    new_filenames = [new_filenames, lris_dir + file_basename(fuel.util.illumflat.raw_files)]
    fuel.util.illumflat.raw_files =  lris_dir + file_basename(fuel.util.illumflat.raw_files)
  endif

  ; same thing for slitflat frames, if provided
  if fuel.util.slitflat.n_frames gt 0 then begin
    old_filenames = [old_filenames, fuel.util.slitflat.raw_files]
    new_filenames = [new_filenames, lris_dir + file_basename(fuel.util.slitflat.raw_files)]
    fuel.util.slitflat.raw_files =  lris_dir + file_basename(fuel.util.slitflat.raw_files)
  endif


  ; now apply conversion to all frames
  print, ''
  print, 'Converting LRIS frames to standard format'
  for i=0, n_elements(old_filenames)-1 do begin

    ; use the readmhdufits.pro routine provided by Keck
    image = readmhdufits(old_filenames[i], header=header, gaindata=gaindata)

    ; make wavelength axis horizontal
    image = transpose(image)

    ; for the blue channel, the exptime is missing; use the elaptime instead
    if fxpar(header, 'EXPTIME', missing=-1.0) eq -1.0 then $
      fxaddpar, header, 'EXPTIME', fxpar(header, 'ELAPTIME', missing=1.0), 'added by flame; see ELAPTIME'

    ; write the new FITS file
    writefits, new_filenames[i], image, header
    print, new_filenames[i]

  endfor


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
  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:instrument, diagnostics:fuel.diagnostics, slits:slits }
  fuel=new_fuel

  ; return the fuel structure
  return, new_fuel


END
