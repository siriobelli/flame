
;
; This file is a template that can be used to develop new
; initialization modules for instruments that are currently not supported
; The name "template" has been used throughout to indicate the unknown instrument
; Please replace all occurences of "template" with the instrument name (including
; in the name of this file)
;

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_template_resolution, instrument
;
; Return the estimated spectral resolution given the instrument settings contained
; in the "instrument" structure, and assuming a 1" slit width.
; Typically, the resolution depends on camera, grating, filter, etc.
;

; for example:

; get the band in which the observations were taken
filter = instrument.filter

; look up the resolution for this band, assuming a 1" slit width ---------------
  case filter of
    'K': R = 3700.
    'H': R = 3600.
    'J': R = 3500.
    'Y': R = 3400.
    else: message, filter + ': filter not supported'
  endcase

  ; return the spectral resolution
  return, R


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_template_waverange, instrument, slit_xlocation
;
; Return the estimated wavelength range for a slit (in micron).
; It usually depends on the instrument settings but also the horizontal location
; of the slit on the mask.
; If the information on the slit location is not available, then assume the slit
; is centered, and use a wider range for the acceptable values
;

; for example:

  ; get the band in which the observations were taken
  filter = instrument.filter

  ; 1 - look up the central wavelength (in um) for a centered slit -------------
  case filter of
    'K': central_wavelength = 2.2
    'H': central_wavelength = 1.6
    'J': central_wavelength = 1.2
    'Y': central_wavelength = 1.0
    else: message, filter + ': filter not supported'
  endcase

  ; 2 - look up the dispersion, in A/pixel -------------------------------------
  case filter of
    'K': dispersion = 1.5
    'H': dispersion = 1.5
    'J': dispersion = 1.5
    'Y': dispersion = 1.5
    else: message, filter + ': filter not supported'
  endcase

  ; 3 - look up the geometric shift, in A/mm -----------------------------------
  case filter of
    'K': shift = 7.0
    'H': shift = 7.0
    'J': shift = 7.0
    'Y': shift = 7.0
    else: message, filter + ': filter not supported'
  endcase


  ; 4 - rough wavelength range given the known geometry of the instrument ------
  lambda_center = central_wavelength - shift * slit_xlocation
  lambda_min = lambda_center - 1024.0*dispersion
  lambda_max = lambda_center + 1024.0*dispersion

  ; return the wavelength range for this slit
  return, [lambda_min, lambda_max]

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_template_instrument, science_header
;
; Read the instrument settings from the FITS header of a science frame
; and create the instrument structure
;

  ; for example:

  ; read things from header - - - - - - - - - - - - - - - - - - -

  ; read instrument name, just to be sure
  ; in some cases this can be useful to distinguish
  ; between two arms or two components
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
  ; (these should all stored in the flame data directory, fuel.util.flame_data_dir)
  default_badpixel_mask = 'TEMPLATE/default_badpixel.fits'
  default_dark = 'none'
  default_pixelflat = 'TEMPLATE/default_pixelflat.fits'
  default_illumflat = 'none'
  default_arc = 'none'

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
    trim_edges: 4, $
    default_badpixel_mask: default_badpixel_mask, $
    default_dark: default_dark, $
    default_pixelflat: default_pixelflat, $
    default_illumflat: default_illumflat, $
    default_arc: default_arc $
    }

  ; now use the instrument structure to calculate the spectral resolution
  instrument.resolution_slit1arcsec = flame_initialize_template_resolution(instrument)

  ; return the instrument structure
  return, instrument

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_initialize_template_arcs, fuel
  ;
  ; Handle the arcs calibrations: identify the lamps used, from the FITS header
  ; and then produce the appropriate line list and save it in the intermediate directory
  ;

  ; 1 - identify what lamps were used

  ; read first arc frame
  arc_hdr = headfits(fuel.util.arc.raw_files[0])

  ; here you need to know exactly how the instrument records the lamps in the
  ; FITS header. This is just an example
  status = strlowcase(strtrim(sxpar(arc_hdr, 'ARCLAMP'), 2))
  lamp_name = ['Neon', 'Argon']
  lamp_status = [0, 0]
  if status eq 1 then lamp_status[0] = 1
  if status eq 2 then lamp_status[1] = 1

  print, ''
  print, 'Lamps used for the arc calibrations:'
  forprint, lamp_name, lamp_status, format='(A10,A7)'

  ; check that some are on
  w_on = where(lamp_status eq 1, /null)
  if w_on eq !NULL then message, 'arc frame was taken with no arc lamp on!'


  ; 2 - produce the appropriate line list and save it in the intermediate directory

  ; filenames with the line lists
  lamp_file = 'linelist_arc_' + ['Ne', 'Ar', 'Xe'] + '.dat'

  ; load line lists for the lamps that were used
  all_lines = []
  for i=0, n_elements(w_on)-1 do begin
    print, 'Loading line list ' + fuel.util.flame_data_dir + lamp_file[w_on[i]]
    readcol, fuel.util.flame_data_dir + lamp_file[w_on[i]], arc_linelist
    all_lines = [all_lines, arc_linelist]
  endfor

	; sort them by wavelength
	all_lines = all_lines[sort(all_lines)]

  ; set the name of local linelist file
  fuel.settings.linelist_arcs_filename = fuel.util.intermediate_dir + 'linelist_arcs.dat'

	; write out the linelist
  forprint, all_lines, replicate(1, n_elements(all_lines)), $
	 	textout=fuel.settings.linelist_arcs_filename, comment='#  arc_lines  trust'

  print, ''
  print, 'Arc line list written to ', fuel.settings.linelist_arcs_filename


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



FUNCTION flame_initialize_template_longslit, header, instrument=instrument, input=input
  ;
  ; Create and return one slit structure for the long slit
  ;

  ; rough wavelength range (slit should be central: x=0.0)
  lambda_range = $
   flame_initialize_template_waverange(instrument, 0.0)

  ; range in lambda0 to be realistically considered
  lambda_width = lambda_range[1] - lambda_range[0]
  range_lambda0 = lambda_range[0] + [-0.3*lambda_width, 0.3*lambda_width]

  ; calculate pixel scale and its possible variation
  pixel_scale = (lambda_range[1]-lambda_range[0])/2048.0
  range_delta_lambda = pixel_scale*[0.5,1.5]

  ; vertical range to be considered. Default is to cut 10% of pixels on each side
  if array_equal( input.longslit_edge, [0,0]) then $
    yrange = [205, 1843] else $   ; for a 2048x2048 detector
    yrange = input.longslit_edge

  ; create the slit structure
  longslit = flame_util_create_slitstructure( $
    number = 1, $
    name = 'longslit', $
    PA = !values.d_nan, $
    approx_bottom = yrange[0], $
    approx_top = yrange[1], $
    approx_target = mean(yrange), $
    width_arcsec = !values.d_nan, $
    approx_R = instrument.resolution_slit1arcsec, $
    range_lambda0 = range_lambda0, $
    range_delta_lambda = range_delta_lambda )

  return, longslit

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



FUNCTION flame_initialize_template_slits, header, instrument=instrument, input=input

;
; read the header of a science frame
; and for each slit finds or calculate the slit number, name, slit PA,
; bottom pixel position, top pixel position, target pixel position,
; and wavelength at the center of the slit. Return the array of slit structures.
;

  ; array of structures that will contain the info for each slit
  slit_hdr = []

  ; need to know how the information on the slitmask is stored in the FITS header,
  ; and extract the relevant information:

  ; here is a naive example:
  slit_name = strarr(3)
  slit_xlocation = fltarr(3)
  for i_slit=0,2 do begin
    slit_name[i_slit] = strtrim( fxpar( header, 'TGT' + i_slit + 'NAM' ), 2)
    slit_xlocation[i_slit] = float(fxpar(header, 'MOS' + i_slit + 'XPO'))
    ; slit_bottom[i_slit] = ... and so on
  endfor

  ;
  ; if needed, identify and exclude the alignment boxes here
  ;

  ; create array of slit structures
  slits = []

  for i_slit=0, n_elements(bottom)-1 do begin

    ; calculate approximate wavelength range
    lambda_range = flame_initialize_template_waverange(instrument, slit_xlocation[i_slit])

    ; range in lambda0 (wavelength of first pixel) to be realistically considered
    lambda_width = lambda_range[1] - lambda_range[0]
    range_lambda0 = lambda_range[0] + [-0.3*lambda_width, 0.3*lambda_width]

    ; calculate pixel scale and its possible variation
    pixel_scale = (lambda_range[1]-lambda_range[0])/2048.0
    range_delta_lambda = pixel_scale*[0.5,1.5]

    ; create one slit structure
    this_slit = flame_util_create_slitstructure( $
      number = i_slit, $
      name = slit_name[i_slit], $
      PA = 0.0, $
      approx_bottom = slit_bottom[i_slit], $
      approx_top = slit_top[i_slit], $
      approx_target = slit_target[i_slit], $
      width_arcsec = slit_width[i_slit], $
      approx_R = instrument.resolution_slit1arcsec / slit_width[i_slit], $
      range_lambda0 = range_lambda0, $
      range_delta_lambda = range_delta_lambda )

    ; stack this slit structure with the other slits
    slits = [slits, this_slit]


  endfor


  ; if manual slit positions are provided, then use them
  ;---------------------------------------------------------
  if input.slit_position_file ne 'none' then begin

    print, 'using slit position file ', input.slit_position_file

    ; read in region file
    readcol, input.slit_position_file, slit_x, slit_y

    ; sort slit_y
    slit_y = slit_y[ sort(slit_y) ]

    ; check that the number of slit edges is correct
    if n_elements(slit_y) NE 2*n_elements(slits) then $
      message, input.slit_position_file + ' must contain exactly ' + $
      strtrim(2*n_elements(slits), 2) + ' slit edges'

    ; now sort the slits from the header by target position
    s = sort(slits.approx_target)

    for i_slit=0, n_elements(slits)-1 do begin
      slits[s[i_slit]].approx_bottom = slit_y[2*i_slit]
      slits[s[i_slit]].approx_top = slit_y[2*i_slit+1]
      slits[s[i_slit]].approx_target = 0.5*(slit_y[2*i_slit] + slit_y[2*i_slit+1])
    endfor

  endif

return, slits

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_template, input


  print, ''
  print, 'Initializing TEMPLATE data reduction'
  print, ''


  ; first, create the fuel structure
  fuel = flame_util_create_fuel(input)


  ; ---------------------   INSTRUMENT structure   -----------------------------

  ; read FITS header of first science frame
  science_header = headfits(fuel.util.science.raw_files[0])

  ; read the instrument settings from the header
  instrument = flame_initialize_template_instrument(science_header)


  ; now read in the slits parameters from the FITS header

  ; LONGSLIT -------------------------------------------------------------------
  if fuel.input.longslit then begin

    ; get the slit edges from the user
    slits = flame_initialize_template_longslit( science_header, instrument=instrument, input=fuel.input)


  ; MOS      -------------------------------------------------------------------
  endif else begin

    ; get all the info from the header
    slits = flame_initialize_template_slits( science_header, instrument=instrument, input=fuel.input)

  endelse


  ; ------------------------   SETTINGS   --------------------------------------
  ; here the default settings are chosen; the user can change any of these values

  ; choose the correct linelist
  linelist = fuel.util.flame_data_dir + 'sky_line_list_R3000.dat'

  ; make a local copy of the line list
  file_copy, linelist, fuel.util.intermediate_dir, /overwrite

  ; save the file name in the settings
  fuel.settings.linelist_sky_filename = fuel.util.intermediate_dir + file_basename(linelist)

  ; is it needed to run L.A.Cosmic on individual frames?
  fuel.settings.clean_individual_frames = 0

  ; are the slits traced with emission lines rather than the continuum?
  fuel.settings.trace_slit_with_emlines = 1

  ; apply illumination correction using the OH lines?
  fuel.settings.illumination_correction = 1

  ; split the spectrum into two when doing the rough wavecal?
  fuel.settings.roughwavecal_split = 0

  ; set the degree of the polynomials for the 2D wavelength solution
  fuel.settings.wavesolution_order_x = 3
  fuel.settings.wavesolution_order_y = 2


  ; ---------------------------------   ARCS    --------------------------------

  if fuel.util.arc.n_frames gt 0 then flame_initialize_template_arcs, fuel


  ; ----------------------------------------------------------------------------
  ; do not change this part

  ; save both instrument and slits in the fuel structure
  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:instrument, diagnostics:fuel.diagnostics, slits:slits }

  ; return the fuel structure
  return, new_fuel

END
