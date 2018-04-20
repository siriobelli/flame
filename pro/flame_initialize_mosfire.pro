
;
; Initialization module for Keck MOSFIRE
;
; All tabulated values are taken from McLean et al. 2012 (SPIE, 8446, 0J)
;


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_mosfire_resolution, instrument
;
; Return the estimated spectral resolution  assuming a 1" slit width.
; For MOSFIRE, the spectral resolution is roughly fixed, with a slight
; dependence on wavelength
;

; get the band in which the observations were taken
filter = instrument.filter

; look up the resolution for this band, assuming a 0.7" slit width -------------
  case filter of
    'K': R_07 = 3620.
    'H': R_07 = 3660.
    'J': R_07 = 3310.
    'Y': R_07 = 3380.
    else: message, filter + ': filter not supported'
  endcase

  ; convert resolution to a standard slit width of 1"
  R_1arcsec = R_07 * 0.7/1.0

  ; return the spectral resolution
  return, R_1arcsec

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_mosfire_waverange, instrument, slit_xlocation
;
; Return the estimated wavelength range for a slit (in micron).
; The range depends on the filter and on the horizontal location
; of the slit on the CSU, in mm (slit_xlocation).
;

  ; get the band in which the observations were taken
  filter = instrument.filter

  ; 1 - look up the central wavelength (in um) for a centered slit -------------
  case filter of
    'K': central_wavelength = 2.1760
    'H': central_wavelength = 1.6321
    'J': central_wavelength = 1.2450
    'Y': central_wavelength = 1.0373
    else: message, filter + ': filter not supported'
  endcase

  ; 2 - look up the dispersion, in A/pixel -------------------------------------
  case filter of
    'K': dispersion = 2.1691
    'H': dispersion = 1.6269
    'J': dispersion = 1.3028
    'Y': dispersion = 1.0855
    else: message, filter + ': filter not supported'
  endcase

  ; 3 - look up the geometric shift, in A/mm -----------------------------------
  case filter of
    'K': shift = 12.209
    'H': shift = 9.157
    'J': shift = 7.460
    'Y': shift = 6.216
    else: message, filter + ': filter not supported'
  endcase


  ; 4 - rough wavelength range given the known geometry of the instrument ------
  lambda_center = central_wavelength + shift*1d-4 * slit_xlocation
  lambda_min = lambda_center - 1024.0*dispersion*1d-4
  lambda_max = lambda_center + 1024.0*dispersion*1d-4

  ; return the wavelength range for this slit
  return, [lambda_min, lambda_max]

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_mosfire_instrument, science_frame
;
; read the MOSFIRE settings from the FITS header of a science frame
; and create the instrument structure
;

  ; read the FITS header of the science frame
  science_header = headfits(science_frame)

  ; read things from FITS header -----------------------------------------------

  ; read instrument name
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)

  ; read filter
  filter = strtrim(fxpar(science_header, 'FILTER'), 2)

  ; read filter and grating
  grating = strtrim(fxpar(science_header, 'OBSMODE'), 2)

  ; extract the band from the grating mode
  grating_band = (strsplit(grating, '-', /extract))[0]

  ; check that the filter and the grating mode match
  if strlowcase(grating_band) NE strlowcase(filter) then $
    message, 'The filter and the grating mode do not match! Filter: ' + filter + ', grating mode: ' + grating

  ; read gain
  gain = fxpar(science_header, 'SYSGAIN')   ; e-/adu

  ; sampling mode
  sampling_mode_number = fxpar(science_header, 'SAMPMODE') ; (1:Single, 2:CDS, 3:MCDS, 4:UTR)
  sampling_mode = (['Single', 'CDS', 'MCDS', 'UTR'])[sampling_mode_number-1]

  ; number of reads
  sampling_reads = fxpar(science_header, 'NUMREADS')


  ; calculate things from hard-coded numbers -----------------------------------

  ; pixel scale in arcsec / pixel
  pixel_scale = 0.1798

  ; linearity correction: the polynomial coefficients describing the transformation
  linearity_correction = [0.0d, 1.0d]   ; do not have linearity correction available

  ; read-out noise (from the MOSFIRE webpage), in units of e-/read
  if sampling_mode eq 'CDS' then readnoise = 21.0
  if sampling_mode eq 'MCDS' then case sampling_reads of
    4: readnoise = 10.8
    8: readnoise = 7.7
    16: readnoise = 5.8
    32: readnoise = 4.2
    64: readnoise = 3.5
    128: readnoise = 3.0
    else: message, 'MCDS-' + strtrim(sampling_reads, 2) + ' not supported!'
  endcase

  if sampling_mode ne 'CDS' and sampling_mode ne 'MCDS' then $
    message, sampling_mode + ': sampling mode not supported!'


  ; calibration files for when the user doesn't have them
  default_badpixel_mask = 'none'
  default_dark = 'none'
  default_pixelflat = 'none'
  default_illumflat = 'none'
  default_arc = 'none'


  ; create the instrument structure --------------------------------------------
  instrument = { $
    instrument_name: instrument_name, $
    filter: filter, $
    grating: grating, $
    pixel_scale: pixel_scale, $
    sampling_mode: sampling_mode, $
    sampling_reads: sampling_reads, $
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
  instrument.resolution_slit1arcsec = flame_initialize_mosfire_resolution(instrument)

  return, instrument


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************




FUNCTION flame_initialize_mosfire_slits, science_frame, instrument, input
;
; read the header of a science frame
; and for each slit finds or calculate the slit number, name, slit PA,
; bottom pixel position, top pixel position, target pixel position,
; and wavelength at the center of the slit. Return the array of slit structures.
;

  ; read in the two HDUs with the relevant information
  hdu_slits = mrdfits(science_frame, 2, /silent) ; information on the science slits
  hdu_bars = mrdfits(science_frame, 3, /silent)  ; information on the individual bars of the CSU

  ; among the science slits there are empty structures corresponding to
  ; the alignment boxes; remove them
  hdu_slits = hdu_slits[where(hdu_slits.slit_number NE '', /null)]

  ; create array of slit structures
  slits = []

  ; loop through the science slits and for each of them make a slit structure
  for i_slit=0, n_elements(hdu_slits)-1 do begin

    ; name of this target
    target_name = hdu_slits[i_slit].target_name

    ; select all the CSU bars corresponding to this target
    w_thistarget = where(hdu_bars.target_in_slit eq target_name, /null)
    if w_thistarget EQ !NULL then $
      message, target_name + ': could not identify any slit for this target!'

    ; x position of the slit on the CSU, in mm (take the average of the bar positions)
    slit_xlocation = mean(float(hdu_bars[w_thistarget].position_of_slit), /nan)

    ; calculate the approximate wavelength range
    waverange = flame_initialize_mosfire_waverange( instrument, slit_xlocation )

    ; range in lambda0 (wavelength of first pixel) to be realistically considered
    lambda_width = waverange[1] - waverange[0]
    range_lambda0 = waverange[0] + [-0.2*lambda_width, 0.2*lambda_width]

    ; calculate pixel scale and its possible variation
    pixel_scale = lambda_width/2048.0
    range_delta_lambda = pixel_scale*[0.7,1.3]

    ; indices of the first and last bars for this slit (first bar is at the top)
    index_first_bar = min(fix(hdu_bars[w_thistarget].slit_number))
    index_last_bar = max(fix(hdu_bars[w_thistarget].slit_number))

    ; simple geometric model of the CSU, estimated from real data
    space_top = 15.5      ; empty space at the top of the detector, in pixels
    bar_height = 44.28    ; height of a bar in pixel
    space_between_slits = 4.0    ; typical empty space between consecutive slits, in pixels

    ; approximate vertical position of the slit on the detector
    slit_top = 2048.0 - space_top - bar_height * (index_first_bar-1.0) - 0.5*space_between_slits
    slit_bottom = 2048.0 - space_top -bar_height * index_last_bar + 0.5*space_between_slits

    ; estimate the expected position of the target
    relative_delta_target = float(hdu_slits[i_slit].target_to_center_of_slit_distance) / float(hdu_slits[i_slit].slit_length)
    target_position = 0.5*(slit_bottom+slit_top) + relative_delta_target * (slit_top-slit_bottom)

    ; the bottom part of the detector is typically not illuminated
    if slit_bottom LT 20.0 then slit_bottom = 20.0

    ; create one slit structure
    this_slit = flame_util_create_slitstructure( $
      number = fix(hdu_slits[i_slit].slit_number), $
      name = target_name, $
      PA = 3.78, $
      approx_bottom = slit_bottom, $
      approx_top = slit_top, $
      approx_target = target_position, $
      width_arcsec = float(hdu_slits[i_slit].slit_width), $
      approx_R = instrument.resolution_slit1arcsec / float(hdu_slits[i_slit].slit_width), $
      range_lambda0 = range_lambda0, $
      range_delta_lambda = range_delta_lambda )

      ; stack this slit structure with the other slits
      slits = [slits, this_slit]

    endfor


  return, slits


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************




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

  ; get a representative science frame (the first one)
  science_frame = fuel.util.science.raw_files[0]

  ; read the instrument settings from the header
  instrument = flame_initialize_mosfire_instrument(science_frame)

  ; now read in the slits parameters from the FITS header
  slits = flame_initialize_mosfire_slits( science_frame, instrument, fuel.input)


  ; ---------------------   SETTINGS   --------------------------------------

  ; no need to run L.A.Cosmic on individual frames
  fuel.settings.clean_individual_frames = 0

  ; use the OH sky emission lines to trace the slit edges
  fuel.settings.trace_slit_with_emlines = 1

  ; apply illumination correction using the OH lines
  fuel.settings.illumination_correction = 1

  ; do not split the spectrum into two when doing the rough wavecal
  fuel.settings.roughwavecal_split = 0

  ; usually need a higher resolution for the bspline model of the sky
  fuel.settings.skysub_bspline_oversample = 4.0

  ; Use the R3000 line list
  linelist = fuel.util.flame_data_dir + 'linelist_sky_R3000.dat'

  ; make a local copy of the line list
  file_copy, linelist, fuel.util.intermediate_dir, /overwrite

  ; save the file name in the settings
  fuel.settings.linelist_sky_filename = fuel.util.intermediate_dir + file_basename(linelist)



  ; -------------------------------------------------------------------------------


  ; save both instrument and slits in the fuel structure
  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:instrument, diagnostics:fuel.diagnostics, slits:slits }

  ; return the fuel structure
  return, new_fuel

END
