
;
; LUCI-specific routine that initializes the fuel structure
;

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_luci_resolution, instrument
;
; return the estimated spectral resolution *for a slit width of 1.0 arcsec*.
; the resolution depends only on the grating (and the slit witdth).
;

; get the first part of the grating name
grating = (strsplit(instrument.grating, /extract))[0]

; look up the resolution for this grating
; Note that the tabulated values for R assume a slit width of 0.5 arcsec

  ; grating G210
  if grating eq 'G210' then case instrument.grating_order of
    '2': R = 5000. ; K
    '3': R = 5900. ; H
    '4': R = 5800. ; J
    '5': R = 5400. ; z
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase

  ; grating G200
  if grating eq 'G200' then case instrument.grating_order of
    '1': R = 2250. ; HK
    '2': R = 2250. ; zJ
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase

  ; grating G150
  if grating eq 'G150' then case instrument.grating_order of
    '2': R = 4150. ; Ks
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase

  ; check whether no grating has been found
  if R EQ !NULL then message, 'grating ' + instrument.grating + ' not supported'

  ; convert to the resolution for a 1"-wide slit:
  R_1arcsec = R / 2.0

  return, R_1arcsec

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_luci_waverange, instrument, slit_xmm
;
; return the estimated wavelength range for a slit
; it depends on grating and camera
;

; get the first part of the grating name
grating = (strsplit(instrument.grating, /extract))[0]

; get the first part of the camera name
camera = (strsplit(instrument.camera, /extract))[0]


  ;1 - look up the wavelength range (in um) covered by the detector, assuming grating N1.8

  ; grating G210
  if grating eq 'G210' then case instrument.grating_order of
    '2': wavelength_range = 0.328 ; K
    '3': wavelength_range = 0.202 ; H
    '4': wavelength_range = 0.150 ; J
    '5': wavelength_range = 0.124 ; z
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase

  ; grating G200
  if grating eq 'G200' then case instrument.grating_order of
    '1': wavelength_range = 0.880 ; HK
    '2': wavelength_range = 0.440 ; zJ
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase

  ; grating G150
  if grating eq 'G150' then case instrument.grating_order of
    '2': wavelength_range = 0.533 ; Ks
    else: message, instrument.grating + ': order ' + instrument.grating_order + ' not supported'
  endcase

  ; check whether no grating has been found
  if wavelength_range EQ !NULL then message, 'grating ' + instrument.grating + ' not supported'


  ; 2 - correct for the camera (scale by the ratio of pixel scales)
  if camera NE 'N1.8' then $
    if camera eq 'N3.75' then wavelength_range *= 0.47 else $
      message, 'camera ' + instrument.camera + ' not supported'


  ; 3 - rough wavelength range; mask is roughly 2*162mm across; x=0mm means centered slit.
  lambda_min = instrument.central_wavelength - (162.0 + slit_xmm) / 162.0 * wavelength_range / 2.0
  lambda_max = lambda_min + wavelength_range


  return, [lambda_min, lambda_max]

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_initialize_luci_instrument, science_header
;
; read the LUCI settings from the FITS header of a science frame
; and create the instrument structure
;

  ; read things from header - - - - - - - - - - - - - - - - - - -

  ; read instrument name - useful to discriminate between LUCI1 and LUCI2
  instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)

  ; read grating name
  grating = strtrim(fxpar(science_header, 'GRATNAME'), 2)

  ; for old data, the grating names were slightly different
  grating_name = (strsplit(grating, /extract))[0]
  if (grating_name eq '150') then grating = 'G150 X'
  if (grating_name eq '200') then grating = 'G200 X'
  if (grating_name eq '210') then grating = 'G210 X'

  ; read in grating order
  grating_order = strtrim(fxpar(science_header, 'GRATORDE'), 2)

  ; read central wavelength
  central_wavelength = fxpar(science_header, 'GRATWLEN')

  ; read camera
  camera = fxpar(science_header, 'CAMERA', missing='NONE')
  if camera eq 'NONE' then camera = fxpar(science_header, 'CAMNAME') ; for old data

  ; read pixel scale
  pixel_scale = fxpar(science_header, 'PIXSCALE')  ; arcsec/pixel

  ; read filters
  filter1 = strtrim(fxpar(science_header, 'FILTER1'), 2)
  filter2 = strtrim(fxpar(science_header, 'FILTER2'), 2)

  ; read in read-out noise
  readnoise = fxpar(science_header, 'RDNOISE', missing=-1)   ; e-/read
  if readnoise eq -1 then readnoise = fxpar(science_header, 'ENOISE') ; for old data

  ; read in gain
  gain = fxpar(science_header, 'GAIN', missing=-1)   ; e-/adu
  if gain eq -1 then gain = fxpar(science_header, 'ELECGAIN') ; for old data

  ; calculate things from hard-coded numbers - - - - - - - - - - - - - - - - - - -

  ; linearity correction: the polynomial coefficients describing the transformation
  linearity_correction = [0.0d, 1.0d, 4.155d-6]

  ; calibration files for when the user doesn't have them - - - - - - - - - - - - - - - - - - -
  ; (these are all stored in the flame data directory, fuel.util.flame_data_dir)
  default_badpixel_mask = 'LUCI/default_badpixel_' + instrument_name + '.fits'
  default_dark = 'none'
  default_pixelflat = 'LUCI/default_pixelflat_' + instrument_name + '_allbands.fits'
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
  instrument.resolution_slit1arcsec = flame_initialize_luci_resolution(instrument)

  return, instrument


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_initialize_luci_arcs, fuel

  ; read first frame
  arc_hdr = headfits(fuel.util.arc.raw_files[0])

  lamp_name = ['NEON', 'ARGON', 'XENON']
  lamp_file = 'linelist_arc_' + ['Ne', 'Ar', 'Xe'] + '.dat'

  ; read lamp status from header
  status_Ne = strlowcase(strtrim(sxpar(arc_hdr, 'STATLMP1'), 2))
  status_Ar = strlowcase(strtrim(sxpar(arc_hdr, 'STATLMP2'), 2))
  status_Xe = strlowcase(strtrim(sxpar(arc_hdr, 'STATLMP3'), 2))

  ; put them together
  lamp_status_string = [status_Ne, status_Ar, status_Xe]

  ; from string of 'on', 'off' to integer
  lamp_status = intarr(n_elements(lamp_status_string))
  lamp_status[where(lamp_status_string eq 'on', /null)] = 1

  print, ''
  print, 'Lamps used for the arc calibrations:'
  forprint, lamp_name, lamp_status, format='(A10,A7)'

  ; check that some are on
  w_on = where(lamp_status eq 1, /null)
  if w_on eq !NULL then message, 'arc frame was taken with no arc lamp on!'

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



FUNCTION flame_initialize_luci_longslit, header, instrument=instrument, input=input

  ; rough wavelength range (slit should be central, x ~ 0 mm)
  lambda_range = $
   flame_initialize_luci_waverange(instrument, 0.0)

  ; range in lambda0 to be realistically considered
  lambda_width = lambda_range[1] - lambda_range[0]
  range_lambda0 = lambda_range[0] + [-0.3*lambda_width, 0.3*lambda_width]

  ; calculate pixel scale and its possible variation
  pixel_scale = (lambda_range[1]-lambda_range[0])/2048.0
  range_delta_lambda = pixel_scale*[0.5,1.5]

  ; vertical range to be considered. Default is to cut 10% of pixels on each side
  if array_equal( input.longslit_edge, [0,0]) then $
    yrange = [205, 1843] else $
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



FUNCTION flame_initialize_luci_slits, header, instrument=instrument, input=input

;
; read the header of a LUCI science frame
; and for each slit finds or calculate the slit number, name, slit PA,
; bottom pixel position, top pixel position, target pixel position,
; and wavelength at the center of the slit. It return the slits structure.
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

  ; identify and exclude the alignment boxes
  ;---------------------------------------------------------

  ; get the slit widths in arcsec
  widths = slit_hdr.width_arcsec

  ; sort them in decreasing order
  sorted_widths = widths[sort(-widths)]

  ; if they all are of the same width, then assume there are no alignment boxes
  if sorted_widths[0] eq sorted_widths[-1] then $
    print, 'All slits have the same width of ', sorted_widths[0], ' arcsec. No alignment boxes found.' $
  else begin

    ; if a maximum width for the science slits was specified, then use that to select the alignment boxes
    if input.max_slitwidth_arcsec NE 0.0 then begin

      print, n_elements(where(slit_hdr.width_arcsec GT input.max_slitwidth_arcsec, /null)), $
        ' alignment boxes ( wider than ', cgnumber_formatter(input.max_slitwidth_arcsec, decimals=2), ' arcsec) found.'

      ; exclude alignment boxes
      slit_hdr = slit_hdr[where( slit_hdr.width_arcsec LE input.max_slitwidth_arcsec, /null) ]

    ; otherwise, select all the slits that have the same width as the wider one
    endif else begin

      maxwidth = sorted_widths[0]

      print, n_elements(where(sorted_widths eq maxwidth, /null)), $
        ' alignment boxes (', cgnumber_formatter(maxwidth, decimals=2), ' arcsec) found.'

      ; exclude alignment boxes
      slit_hdr = slit_hdr[where( slit_hdr.width_arcsec LT maxwidth, /null) ]

    endelse

  endelse

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
  slit_width = slit_hdr.width_arcsec

  ; create array of slit structures
  slits = []

  for i_slit=0, n_elements(bottom)-1 do begin

    ; calculate approximate wavelength range
    lambda_range = flame_initialize_luci_waverange(instrument, slit_hdr[i_slit].x_mm)

    ; range in lambda0 (wavelength of first pixel) to be realistically considered
    lambda_width = lambda_range[1] - lambda_range[0]
    range_lambda0 = lambda_range[0] + [-0.3*lambda_width, 0.3*lambda_width]

    ; calculate pixel scale and its possible variation
    pixel_scale = (lambda_range[1]-lambda_range[0])/2048.0
    range_delta_lambda = pixel_scale*[0.5,1.5]

    ; create one slit structure
    this_slit = flame_util_create_slitstructure( $
      number = slit_num[i_slit], $
      name = slit_name[i_slit], $
      PA = slit_PA[i_slit], $
      approx_bottom = bottom[i_slit], $
      approx_top = top[i_slit], $
      approx_target = target[i_slit], $
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


FUNCTION flame_initialize_luci, input


  print, ''
  print, 'Initializing LUCI data reduction'
  print, ''


  ; first, create the fuel structure
  fuel = flame_util_create_fuel(input)


  ; ---------------------   INSTRUMENT structure   --------------------------------------

  ; read FITS header of first science frame
  science_header = headfits(fuel.util.science.raw_files[0])

  ; read the instrument settings from the header
  instrument = flame_initialize_luci_instrument(science_header)


  ; now read in the slits parameters from the FITS header

  ; LONGSLIT ---------------------------------------------------------------------
  if fuel.input.longslit then begin

    ; get the slit edges from the user
    slits = flame_initialize_luci_longslit( science_header, instrument=instrument, input=fuel.input)


  ; MOS      ---------------------------------------------------------------------
  endif else begin

    ; get all the info from the header
    slits = flame_initialize_luci_slits( science_header, instrument=instrument, input=fuel.input)

  endelse


  ; ---------------------   SETTINGS   --------------------------------------

  ; calculate the typical spectral resolution of the data
  spectral_res = median([slits.approx_r])

  ; pick the most appropriate sky line list
  if spectral_res LT 2000.0 then $
    linelist = fuel.util.flame_data_dir + 'linelist_sky_R1000.dat'
  if spectral_res GE 2000.0 and spectral_res LT 4500.0 then $
    linelist = fuel.util.flame_data_dir + 'linelist_sky_R3000.dat'
  if spectral_res GE 4500.0 then $
    linelist = fuel.util.flame_data_dir + 'linelist_sky_R6000.dat'

  ; make a local copy of the line list
  file_copy, linelist, fuel.util.intermediate_dir, /overwrite

  ; save the file name in the settings
  fuel.settings.linelist_sky_filename = fuel.util.intermediate_dir + file_basename(linelist)

  ; no need to run L.A.Cosmic on individual frames
  fuel.settings.clean_individual_frames = 0

  ; usually AB nodding is used, which means that the dark is automatically subtracted
  fuel.settings.darksub_data = 0

  ; use the OH sky emission lines to trace the slit edges
  fuel.settings.trace_slit_with_emlines = 1

  ; scale the roughwavecal_R to the spectral resolution - useful for argos
  fuel.settings.roughwavecal_R = [ 0.05*slits[0].approx_R > 500.0 , $
		0.15*slits[0].approx_R > 1000.0 , slits[0].approx_R < 5000.0 ]

  ; apply illumination correction using the OH lines
  fuel.settings.illumination_correction = 1

  ; do not split the spectrum into two when doing the rough wavecal
  fuel.settings.roughwavecal_split = 0

  ; set the degree of the polynomials for the 2D wavelength solution
  fuel.settings.wavesolution_order_x = 3
  fuel.settings.wavesolution_order_y = 2


  ; ---------------------   ARC LINE LIST    -----------------------------------

  if fuel.util.arc.n_frames gt 0 then flame_initialize_luci_arcs, fuel


  ; -------------------------------------------------------------------------------

  ; save both instrument and slits in the fuel structure
  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:instrument, diagnostics:fuel.diagnostics, slits:slits }

  ; return the fuel structure
  return, new_fuel

END
