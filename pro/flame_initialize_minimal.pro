;
; Minimal initialization for a generic instrument - its only purpose is to
; run the first two steps of flame (flame_diagnostics and flame_quickstack)
; for on-the-fly diagnostics
;
; The only parameter that needs to be specified is the spatial pixel scale,
; in units of arcsec per pixel
;


FUNCTION flame_initialize_minimal, input, pixel_scale=pixel_scale

  print, ''
  print, 'Initializing minimal data reduction'
  print, ''


  ; create the fuel structure
  fuel = flame_util_create_fuel(input)

  ; check that pixel scale is set
  if ~keyword_set(pixel_scale) then begin
    print, 'WARNING: pixel scale has not been specified; assuming 0.25 arcsec/pixel'
    pixel_scale=0.25
  endif


  ; ---------------------  make INSTRUMENT structure   --------------------------------------


    ; read FITS header of first science frame
    science_header = headfits(fuel.util.science.raw_files[0])

    ; read instrument name and detector parameters
    instrument_name = strtrim(fxpar(science_header, 'INSTRUME'), 2)
    gain = strtrim(fxpar(science_header, 'GAIN'), 2)
    rdnoise = strtrim(fxpar(science_header, 'RDNOISE'), 2)

    ; create the instrument structure
    instrument = { $
      instrument_name: instrument_name, $
      pixel_scale: pixel_scale, $
      gain: gain, $
      readnoise: rdnoise, $
      resolution_slit1arcsec: 3000.0, $
      linearity_correction: [0.0d, 1.0d], $   ; no linearity correction
      trim_edges: 4, $
      default_badpixel_mask: 'none', $
      default_dark: 'none', $
      default_pixelflat: 'none', $
      default_illumflat: 'none', $
      default_arc: 'none' $
      }


  ; -------------------------------------------------------------------------------


  ; save both instrument and slits in the fuel structure
  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:instrument, diagnostics:fuel.diagnostics, slits:ptr_new() }

  ; return the fuel structure
  return, new_fuel

END
