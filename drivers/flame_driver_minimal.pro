;
; Driver file for monitoring data quality using the
; Flame pipeline (https://github.com/siriobelli/flame)
;
; This will work with nearly any instrument as long as the data are
; saved as FITS files with just one extension, and the wavelength axis
; runs along the horizontal direction. The only instrument parameter needed
; is the spatial pixel scale in arcsec per pixel:
pixel_scale = 0.18
;
; Copy this file in a new directory,
; and then run it in IDL. You can do this either by copying each line
; of code on an interactive IDL session, or by running in IDL:
; .run flame_driver_minimal
;

  ; suppress harmless messages (e.g., module compilation)
  !quiet = 1

  ; create the input structure
  input = flame_create_input()

  ; text file containing the list of science FITS files that need to be reduced
  input.science_filelist = 'science.txt'

  ; do you want to apply A-B sky subtraction?
  input.AB_subtraction = 1

  ; array with y-pixel positions for the traces of the reference star. 0 if there is no reference star
  input.star_y_A = 1281
  input.star_y_B = 1300

  ; initialize and create the fuel structure
  fuel = flame_initialize_minimal(input, pixel_scale=pixel_scale)

  ; data reduction: only the first two steps

  flame_diagnostics, fuel

  flame_quickstack, fuel

END
