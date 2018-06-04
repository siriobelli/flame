;
; Driver file for reducing data from LBT/LUCI using the
; Flame pipeline (https://github.com/siriobelli/flame)
;
; Copy this file in a new directory, edit the input and settings,
; and then run it in IDL. You can do this either by copying each line
; of code on an interactive IDL session, or by running in IDL:
; .run flame_driver_luci
;

  ;****************************************************
  ;****************************************************
  ;                       INPUT
  ;****************************************************
  ;****************************************************

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

  ; if 0, then reduce all slits. If n, then reduce slit number n (starting from 1).
  input.reduce_only_oneslit = 0

  ; for longslit, set this to 1
  input.longslit = 0

  ; and specify the y-range in pixel that you want to reduce
  input.longslit_edge = [0, 0]

  ; see the Flame manual for additional input fields
  ; you can also get a complete list with:
  ; help, input


  ;****************************************************
  ;****************************************************
  ;                 INITIALIZATION
  ;****************************************************
  ;****************************************************


  ; initialize and create the fuel structure
  fuel = flame_initialize_luci(input)

  ; here you can customize the settings
  ; for example, to change the range in x-coordinates
  ; used to extract the star traces:
  ; fuel.settings.star_x_range = [100, 500]
  ;
  ; see the Flame manual for details
  ; you can also get a complete list of settings with:
  ; help, fuel.settings


  ;****************************************************
  ;****************************************************
  ;                   DATA REDUCTION
  ;****************************************************
  ;****************************************************


  flame_diagnostics, fuel

  flame_quickstack, fuel

  flame_calibrations, fuel

  flame_slitid, fuel

  flame_cutouts, fuel

  flame_spatialcal, fuel

  flame_roughwavecal, fuel

  flame_findlines, fuel

  flame_wavecal, fuel

  flame_illumcorr, fuel

  flame_skysub, fuel

  flame_rectify, fuel

  flame_combine, fuel

  flame_checkdata, fuel

  flame_extract, fuel

END
