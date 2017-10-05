
  ;****************************************************
  ;****************************************************
  ;              PART 1: INPUT
  ;****************************************************
  ;****************************************************

  ; suppress harmless messages (e.g., module compilation)
  !quiet = 1

  ; create the input structure
  input = flame_create_input()

  ; text file containing the list of science FITS files that need to be reduced
  input.science_filelist = 'science.txt'

  ; do you want to apply A-B sky subtraction?
  input.AB_subtraction = 0

  ; array with y-pixel positions for the traces of the reference star. 0 if there is no reference star
  input.star_y_A = 1281
  input.star_y_B = 0

  ; if 0, then reduce all slits. If n, then reduce slit number n (starting from 1).
  input.reduce_only_oneslit = 0

  ; for longslit, set this to 1
  input.longslit = 0

  ; and specify the y-range in pixel that you want to reduce
  input.longslit_edge = [0, 0]

  ; manual slit positions
  input.slit_position_file = 'slit_edges.reg'

  ; text file containing the list of FITS files for slit-flat field
  input.slitflat_filelist = 'slitflat.txt'
  input.slitflat_pixelshift = 0



  ;****************************************************
  ;****************************************************
  ;              PART 2: INITIALIZATION
  ;****************************************************
  ;****************************************************


  ; initialize and create the fuel structure
  fuel = flame_initialize_lris(input)


  ;****************************************************
  ;****************************************************
  ;               PART 3: DATA REDUCTION
  ;****************************************************
  ;****************************************************


  flame_diagnostics, fuel

  flame_quickstack, fuel

  flame_calibrations, fuel

  flame_getslits, fuel

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
