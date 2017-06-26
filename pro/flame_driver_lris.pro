
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

  ; vertical offset, in pixel, between the slit-flat field and the science frames
  input.slitflat_offset = 0



  ;****************************************************
  ;****************************************************
  ;              PART 2: INITIALIZATION
  ;****************************************************
  ;****************************************************


  ; initialize and create the fuel structure
  fuel = flame_initialize_lris(input)


  ; set this to zero if you want to use the sky background to trace the slit edges
  ; instead of the OH lines (e.g. in the K band or in the optical or with slit flats)
  fuel.util.trace_slit_with_skylines = 0

  ; identify cosmic rays using L.A.Cosmic in each science frame
  fuel.util.clean_individual_frames = 1



  ;****************************************************
  ;****************************************************
  ;               PART 3: DATA REDUCTION
  ;****************************************************
  ;****************************************************


  flame_diagnostics, fuel

  flame_quickstack, fuel

  flame_correct, fuel

  flame_getslits, fuel

  flame_cutout_slits, fuel

  flame_wavecal_rough, fuel

  flame_identify_lines, fuel

  flame_wavecal_accurate, fuel

  flame_skysub, fuel

  flame_rectify, fuel

  flame_combine, fuel
  
  flame_checkdata, fuel
