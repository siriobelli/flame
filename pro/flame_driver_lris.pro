
  ;****************************************************
  ;****************************************************
  ;              PART 1: INITIALIZATION
  ;****************************************************
  ;****************************************************

  ; suppress harmless messages (e.g., module compilation)
  !quiet = 1

  ; create the input structure
  input = flame_create_input()


  ; BASIC INPUT
  ;**********************************

  ; text file containing the list of science FITS files that need to be reduced
  input.science_filelist = 'science.txt'

  ; text file containing the list of FITS files with dark frames
  input.dark_filelist = 'none'

  ; text file containing the list of FITS files for pixel-flat field
  input.pixelflat_filelist = 'none'

  ; text file containing the list of FITS files for illumination-flat field
  input.illumflat_filelist = 'none'

  ; text file containing the list of FITS files with arcs for wavelength calibration
  input.arc_filelist = 'none'

  ; do you want to apply A-B sky subtraction?
  input.AB_subtraction = 0

  ; array with y-pixel positions for the traces of the reference star. 0 if there is no reference star
  input.star_y_A = 1575
  input.star_y_B = 0


  ; ADVANCED OPTIONS
  ;**********************************

  ; if 0, then reduce all slits. If n, then reduce slit number n (starting from 1).
  ;input.reduce_only_oneslit = 5

  ; if you want to change the range in x-coordinates used to extract the star traces:
  ;input.star_x_range = [100, 500]

  ; if we don't have a star on the slit then we have to specify the dithering
  ;input.dither_filelist = 'input/dither.txt'

  ; set to zero if you want to use the sky background to trace the slit edges
  ; for when OH lines are not enough (e.g. in the K band or in the optical)
  input.trace_slit_with_skylines = 0

  ; for longslit
  ;input.longslit = 1
  ;input.longslit_edge = [1133, 1179]

  ; frames to be used for slit identification
  ;input.slitflat_filelist = 'slitflat.txt'
  ;input.slitflat_offset = 4

  ; manual slit positions
  input.slit_position_file = 'slit_edges.reg'

  ; there are lots of cosmic rays on each individual frame
  input.clean_individual_frames = 1

  ;**********************************
  ;**********************************


  ; initialize and create the fuel structure
  fuel = flame_initialize_lris(input)


  ; ; check that everything is good
  ; help, fuel.input
  ; help, fuel.util
  ; help, fuel.instrument
  ; help, fuel.slits



  ;****************************************************
  ;****************************************************
  ;               PART 2: DATA REDUCTION
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


  save, fuel, filename='fuel.sav'


;END
