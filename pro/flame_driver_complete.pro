
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

  ; text file containing the list of FITS files with dark frames
  input.dark_filelist = 'none'

  ; text file containing the list of FITS files for pixel-flat field
  input.pixelflat_filelist = 'none'

  ; text file containing the list of FITS files for illumination-flat field
  input.illumflat_filelist = 'none'

  ; text file containing the list of FITS files with arcs for wavelength calibration
  input.arc_filelist = 'none'

  ; text file containing the list of FITS files for slit-flat field
  input.slitflat_filelist = 'none'

  ; vertical offset, in pixel, between the slit-flat field and the science frames
  input.slitflat_offset = 0

  ; manual slit positions
  input.slit_position_file = 'none'

  ; if we don't have a star on the slit then we have to specify the dithering
  ; in this case, you must set input.star_y_A and input.star_y_B to 0
  input.dither_file = 'none'

  ; for when the alignment boxes have varying width (otherwise leave to 0)
  input.max_slitwidth_arcsec = 0.0

  ; specify the directory for the intermediate products
  input.intermediate_dir = 'intermediate/'

  ; specify the directory for the final output
  input.output_dir = 'output/'



  ;****************************************************
  ;****************************************************
  ;              PART 2: INITIALIZATION
  ;****************************************************
  ;****************************************************


  ; initialize and create the fuel structure
  ; need to use the initialization routine of the correct instrument
  fuel = flame_initialize_luci(input)
  ;fuel = flame_initialize_lris(input)

  ; Here are further options for special cases or troubleshooting.
  ; Most times the settings in fuel.settings should be left to their default values

  ; change the range in x-coordinates used to extract the star traces:
  ; fuel.settings.star_x_range = [100, 500]

  ; change the extension of the range in y-coordinates used to fit the star traces:
  ; fuel.settings.star_y_window = 40

  ; set this to zero if you want to use the sky background to trace the slit edges
  ; instead of the OH lines (e.g. in the K band or in the optical or with slit flats)
  ; fuel.settings.trace_slit_with_skylines = 0

  ; change the spectral resolution used in the successive steps during wavecal_rough
  ; fuel.settings.wavecal_rough_R = [500, 1000, 3000]

  ; identify cosmic rays using L.A.Cosmic in each science frame
  ; fuel.settings.clean_individual_frames = 1

  ; subtract sky continuum before wavecal_rough. This specifies the smoothing window
  ; fuel.settings.wavecal_rough_smooth_window = 20

  ; if the sky lines are weak, stack more than one pixel row when identifying lines
  ;fuel.settings.identify_lines_stack_rows = 2

  ; change the window used to look for a spectral line, in units of expected line width
  ;fuel.settings.identify_lines_linefit_window = 6.0

  ; change the minimum number of emission lines required for a valid wavelength
  ; solution in one pixel row during identify_lines
  ; fuel.settings.identify_lines_Nmin_lines = 4

  ; change the polynomial degree of the wavelength solution for one pixel row
  ; during identify_lines
  ; fuel.settings.identify_lines_poly_degree = 2

  ; change the degree of the 2D wavelength solution describing one entire cutout
  ; fuel.settings.wavesolution_order_x = 3
  ; fuel.settings.wavesolution_order_y = 2

  ; change the threshold for sigma-clipping of each pixel when combining all the frames
  ; fuel.settings.combine_sigma_clip = 2.0



  ; check that everything is good
  ;help, fuel.input
  ;help, fuel.util
  ;help, fuel.settings
  ;help, fuel.instrument
  ;help, fuel.slits



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
