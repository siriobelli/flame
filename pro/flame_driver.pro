; driver for May run data -  mask

  ;****************************************************
  ;                 INPUT
  ;****************************************************

  fuel = flame_create_fuel()

  ; text file containing the list of science FITS files that need to be reduced
  fuel.science_filelist = 'input/science.txt'

  ; text file containing the list of FITS files with dark frames (used for bad pixel mask)
  fuel.darks_filelist = 'input/darks.txt'

  ; text file containing the list of FITS files with flat field
  fuel.flats_filelist = 'none'
  
  ; name of the directory where intermediate data products will be saved
  fuel.intermediate_dir = 'intermediate/'

  ; name of the directory where the final output files will be saved
  fuel.output_dir = 'output/'

  ; if 0, then reduce all slits. If n, then reduce slit number n (starting from 1).
  fuel.reduce_only_oneslit = 0

  ; array with y-pixel positions for the traces of the reference star. [0,0] if there is no reference star
  fuel.startrace_y_pos = [0, 0]

  ; if we don't have a star on the slit then we have to specify the dithering
  fuel.dither_filelist = 'input/dither.txt'

  ; for longslit
  fuel.longslit = 1
  fuel.longslit_edge = [960, 1090]

  ; for this particular dataset
  fuel.OUTPUT_LAMBDA_0 = 1.12
  fuel.OUTPUT_LAMBDA_DELTA = 7.5d-5
  fuel.OUTPUT_LAMBDA_NPIX = 450
  
  ; create the fuel structure
  flame_initialize_luci, fuel=fuel

;  help, fuel


  ;****************************************************
  ;                 MONITOR STAR
  ;****************************************************
  ; using the star on the reference slit, get seeing, flux, vertical shift,
  ; and dither position for each frame. Also outputs a ps file with the plots.
  
  
  flame_monitor_star, fuel=fuel
  

  ;****************************************************
  ;                 QUICK LOOK
  ;****************************************************
  ; in order to have a quick look at the data, create the simple A-B stack
  

  flame_simple_stack, fuel=fuel


  ;****************************************************
  ;                 BAD PIXEL MASK
  ;****************************************************
  ; generate a bad pixel mask from the dark frames. Should be able to skip and specify a bad pixel mask.
  

  flame_make_badpix, fuel=fuel
  

  ;****************************************************
  ;                 DATA CORRECTION
  ;****************************************************
  ; this step corrects the science frames for linearization and bad pixels
  ; it will output corrected science frames in the intermediate directory
  
  
  flame_correct_data, fuel=fuel
  
  
  ;****************************************************
  ;                 CUTOUT SLITS
  ;****************************************************
  ; create a fits file for each slit in each frame. They will be in intermediate/slitXX/
  
  
  ; first identify slits and create slitim.fits
  flame_do_slitim, fuel=fuel

  ; then cut out the slits
  flame_cutout_slits, fuel=fuel
  
  
  ;****************************************************
  ;                 WAVELENGTH CALIBRATION
  ;****************************************************
  

  flame_wavelength_calibration, fuel=fuel


  ;****************************************************
  ;                 SKY SUBTRACTION
  ;****************************************************


  flame_sky_subtraction, fuel=fuel
  
  

  ;****************************************************
  ;                 RECTIFICATION
  ;****************************************************
  

  flame_rectify, fuel=fuel

  
  ;****************************************************
  ;                 COMBINE FRAMES
  ;****************************************************
  
  
  flame_combine_frames, fuel=fuel
  

  ;****************************************************
  ;                 END: save fuel structure
  ;****************************************************
  

  save, fuel, filename='fuel.sav'

  
;END