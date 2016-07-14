
  ;****************************************************
  ;****************************************************
  ;                 PART 1: INPUT
  ;****************************************************
  ;****************************************************

  fuel = flame_create_fuel()

  ; text file containing the list of science FITS files that need to be reduced
  fuel.science_filelist = 'input/science.txt'

  ; text file containing the list of FITS files with dark frames (used for bad pixel mask)
  ; if 'none', the default dark will be used
  fuel.darks_filelist = 'none'

  ; text file containing the list of FITS files with flat field
  ; if 'none', the default flat field will be used
  fuel.flats_filelist = 'none'
  
  ; name of the directory where intermediate data products will be saved
  fuel.intermediate_dir = 'intermediate/'

  ; name of the directory where the final output files will be saved
  fuel.output_dir = 'output/'

  ; array with y-pixel positions for the traces of the reference star. [0,0] if there is no reference star
  fuel.startrace_y_pos = [1283, 1303]

  
  ; ADVANCED OPTIONS
  ;**********************************

  ; if 0, then reduce all slits. If n, then reduce slit number n (starting from 1).
  fuel.reduce_only_oneslit = 4

  ; if you want to change the range in x-coordinates used to extract the star traces:
  ;fuel.xrange_star = [100, 500]

  ; if we don't have a star on the slit then we have to specify the dithering
  ;fuel.dither_filelist = 'input/dither.txt'

  ; for longslit
  ;fuel.longslit = 1
  ;fuel.longslit_edge = [1133, 1179]

  ;**********************************


  ; create the fuel structure
  flame_initialize_luci, fuel=fuel

  ; check that everything makes sense
  ;help, fuel



  ;****************************************************
  ;****************************************************
  ;                 PART 2: DATA REDUCTION
  ;****************************************************
  ;****************************************************


  ;****************************************************
  ;                 MONITOR STAR
  ;****************************************************
  
  
  flame_diagnostics, fuel=fuel
  

  ;****************************************************
  ;                 QUICK LOOK
  ;****************************************************
  

  flame_quickstack, fuel=fuel


  ;****************************************************
  ;                 DATA CORRECTION
  ;****************************************************
  

  flame_correct, fuel=fuel

  
  ;****************************************************
  ;               IDENTIFY AND CUTOUT SLITS
  ;****************************************************
  

  flame_getslits, fuel=fuel

  
  ;****************************************************
  ;                 WAVELENGTH CALIBRATION
  ;****************************************************
  

  flame_wavecal, fuel=fuel


  ;****************************************************
  ;                 SKY SUBTRACTION
  ;****************************************************


  flame_skysub, fuel=fuel
  

  ;****************************************************
  ;                 RECTIFICATION
  ;****************************************************
  

  flame_rectify, fuel=fuel

  
  ;****************************************************
  ;                 COMBINE FRAMES
  ;****************************************************
  
  
  flame_combine, fuel=fuel
  

  ;****************************************************
  ;                 END: save fuel structure
  ;****************************************************
  

  save, fuel, filename='fuel.sav'

  
;END