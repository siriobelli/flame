;
; reads the raw science frames and writes out the "corrected" science frames, which are:
; - cleaned from cosmic rays (if needed)
; - corrected for non-linearity of the detector
; - corrected for variations in the pixel flat field
; - corrected for bad pixels (replaced by NaNs)
;  (generate the bad pixel mask from the darks and flats, if needed)
; - converted from ADU to electrons
;
; TO-DO: illumflat (for when OH lines are not enough) and arcs
;

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_correct_median_combine, filenames, outfilename

  ; utility to median-combine frames

  ; read the filenames
  print, 'median combining the following frames:'
  forprint, filenames

  ; read in first frame
  first_frame = readfits(filenames[0], hdr)

  ; read the type of data (long, float, etc)
  data_type = size(first_frame, /type)

  ; make 3D array containing all frames at once
  cube = make_array( (size(first_frame))[1], (size(first_frame))[2], n_elements(filenames), type=data_type )

  ; read in all frames into the cube
  for i=0, n_elements(filenames)-1 do cube[*,*,i] = readfits(filenames[i])

  ; take the median
  master = median(cube, dimension=3)

  ; write out the master file
  writefits, outfilename, master, hdr

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_correct_master_dark, fuel

  ;
  ; Make the master dark by combining the dark frames provided by the user
  ; or, if requested, by copying the default master dark file.
  ; Return the frame or !NULL
  ;

  print, ''

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.dark_filelist, 2) )

  ; this will be the output file name of the master dark
  master_file = fuel.util.master_dark

  ; case 1: do not use darks ---------------------------------------------
  if filelist eq '' or filelist eq 'none' then begin
    print, 'dark not used'
    return, !NULL
  endif

  ; case 2: use the default master file ---------------------------------------------
  if filelist eq 'default' then begin

    ; default master file
    default_master = fuel.util.flame_data_dir + fuel.instrument.default_dark

    ; check that the default dark frame is defined
    if ~file_test(default_master) then $
    message, 'default master dark file not found!'

    ; copy the default master frame to the local (intermediate) directory
    file_copy, default_master, master_file, /overwrite
    print, 'using default master dark ', default_master

    ; return the master dark
    master = readfits(master_file, hdr)
    return, master

  endif


  ; case 3: use the frames provided by the user -------------------------------------

  ; median combine the dark frames
  flame_correct_median_combine, fuel.util.filenames_dark, master_file
  print, 'master dark file created: ', master_file

  ; return the master dark
  master = readfits(master_file, hdr)
  return, master

END




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************

FUNCTION flame_correct_master_pixelflat, fuel

  ;
  ; Make the master pixel flat by combining the flat frames provided by the user
  ; and then removing the large-scale features
  ; or, if requested, by copying the default pixel flat file.
  ; Return the file name of the master pixel flat
  ;

  print, ''

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.pixelflat_filelist, 2) )

  ; this will be the output file name of the master dark
  master_file = fuel.util.master_pixelflat

  ; case 1: do not use pixel flat ---------------------------------------------
  if filelist eq '' or filelist eq 'none' then begin
    print, 'pixel flat field not used'
    return, !NULL
  endif

  ; case 2: use the default master file ---------------------------------------------
  ; note: this is the normalized pixel flat field, not the median of the flat field frames
  if filelist eq 'default' then begin

    ; default master file
    default_master = fuel.util.flame_data_dir + fuel.instrument.default_pixelflat

    ; check that the default pixel flat field is defined
    if ~file_test(default_master) then $
    message, 'default master pixel flat field file not found!'

    ; copy the default master frame to the local (intermediate) directory
    file_copy, default_master, master_file, /overwrite
    print, 'using default pixel flat field ', default_master

    ; return the master pixel flat
    master = readfits(master_file, hdr)
    return, master

  endif


  ; case 3: use the frames provided by the user -------------------------------------

  ; first we need to median combine the pixel flat frames
  median_pixelflat_file = fuel.input.intermediate_dir + 'median_pixelflat.fits'

  ; median combine the pixel flat field frames
  flame_correct_median_combine, fuel.util.filenames_pixelflat, median_pixelflat_file

  ; read the median combined pixel flat
  median_pixelflat = readfits(median_pixelflat_file, hdr)

  ; create horizontal and vertical profiles (the trick is to normalize each row first)
  norm_y = median(median_pixelflat, dimension=1 ) ## replicate(1, 2048)
  x_profile = median(median_pixelflat / norm_y, dimension=2)

  norm_x = transpose( replicate(1, 2048) # median(median_pixelflat, dimension=2 ) )
  y_profile = median(median_pixelflat / norm_x, dimension=1)

  ; make a model for the illumination using the smoothed profiles
  model_pixelflat = median(x_profile, 101) # median(y_profile, 101)

  ; divide the model out
  pixelflat = median_pixelflat / model_pixelflat

  ; normalize
  pixelflat /= median(pixelflat)

  ; write the pixel flat
  writefits, master_file, pixelflat, hdr
  print, 'master pixel flat field file created: ', master_file

  ; return the master pixel flat
  master = readfits(master_file, hdr)
  return, master

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_correct_badpixel, fuel, master_dark, master_pixelflat

  ;
  ; make badpixel mask using the master dark and/or pixel flat field
  ;

  print, ''

  ; pixels in the master dark that are outliers
  ; by more than this value are considered bad pixels
  sig_clip = 7.0

  ; output file name
  badpix_filename = fuel.input.intermediate_dir + 'badpixel_mask.fits'

  ; if there are no inputs, then use the default bad pixel mask ----------------
  if master_dark eq !NULL and master_pixelflat eq !NULL then begin

    if fuel.instrument.default_badpixel_mask eq 'none' then begin
      print, 'no bad pixel mask found'
      return, !NULL
    endif

    ; default badpix mask
    default_badpix_file = fuel.util.flame_data_dir + fuel.instrument.default_badpixel_mask

    ; check that the default mask is defined
    if ~file_test(default_badpix_file) then $
      message, 'default bad pixel mask file not found!'

    ; copy the default mask to the local (intermediate) directory
    file_copy, default_badpix_file, badpix_filename, /overwrite
    print, 'using the default badpixel mask ', default_badpix_file

    ; return the default mask
    badpix = readfits(badpix_filename, hdr)
    return, badpix
  endif


  ; if it is defined, the use the master dark to identify bad pixels ----------------
  if master_dark ne !NULL then begin

      print, ''
      print, 'using master dark frame to identify bad pixels'

      ; calculate typical value and dispersion for pixel values in a robust way
      mmm, master_dark, dark_bias, dark_sigma

      ; cut everything outside the central +/- sig_clip sigmas
      low_cut = dark_bias - sig_clip * dark_sigma
      high_cut = dark_bias + sig_clip * dark_sigma
      w_badpixels = where(master_dark LT low_cut or master_dark GT high_cut, /null)

      ; calculate the fraction of bad pixels
      badpix_fraction = float(n_elements(w_badpixels))/float(n_elements(master_dark))
      print, 'Fraction of bad pixels: ' + cgnumber_formatter(badpix_fraction*100.0, decimals=4) + ' %'

      ; plot distribution
      cgPS_open, fuel.input.intermediate_dir + 'master_dark_histogram.ps', /nomatch
      cghistoplot, master_dark, /freq, binsize=max([0.1*dark_sigma, 1.0]), $
        xra=dark_bias+[-10.0, 10.0]*dark_sigma, /fillpoly, $
        xtit='pixel value', ytit='frequency', charsize=1.0, xthick=4, ythick=4, $
        title = strtrim(n_elements(w_badpixels),2) + ' bad pixels (' + $
        cgnumber_formatter(badpix_fraction*100.0, decimals=4) + ' % of the total)'
      cgplot, low_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
      cgplot, high_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
      cgplot, dark_bias +[0,0], [0,1.0], /overplot, thick=4
      cgPS_close

      ; create bad pixel mask
      badpix = byte(master_dark*0.0)

      ; add bad pixels found in the master dark
      badpix[w_badpixels] = 1

  endif $
  else badpix = !NULL


  ; if it is defined, the use the master pixel flat to identify bad pixels ----------------
  if master_pixelflat ne !NULL then begin

      print, ''
      print, 'using master pixel flat field to identify bad pixels'

      ; calculate typical value and dispersion for pixel values in a robust way
      mmm, master_pixelflat, flat_bias, flat_sigma

      ; we trust corrections up to 20%. Beyond that, we consider it a bad pixel,
      max_correction = 0.20

      ; however if the standard deviation is very large, then use that as threshold
      if sig_clip*flat_sigma GT max_correction then max_correction = sig_clip*flat_sigma

      ; identify bad pixels
      low_cut = flat_bias - max_correction
      high_cut = flat_bias + max_correction
      w_badpixels = where(master_pixelflat LT low_cut or master_pixelflat GT high_cut, /null)

      ; calculate the fraction of bad pixels
      badpix_fraction = float(n_elements(w_badpixels))/float(n_elements(master_pixelflat))
      print, 'Fraction of bad pixels: ' + cgnumber_formatter(badpix_fraction*100.0, decimals=4) + ' %'

      ; plot distribution
      cgPS_open, fuel.input.intermediate_dir + 'master_pixelflat_histogram.ps', /nomatch
      cghistoplot, master_pixelflat, /freq, binsize=max([flat_sigma, 0.005]), $
        xra=flat_bias+[-10.0, 10.0]*flat_sigma, /fillpoly, $
        xtit='pixel value', ytit='frequency', charsize=1.0, xthick=4, ythick=4, $
        title = strtrim(n_elements(w_badpixels),2) + ' bad pixels (' + $
        cgnumber_formatter(badpix_fraction*100.0, decimals=4) + ' % of the total)'
      cgplot, low_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
      cgplot, high_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
      cgplot, flat_bias +[0,0], [0,1.0], /overplot, thick=4
      cgPS_close

      ; create bad pixel mask, if needed
      if badpix EQ !NULL then badpix = byte(master_pixelflat*0.0)

      ; add the bad pixels found in the pixel flat
      badpix[w_badpixels] = 1

  endif

  ; write bad pixel mask
  writefits, badpix_filename, badpix

  ; return bad pixel mask
  return, badpix

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_correct, fuel

	flame_util_module_start, fuel, 'flame_correct'


  ; create the master dark
  master_dark = flame_correct_master_dark( fuel )

  ; create the master pixel flat
  master_pixelflat = flame_correct_master_pixelflat( fuel )

  ; make bad pixel mask using darks and/or flats
  badpixel_mask = flame_correct_badpixel( fuel, master_dark, master_pixelflat )

  ; create the master slit flat
  if fuel.util.filenames_slitflat NE !NULL then $
    flame_correct_median_combine, fuel.util.filenames_slitflat, fuel.util.master_getslit


  ; apply corrections ----------------------------------------------------------

  print, ''
  for i_frame=0,fuel.util.N_frames-1 do begin

    print, 'Correcting frame ', fuel.util.science_filenames[i_frame]

    ; name for the corrected frame (i.e. the output)
    filename_corr = (fuel.util.corrscience_filenames)[i_frame]

    ; read in science frame
    frame = readfits(fuel.util.science_filenames[i_frame], header)

    ; CORRECTION 0: cosmic rays
    if fuel.util.clean_individual_frames then begin

      ; identify cosmic rays using L.A. Cosmic
      la_cosmic, fuel.util.science_filenames[i_frame], gain=fuel.instrument.gain, readn=fuel.instrument.readnoise, $
      masklist = flame_util_replace_string( filename_corr, '_corr.fits', '_mask.fits'), $
      outlist = flame_util_replace_string( filename_corr, '_corr.fits', '_cleaned.fits')

      ; read in the cosmic ray mask
      cr_mask = readfits( flame_util_replace_string( filename_corr, '_corr.fits', '_mask.fits') )

      ; set the CRs to NaNs
      frame[where(cr_mask, /null)] = !values.d_nan

    endif

    ; CORRECTION 1: non-linearity
    frame_corr1 = poly(frame, fuel.instrument.linearity_correction )

    ; CORRECTION 2: pixel flat field
    if master_pixelflat NE !NULL then $
      frame_corr2 = frame_corr1 / master_pixelflat $
      else frame_corr2 = frame_corr1

    ; CORRECTION 3: bad pixels
    frame_corr3 = frame_corr2
    if badpixel_mask NE !NULL then $
      frame_corr3[where(badpixel_mask, /null)] = !values.d_nan

    ; CORRECTION 4: convert to electrons per second

    ; first convert from ADU to electrons
    frame_electrons = frame_corr3 * fuel.instrument.gain

    ; find the exposure time
    exptime = fxpar(header, 'EXPTIME', missing=-1.0)
    if exptime EQ -1.0 then begin
      print, 'EXPTIME not found; setting to 1 second.'
      exptime = 1.0
    endif

    ; now convert electrons to electrons per second
    frame_corr4 = frame_electrons / exptime

    ; change the flux units in the header
    fxaddpar, header, 'BUNIT', 'electrons per second', ' '

    ; ------------------------------------
    ; error spectrum

    ; make the error image in units of electrons (Poisson + readnoise )
    frame_sigma_electrons = fuel.instrument.readnoise + sqrt(frame_electrons)

    ; convert to electrons per second
    frame_sigma = frame_sigma_electrons / exptime

    ; ------------------------------------
    ; write output

    ; corrected frame in the first HDU, error frame in the second one
    writefits, filename_corr, frame_corr4, header
    writefits, filename_corr, frame_sigma, /append


  endfor


  flame_util_module_end, fuel

END
