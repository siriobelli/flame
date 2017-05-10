;
; reads the raw science frames and writes out the "corrected" science frames, which are:
; - cleaned from cosmic rays (if needed)
; - corrected for non-linearity of the detector
; - corrected for bad pixels (replaced by NaNs)
;  (generate the bad pixel mask from the dark frames, if needed)
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


FUNCTION flame_correct_master_dark, fuel=fuel

  ;
  ; Make the master dark by combining the dark frames provided by the user
  ; or, if requested, by copying the default master dark file.
  ; Return the file name of the master dark
  ;

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.dark_filelist, 2) )

  ; this will be the output file name of the master dark
  master_file = fuel.util.master_dark

  ; case 1: do not use darks ---------------------------------------------
  if filelist eq '' or filelist eq 'none' then begin
    print, 'dark not used'
    return, ''
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

    ; return the filename
    return, master_file

  endif


  ; case 3: use the frames provided by the user -------------------------------------

  ; median combine the dark frames
  flame_correct_median_combine, fuel.util.filenames_dark, master_file
  print, 'master dark file created: ', master_file

  ; return the filename
  return, master_file


END




;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_correct_master_pixelflat, fuel=fuel

  ;
  ; Make the master pixel flat by combining the flat frames provided by the user
  ; and then removing the large-scale features
  ; or, if requested, by copying the default pixel flat file.
  ; Return the file name of the master pixel flat
  ;

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.pixelflat_filelist, 2) )

  ; this will be the output file name of the master dark
  master_file = fuel.util.master_pixelflat

  ; case 1: do not use darks ---------------------------------------------
  if filelist eq '' or filelist eq 'none' then begin
    print, 'pixel flat field not used'
    return, ''
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

    ; return the filename
    return, master_file

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

  ; return the filename
  return, master_file


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_correct_makeflat, fuel=fuel


  ; read in all flat frames
  flats = []
  for i_frame=0, n_elements(fuel.util.flats_filenames)-1 do $
    flats = [ [[flats]], [[ readfits(fuel.util.flats_filenames[i_frame], header) ]] ]

  ; median combine all the flat frames
  master_flat = median(flats, dimension=3)

  ; heavily smoothed version of the flat - this is to ignore the hot pixels
  smoothed_flat = gauss_smooth(master_flat, 5, /edge_truncate)

  ; determine the brightest spot on the master flat - excluding hot pixels
  flat_peak = max(smoothed_flat)

  ; assume that pixels with less than 15% of the peak flux are not being illuminated
  w_notilluminated = where(smoothed_flat LT 0.15 * flat_peak )
  print, 100.0*double( n_elements(w_notilluminated) ) / double( n_elements(smoothed_flat) ), $
    '% of the flat pixels are not illuminated'

  ; set NaNs for pixels that are not illuminated
  master_flat[w_notilluminated] = !values.d_nan
  smoothed_flat[w_notilluminated] = !values.d_nan

  ; create horizontal and vertical profiles (the trick is to normalize each row first)
  norm_y = median( master_flat, dimension=1 ) ## replicate(1, 2048)
  x_profile = median(master_flat / norm_y, dimension=2)

  norm_x = transpose( replicate(1, 2048) # median( master_flat, dimension=2 ) )
  y_profile = median(master_flat / norm_x, dimension=1)

  ; make a model for the illumination using the smoothed profiles
  flat_model = median(x_profile, 101) # median(y_profile, 101)

  ; divide the model out
  smallscale_field = master_flat / flat_model

  ; normalize
  smallscale_field /= median(smallscale_field)

  ; save the master flat
  writefits, fuel.input.intermediate_dir + 'master_flat.fits', smallscale_field, header

  ; now identify the bad pixels, truncate edges, etc

  ;
  ; rename variables
  ; consider outputting a PS file and some info on bad pixels and flat fielding
  ;

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_correct_badpix, fuel=fuel, master_dark_file=master_dark_file, master_pixelflat_file=master_pixelflat_file

  ; if there is input, then return !NULL
  if master_dark_file eq '' then return, !NULL

  ; read in the master_dark
  master_dark = readfits(master_dark_file, hdr)

  ; calculate typical value and dispersion for pixel values in a robust way
  mmm, master_dark, dark_bias, dark_sigma

  ; cut everything outside the central +/- 7 sigmas
  low_cut = dark_bias - 7.0 * dark_sigma
  high_cut = dark_bias + 7.0 * dark_sigma
  w_badpixels = where(master_dark LT low_cut or master_dark GT high_cut, /null)

  ; calculate the fraction of bad pixels
  badpix_fraction = float(n_elements(w_badpixels))/float(n_elements(master_dark))
  print, 'Fraction of bad pixels: ' + number_formatter(badpix_fraction*100.0, decimals=4) + ' %'

  ; plot distribution
  cgPS_open, fuel.input.intermediate_dir + 'master_dark_histogram.ps', /nomatch
  cghistoplot, master_dark, /freq, binsize=max([0.1*dark_sigma, 1.0]), $
    xra=dark_bias+[-10.0, 10.0]*dark_sigma, /fillpoly, $
    xtit='pixel value', ytit='frequency', charsize=1.0, xthick=4, ythick=4, $
    title = strtrim(n_elements(w_badpixels),2) + ' bad pixels (' + $
    number_formatter(badpix_fraction*100.0, decimals=4) + ' % of the total)'
  cgplot, low_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
  cgplot, high_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
  cgplot, dark_bias +[0,0], [0,1.0], /overplot, thick=4
  cgPS_close

  ; create bad pixel mask
  badpixel_mask = byte(master_dark*0.0)
  badpixel_mask[w_badpixels] = 1

  ; save bad pixel mask
  badpix_filename = fuel.input.intermediate_dir + 'badpixel_mask.fits'
  writefits, badpix_filename, badpixel_mask, header

  ; return bad pixel mask
  return, badpixel_mask

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_correct, fuel=fuel

  ; create the master dark file
  master_dark_file = flame_correct_master_dark( fuel=fuel )

  ; create the master pixel flat file
  master_pixelflat_file = flame_correct_master_pixelflat( fuel=fuel )

  ; make bad pixel mask using darks and/or flats
  ; badpix = flame_correct_badpix( fuel=fuel, master_dark_file=master_dark_file, $
  ;   master_pixelflat_file=master_pixelflat_file)
  badpix = !NULL

  ; apply corrections ----------------------------------------------------------

  for i_frame=0,fuel.util.N_frames-1 do begin

    print, 'Correcting frame ', fuel.util.science_filenames[i_frame]

    ; name for the corrected frame (i.e. the output)
    filename_corr = (fuel.util.corrscience_filenames)[i_frame]

    ; read in science frame
    frame = readfits(fuel.util.science_filenames[i_frame], header)

    ; CORRECTION 0: cosmic rays
    if fuel.input.clean_individual_frames then begin

      ; identify cosmic rays using L.A. Cosmic
      la_cosmic, fuel.util.science_filenames[i_frame], gain=fuel.instrument.gain, readn=fuel.instrument.read_noise, $
      masklist = flame_util_replace_string( filename_corr, '_corr.fits', '_mask.fits'), $
      outlist = flame_util_replace_string( filename_corr, '_corr.fits', '_cleaned.fits')

      ; read in the cosmic ray mask
      cr_mask = readfits( flame_util_replace_string( filename_corr, '_corr.fits', '_mask.fits') )

      ; set the CRs to NaNs
      frame[where(cr_mask, /null)] = !values.d_nan

    endif

    ; CORRECTION 1: non-linearity
    frame_corr1 = poly(frame, fuel.instrument.linearity_correction )

    ; CORRECTION 2: bad pixels
    frame_corr2 = frame_corr1
    if badpix NE !NULL then $
      frame_corr2[where(badpix, /null)] = !values.d_nan

    ; CORRECTION 3: convert to electrons per second
    frame_corr3 = frame_corr2 * fuel.instrument.gain

    ; change the flux units in the header
    fxaddpar, header, 'BUNIT', 'electrons', ' '

    ; save corrected frame
    writefits, filename_corr, frame_corr3, header

  endfor


END
