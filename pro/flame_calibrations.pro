;
; reads the raw frames and writes out the "corrected" frames, which are:
; - cleaned from cosmic rays (if needed)
; - corrected for non-linearity of the detector
; - corrected for variations in the pixel flat field
; - corrected for bad pixels (replaced by NaNs)
;  (generate the bad pixel mask from the darks and flats, if needed)
; - converted from ADU to electrons
;

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************

FUNCTION flame_calibrations_median_filter, im, radius
  ;
  ; Median-filter an image. The pixel (x,y) in the output image
  ; is equal to the median of the pixels in the input image contained in a box
  ; centered on (x,y) and with side equal to 2*radius+1
  ;

      ; kernel size for finding outliers
      krad = radius

      ; dimensions of frame
      Nx = (size(im))[1]
      Ny = (size(im))[2]

      ; divide the frame into chunks for computational reasons
      ; maximum size of chunks (including overlapping region)
      max_chunk_size = 1500

      ; how many chunks do we need along x and along y
      N_chunks_x = ceil( float(Nx) / float(max_chunk_size-2*krad) )
      N_chunks_y = ceil( float(Ny) / float(max_chunk_size-2*krad) )

      ; x and y sizes of chunks
      chunk_size_x = Nx / N_chunks_x
      chunk_size_y = Ny / N_chunks_y

      ; empty array that will contain the local median values
      median_im = im*0.0

      print, 'Applying median filter...'

      ; loop through all the chunks
      for ix=0, N_chunks_x-1 do $
        for iy=0, N_chunks_y-1 do begin

          print, 'working on chunk ' + strtrim(N_chunks_y*ix + iy + 1, 2) + ' of ' + strtrim(N_chunks_x*N_chunks_y, 2)

          ; calculate limits of this chunk
          x1 = ix*chunk_size_x
          x2 = (x1+chunk_size_x-1) < (Nx-1)
          y1 = iy*chunk_size_y
          y2 = (y1+chunk_size_y-1) < (Nx-1)

          ; to avoid edge issues, expand the chunk
          x1_ext = (x1-krad) > 0
          x2_ext = (x2+krad) < (Nx-1)
          y1_ext = (y1-krad) > 0
          y2_ext = (y2+krad) < (Ny-1)

          ; extract chunk
          chunk = im[x1_ext:x2_ext, y1_ext:y2_ext]

          ; for each pixel in the image, calculate the median of its neighbors
          chunk_median = estimator_filter(chunk, 2*krad+1, /median)
          median_im[x1:x2, y1:y2] = chunk_median[x1-x1_ext:x2-x2_ext-1, y1-y1_ext:y2-y2_ext-1]

        endfor

        return, median_im

END

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_calibrations_median_combine, filenames, outfilename, vertical_shift=vertical_shift

  ; utility to median-combine frames
  ; optionally, apply a vertical shift
  ; NB: first apply a scaling to each frame to match the median values; this is
  ; useful when using sky flats with a quickly varying sky brightness

  ; read the filenames
  print, 'median combining the following frames:'
  forprint, filenames

  ; read in first frame
  first_frame = readfits(filenames[0], hdr)

  ; if there is only one frame, then it's easy
  if n_elements(filenames) eq 1 then master = first_frame else begin

    ; read the type of data (long, float, etc)
    data_type = size(first_frame, /type)

    ; make 3D array containing all frames at once
    cube = make_array( (size(first_frame))[1], (size(first_frame))[2], n_elements(filenames), type=data_type )

    ; read in all frames into the cube
    for i=0, n_elements(filenames)-1 do cube[*,*,i] = readfits(filenames[i])

    ; calculate the sky value (here taken as the median pixel value) of each frame
    sky_value = make_array( n_elements(filenames), type=data_type )
    for i=0, n_elements(filenames)-1 do sky_value[i] = median(cube[*,*,i])

    ; take as reference the median sky value
    ref_value = median(sky_value)

    ; scale every frame so that the combination is not biased
    for i=0, n_elements(filenames)-1 do cube[*,*,i] *= ref_value / sky_value[i]

    ; median-combine the frames
    master = median(cube, dimension=3, /even)

  endelse

  ; if needed, apply vertical shift (without circular boundary conditions)
  if keyword_set(vertical_shift) then begin
    master = shift(master, 0, vertical_shift)
    if vertical_shift GT 0 then master[*,0:vertical_shift-1] = !values.d_nan
    if vertical_shift LT 0 then master[*,-abs(vertical_shift):-1] = !values.d_nan
  endif

  ; write out the master file
  writefits, outfilename, master, hdr

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_calibrations_master_dark, fuel

  ;
  ; Make the master dark by combining the dark frames provided by the user
  ; or, if requested, by copying the default master dark file.
  ; Return the master frame or !NULL
  ;

  print, ''
  print, ''
  print, 'dark frame'
  print, '----------'
  print, ''

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.dark_filelist, 2) )

  ; this will be the output file name of the master file
  master_file = fuel.util.dark.master_file

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
  flame_calibrations_median_combine, fuel.util.dark.raw_files, master_file
  print, 'master dark file created: ', master_file

  ; return the master dark
  master = readfits(master_file, hdr)
  return, master

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_calibrations_master_arc, fuel

  ;
  ; Make the master arcs by combining the arc frames provided by the user
  ; and, if needed, apply a vertical pixel shift
  ; Return the master frame or !NULL
  ;

  print, ''
  print, ''
  print, 'arcs'
  print, '----------'
  print, ''

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.arc_filelist, 2) )

  ; this will be the output file name of the master file
  master_file = fuel.util.arc.master_file

  ; case 1: do not use arcs ---------------------------------------------
  if filelist eq '' or filelist eq 'none' or filelist eq 'default' then begin
    print, 'arcs not used'
    return, !NULL
  endif

  ; case 2: use the frames provided by the user -------------------------------------

  ; median combine the arc frames
  flame_calibrations_median_combine, fuel.util.arc.raw_files, master_file, vertical_shift=fuel.input.arc_pixelshift
  print, 'master arc file created: ', master_file

  ; return the master arcs
  master = readfits(master_file, hdr)
  return, master

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************

FUNCTION flame_calibrations_master_pixelflat, fuel

  ;
  ; Make the master pixel flat by combining the flat frames provided by the user
  ; and then removing the large-scale features
  ; or, if requested, by copying the default pixel flat file.
  ; Return the file name of the master pixel flat
  ;

  print, ''
  print, ''
  print, 'pixel flat'
  print, '-----------------'
  print, ''

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.pixelflat_filelist, 2) )

  ; this will be the output file name of the master flat
  master_file = fuel.util.pixelflat.master_file

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
  median_pixelflat_file = fuel.util.intermediate_dir + 'median_pixelflat.fits'

  ; median combine the pixel flat field frames
  flame_calibrations_median_combine, fuel.util.pixelflat.raw_files, median_pixelflat_file

  ; read the median combined pixel flat
  median_pixelflat = readfits(median_pixelflat_file, hdr)

  ; dimensions of the pixel flat
  Nx = (size(median_pixelflat))[1]
  Ny = (size(median_pixelflat))[2]

  ; create horizontal and vertical profiles (the trick is to normalize each row first)
  norm_y = median(median_pixelflat, dimension=1 ) ## replicate(1, Nx)
  x_profile = median(median_pixelflat / norm_y, dimension=2)

  norm_x = transpose( replicate(1, Ny) # median(median_pixelflat, dimension=2 ) )
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


FUNCTION flame_calibrations_master_illumflat, fuel

  ;
  ; Make the master illumination flat by combining the illumflat frames provided by the user
  ; and, if needed, apply a vertical pixel shift
  ; Return the master frame or !NULL
  ;

  print, ''
  print, ''
  print, 'illumination flat'
  print, '-----------------'
  print, ''

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.illumflat_filelist, 2) )

  ; this will be the output file name of the master file
  master_file = fuel.util.illumflat.master_file

  ; case 1: do not use illumflats ---------------------------------------------
  if filelist eq '' or filelist eq 'none' or filelist eq 'default' then begin
    print, 'illumination flat not used'
    return, !NULL
  endif

  ; case 2: use the frames provided by the user -------------------------------------

  ; median combine the arc frames
  flame_calibrations_median_combine, fuel.util.illumflat.raw_files, master_file, vertical_shift=fuel.input.illumflat_pixelshift
  print, 'master illumination flat file created: ', master_file

  ; return the master illumination flat
  master = readfits(master_file, hdr)
  return, master

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_calibrations_master_slitflat, fuel

  ;
  ; Make the master slit flat by combining the slitflat frames provided by the user
  ; or using the science frames
  ; and, if needed, apply a vertical pixel shift
  ; Return the master frame
  ;

  print, ''
  print, ''
  print, 'slit flat'
  print, '-----------------'
  print, ''

  ; read the file list specified by the user
  filelist = strlowcase( strtrim(fuel.input.slitflat_filelist, 2) )

  ; this will be the output file name of the master file
  master_file = fuel.util.slitflat.master_file

  ; case 1: slitflats not specified ---------------------------------------------
  if filelist eq '' or filelist eq 'none' or filelist eq 'default' then begin

    print, 'slitflat not specified; using science frames'

    ; number of science frames
    N_frames = fuel.util.science.n_frames

    ; if there are three or more science frames, then use the central three as slit flats
    if N_frames GE 3 then $
      files = fuel.util.science.raw_files[fix(N_frames/2)-1:fix(N_frames/2)+1] $
    ; otherwise use all of them
    else files = fuel.util.science.raw_files

    ; case 2: use the frames provided by the user -------------------------------------
  endif else files =fuel.util.slitflat.raw_files

  ; median combine the frames
  flame_calibrations_median_combine, files, master_file, vertical_shift=fuel.input.slitflat_pixelshift
  print, 'master slit flat file created: ', master_file

  ; return the master slit flat
  master = readfits(master_file, hdr)
  return, master

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_calibrations_badpixel, fuel, master_dark, master_pixelflat

  ;
  ; make badpixel mask using the master dark and/or pixel flat field
  ;

  print, ''

  ; pixels in the master dark that are outliers
  ; by more than this value are considered bad pixels
  sig_clip = fuel.settings.badpix_sigma

  ; pixel radius for the median filtering used to remove large-scale structures
  kernel_radius = 10

  ; output file name
  badpix_filename = fuel.util.intermediate_dir + 'badpixel_mask.fits'

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
  if master_dark ne !NULL and fuel.settings.badpix_usedark then begin

      print, ''
      print, 'using master dark frame to identify bad pixels'

      ; construct the median-filtered version of the master dark
      master_dark_smooth = flame_calibrations_median_filter(master_dark, kernel_radius)

      ; subtract the smooth version and get the pixel map
      pixel_map = master_dark - master_dark_smooth

      ; write pixel map
      writefits, fuel.util.intermediate_dir + 'master_dark_filtered.fits', pixel_map

      ; calculate typical value and dispersion for pixel values in a robust way
      mmm, pixel_map, dark_bias, dark_sigma

      ; cut everything outside the central +/- sig_clip sigmas
      low_cut = dark_bias - sig_clip * dark_sigma
      high_cut = dark_bias + sig_clip * dark_sigma
      w_badpixels = where(pixel_map LT low_cut or pixel_map GT high_cut, /null)

      ; calculate the fraction of bad pixels
      badpix_fraction = float(n_elements(w_badpixels))/float(n_elements(pixel_map))
      print, 'Fraction of bad pixels: ' + cgnumber_formatter(badpix_fraction*100.0, decimals=4) + ' %'

      ; remove bad pixels from the master dark for plotting purposes
      pixel_map[w_badpixels] = !values.d_nan

      ; plot distribution
      cgPS_open, fuel.util.intermediate_dir + 'master_dark_histogram.ps', /nomatch
      cghistoplot, pixel_map, /freq, binsize=max([dark_sigma/5.0, 1.0]), $
        xra=dark_bias+[-10.0, 10.0]*dark_sigma, /fillpoly, $
        xtit='pixel value', ytit='frequency', charsize=1.0, xthick=4, ythick=4, $
        title = strtrim(n_elements(w_badpixels),2) + ' bad pixels (' + $
        cgnumber_formatter(badpix_fraction*100.0, decimals=4) + ' % of the total)'
      cgplot, low_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
      cgplot, high_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
      cgplot, dark_bias +[0,0], [0,1.0], /overplot, thick=4
      cgPS_close

      ; create bad pixel mask
      badpix = byte(pixel_map*0.0)

      ; add bad pixels found in the master dark
      badpix[w_badpixels] = 1

  endif $
  else badpix = !NULL


  ; if it is defined, the use the master pixel flat to identify bad pixels ----------------
  if master_pixelflat ne !NULL and fuel.settings.badpix_useflat then begin

      print, ''
      print, 'using master pixel flat field to identify bad pixels'

      ; construct the median-filtered version of the master flat
      master_pixelflat_smooth = flame_calibrations_median_filter(master_pixelflat, kernel_radius)

      ; subtract the smooth version and get the pixel map
      pixelflat_map = (master_pixelflat - master_pixelflat_smooth) / master_pixelflat_smooth

      ; write pixel map
      writefits, fuel.util.intermediate_dir + 'median_pixelflat_filtered.fits', pixelflat_map

      ; calculate typical value and dispersion for pixel values in a robust way
      mmm, pixelflat_map, flat_bias, flat_sigma

      ; we trust corrections up to a point. Beyond that, we consider it a bad pixel
      max_correction = fuel.settings.badpix_flatcorrection

      ; identify bad pixels by either sigma clipping or max correction criteria
      low_cut = (flat_bias-max_correction) > (flat_bias-sig_clip*flat_sigma)
      high_cut = (flat_bias+max_correction) < (flat_bias+sig_clip*flat_sigma)
      w_badpixels = where(pixelflat_map LT low_cut or pixelflat_map GT high_cut, /null)

      ; calculate the fraction of bad pixels
      badpix_fraction = float(n_elements(w_badpixels))/float(n_elements(pixelflat_map))
      print, 'Fraction of bad pixels: ' + cgnumber_formatter(badpix_fraction*100.0, decimals=4) + ' %'

      ; remove bad pixels from the flat field for plotting purposes
      pixelflat_map[w_badpixels] = !values.d_nan

      ; plot distribution
      cgPS_open, fuel.util.intermediate_dir + 'master_pixelflat_histogram.ps', /nomatch
      cghistoplot, pixelflat_map, /freq, binsize=max([flat_sigma/5.0, 0.005]), $
        xra=flat_bias+[-10.0, 10.0]*flat_sigma, /fillpoly, /nan, $
        xtit='pixel value', ytit='frequency', charsize=1.0, xthick=4, ythick=4, $
        title = strtrim(n_elements(w_badpixels),2) + ' bad pixels (' + $
        cgnumber_formatter(badpix_fraction*100.0, decimals=4) + ' % of the total)'
      cgplot, low_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
      cgplot, high_cut +[0,0], [0,1.0], /overplot, thick=3, linestyle=2
      cgplot, flat_bias +[0,0], [0,1.0], /overplot, thick=4
      cgPS_close

      ; create bad pixel mask, if needed
      if badpix EQ !NULL then badpix = byte(pixelflat_map*0.0)

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

PRO flame_calibrations_oneframe, fuel, filename_raw, filename_corr, $
    master_pixelflat=master_pixelflat, master_dark=master_dark, $
    badpixel_mask=badpixel_mask, lacosmic=lacosmic

  ; read in raw frame and convert to double
  frame = double(readfits(filename_raw, header))
  size_science = size(frame)


  ; CORRECTION 0: cosmic rays
  if keyword_set(lacosmic) then begin

    ; identify cosmic rays using L.A. Cosmic
    la_cosmic, filename_raw, gain=fuel.instrument.gain, readn=fuel.instrument.readnoise, niter=5, $
    masklist = flame_util_replace_string( filename_corr, '_corr.fits', '_mask.fits'), $
    outlist = flame_util_replace_string( filename_corr, '_corr.fits', '_cleaned.fits')

    ; read in the cosmic ray mask
    cr_mask = readfits( flame_util_replace_string( filename_corr, '_corr.fits', '_mask.fits') )

    ; set the CRs to NaNs
    frame[where(cr_mask, /null)] = !values.d_nan

  endif


  ; CORRECTION 1: bad pixels
  if badpixel_mask NE !NULL then begin

    ; check that dimensions are right
    if (size(badpixel_mask))[1] NE size_science[1] OR (size(badpixel_mask))[1] NE size_science[1] then $
    message, 'Dimensions of bad pixel mask do not match dimensions of science frame!'

    ; mask bad pixels
    frame[where(badpixel_mask, /null)] = !values.d_nan

  endif


  ; CORRECTION 2: dark frame
  if master_dark NE !NULL and fuel.settings.darksub_data then begin

    ; check that dimensions are right
    if (size(master_dark))[1] NE size_science[1] OR (size(master_dark))[1] NE size_science[1] then $
    message, 'Dimensions of dark frame do not match dimensions of science frame!'

    ; subtract master dark
    frame -= master_dark

  endif


  ; CORRECTION 3: pixel flat field
  if master_pixelflat NE !NULL and fuel.settings.flatfield_data then begin

    ; check that dimensions are right
    if (size(master_pixelflat))[1] NE size_science[1] OR (size(master_pixelflat))[1] NE size_science[1] then $
    message, 'Dimensions of pixel flat do not match dimensions of science frame!'

    ; divide by pixel flat
    frame /= master_pixelflat

  endif


  ; CORRECTION 4: non-linearity
  frame = poly(frame, fuel.instrument.linearity_correction )


  ; CORRECTION 5: trim the edges of the detector

  ; set the size of the margin to trim (default is 1)
  if tag_exist(fuel.instrument, 'trim_edges') then $
    trim = fuel.instrument.trim_edges else $
    trim = 1

  ; set to NaN the pixels at the edge of the detector
  frame[0:trim-1,*] = !values.d_nan
  frame[-trim:-1,*] = !values.d_nan
  frame[*,0:trim-1] = !values.d_nan
  frame[*,-trim:-1] = !values.d_nan


  ; CORRECTION 6: convert to electrons per second

  ; first convert from ADU to electrons
  frame_electrons = frame * fuel.instrument.gain

  ; find the exposure time: try EXPTIME first, then TRUITIME
  exptime = fxpar(header, 'EXPTIME', missing=-1.0)
  if exptime EQ -1.0 then exptime = fxpar(header, 'TRUITIME', missing=-1.0)

  ; if no valid exptime is found, set it to 1 second
  if exptime EQ -1.0 then begin
    print, 'EXPTIME not found; setting to 1 second.'
    exptime = 1.0
  endif

  ; now convert electrons to electrons per second
  frame_eps = frame_electrons / exptime

  ; change the flux units in the header
  fxaddpar, header, 'BUNIT', 'electrons per second', ' '


  ; ------------------------------------
  ; error spectrum

  ; make the error image in units of electrons (Poisson + readnoise )
  frame_sigma_electrons = sqrt( (fuel.instrument.readnoise)^2 + frame_electrons )

  ; convert to electrons per second
  frame_sigma_eps = frame_sigma_electrons / exptime

  ; ------------------------------------
  ; write output

  ; corrected frame in the first HDU, error frame in the second one
  writefits, filename_corr, frame_eps, header
  writefits, filename_corr, frame_sigma_eps, /append


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_calibrations, fuel

	flame_util_module_start, fuel, 'flame_calibrations'

  ; number of science frames
  N_frames = fuel.util.science.n_frames


  ; --------------------------
  ; --- make master frames ---
  ; --------------------------

  ; create the master dark frame
  master_dark = flame_calibrations_master_dark( fuel )

  ; create the master arc frame (!NULL because we don't need the master arc right now)
  !NULL = flame_calibrations_master_arc( fuel )

  ; create the master pixel flat
  master_pixelflat = flame_calibrations_master_pixelflat( fuel )

  ; create the master illumination flat
  !NULL = flame_calibrations_master_illumflat( fuel )

  ; create the master slit flat
  !NULL = flame_calibrations_master_slitflat( fuel )


  ; make bad pixel mask using darks and/or flats
  badpixel_mask = flame_calibrations_badpixel( fuel, master_dark, master_pixelflat )


  ; -------------------------------------------
  ; --- apply corrections to science frames ---
  ; -------------------------------------------

  ; L.A. Cosmic notice
  if fuel.settings.clean_individual_frames then begin
    print, ''
    print, '----------------------------------------------------'
    print, ''
    print, ' L.A. Cosmic: Laplacian cosmic ray removal'
    print, ''
    print, '      by Pieter van Dokkum'
    print, '    IDL version by Josh Bloom'
    print, ''
    print, ' see van Dokkum P. G., 2001, PASP, 113, 1420  '
    print, '----------------------------------------------------'
  endif


  ; filenames for science frames
  raw_filenames = fuel.util.science.raw_files
  corr_filenames = fuel.util.science.corr_files

  print, ''
  for i_frame=0, n_elements(raw_filenames)-1 do begin

    print, ''
    print, 'Correcting frame ', raw_filenames[i_frame]

    ; if running LACosmic, then check if the corrected frame already exist to save time
    if fuel.settings.clean_individual_frames AND file_test(corr_filenames[i_frame]) then begin
      print, 'file already exists; skipping frame correction'
      continue
    endif

    ; apply corrections and create "corrected" file
    flame_calibrations_oneframe, fuel, raw_filenames[i_frame], corr_filenames[i_frame], $
      master_pixelflat=master_pixelflat, master_dark=master_dark, $
      badpixel_mask=badpixel_mask, lacosmic=fuel.settings.clean_individual_frames

  endfor


  flame_util_module_end, fuel

END
