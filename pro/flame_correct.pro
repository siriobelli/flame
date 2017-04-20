;
; reads the raw science frames and writes out the "corrected" science frames, which are:
; - cleaned from cosmic rays (if needed)
; - corrected for non-linearity of the detector
; - corrected for bad pixels (replaced by NaNs)
;  (generate the bad pixel mask from the dark frames, if needed)
; - converted from ADU to electrons
;
; TO-DO: flat field correction
;


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_correct_makemaster, fuel=fuel, calib=calib
  ;
  ; calib (input): 'dark', 'pixelflat', 'illumflat', 'arc'
  ; output the name of the master file,
  ; or an empty string if no calibration will be applied
  ;

  ; select the type of calibration ---------------------------------------------
  case calib of
    'dark': begin
      filelist = fuel.input.dark_filelist
      default_master = fuel.instrument.default_dark
      filenames = fuel.util.filenames_dark
      master_file = fuel.util.master_dark
    end
    'pixelflat': begin
      filelist = fuel.input.pixelflat_filelist
      default_master = fuel.instrument.default_pixelflat
      filenames = fuel.util.filenames_pixelflat
      master_file = fuel.util.master_pixelflat
    end
    'illumflat': begin
      filelist = fuel.input.illumflat_filelist
      default_master = fuel.instrument.default_illumflat
      filenames = fuel.util.filenames_illumflat
      master_file = fuel.util.master_illumflat
    end
    'arc': begin
      filelist = fuel.input.arc_filelist
      default_master = fuel.instrument.default_arc
      filenames = fuel.util.filenames_arc
      master_file = fuel.util.master_arc
    end
    else: message, 'calib keyword not valid: ' + calib
  endcase

  ; get input from user
  filelist_norm = strlowcase( strtrim(filelist, 2) )

  ; 1) do not use this calibration ---------------------------------------------
  if filelist_norm eq '' or filelist_norm eq 'none' then begin
    print, calib + ' not used'
    return, ''
  endif

  ; 2) use the default master file ---------------------------------------------
  if filelist_norm eq 'default' then begin

    ; check that the default dark frame is defined
    if ~file_test(default_master) then $
      message, 'default calibration ' + default_master + ' not found!'

    ; copy the default master frame to the local (intermediate) directory
    file_copy, default_master, master_file, /overwrite

    return, master_file

  endif


  ; 3) use the frames provided by the user -------------------------------------

  print, 'making master out of', filenames



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


PRO flame_correct_makemask, fuel=fuel

  ; identify the darks files
  print, n_elements(fuel.util.darks_filenames), ' dark frames'

  ; read in all dark frames
  darks = []
  for i_frame=0, n_elements(fuel.util.darks_filenames)-1 do darks = [ [[darks]], [[ readfits(fuel.util.darks_filenames[i_frame], header) ]] ]

  ; median combine all the dark frames
  master_dark = median(darks, dimension=3)

  ; save the master dark
  writefits, fuel.input.intermediate_dir + 'master_dark.fits', master_dark, header

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
  writefits, fuel.input.intermediate_dir + 'badpixel_mask.fits', badpixel_mask, header

END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_correct, fuel=fuel

  ; create master files ('' for those that are not to be used)
  master_dark = flame_correct_makemaster( fuel=fuel, calib='dark')
  master_pixelflat = flame_correct_makemaster( fuel=fuel, calib='pixelflat')
  master_illumflat = flame_correct_makemaster( fuel=fuel, calib='illumflat')
  master_arc = flame_correct_makemaster( fuel=fuel, calib='arc')

  ; also need to make bad pixel mask here
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
