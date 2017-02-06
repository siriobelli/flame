;
; reads the raw science frames and writes out the "corrected" science frames, which are:
; - corrected for non-linearity of the detector
; - corrected for bad pixels (replaced by NaNs) 
;  (generate the bad pixel mask from the dark frames, if needed)
; - converted from ADU to electrons
;
; TO-DO: flat field correction
;




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



;******************************************************************



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



;******************************************************************


PRO flame_correct, fuel=fuel


  ; bad pixel mask - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  ; decide whether we are going to make a bad pixel mask 
  if fuel.util.darks_filenames[0] eq '' then begin

    ; copy the default bad pixel mask to the intermediate directory
    file_copy, fuel.util.flame_data_dir + fuel.instrument.default_badpixel_mask, $
      fuel.input.intermediate_dir + 'badpixel_mask.fits', /overwrite

  endif else begin

    ; make bad pixel mask 
    flame_correct_makemask, fuel=fuel

  endelse

  ; read in the bad pixel mask
  badpix = readfits(fuel.input.intermediate_dir + 'badpixel_mask.fits')


  ; master flat field - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  ; decide whether we are going to make a master flat field 
  if fuel.util.flats_filenames[0] eq '' then begin

    print, 'no flat field'

  endif else begin 

    ; make master flat field
    flame_correct_makeflat, fuel=fuel

  endelse


  ; apply corrections - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    
  for i_frame=0,fuel.util.N_frames-1 do begin
  
    print, 'Correcting frame ', fuel.util.science_filenames[i_frame]

    ; read in science frame
    frame = readfits(fuel.util.science_filenames[i_frame], header)

    ; CORRECTION 1: non-linearity
    frame_corr1 = poly(frame, fuel.instrument.linearity_correction )

    ; CORRECTION 2: bad pixels
    frame_corr2 = frame_corr1
    frame_corr2[where(badpix, /null)] = !values.d_nan

    ; CORRECTION 3: convert to electrons per second
    frame_corr3 = frame_corr2 * fuel.instrument.gain

    ; change the flux units in the header
    fxaddpar, header, 'BUNIT', 'electrons', ' '

    ; save corrected frame
    writefits, (fuel.util.corrscience_filenames)[i_frame], frame_corr3, header
    
  endfor
  
  
END