;
; reads the raw science frames and writes out the "corrected" science frames, which are:
; - corrected for non-linearity of the detector
; - corrected for bad pixels (replaced by NaNs) 
;  (generate the bad pixel mask from the dark frames, if needed)
; - converted from ADU to electrons
;
; TO-DO: flat field correction
;


PRO flame_correct_makemask, fuel=fuel

  ; identify the darks files
  print, n_elements(fuel.util.darks_filelist), ' dark frames'

  ; read in all dark frames
  darks = []
  for i_frame=0, n_elements(fuel.util.darks_filelist)-1 do darks = [ [[darks]], [[ readfits(fuel.util.darks_filelist[i_frame], header) ]] ]

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

  ; decide whether we are going to make a bad pixel mask 
  if fuel.util.darks_filenames[0] eq '' then begin

    ; copy the default bad pixel mask to the intermediate directory
    spawn, 'cp ' + fuel.util.flame_data_dir + 'default_badpixel_mask.fits ' + $
      fuel.input.intermediate_dir + 'badpixel_mask.fits'
  
  endif else begin

    ; make bad pixel mask 
    flame_correct_makemask, fuel=fuel

  endelse

  ; read in the bad pixel mask
  badpix = readfits(fuel.input.intermediate_dir + 'badpixel_mask.fits')
    
  for i_frame=0,fuel.util.N_frames-1 do begin
  
    print, 'Correcting frame ', fuel.util.science_filenames[i_frame]

    ; read in science frame
    frame = readfits(fuel.util.science_filenames[i_frame], header, /silent)

    ; CORRECTION 1: non-linearity
    frame_corr1 = frame + 4.155d-6 * frame^2    ; from http://scienceops.lbto.org/sciops_cookbook/luci2-vs-luci1/  (Jan 2016)
    
    ; CORRECTION 2: bad pixels
    frame_corr2 = frame_corr1
    frame_corr2[where(badpix, /null)] = !values.d_nan

    ; CORRECTION 3: convert to electrons per second
    frame_corr3 = frame_corr2 * (*fuel.instrument).gain

    ; change the flux units in the header
    fxaddpar, header, 'BUNIT', 'electrons', ' '
    
    ; save corrected frame
    writefits, (fuel.util.corrscience_filenames)[i_frame], frame_corr3, header
    
  endfor
  
  
END