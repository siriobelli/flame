;
; reads the raw science frames and writes out the "corrected" science frames, which are:
; - corrected for non-linearity of the detector
; - corrected for bad pixels (replaced by NaNs)
; - converted from ADU to electrons
;


PRO flame_correct_data, fuel=fuel

  ; read in the bad pixel mask
  badpix = readfits(fuel.intermediate_dir + fuel.badpix_filename)
  
  ; identify the science files
  readcol, fuel.science_filelist, science_filenames, format='A'
  
  for i_frame=0,fuel.N_frames-1 do begin
  
    print, 'Correcting frame ', science_filenames[i_frame]

    ; read in science frame
    frame = readfits(science_filenames[i_frame], header, /silent)

    ; CORRECTION 1: non-linearity
    frame_corr1 = frame + 4.155d-6 * frame^2    ; from http://scienceops.lbto.org/sciops_cookbook/luci2-vs-luci1/  (Jan 2016)
    
    ; CORRECTION 2: bad pixels
    frame_corr2 = frame_corr1
    frame_corr2[where(badpix, /null)] = !values.d_nan

    ; CORRECTION 3: convert to electrons per second
    frame_corr3 = frame_corr2 * fuel.gain

    ; change the flux units in the header
    fxaddpar, header, 'BUNIT', 'electrons', ' '
    
    ; save corrected frame
    writefits, (*fuel.corrscience_files)[i_frame], frame_corr3, header
    
  endfor
  
  
END