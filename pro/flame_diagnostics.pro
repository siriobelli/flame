

FUNCTION flame_monitor_fit_trace, frame, xrange=xrange, yrange=yrange, est_width=est_width, $
    plot_extra=plot_extra
;
; this function fits a Gaussian to the spatial profile and returns 
; the results of the fit in a structure 
;

  ; make cutout where the trace is
  cutout = frame[ xrange[0]:xrange[1] , yrange[0]:yrange[1] ]

  ; median profile (avoid outliers)
  profile = median(cutout, dimension=1)

  ; make vertical axis with the pixel position
  ypixel = yrange[0] + indgen(n_elements(profile))
  
  ; initial estimates for the Gaussian curve
  est_peak = max(profile, i_max, /nan)-median(profile)
  est_center = ypixel[i_max] 
  est_width = est_width
  est_continuum = median(profile)

  ; fit a Gaussian
  fit_result = gaussfit(ypixel, profile, coeff, nterms=4, $
    estimates=[ est_peak, est_center, est_width, est_continuum], $
    chisq=chisq, sigma=coeff_err)

  cgplot, ypixel, profile, /ynozero, $
    xtit='y pixel coordinate', ytit='spatial profile', charsize=1, _extra=plot_extra
  cgplot, ypixel, fit_result, /overplot, color='red'

  ; put the results in a nice structure
  fit_structure = { $
  peak: coeff[0], $
  peak_err: coeff_err[0], $
  center: coeff[1], $
  center_err: coeff_err[1], $
  width: coeff[2], $
  width_err: coeff_err[2], $
  continuum: coeff[3], $
  continuum_err: coeff_err[3], $
  chisq: chisq }

  ; return the structure
  return, fit_structure

END




;****************************************************************



FUNCTION flame_diagnostics_AorB, frame_filename, fuel=fuel
  ;
  ; Takes a raw frame and identifies whether the star trace is in the A 
  ; or in the B position, by trying to fit a Gaussian to the spatial profile.
  ; It outputs 'A', 'B', or 'X' when cannot find a trace
  ; It uses flame_monitor_fit_trace
  ;

  ; read in the frame
  frame = readfits(frame_filename, /silent)

  ; x range in which to fit a Gaussian to the star trace
  xrange=fuel.xrange_star

  ; estimate width of trace in pixels: 0.8 arcsec (compromise between seeing-limited and ARGOS)
  est_seeing = 0.8  ; FWHM in arcsec
  est_width = (est_seeing / 2.355) / fuel.pixel_scale   ; sigma in pixels

  ; determine the vertical range to consider for finding the star trace
  half_range = 10.0 * est_width

  ; A frames (always the one on top)
  A_yrange = max(fuel.startrace_y_pos) + [-half_range, half_range]

  ; B frames
  B_yrange = min(fuel.startrace_y_pos) + [-half_range, half_range]

  ; fit a Gaussian at the A position
  fit_A = flame_monitor_fit_trace( frame, xrange=xrange, yrange=A_yrange, est_width=est_width, $
      plot_extra = {title:frame_filename, layout:[1,2,1]})

  ; fit a Gaussian at the B position
  fit_B = flame_monitor_fit_trace( frame, xrange=xrange, yrange=B_yrange, est_width=est_width, $
      plot_extra = {layout:[1,2,2]})

  ; label the plots
  cgtext, 'A', 0.18, 0.85, /normal
  cgtext, 'B', 0.18, 0.38, /normal

  ; new page for the next plots
  erase

  ; let's define a "good" fit one that satisfies these criteria:
  ; - Gaussian width not larger than 5 times the estimated value
  ; - peak must be positive
  ; - signal/noise ratio of the peak must be larger than 5

  ; check if the A fit makes sense 
  if fit_A.width LT 5.0*est_width and fit_A.peak GT 0.0 and fit_A.peak/fit_A.peak_err GT 5.0 $
    then A_ok = 1 else A_ok = 0

  ; check if the B fit makes sense 
  if fit_B.width LT 5.0*est_width and fit_B.peak GT 0.0 and fit_B.peak/fit_B.peak_err GT 5.0 $
    then B_ok = 1 else B_ok = 0

  ; decide if A or B contains the trace 

  ; if only one of the two makes sense, then it's that one
  if A_ok eq 1 and B_ok eq 0 then return, 'A'
  if A_ok eq 0 and B_ok eq 1 then return, 'B'

  ; if both make sense, take the one with the larger peak
  if A_ok eq 1 and B_ok eq 1 then $
    if fit_A.peak GT fit_B.peak then return, 'A' $
      else return, 'B'

  ; if both are wrong fits, then return X and issue a warning
  if A_ok eq 0 and B_ok eq 0 then begin
    message, frame_filename + ': did not find a valid star trace!', /informational
    return, 'X'
  endif

END



;****************************************************************



FUNCTION flame_diagnostics_fit, frame_filename, sky_filename, offset_pos=offset_pos, fuel=fuel
  ;
  ; Takes a frame, subtracts the sky, and fits a Gaussian 
  ; to the spatial profile at the A or B position.
  ; It outputs a diagnostics structure
  ; It uses flame_monitor_fit_trace
  ;

  ; read frame number from file name
  frame_num = (strsplit(frame_filename, '.', /extract))[-2]

  ; read in the frame
  frame_star = readfits(frame_filename, header, /silent)

  ; read the airmass from the header
  airmass = sxpar(header, 'AIRMASS')

  ; if no trace was identified, then return all NaNs
  if offset_pos eq 'X' then begin
    diagnostics = { $
      frame_num: frame_num, $
      offset_pos:offset_pos, $
      seeing: !values.f_NaN, $
      flux: !values.f_NaN, $
      position: !values.f_NaN, $
      airmass: airmass}
    return, diagnostics
  endif

  ; read in the frame to be used a sky
  frame_sky = readfits(sky_filename, /silent)
  
  ; subtract sky
  frame = frame_star - frame_sky

  ; x range in which to fit a Gaussian to the star trace
  xrange=fuel.xrange_star

  ; estimate width of trace in pixels: 0.8 arcsec (compromise between seeing-limited and ARGOS)
  est_seeing = 0.8  ; FWHM in arcsec
  est_width = (est_seeing / 2.355) / fuel.pixel_scale   ; sigma in pixels

  ; determine the vertical range to consider for finding the star trace
  half_range = 5.0 * est_width

  ; find expected position of star trace
  if offset_pos eq 'A' then ycenter = max(fuel.startrace_y_pos) $
    else ycenter = min(fuel.startrace_y_pos) 

  ; vertical range to be considered for the fit
  yrange = ycenter + [-half_range, half_range]

  ; fit a Gaussian
  fit_result = flame_monitor_fit_trace(frame, xrange=xrange, yrange=yrange, est_width=est_width, $
      plot_extra = {title:frame_num})

  ; make the diagnostics structure for this frame
  diagnostics = { $
    frame_num: frame_num, $
    offset_pos: offset_pos, $
    seeing: 2.355 * fit_result.width * fuel.pixel_scale, $
    flux: sqrt(2.0*3.14) * fit_result.width * fit_result.peak, $
    position: fit_result.center, $
    airmass: airmass}

    return, diagnostics


END


;****************************************************************

FUNCTION flame_diagnostics_fromdata, fuel

;
; used when there is a reference star on the slit 
; the trace is fit with a Gaussian profile in each frame
; and offset, transmission, and seeing are calculated
;


  ; identify the science files
  readcol, fuel.science_filelist, science_filenames, format='A'

  ; 1)  For each frame, decide whether it is an A or B position 
  ;     by looking at which one of the two "windows" is more likely
  ;     to contain a star trace

  ; *********A is always the one on top***************

  ; this string array will keep track of the 'A' and 'B' position. 'X' means undecided
  offset_pos = strarr(fuel.N_frames)
  
  ; identify A and B frames  
  cgPs_open, fuel.intermediate_dir + 'startrace_identify_AB.ps', /nomatch
    for i_frame=0,fuel.N_frames-1 do begin
       offset_pos[i_frame] = flame_diagnostics_AorB( science_filenames[i_frame], fuel=fuel )
       print, i_frame, ' ', science_filenames[i_frame], ' offset position: ', offset_pos[i_frame]
    endfor
  cgPS_close


  ; 2)  Now do a proper fit from the sky-subtracted frame. 
  ;     Need to identify the closest frame with a different offset position
  ;     to use as sky frame

  cgPS_open, fuel.intermediate_dir + 'startraces.ps', /nomatch

  ; array that will contain the diagnostics
  diagnostics = []

  print, 'Fitting the star trace for each frame...'
  
  for i_frame=0, fuel.N_frames-1 do begin
 
    print, 'fitting star trace for ', science_filenames[i_frame] + ' at position ' + offset_pos[i_frame]
    
    ; need to find the closest available frame (starting with the next one)
    distance = abs( i_frame+0.01 - indgen(fuel.N_frames) )
    closest_frames = sort(distance) ; this array contains the frame numbers from the closest to the farthest
    ii = 1
    while offset_pos[closest_frames[ii]] EQ offset_pos[i_frame] do ii++
    i_frame_background = closest_frames[ii]   ; this is the frame to use a background
    
    ; fit a Gaussian to the sky-subtracted frame and obtain diagnostics
    diagnostics_thisframe = flame_diagnostics_fit( science_filenames[i_frame], $
      science_filenames[i_frame_background], offset_pos=offset_pos[i_frame], fuel=fuel )

    ; add these to the total diagnostics
    diagnostics = [ diagnostics , diagnostics_thisframe ]

  endfor
  
  cgPS_close

  return, diagnostics

END


;****************************************************************


FUNCTION flame_diagnostics_blind, fuel

;
; used when there is no reference star on the slit 
; a dither file is needed, where the dither positions
; are given for all frames
;

  ; read pixel scale
  pixel_scale = fuel.pixel_scale    ; arcsec/pixel
  
  ; identify the science files
  readcol, fuel.science_filelist, science_filenames, format='A'

  ; identify the dither file
  if fuel.dither_filelist eq 'none' then message, 'dither_filelist needs to be specified'
  readcol, fuel.dither_filelist, dither_position, format='D'

  ; create a template of the diagnostics structure
  diagnostics_tmp = { $
      frame_num: 0, $
      offset_pos:'', $
      seeing: !values.f_NaN, $
      flux: !values.f_NaN, $
      position: !values.f_NaN, $
      airmass: !values.f_NaN}

  ; array of diagnostics 
  diagnostics = replicate(diagnostics_tmp, fuel.N_frames)

  ; read frame numbers from file names
  frame_num = intarr(fuel.N_frames)
  for i_frame=0, fuel.N_frames-1 do $
    frame_num[i_frame] = (strsplit(science_filenames[i_frame], '.', /extract))[-2]

  ; fill the frame_num field
  diagnostics.frame_num = frame_num

  ; fill the dither position field (in pixels, not arcsec)
  diagnostics.position = dither_position / pixel_scale

  ; detect the As and the Bs
  w_A = where( dither_position GT mean(dither_position, /nan), complement=w_B )
  diagnostics[w_A].offset_pos = 'A'
  diagnostics[w_B].offset_pos = 'B'

  ; read the airmass from the headers
  for i_frame=0, fuel.N_frames-1 do $
    diagnostics[i_frame].airmass = sxpar(headfits(science_filenames[i_frame]), 'AIRMASS')

  return, diagnostics

END


;****************************************************************


PRO flame_diagnostics_plot, diagnostics

;
; make a plot with five panels and show the trend of flux, FWHM,
; vertical position in A frames, vertical position in B frames, 
; and airmass for the star trace as a function of frame number
;

  ; plot parameters
  extra_structure = {noerase:1, xtickformat:'(A1)', charsize:0.7, xsty:1, ynozero:1, psym:-16}
  x0 = 0.1
  x1 = 0.95
  y0 = 0.10
  y1 = 0.95
  delta_y = (y1-y0)/5.0

  ; check if a reference star was measured
  if where( finite(diagnostics.flux), /null) NE !NULL then begin   
    
    ; plot flux
    cgplot, diagnostics.frame_num, diagnostics.flux/median(diagnostics.flux), $
      _extra = extra_structure, xra=xra, $
      ytit='flux / median', position=[x0,y1-delta_y,x1,y1]
    cgplot, [0, 10000], 1+[0,0], /overplot

    ; plot seeing
    cgplot, diagnostics.frame_num, diagnostics.seeing, $
      _extra = extra_structure, xra=xra, $
      ytit='FWHM (arcsec)', position=[x0,y1-2.0*delta_y,x1,y1-delta_y]

  endif

  ; plot position of A frames
  if where(diagnostics.offset_pos eq 'A', /null) ne !NULL then $
    cgplot, diagnostics[where(diagnostics.offset_pos eq 'A', /null)].frame_num, $
      diagnostics[where(diagnostics.offset_pos eq 'A', /null)].position, $
      _extra = extra_structure, xra=xra, $
      ytit='A $\Delta$ y (pixels)', position=[x0,y1-3.0*delta_y,x1,y1-2.0*delta_y]

  ; plot position of B frames
  if where(diagnostics.offset_pos eq 'B', /null) ne !NULL then $
    cgplot, diagnostics[where(diagnostics.offset_pos eq 'B', /null)].frame_num, $
      diagnostics[where(diagnostics.offset_pos eq 'B', /null)].position, $
      _extra = extra_structure, xra=xra, $
      ytit='B $\Delta$ y (pixels)', position=[x0,y1-4.0*delta_y,x1,y1-3.0*delta_y]

  ; plot airmass
  extra_structure.xtickformat = ''
  cgplot, diagnostics.frame_num, diagnostics.airmass, $
    _extra = extra_structure, xra=xra, $
    ytit='airmass', xtit='frame number', position=[x0,y1-5.0*delta_y,x1,y1-4.0*delta_y]
  
  

END


;****************************************************************


PRO flame_diagnostics, fuel=fuel

  ; 1 - create the diagnostics structure
  ;----------------------------------------

  ; check if a reference star position has been set
  if fuel.startrace_y_pos[0] GT 0.0 and fuel.startrace_y_pos[1] GT 0.0 then $

    ; if a valid reference star position is given, then monitor the star
    diagnostics = flame_diagnostics_fromdata(fuel) $
  
  else $

    ; otherwise use the dither file to get the offset position for each frame
    diagnostics = flame_diagnostics_blind(fuel)


  ; 2 - plot diagnostics
  ;----------------------------------------

  cgPS_open, fuel.intermediate_dir + 'diagnostics.ps', /nomatch
    flame_diagnostics_plot, diagnostics
  cgPS_close
  

  ; 3 - write text file with diagnostics
  ;----------------------------------------
  
  forprint, diagnostics.frame_num, '    ' + '    ' + cgnumber_formatter(diagnostics.offset_pos, decimals=2), $
    '    ' + cgnumber_formatter(diagnostics.seeing, decimals=2), '    ' + cgnumber_formatter(diagnostics.flux, decimals=1), $
    '    ' + cgnumber_formatter(diagnostics.flux/median(diagnostics.flux), decimals=2), '    ' + cgnumber_formatter(diagnostics.position, decimals=1), $
    textout=fuel.intermediate_dir + 'diagnostics.txt', $
    comment = '# frame number   offset    seeing    flux    normalized flux    position '


  ; 4 - add diagnostics to the fuel structure
  ;----------------------------------------

  *fuel.diagnostics = diagnostics
    
    
END