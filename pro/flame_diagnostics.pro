
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

  ; let's define a "good" fit one that satisfies these criteria:
  ; - Gaussian width not larger than 5 times the estimated value
  ; - peak must be positive
  ; - peak is within the fitted range
  ; - signal/noise ratio of the peak must be larger than 5

  ; read in the frame
  frame = readfits(frame_filename, /silent)

  ; x range in which to fit a Gaussian to the star trace
  xrange=fuel.settings.star_x_range

  ; estimate width of trace in pixels: starting guess is 0.8 arcsec
  est_seeing = 0.8  ; FWHM in arcsec
  est_width = (est_seeing / 2.355) / fuel.instrument.pixel_scale   ; sigma in pixels

  ; determine the vertical range to consider for finding the star trace
  if fuel.settings.star_y_window eq 0 then $
    half_range = 6.0 * est_width else $
    half_range = fuel.settings.star_y_window / 2

  ; let's see if there is a star in the A position
  if fuel.input.star_y_A ne 0 then begin

    ; fit a Gaussian at the A position
    A_yrange = fuel.input.star_y_A + [-half_range, half_range]
    fit_A = flame_monitor_fit_trace( frame, xrange=xrange, yrange=A_yrange, est_width=est_width, $
        plot_extra = {title:frame_filename, layout:[1,2,1]})

    ; label the plot
    cgtext, 'A', 0.18, 0.85, /normal

    ; check if the A fit makes sense
    if fit_A.width LT 5.0*est_width and fit_A.peak GT 0.0 and $
      fit_A.center GE A_yrange[0] and fit_A.center LE A_yrange[1] and $
      fit_A.peak/fit_A.peak_err GT 5.0 $
      then A_ok = 1 else A_ok = 0

  endif else A_ok = 0

  ; let's see if there is a star in the B position
  if fuel.input.star_y_B ne 0 then begin

    ; fit a Gaussian at the B position
    B_yrange = fuel.input.star_y_B + [-half_range, half_range]
    fit_B = flame_monitor_fit_trace( frame, xrange=xrange, yrange=B_yrange, est_width=est_width, $
        plot_extra = {layout:[1,2,2]})

    ; label the plot
    cgtext, 'B', 0.18, 0.38, /normal

    ; check if the B fit makes sense
    if fit_B.width LT 5.0*est_width and fit_B.peak GT 0.0 and $
     fit_B.center GE B_yrange[0] and fit_B.center LE B_yrange[1] and $
     fit_B.peak/fit_B.peak_err GT 5.0 $
     then B_ok = 1 else B_ok = 0

  endif else B_ok = 0


  ; new page for the next plots
  erase

  ; decide if A or B contains the trace

  ; if only one of the two makes sense, then it's that one
  if A_ok eq 1 and B_ok eq 0 then return, 'A'
  if A_ok eq 0 and B_ok eq 1 then return, 'B'

  ; if both make sense, take the one with the largest peak
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
  ; Takes a frame, subtracts the sky if needed, and fits a Gaussian
  ; to the spatial profile at the A or B position.
  ; It outputs a diagnostics structure
  ; It uses flame_monitor_fit_trace
  ;

  ; read frame number from file name
  if strmatch(frame_filename, '*.gz') then $
    frame_num = (strsplit(frame_filename, '._-', /extract))[-3] else $
    frame_num = (strsplit(frame_filename, '._-', /extract))[-2]

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

  ; do we have to subtract the sky?
  if fuel.input.AB_subtraction then begin

    ; read in the frame to be used a sky
    frame_sky = readfits(sky_filename, /silent)

    ; subtract sky
    frame = frame_star - frame_sky

  endif else frame = frame_star

  ; x range in which to fit a Gaussian to the star trace
  xrange=fuel.settings.star_x_range

  ; estimate width of trace in pixels: 0.8 arcsec (compromise between seeing-limited and ARGOS)
  est_seeing = 0.8  ; FWHM in arcsec
  est_width = (est_seeing / 2.355) / fuel.instrument.pixel_scale   ; sigma in pixels

  ; determine the vertical range to consider for finding the star trace
  half_range = 5.0 * est_width

  ; find expected position of star trace
  if offset_pos eq 'A' then ycenter = fuel.input.star_y_A $
    else ycenter = fuel.input.star_y_B

  ; vertical range to be considered for the fit
  yrange = ycenter + [-half_range, half_range]

  ; fit a Gaussian
  fit_result = flame_monitor_fit_trace(frame, xrange=xrange, yrange=yrange, est_width=est_width, $
      plot_extra = {title:frame_num})

  ; make the diagnostics structure for this frame
  diagnostics = { $
    frame_num: frame_num, $
    offset_pos: offset_pos, $
    seeing: 2.355 * fit_result.width * fuel.instrument.pixel_scale, $
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


  ; 1)  For each frame, decide whether it is an A or B position
  ;     by looking at which one of the two "windows" is more likely
  ;     to contain a star trace

  ; this string array will keep track of the 'A' and 'B' position. 'X' means undecided (i.e. sky)
  offset_pos = strarr(fuel.util.science.n_frames)

  ; identify A and B frames
  cgPs_open, fuel.util.intermediate_dir + 'startrace_identify_AB.ps', /nomatch
    for i_frame=0,fuel.util.science.n_frames-1 do begin
       offset_pos[i_frame] = flame_diagnostics_AorB( fuel.util.science.raw_files[i_frame], fuel=fuel )
       print, i_frame, ' ', fuel.util.science.raw_files[i_frame], ' offset position: ', offset_pos[i_frame]
    endfor
  cgPS_close


  ; 2)  Now do a proper fit, with sky-subtraction if needed.
  ;     Need to identify the closest frame with a different offset position
  ;     to use as sky frame

  cgPS_open, fuel.util.intermediate_dir + 'startraces.ps', /nomatch

  ; array that will contain the diagnostics
  diagnostics = []

  print, 'Fitting the star trace for each frame...'

  for i_frame=0, fuel.util.science.n_frames-1 do begin

    ; if doing an A-B subtraction, then needs to find the sky frame
    if fuel.input.AB_subtraction then begin

      ; for AB subtraction, you need at least two types of offset positions
      if array_equal(offset_pos, offset_pos[0]) then $
        message, 'input.AB_subtraction is set, but all frames are taken in the same offset position!'

      ; need to find the closest available frame (starting with the next one)
      distance = abs( i_frame+0.01 - indgen(fuel.util.science.n_frames) )
      closest_frames = sort(distance) ; this array contains the frame numbers from the closest to the farthest
      ii = 1
      while offset_pos[closest_frames[ii]] EQ offset_pos[i_frame] do ii++
      i_frame_background = closest_frames[ii]   ; this is the frame to use a background

    ; if not doing an A-B subtraction, then the background frame is not meaningful
    endif else i_frame_background = i_frame

    ; fit a Gaussian and obtain diagnostics
    diagnostics_thisframe = flame_diagnostics_fit( fuel.util.science.raw_files[i_frame], $
      fuel.util.science.raw_files[i_frame_background], offset_pos=offset_pos[i_frame], fuel=fuel )

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
  pixel_scale = fuel.instrument.pixel_scale    ; arcsec/pixel

  ; check that dither blind positions were specified
  if ~finite(fuel.util.dither_blind_positions[0]) then message, 'dither_file needs to be specified'

  ; create a template of the diagnostics structure
  diagnostics_tmp = { $
      frame_num: 0, $
      offset_pos:'', $
      seeing: !values.f_NaN, $
      flux: !values.f_NaN, $
      position: !values.f_NaN, $
      airmass: !values.f_NaN}

  ; array of diagnostics
  diagnostics = replicate(diagnostics_tmp, fuel.util.science.n_frames)

  ; read frame numbers from file names
  frame_num = intarr(fuel.util.science.n_frames)
  for i_frame=0, fuel.util.science.n_frames-1 do $
    frame_num[i_frame] = (strsplit(fuel.util.science.raw_files[i_frame], '.', /extract))[-2]

  ; fill the frame_num field
  diagnostics.frame_num = frame_num

  ; fill the dither position field (in pixels, not arcsec)
  diagnostics.position = fuel.util.dither_blind_positions / pixel_scale

  ; detect the As and the Bs
  w_A = where( fuel.util.dither_blind_positions GT mean(fuel.util.dither_blind_positions, /nan), complement=w_B )
  diagnostics[w_A].offset_pos = 'A'
  diagnostics[w_B].offset_pos = 'B'

  ; read the airmass from the headers
  for i_frame=0, fuel.util.science.n_frames-1 do $
    diagnostics[i_frame].airmass = sxpar(headfits(fuel.util.science.raw_files[i_frame]), 'AIRMASS')

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
  extra_structure = {noerase:1, xtickformat:'(A1)', charsize:0.9, xsty:1, ynozero:1, psym:-16}
  x0 = 0.15
  x1 = 0.95
  y0 = 0.07
  y1 = 0.95
  delta_y = (y1-y0)/5.0

  ; number of frames
  Nfr = n_elements(diagnostics)

  ; sequential number for frames: 0, 1, 2....
  frame_seqnum = indgen(Nfr)

  ; set the range for the x axis
  xra=[0, Nfr-1]

  ; these are the values used to label the x axis
  xtickname = strtrim(fix(diagnostics.frame_num), 2)
  xtickv = frame_seqnum
  xminor = 1

  ; if there are too many frames, only label some of them
  while n_elements(xtickv) GT 18 do begin

    ; make the array of subindices used to select every other element
    if n_elements(xtickname) mod 2 eq 0 then $
      subset = 2*indgen(n_elements(xtickname)/2) else $   ; if odd
      subset = 2*indgen((1+n_elements(xtickname))/2)      ; if even, select also the last element

    ; keep only every other element for the labeling of the axis
    xtickname = xtickname[ subset ]
    xtickv = xtickv[ subset ]

    ; therefore, need to double the number of minor tick marks between two major marks
    xminor *= 2

  endwhile


  ; check if a reference star was measured
  if where( finite(diagnostics.flux), /null) NE !NULL then begin

    ; plot flux
    cgplot, frame_seqnum, diagnostics.flux/median(diagnostics.flux), $
      _extra = extra_structure, xra=xra, /xsty, $
      ytit='flux / median', position=[x0,y1-delta_y,x1,y1], $
      xtickv = xtickv, xticks=n_elements(xtickv)-1, xminor=xminor, xticklen=0.04

    cgplot, [0, 10000], 1+[0,0], /overplot

    ; plot seeing
    cgplot, frame_seqnum, diagnostics.seeing, $
      _extra = extra_structure, xra=xra, /xsty, $
      ytit='FWHM (arcsec)', position=[x0,y1-2.0*delta_y,x1,y1-delta_y], $
      xtickv = xtickv, xticks=n_elements(xtickv)-1, xminor=xminor, xticklen=0.04

  endif

  ; plot position of A frames
  if where(diagnostics.offset_pos eq 'A', /null) ne !NULL then $
    cgplot, frame_seqnum[where(diagnostics.offset_pos eq 'A', /null)], $
      diagnostics[where(diagnostics.offset_pos eq 'A', /null)].position, $
      _extra = extra_structure, xra=xra, /xsty, $
      ytit='A y-position (pixels)', position=[x0,y1-3.0*delta_y,x1,y1-2.0*delta_y], $
      xtickv = xtickv, xticks=n_elements(xtickv)-1, xminor=xminor, xticklen=0.04

  ; plot position of B frames
  if where(diagnostics.offset_pos eq 'B', /null) ne !NULL then $
    cgplot, frame_seqnum[where(diagnostics.offset_pos eq 'B', /null)], $
      diagnostics[where(diagnostics.offset_pos eq 'B', /null)].position, $
      _extra = extra_structure, xra=xra, /xsty, $
      ytit='B y-position (pixels)', position=[x0,y1-4.0*delta_y,x1,y1-3.0*delta_y], $
      xtickv = xtickv, xticks=n_elements(xtickv)-1, xminor=xminor, xticklen=0.04

  ; plot airmass
  extra_structure.xtickformat = ''
  cgplot, frame_seqnum, diagnostics.airmass, $
    _extra = extra_structure, xra=xra, /xsty, $
    ytit='airmass', xtit='frame number', position=[x0,y1-5.0*delta_y,x1,y1-4.0*delta_y], $
    xtickv = xtickv, xticks=n_elements(xtickv)-1, xtickname=xtickname, xminor=xminor, xticklen=0.04

END


;****************************************************************


PRO flame_diagnostics, fuel

  flame_util_module_start, fuel, 'flame_diagnostics'


  ; 1 - create the diagnostics structure
  ;----------------------------------------

  ; check if a reference star position has been set
  if fuel.input.star_y_A GT 0.0 or fuel.input.star_y_B GT 0.0 then $

    ; if a valid reference star position is given, then monitor the star
    diagnostics = flame_diagnostics_fromdata(fuel) $

  else $

    ; otherwise use the dither file to get the offset position for each frame
    diagnostics = flame_diagnostics_blind(fuel)


  ; 2 - plot diagnostics
  ;----------------------------------------

  cgPS_open, fuel.util.intermediate_dir + 'diagnostics.ps', /nomatch, ysize=10, xsize=6, /cm
    flame_diagnostics_plot, diagnostics
  cgPS_close


  ; 3 - write text file with diagnostics
  ;----------------------------------------

  forprint, diagnostics.frame_num, '    ' + '    ' + cgnumber_formatter(diagnostics.offset_pos, decimals=2), $
    '    ' + cgnumber_formatter(diagnostics.seeing, decimals=2), '    ' + cgnumber_formatter(diagnostics.flux, decimals=1), $
    '    ' + cgnumber_formatter(diagnostics.flux/median(diagnostics.flux), decimals=2), '    ' + cgnumber_formatter(diagnostics.position, decimals=1), $
    textout=fuel.util.intermediate_dir + 'diagnostics.txt', $
    comment = '# frame number   offset    seeing    flux    normalized flux    position '


  ; 4 - add diagnostics to the fuel structure
  ;----------------------------------------

  new_fuel = { input:fuel.input, settings:fuel.settings, util:fuel.util, instrument:fuel.instrument, diagnostics:diagnostics, slits:fuel.slits }
  fuel=new_fuel


  ; 5 - print out info
  ;----------------------------------------

  N_A = n_elements( where(diagnostics.offset_pos eq 'A', /null) )
  N_B = n_elements( where(diagnostics.offset_pos eq 'B', /null) )
  N_X = n_elements( where(diagnostics.offset_pos eq 'X', /null) )

  print, ''
  if N_A ne 0 then print, strtrim(N_A, 2) + ' frames in the A offset position'
  if N_B ne 0 then print, strtrim(N_B, 2) + ' frames in the B offset position'
  if N_X ne 0 then print, strtrim(N_X, 2) + ' frames in the X offset position (sky)'

  ; for AB subtraction, you need at least two types of offset positions
  if fuel.input.AB_subtraction then $
    if N_A*N_B eq 0 and N_A*N_X eq 0 and N_B*N_X eq 0 then $
      message, 'input.AB_subtraction is set, but all frames are taken in the same offset position!'

  ; if you have A and B then you should not have X
  if N_A ne 0 and N_B ne 0 and N_X ne 0 then begin
    print, ''
    print, '**************************'
    print, '***       WARNING      ***'
    print, '**************************'
    print, 'A and B frames are present.'
    print, 'The frames marked X should not be present.'
    print, 'Please check the offset positions and the input frames.'
  endif

  flame_util_module_end, fuel

END
