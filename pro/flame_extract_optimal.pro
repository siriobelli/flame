

FUNCTION flame_util_detect_trace, spec2d

  ; integrate 2d spectrum
  profile = total(spec2d, 1, /nan)
  x_profile = indgen(n_elements(profile))
  
  cgplot, x_profile, median(profile,3)

  ; find peak of the profile around the expected location
  peak_profile = max(median(profile,3), ind_peak)
  
  ; cut away non positive part right and left of the peak
  w_negative_plus = where(x_profile GT ind_peak and profile LE 0.)
  w_negative_minus = where(x_profile LT ind_peak and profile LE 0.)
  
  ; calculate width of the positive peak
  width = w_negative_plus[0] - w_negative_minus[-1]
  
  ; this parameter sets how many pixels beyond the edge we will include in the trace
  extra_edge = 0
  
  ; define the range of pixels within the trace, with generous space at the edges (set by extra_edge)
  w_trace = where( x_profile GT w_negative_minus[-1] - extra_edge AND x_profile LT w_negative_plus[0] + extra_edge , /null)
  
  ; if nothing was found, take a fixed range
;  if w_trace eq !NULL then w_trace = where( x_profile GT spatial_position - 8 and x_profile LT spatial_position )
  
  ; return the spatial pixels that are part of the trace
  return, w_trace
  
END


; **************************************************

FUNCTION flame_util_calculate_weights, im, w_trace
  ;
  ; since the mosfire frames are rectified, I'm going to use just the profile as weight (i.e. weight equals the spectrum collapsed on the lambda axis)
  ; note that the normalization is important
  ;

  ; number of spatial pixels
  Nspace = n_elements(w_trace)
  
  ; number of wavelength pixels
  Nlambda = n_elements(im[*,0])
  
  ; create 2D weight map (only in the trace region)
  weights = dblarr(Nlambda,Nspace)
  
  ; cut out the part of the image where the trace of the object is
  profile = im[*,w_trace]
  profile_ivar = profile*0.0 + 1.0
  profile[where(profile_ivar eq 0., /null)] = !values.f_nan
  
  ; calculate the median profile that will be used as weight for each wavelength
  median_profile = median( profile, dimension=1 )
  
  ; weights cannot be negative!
  median_profile[ where(median_profile LT 0.) ] = 0.
  
  ; normalize the profile so that the at each wavelength the sum over all spatial pixels gives 1
  median_profile /= total(median_profile, /nan)
  
  ; weight map is simply proportional to the spatial profile (actually median profile, not average, to discard outliers)
  for ilambda=0,Nlambda-1 do weights[ilambda,*] = median_profile
  
  ; now extend the weights image to match the size of the whole slit by adding lots of zeros
  w2d = im
  w2d[*,*]=0.
  w2d[*,w_trace] = weights
  
  ; output the 2s weight map
  return, w2d
  
END

; **************************************************





PRO flame_extract_optimal
  
    filename = '/data/luci/pipeline/output/slit07-11494_ABcombined.fits'
  
    ; read in 2d spectrum
    spec2d = readfits(filename, header_eps)
    
    ; ; read in 2d error spectrum (in e-/second)
    ; err2d = readfits(spec2d_sig_filename, header_sig)
    
    ; ; create ivar image
    ; ivar2d = 1d/err2d^2
    
    ; read in wavelength axis ************************** need header_lambda.pro
    lambda = header_lambda(filename, crval2=crval2)
    
    spectrum = total(spec2d[*,18:24], 2)
    lambda_rf = lambda / 3.09
  
  cgPS_open, 'monster_luci.ps', /nomatch
	cgplot, lambda_rf*1d4, median(spectrum,25), yra=[0,10], $
		xtitle = 'rest-frame lambda (angstrom)'
	cgPS_close
  
    ; detect the positive trace
    w_trace = flame_util_detect_trace(spec2d)
  
  
  
    ; calculate the weights for optimal extraction
    weights = flame_util_calculate_weights(spec2d, w_trace)
    
    ; calculate weight profile
    spatial_axis = indgen(n_elements(weights[0,*]))
    weight_profile = total(weights,1) / max(total(weights,1))
    
    ; plot profile
    cgPS_open, outdir + 'profile.' + mask + '.0' + STRING(slit, FORMAT='(I2.2)') + '.' + name + '.ps', /nomatch
    mosred_util_plot_profile, spec2d, w_trace, spatial_position, spatial_axis, weight_profile
    cgPS_close
    
    ; perform optimal extraction
    ;spec1d_opt = mosred_util_optimal_extract(spec2d, ivar2d, weights, ivar1d=ivar1d_opt)
    
    
    ;-------------------------
    ; OUTPUT
    
    ;++++++++++
    sky_box = 1./sqrt(ivar1d_opt)
    sky_box[ where(~finite(sky_box)) ] = !values.d_NaN
    ;++++++++++
    
    ; write specpro files
    
    str1d = {flux:spec1d_opt, ivar:ivar1d_opt, lambda:lambda, sky:sky_box}
    mwrfits, str1d, outdir + 'spec1d.' + mask + '.0' + STRING(slit, FORMAT='(I2.2)') + '.' + name + '.fits', /create;, head
    
    str2d = { flux:spec2d,  lambda:lambda#replicate(1,n_elements(spec2d[0,*])) , ivar:ivar2d }
    mwrfits, str2d, outdir + 'spec2d.' + mask + '.0' + STRING(slit, FORMAT='(I2.2)') + '.' + name + '.fits', /create;, head
    
    ; write the INFO file that will help in the visualization of the 2d spectrum in specpro
    weight_profile /= total(weight_profile)
    y_mean = total( spatial_axis * weight_profile )
    y_width = sqrt( total( ( spatial_axis-y_mean)^2 * weight_profile ) )
    writecol, outdir + 'info.' + mask + '.0' +  STRING(slit, FORMAT='(I2.2)') + '.' + name + '.dat', [ 'ID' , 'extractpos', 'extractwidth' ] , [ 0 , y_mean, 3*y_width ]
    
    ; write also fits file with weights, just to visualize the extraction
    mwrfits, weights, outdir + 'weights2d.' + mask + '.0' + STRING(slit, FORMAT='(I2.2)') + '.' + name + '.fits', /create
  
  END