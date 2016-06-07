PRO flame_initialize, settings=settings, fuel=fuel

  ; create the fuel structure
  fuel = {}

  ; add the user settings
  fuel = create_struct(fuel, settings)

  ; find the Flame directory
  path_to_thisfile = file_which('flame_initialize.pro', /include_current_dir)
  data_dir = flame_util_replace_string(path_to_thisfile, 'pro/flame_initialize.pro', 'data/')

  ; these are additional fields
  init = { $
    N_frames: 0, $
    darktot_filename : settings.intermediate_dir + 'darktot.fits', $
    badpix_filename : settings.intermediate_dir + 'badpix.fits', $
    diagnostics_ps_filename : settings.intermediate_dir + 'diagnostics.ps', $
    diagnostics_txt_filename : settings.intermediate_dir + 'diagnostics.txt', $
    startraces_filename : settings.intermediate_dir + 'startraces.ps', $
    slitim_filename : settings.intermediate_dir + 'slitim.fits', $
    sky_emission_filename : data_dir + 'sky_emission_model.dat', $
    corrscience_files : strarr(file_lines(settings.science_filelist)), $
    slit_coarse_wavecal : 0, $     
    band : '', $
    linelist_filename: '', $
    pixel_scale : 0.0, $
    xrange_star:[1000, 1200], $
    output_lambda_0 : 0d, $
    output_lambda_delta : 0d, $
    output_lambda_Npix : 0 $
    }

  ; add them to fuel
  fuel = create_struct(fuel, init)

  ; if we are reducing only one slit, then that slit will be used for the coarse wavelength calibration
  if settings.reduce_only_oneslit ne 0 then fuel.slit_coarse_wavecal = settings.reduce_only_oneslit-1

  ; check and setup directory structure
  if file_test(fuel.intermediate_dir) eq 0 then spawn, 'mkdir ' + fuel.intermediate_dir
  if file_test(fuel.output_dir) eq 0 then spawn, 'mkdir ' + fuel.output_dir

  ; set the number of science frames
  fuel.N_frames = file_lines(fuel.science_filelist)

  ; create file name for corrected science frames
  readcol, fuel.science_filelist, science_filenames, format='A'
  for i_frame=0, fuel.N_frames-1 do begin
    components = strsplit( science_filenames[i_frame], '/', /extract)
    out_filename = components[-1] ; keep just the filename
    out_filename = flame_util_replace_string( out_filename, '.fits', '_corr.fits' )
    print, i_frame
    fuel.corrscience_files[i_frame] = fuel.intermediate_dir + out_filename
  endfor

  ; determine what band we are in
  hdr = headfits(science_filenames[0])
  fuel.band =  strtrim(fxpar(hdr, 'FILTER2'),2)

  ; choose the file with the line list
  fuel.linelist_filename = data_dir + 'line_list_' + fuel.band + '.dat'

  ; read in the pixel scale 
  fuel.pixel_scale = fxpar(hdr, 'PIXSCALE')  ; arcsec/pixel

  ; set the wavelength grid for the final, resampled spectra
  case fuel.band of

    'J': begin
      fuel.output_lambda_0 = 1.13d           ; wavelength for the first pixel, in micron
      fuel.output_lambda_delta = 0.75d-4     ; micron per pixel
      fuel.output_lambda_Npix =  3200        ; number of pixels (determines wavelength range)
    end

    'K': begin
      fuel.output_lambda_0 = 1.90d           ; wavelength for the first pixel, in micron
      fuel.output_lambda_delta = 1.5d-4      ; micron per pixel
      fuel.output_lambda_Npix =  4200        ; number of pixels (determines wavelength range)
    end

    else: message, fuel.band + ' band not supported yet'

  endcase



  ; add empty pointers for future sub-structures
  fuel = create_struct(fuel, $
    'diagnostics', ptr_new(/allocate_heap), $
    'slits', ptr_new(/allocate_heap) )

END