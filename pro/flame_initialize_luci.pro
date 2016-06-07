PRO flame_initialize_luci, fuel=fuel
  ;
  ; LUCI-specific routine that initializes the fuel structure after
  ; the user has changed the default values.
  ;

  ; check and setup directory structure
  if file_test(fuel.intermediate_dir) eq 0 then spawn, 'mkdir ' + fuel.intermediate_dir
  if file_test(fuel.output_dir) eq 0 then spawn, 'mkdir ' + fuel.output_dir

  ; set the number of science frames
  fuel.N_frames = file_lines(fuel.science_filelist)

  ; create file name for corrected science frames
  readcol, fuel.science_filelist, science_filenames, format='A'
  corrscience_files = strarr(fuel.N_frames)
  for i_frame=0, fuel.N_frames-1 do begin
    components = strsplit( science_filenames[i_frame], '/', /extract)
    out_filename = components[-1] ; keep just the filename
    out_filename = flame_util_replace_string( out_filename, '.fits', '_corr.fits' )
    print, i_frame
    corrscience_files[i_frame] = fuel.intermediate_dir + out_filename
  endfor

  ; save it to fuel
  *fuel.corrscience_files = corrscience_files

  ; determine what band we are in
  hdr = headfits(science_filenames[0])
  fuel.band =  strtrim(fxpar(hdr, 'FILTER2'),2)

  ; choose the file with the line list
  fuel.linelist_filename = fuel.flame_data_dir + 'line_list_' + fuel.band + '.dat'

  ; read in the pixel scale 
  fuel.pixel_scale = fxpar(hdr, 'PIXSCALE')  ; arcsec/pixel

  ; set the wavelength grid for the final, resampled spectra (unless it's already been specified)
  if fuel.output_lambda_0 eq 0d then $
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


END