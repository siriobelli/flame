FUNCTION flame_create_fuel, input
;
; initialize flame and output the fuel structure. The only input is the input structure.
;

  ; find the Flame data directory and check that it exists
  path_to_thisfile = file_which('flame_create_fuel.pro', /include_current_dir)
  data_dir = flame_util_replace_string(path_to_thisfile, 'pro/flame_create_fuel.pro', 'data/')
  if ~file_test(data_dir, /directory) then message, 'data directory not found. Check the flame directory structure.'

  ; check and setup directory structure
  if file_test(input.intermediate_dir) eq 0 then spawn, 'mkdir ' + input.intermediate_dir
  if file_test(input.output_dir) eq 0 then spawn, 'mkdir ' + input.output_dir

  ; check that the science file exists
  if ~file_test(input.science_filelist) then message, 'file ' + input.science_filelist + ' not found'

  ; read the science file
  readcol, input.science_filelist, science_filenames, format='A'

  ; read the number of science frames
  N_frames = n_elements(science_filenames)

  ; check the all the science frames exist
  if ~array_equal(file_test(science_filenames), 1+intarr(N_frames)) then message, 'not all science frames exist'

  ; create file names for corrected science frames
  corrscience_filenames = strarr(N_frames)
  for i_frame=0, N_frames-1 do begin
    components = strsplit( science_filenames[i_frame], '/', /extract)
    out_filename = components[-1] ; keep just the filename
    out_filename = flame_util_replace_string( out_filename, '.fits', '_corr.fits' )
    print, i_frame
    corrscience_filenames[i_frame] = input.intermediate_dir + out_filename
  endfor

  ; check that the darks file exists and read it
  if strlowcase(input.darks_filelist) EQ 'none' then darks_filenames = '' else begin
    if ~file_test(input.darks_filelist) then message, 'file ' + input.darks_filelist + ' not found'
    readcol, input.darks_filelist, darks_filenames, format='A'
  endelse 

  ; check that the flats file exists and read it
  if strlowcase(input.flats_filelist) EQ 'none' then flats_filenames = '' else begin
    if ~file_test(input.flats_filelist) then message, 'file ' + input.flats_filelist + ' not found'
    readcol, input.flats_filelist, flats_filenames, format='A'
  endelse

  ; check that the dither file exists and read it
  if strlowcase(input.dither_filelist) EQ 'none' then dither_blind_positions = !values.d_NaN else begin
    if ~file_test(input.dither_filelist) then message, 'file ' + input.dither_filelist + ' not found'
    readcol, input.dither_filelist, dither_blind_positions, format='D'
  endelse

  ; create the util substructure
  util = { $
    N_frames: N_frames, $
    science_filenames: science_filenames, $
    corrscience_filenames: corrscience_filenames, $
    darks_filenames: darks_filenames, $
    flats_filenames: flats_filenames, $
    dither_blind_positions: dither_blind_positions, $	
    slitim_filename: 'slitim.fits', $
    flame_data_dir : data_dir, $
    sky_emission_filename : data_dir + 'sky_emission_model.dat', $
    linelist_filename: data_dir + 'line_list.dat' $    
   }

  ; create the fuel structure
  fuel = { $
    input: input, $
    util: util, $
    instrument: ptr_new(/allocate_heap), $
    diagnostics : ptr_new(/allocate_heap), $
    slits_fromheader : ptr_new(/allocate_heap), $
    slits : ptr_new(/allocate_heap) $
    }
      
return, fuel

END