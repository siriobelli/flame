
FUNCTION flame_create_fuel_loadfiles, filelist
;
; Load file names from the text file "filelist"
; Check for the existence of all files and expand the paths
;

  ; if the input is a meaningful string, then return empty string
  filelist_norm = strlowcase(strtrim( filelist, 2 ))
  if filelist_norm eq '' or filelist_norm eq 'none' or filelist_norm eq 'default' then $
    return, ptr_new()

  ; check that the input file exists
  if ~file_test(filelist) then message, 'file ' + filelist + ' not found!''

  ; read the input file
  readcol, filelist, filenames, format='A'

  ; make sure we don't have to work with relative paths
  filenames = file_expand_path(filenames)

  ; check that all the frames exist
  if ~array_equal(file_test(filenames), 1+intarr(n_elements(filenames))) then $
    message, filelist + ': not all frames exist!'

  ; output the list of file names
  return, filenames


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



FUNCTION flame_util_create_fuel, input
;
; initialize flame and output the fuel structure. The only input is the input structure.
;

  ; find the Flame data directory and check that it exists
  path_to_thisfile = file_which('flame_util_create_fuel.pro', /include_current_dir)
  data_dir = flame_util_replace_string(path_to_thisfile, 'pro/flame_util_create_fuel.pro', 'data/')
  if ~file_test(data_dir, /directory) then message, 'data directory not found! Check the flame directory structure.'

  ; setup directory structure
  if file_test(input.intermediate_dir) eq 0 then file_mkdir, input.intermediate_dir
  if file_test(input.intermediate_dir + '/frames/') eq 0 then file_mkdir, input.intermediate_dir + '/frames/'
  if file_test(input.output_dir) eq 0 then file_mkdir, input.output_dir

  ; read the science filenames
  science_filenames = flame_create_fuel_loadfiles(input.science_filelist)
  if ~keyword_set(science_filenames) then message, 'list of science frames is required!'
  forprint, science_filenames

  ; read the number of science frames
  N_frames = n_elements(science_filenames)
  print, strtrim(N_frames, 2) + ' science frames found'

  ; make filenames for the corrected science frames
  if strmatch(science_filenames[0], '*.fits.gz') then $
    corrscience_basenames = file_basename(science_filenames, '.fits.gz') + '_corr.fits' else $
    corrscience_basenames = file_basename(science_filenames, '.fits') + '_corr.fits'
  corrscience_filenames = file_expand_path(input.intermediate_dir + '/frames/' + corrscience_basenames)

  ; read dark filenames
  filenames_dark = flame_create_fuel_loadfiles(input.dark_filelist)
  if keyword_set(filenames_dark) then print, 'Reading dark frames: ', filenames_dark

  ; read pixel-flat filenames
  filenames_pixelflat = flame_create_fuel_loadfiles(input.pixelflat_filelist)
  if keyword_set(filenames_pixelflat) then print, 'Reading pixel flat frames: ', filenames_pixelflat

  ; read illumination-flat filenames
  filenames_illumflat = flame_create_fuel_loadfiles(input.illumflat_filelist)
  if keyword_set(filenames_illumflat) then print, 'Reading illumination flat frames: ', filenames_illumflat

  ; read arc filenames
  filenames_arc = flame_create_fuel_loadfiles(input.arc_filelist)
  if keyword_set(filenames_arc) then print, 'Reading arc frames: ', filenames_arc

  ; read slit-flat filenames
  filenames_slitflat = flame_create_fuel_loadfiles(input.slitflat_filelist)
  if keyword_set(filenames_slitflat) then print, 'Reading slit flat frames: ', filenames_slitflat $
    else $
    if N_frames GE 3 then filenames_slitflat = corrscience_filenames[fix(N_frames/2)-1:fix(N_frames/2)+1] $
      else filenames_slitflat = corrscience_filenames[N_frames/2]

  ; check that the dither file exists and read it
  dither_file_norm = strlowcase( strtrim( input.dither_file, 2 ))
  if dither_file_norm EQ 'none' or dither_file_norm EQ '' then $
    dither_blind_positions = ptr_new() else begin
      if ~file_test(input.dither_file) then message, 'file ' + input.dither_file + ' not found!'
      readcol, input.dither_file, dither_blind_positions, format='D'
    endelse

  ; create the util substructure
  util = { $
    N_frames: N_frames, $
    science_filenames: science_filenames, $
    corrscience_filenames: corrscience_filenames, $
    filenames_dark: filenames_dark, $
    filenames_pixelflat: filenames_pixelflat, $
    filenames_illumflat: filenames_illumflat, $
    filenames_arc: filenames_arc, $
    filenames_slitflat: filenames_slitflat, $
    master_dark: input.intermediate_dir + 'master_dark.fits', $
    master_pixelflat: input.intermediate_dir + 'master_pixelflat.fits', $
    master_illumflat: input.intermediate_dir + 'master_illumflat.fits', $
    master_arc: input.intermediate_dir + 'master_arc.fits', $
    master_getslit: input.intermediate_dir + 'master_slitflat.fits', $
    dither_blind_positions: dither_blind_positions, $
    slitim_filename: 'slitim.fits', $
    flame_data_dir: data_dir, $
    sky_emission_filename: data_dir + 'sky_emission_model_nir.dat', $
    linelist_filename: data_dir + 'line_list_midres.dat', $
    star_x_range: [1000, 1200], $
    star_y_window: 0, $
    trace_slit_with_skylines : 1, $
    wavecal_rough_R : [500, 1000, 3000], $
    clean_individual_frames : 0, $
    illumination_correction : 1, $
    wavecal_rough_continuum_degree: 0, $
    identify_lines_poly_degree: 5, $
    identify_lines_Nmin_lines: 8, $
    identify_lines_linefit_window: 5.0, $       ; in units of expected linewidth
    wavesolution_degree: 3, $
    skysub_plot_range: [0.4, 0.6], $
    combine_sigma_clip : 2.0, $
    start_time: systime(/seconds), $
    last_routine_time: systime(/seconds), $
    last_routine_name: 'flame_util_create_fuel' $
   }

  ; create the fuel structure
  fuel = { $
    input: input, $
    util: util, $
    instrument: ptr_new(), $
    diagnostics : ptr_new(), $
    slits : ptr_new() $
    }

  return, fuel

END
