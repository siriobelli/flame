

FUNCTION flame_create_fuel_loadfiles, filelist, corr_dir=corr_dir
;
; Load file names from the text file "filelist"
; Check for the existence of all files and expand the paths
; Create the filenames for the corrected frames
; output structure with raw filenames, corrected filenames, and master file name
;

  ; if the input is empty or 'none' or 'default', then return 'zero' 'structure
  filelist_norm = strlowcase(strtrim( filelist, 2 ))
  if filelist_norm eq '' or filelist_norm eq 'none' or filelist_norm eq 'default' then $
    return, { $
      n_frames:0, $
      raw_files:'', $
      corr_files:'', $
      master_file:''}

  ; check that the input file exists
  if ~file_test(filelist) then message, 'file ' + filelist + ' not found!''

  ; read the input file
  readcol, filelist, filenames, format='A'

  ; make sure we don't have to work with relative paths
  filenames = file_expand_path(filenames)

  ; check that all the frames exist
  print, 'Testing that files exist: ', filenames
  if ~array_equal(1, file_test(filenames)) then $
    message, filelist + ': not all frames exist!'

  ; make filenames for the corrected frames
  if strmatch(filenames[0], '*.fits.gz') then $
    corr_basenames = file_basename(filenames, '.fits.gz') + '_corr.fits' else $
    corr_basenames = file_basename(filenames, '.fits') + '_corr.fits'

  ; put them in the directory for corrected frames
  if keyword_set(corr_dir) then corr_filenames = corr_dir + corr_basenames $
    else corr_filenames = strarr(n_elements(filenames))

    ; output a structure containing the filenames
  return, { $
    n_frames:n_elements(filenames), $
    raw_files:filenames, $
    corr_files:corr_filenames, $
    master_file:''}


END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


FUNCTION flame_util_create_fuel, input
;
; initialize flame and output the fuel structure. The only input is the input structure.
;


  ; set up directories ---------------------------------------------------------

  ; find the flame data directory and check that it exists
  path_to_thisfile = file_which('flame_util_create_fuel.pro', /include_current_dir)
  data_dir = flame_util_replace_string(path_to_thisfile, 'pro/flame_util_create_fuel.pro', 'data/')
  if ~file_test(data_dir, /directory) then message, 'data directory not found! Check the flame directory structure.'

  ; intermediate directory
  if file_test(input.intermediate_dir) eq 0 then file_mkdir, input.intermediate_dir
  intermediate_dir = file_expand_path(input.intermediate_dir) + '/'

  ; frames directory
  if file_test(intermediate_dir + 'frames/') eq 0 then file_mkdir, intermediate_dir + 'frames/'
  frames_dir = intermediate_dir + 'frames/'

  ; output directory
  if file_test(input.output_dir) eq 0 then file_mkdir, input.output_dir
  output_dir = file_expand_path(input.output_dir) + '/'


  ; read in filenames ----------------------------------------------------------

  ; read science filenames
  print, ''
  print, 'Science frames'
  science = flame_create_fuel_loadfiles(input.science_filelist, corr_dir=frames_dir)
  print, science.n_frames, ' files read.'

  ; read dark filenames
  print, ''
  print, 'Dark frames'
  calib_dark = flame_create_fuel_loadfiles(input.dark_filelist, corr_dir=frames_dir)
  print, calib_dark.n_frames, ' files read.'
  calib_dark.master_file = intermediate_dir + 'master_dark.fits'

  ; read arc filenames
  print, ''
  print, 'Arc frames'
  calib_arc = flame_create_fuel_loadfiles(input.arc_filelist, corr_dir=frames_dir)
  print, calib_arc.n_frames, ' files read.'
  calib_arc.master_file = intermediate_dir + 'master_arc.fits'

  ; read pixelflat filenames
  print, ''
  print, 'Pixel flat frames'
  calib_pixelflat = flame_create_fuel_loadfiles(input.pixelflat_filelist, corr_dir=frames_dir)
  print, calib_pixelflat.n_frames, ' files read.'
  calib_pixelflat.master_file = intermediate_dir + 'master_pixelflat.fits'

  ; read illumflat filenames
  print, ''
  print, 'Illumination flat frames'
  calib_illumflat = flame_create_fuel_loadfiles(input.illumflat_filelist, corr_dir=frames_dir)
  print, calib_illumflat.n_frames, ' files read.'
  calib_illumflat.master_file = intermediate_dir + 'master_illumflat.fits'

  ; read slitflat filenames
  print, ''
  print, 'Slit flat frames'
  calib_slitflat = flame_create_fuel_loadfiles(input.slitflat_filelist, corr_dir=frames_dir)
  print, calib_slitflat.n_frames, ' files read.'
  calib_slitflat.master_file = intermediate_dir + 'master_slitflat.fits'


  ; read in dither file ----------------------------------------------------------

  ; check that the dither file exists and read it
  dither_file_norm = strlowcase( strtrim( input.dither_file, 2 ))
  if dither_file_norm EQ 'none' or dither_file_norm EQ '' then $
    dither_blind_positions = ptr_new() else begin
      if ~file_test(input.dither_file) then message, 'file ' + input.dither_file + ' not found!'
      readcol, input.dither_file, dither_blind_positions, format='D'
    endelse

  ; if we are nodding, then we want to discard pixels with data from few frames.
  ; otherwise the ABcombined frame will be noisy
  if input.AB_subtraction then combine_min_framefrac = 0.49 else combine_min_framefrac = 0.2

  ; create the fuel structure ----------------------------------------------------------

  ; create the util substructure
  util = { $
    science: science, $
    dark: calib_dark, $
    arc: calib_arc, $
    pixelflat: calib_pixelflat, $
    illumflat: calib_illumflat, $
    slitflat: calib_slitflat, $
    dither_blind_positions: dither_blind_positions, $
    flame_data_dir: data_dir, $
    intermediate_dir: intermediate_dir, $
    output_dir: output_dir, $
    start_time: systime(/seconds), $
    last_routine_time: systime(/seconds), $
    last_routine_name: 'flame_util_create_fuel' $
    }


  ; create the settings substructure
  settings = { $
    sky_emission_filename: data_dir + 'sky_emission_model_nir.dat', $
    linelist_filename: data_dir + 'sky_line_list_R3000.dat', $
    star_x_range: [1000, 1200], $
    star_y_window: 0, $
    clean_individual_frames : 0, $
    trace_slit_with_skylines : 0, $
    trace_slit_ymargin: 12, $
    trace_slit_polydegree: 3, $
    wavecal_rough_R : [500, 1000, 3000], $
    illumination_correction : 1, $
    wavecal_rough_smooth_window: 20, $
    wavecal_rough_split: 1, $
    identify_lines_stack_rows: 0, $
    identify_lines_poly_degree: 5, $
    identify_lines_Nmin_lines: 6, $
    identify_lines_linefit_window: 6.0, $       ; in units of expected linewidth
    wavesolution_order_x: 3, $
    wavesolution_order_y: 2, $
    shift_arcs_Nmin_lines: 1, $
    skysub_plot_range: [0.4, 0.6], $
    skysub_bspline_oversample: 1.0, $
    skysub_reject_fraction: 0.5, $
    skysub_reject_window: 2.0, $
    combine_sigma_clip : 2.0, $
    combine_min_framefrac : combine_min_framefrac, $
    debugging:0 $
   }


  ; create the fuel structure
  fuel = { $
    input: input, $
    settings: settings, $
    util: util, $
    instrument: ptr_new(), $
    diagnostics : ptr_new(), $
    slits : ptr_new() $
    }


  return, fuel

END
