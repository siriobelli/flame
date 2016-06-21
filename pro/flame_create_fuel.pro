FUNCTION flame_create_fuel
;
; Creates a default 'fuel' structure that the user can then edit.
; Usage:
; fuel = flame_create_fuel()
;

  ; find the Flame directory
  path_to_thisfile = file_which('flame_create_fuel.pro', /include_current_dir)
  data_dir = flame_util_replace_string(path_to_thisfile, 'pro/flame_create_fuel.pro', 'data/')

  fuel = { $
    science_filelist : 'input/science.txt', $
    darks_filelist: 'input/darks.txt', $
    flats_filelist: 'input/flats.txt', $
    dither_filelist: 'none', $
    intermediate_dir : 'intermediate/', $
    output_dir: 'output/', $
    startrace_y_pos: [0, 0], $
    reduce_only_oneslit : 0, $
    longslit: 0, $
    longslit_edge: [0,0], $
    N_frames: 0, $
    badpix_filename : 'badpix.fits', $
    slitim_filename : 'slitim.fits', $
    flame_data_dir : data_dir, $
    sky_emission_filename : data_dir + 'sky_emission_model.dat', $
    band : '', $
    linelist_filename: '', $
    pixel_scale : 0.0, $
    readnoise : 0.0, $
    gain : 0.0, $
    xrange_star:[1000, 1200], $
    output_lambda_0 : 0d, $
    output_lambda_delta : 0d, $
    output_lambda_Npix : 0, $
    corrscience_files : ptr_new(/allocate_heap), $
    diagnostics : ptr_new(/allocate_heap), $
    slits : ptr_new(/allocate_heap) $
    }
      

    return, fuel



END