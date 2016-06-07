FUNCTION flame_create_settings
;
; Creates a default 'settings' structure that the user can then edit
;

  settings = { $
    science_filelist : './input/science.txt', $
    darks_filelist: './input/darks.txt', $
    flats_filelist: './input/flats.txt', $
    dither_filelist: 'none', $
    intermediate_dir : './intermediate/', $
    output_dir: './output/', $
    startrace_y_pos: [0, 0], $
    reduce_only_oneslit : 0, $
    longslit: 0, $
    longslit_edge: [0,0] $
    }
      

    return, settings



END