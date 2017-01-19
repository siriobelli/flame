FUNCTION flame_create_input
;
; Creates a default 'input' structure that the user can then edit.
; Usage:
; input = flame_create_input()
;

  ; create the input substructure and set the default values
  input = { $
    science_filelist : 'science.txt', $
    darks_filelist: 'none', $
    flats_filelist: 'none', $
    dither_filelist: 'none', $
    intermediate_dir : 'intermediate/', $
    output_dir: 'output/', $
    startrace_y_pos: [0, 0], $
    reduce_only_oneslit : 0, $
    longslit: 0, $
    longslit_edge: [0,0], $
    xrange_star:[1000, 1200], $
    wavecal_approx_smooth : 5, $
    use_sky_edge : 0 $
    }

  return, input

END