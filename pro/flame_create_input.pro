FUNCTION flame_create_input
;
; Creates a default 'input' structure that the user can then edit.
; Usage:
; input = flame_create_input()
;

  ; create the input substructure and set the default values
  input = { $
    science_filelist : 'science.txt', $
    dark_filelist: 'none', $
    pixelflat_filelist: 'none', $
    illumflat_filelist: 'none', $
    arc_filelist: 'none', $
    dither_file: 'none', $
    slitflat_filelist: 'none', $
    slit_position_file: 'none', $
    arc_offset: 0.0, $
    slitflat_offset: 0.0, $
    AB_subtraction: 1, $
    star_y_A: 0, $
    star_y_B: 0, $
    reduce_only_oneslit : 0, $
    longslit: 0, $
    longslit_edge: [0,0], $
    max_slitwidth_arcsec : 0.0, $
    intermediate_dir : 'intermediate/', $
    output_dir: 'output/' $
    }

  return, input

END
