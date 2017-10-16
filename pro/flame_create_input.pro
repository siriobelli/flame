FUNCTION flame_create_input
;
; Creates a default 'input' structure that the user can then edit.
; Usage:
; input = flame_create_input()
;

  ; create the input substructure and set the default values
  input = { $
    science_filelist : 'science.txt', $
    AB_subtraction: 0, $
    star_y_A: 0, $
    star_y_B: 0, $
    reduce_only_oneslit : 0, $
    longslit: 0, $
    longslit_edge: [0,0], $
    dark_filelist: 'none', $
    pixelflat_filelist: 'none', $
    illumflat_filelist: 'none', $
    illumflat_pixelshift: 0.0, $
    slitflat_filelist: 'none', $
    slitflat_pixelshift: 0.0, $
    arc_filelist: 'none', $
    arc_pixelshift: 0.0, $
    slit_position_file: 'none', $
    dither_file: 'none', $
    max_slitwidth_arcsec : 0.0, $
    intermediate_dir : 'intermediate' + path_sep(), $
    output_dir: 'output' + path_sep() $
    }

  return, input

END
