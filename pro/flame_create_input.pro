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
    dither_filelist: 'none', $
    slitflat_filelist: 'none', $
    slit_position_file: 'none', $
    AB_subtraction: 1, $
    star_y_A: 0, $
    star_y_B: 0, $
    star_x_range: [1000, 1200], $
    reduce_only_oneslit : 0, $
    longslit: 0, $
    longslit_edge: [0,0], $
    use_sky_edge : 0, $
    rough_wavecal_R : [500, 1000, 3000], $
    clean_individual_frames : 0, $
    sigma_clip : 2.0, $
    max_slitwidth_arcsec : 0.0, $
    intermediate_dir : 'intermediate/', $
    output_dir: 'output/' $
    }

  return, input

END
