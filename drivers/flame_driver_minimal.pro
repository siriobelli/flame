
  ; suppress harmless messages (e.g., module compilation)
  !quiet = 1

  ; create the input structure
  input = flame_create_input()

  ; text file containing the list of science FITS files that need to be reduced
  input.science_filelist = 'science.txt'

  ; do you want to apply A-B sky subtraction?
  input.AB_subtraction = 1

  ; array with y-pixel positions for the traces of the reference star. 0 if there is no reference star
  input.star_y_A = 1281
  input.star_y_B = 1300

  ; initialize and create the fuel structure
  fuel = flame_initialize_minimal(input, pixel_scale=0.18)

  ; data reduction: only the first two steps

  flame_diagnostics, fuel

  flame_quickstack, fuel

END
