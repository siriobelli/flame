;
; make cutout for each slit and for each frame and save them in
; slitXX/ directories
;
;******************************************************************

PRO flame_cutout_slits_extract, fuel, slit_structure, input_filename, output_filename

  ; if we are not applying the illumination correction, then we need to be more generous in cutting the margin
  if ~fuel.settings.illumination_correction then margin = 3 else margin = 1

  ; read in science frame
  im = mrdfits(input_filename, 0, header, /silent)

  ; construct the coordinates for the pixels in the image
  N_pix_x = (size(im))[1]
  N_pix_y = (size(im))[2]
  x_axis = indgen(N_pix_x)
  y_axis = indgen(N_pix_y)
  pixel_x = x_axis # replicate(1, N_pix_y)
  pixel_y = transpose(y_axis # replicate(1, N_pix_x))

  ; calculate slit edges
  top_y = (poly(x_axis, slit_structure.bottom_poly) + slit_structure.height) # replicate(1, N_pix_x)
  bottom_y = poly(x_axis, slit_structure.bottom_poly) # replicate(1, N_pix_x)

  ; select pixels belonging to this slit
  w_slit = where( pixel_y LT top_y - margin AND pixel_y GT bottom_y + margin, /null, complement=w_outside_slit)
  if w_slit eq !NULL then message, slit_structure.name + ': slit not valid!'

  ; Set to NaN all pixels outside the slit
  im[w_outside_slit] = !values.d_nan

  ; convert indices to 2D
  w_slit2d = array_indices(im, w_slit)

  ; calculate upper and lower limits
  max_y = max(w_slit2d[1,*])
  min_y = min(w_slit2d[1,*])

  ; extract the slit as a rectangle
  this_slit = im[ * , min_y:max_y]

  ; set to NaN the top and bottom edges of the cutout (to avoid extrapolations)
  this_slit[*, 0] = !values.d_nan
  this_slit[*,-1] = !values.d_nan

  ; write this out
  writefits, output_filename, this_slit, header


  ; cutout the error spectrum
  ; -------------------------

	; check whether there is an extension in the input file (containing the error spectrum)
	rdfits_struct, input_filename, struct, /silent, /header_only
	Next = n_tags(struct)


  if Next GT 1 then begin

    ; read in error frame
    im_sigma = mrdfits(input_filename, 1, /silent)

    ; Set to NaN all pixels outside the slit
    im_sigma[w_outside_slit] = !values.d_nan

    ; extract the slit as a rectangle
    this_slit_sigma = im_sigma[ * , min_y:max_y]

    ; set to NaN the top and bottom edges of the cutout (to avoid extrapolations)
    this_slit_sigma[*, 0] = !values.d_nan
    this_slit_sigma[*,-1] = !values.d_nan

    ; write this out
    writefits, output_filename, this_slit_sigma, /append

  endif


END


;******************************************************************


PRO flame_cutout_slits, fuel

	flame_util_module_start, fuel, 'flame_cutout_slits'


  ; extract slits structure
  slits = fuel.slits

  for i_slit=0, n_elements(slits)-1 do begin

    if slits[i_slit].skip then continue

    print,'Working on slit ', slits[i_slit].name, ' - ', slits[i_slit].number

    ; create directory
    slitdir = fuel.util.intermediate_dir + $
      'slit' + string(slits[i_slit].number,format='(I02)') + '/'
    file_delete, slitdir, /allow_nonexistent, /recursive
    file_mkdir, slitdir

    ; file names for the cutouts
    output_filenames = strarr(fuel.util.science.n_frames)
    for i_frame=0, fuel.util.science.n_frames-1 do begin
      naked_filename = ( strsplit((fuel.util.science.corr_files)[i_frame], '/', /extract) )[-1]
      output_filenames[i_frame] = flame_util_replace_string( slitdir + naked_filename, '.fits', '_slit' + string(slits[i_slit].number,format='(I02)') + '.fits'  )
    endfor

    ; extract slit
    print,'*** Cutting out slit ', slits[i_slit].name

    ; loop through science frames
    for i_frame=0, n_elements(output_filenames)-1 do $
      flame_cutout_slits_extract, fuel, slits[i_slit], fuel.util.science.corr_files[i_frame], output_filenames[i_frame]

    ; add the cutout filenames
    fuel.slits[i_slit].cutouts.filename_step1 = output_filenames

    ; cutout the corresponding arc
    if fuel.util.arc.n_frames gt 0 then begin

      print,'*** Cutting out slit from master arc frame'

      output_filename = slitdir + 'arc_slit' + string(slits[i_slit].number,format='(I02)') + '.fits'
      flame_cutout_slits_extract, fuel, slits[i_slit], fuel.util.arc.master_file, output_filename
      fuel.slits[i_slit].arc_cutout.filename_step1 = output_filename

    endif



  endfor


  flame_util_module_end, fuel

END
