;
; make cutout for each and for each frame and save them in
; slitXX/ directories
;
;******************************************************************

PRO flame_cutout_slits_extract, slit_structure, science_filenames, output_filenames

  ; loop through science frames
  for i_frame=0, n_elements(science_filenames)-1 do begin

    ; read in science frame
    im = mrdfits(science_filenames[i_frame], 0, header, /silent)

    ; read in error frame
    im_sigma = mrdfits(science_filenames[i_frame], 1, /silent)

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
    w_slit = where( pixel_y LT top_y AND pixel_y GT bottom_y, /null, complement=w_outside_slit)
    if w_slit eq !NULL then message, slit_structure.name + ': slit not valid!'

    ; Set to NaN all pixels outside the slit
    im[w_outside_slit] = !values.d_nan
    im_sigma[w_outside_slit] = !values.d_nan

    ; convert indices to 2D
    w_slit2d = array_indices(im, w_slit)

    ; calculate upper and lower limits
    max_y = max(w_slit2d[1,*])
    min_y = min(w_slit2d[1,*])

    ; extract the slit as a rectangle
    this_slit = im[ * , min_y:max_y]
    this_slit_sigma = im_sigma[ * , min_y:max_y]

    ; write this out
    writefits, output_filenames[i_frame], this_slit, header
    writefits, output_filenames[i_frame], this_slit_sigma, /append

  endfor

END


;******************************************************************


PRO flame_cutout_slits, fuel

	start_time = systime(/seconds)

  print, ''
  print, '-------------------------------------'
  print, '---      flame_cutout_slits       ---'
  print, '-------------------------------------'
  print, ''


  ; extract slits structure
  slits = fuel.slits

  print,'Slits: ', n_elements(slits)
  for i_slit=0, n_elements(slits)-1 do begin

    print,'Working on slit ', slits[i_slit].name, ' - ', slits[i_slit].number

    ; create directory
    slitdir = file_expand_path(fuel.input.intermediate_dir) + $
      '/slit' + string(slits[i_slit].number,format='(I02)') + '/'
    file_delete, slitdir, /allow_nonexistent, /recursive
    file_mkdir, slitdir

    ; file names for the cutouts
    output_filenames = strarr(fuel.util.n_frames)
    for i_frame=0, fuel.util.n_frames-1 do begin
      naked_filename = ( strsplit((fuel.util.corrscience_filenames)[i_frame], '/', /extract) )[-1]
      output_filenames[i_frame] = flame_util_replace_string( slitdir + naked_filename, '.fits', '_slit' + string(slits[i_slit].number,format='(I02)') + '.fits'  )
    endfor

    ; extract slit
    print,'*** Cutting out slit ', slits[i_slit].name
    flame_cutout_slits_extract, slits[i_slit], (fuel.util.corrscience_filenames), output_filenames

    ; add the cutout filenames
    fuel.slits[i_slit].cutouts.filename = output_filenames

  endfor


	print, ''
  print, 'flame_cutout_slits took ', $
    cgnumber_formatter( systime(/seconds) - start_time, decimals=2), ' seconds'
  print, ''

END
