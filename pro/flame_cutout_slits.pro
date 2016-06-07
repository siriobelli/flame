 
PRO flame_cutout_slits_extract, slitim, slit_structure, science_filenames, output_filenames

  ; loop through science frames
  for i_frame=0, n_elements(science_filenames)-1 do begin

    ; read in science frame 
    im = readfits(science_filenames[i_frame], header)

    ; select pixels belonging to this slit
    w_slit = where(slitim eq slit_structure.number, /null)
    if w_slit eq !NULL then message, slit_structure.name + ': slit not found in slitim!'

    ; create a mask that selects only pixels belonging to the slit 
    slit_mask = im
    slit_mask[*] = !values.d_nan
    slit_mask[w_slit] = 1.0

    ; create a new frame where everything outside the slit is a Nan 
    im_masked = im * slit_mask

    ; convert indices to 2D
    w_slit2d = array_indices(slitim, w_slit)

    ; calculate upper and lower limits
    max_y = max(w_slit2d[1,*])
    min_y = min(w_slit2d[1,*])

    ; extract the slit as a rectangle
    this_slit = im_masked[ * , min_y:max_y]

    ; write this out
    writefits, output_filenames[i_frame], this_slit, header

  endfor

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 

PRO flame_cutout_slits, fuel=fuel

  ; extract slits structure
  slits = *fuel.slits

  ; read in the slitim image
  slitim = readfits(fuel.intermediate_dir + fuel.slitim_filename)

  print,'Slits: ', n_elements(slits)
  for i_slit=0, n_elements(slits)-1 do begin
  
    print,'Working on slit ', slits[i_slit].name, ' - ', slits[i_slit].number
    
    ; create directory
    slitdir = fuel.intermediate_dir + 'slit' + string(slits[i_slit].number,format='(I02)') + '/'
    spawn,'rm -rf ' + slitdir
    if file_test(slitdir) eq 0 then spawn, 'mkdir ' + slitdir

    ; file names for the cutouts
    output_filenames = strarr(fuel.n_frames)
    for i_frame=0, fuel.n_frames-1 do begin
      naked_filename = ( strsplit((*fuel.corrscience_files)[i_frame], '/', /extract) )[-1]
      output_filenames[i_frame] = flame_util_replace_string( slitdir + naked_filename, '.fits', '_slit' + string(slits[i_slit].number,format='(I02)') + '.fits'  )
    endfor

    ; extract slit
    print,'*** Cutting out slit ', slits[i_slit].name
    flame_cutout_slits_extract, slitim, slits[i_slit], (*fuel.corrscience_files), output_filenames

    ; add filenames to the slit structure
    *slits[i_slit].filenames = output_filenames

  endfor

END
  