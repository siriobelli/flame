;
; make cutout for each slit and for each frame and save them in
; slitXX/ directories
;
;******************************************************************

PRO flame_cutout_extract, fuel, slit_structure, input_filename, output_filename, yref

  ; how much margin to leave beyond the slit edge, in pixels
  margin = 2

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

  ; extract the slit as a rectangle
  this_cutout = im[ * , slit_structure.yrange_cutout[0] : slit_structure.yrange_cutout[1] ]

  ; set to NaN the top and bottom edges of the cutout (to avoid extrapolations)
  this_cutout[*, 0] = !values.d_nan
  this_cutout[*,-1] = !values.d_nan

  ; add the horizontal trivial coordinate to the FITS header
	SXADDPAR, Header, 'CTYPE1', 'LINEAR'
	SXADDPAR, Header, 'CUNIT1', 'PIXEL'
	SXADDPAR, Header, 'CRPIX1', 1
	SXADDPAR, Header, 'CRVAL1', 1.0
	SXADDPAR, Header, 'CDELT1', 1.0

  ; add the vertical coordinate (relative to the reference star!) to the FITS header
  SXADDPAR, Header, 'CTYPE2', 'LINEAR'
	SXADDPAR, Header, 'CUNIT2', 'PIXEL'
	SXADDPAR, Header, 'CRPIX2', 1
	SXADDPAR, Header, 'CRVAL2', slit_structure.yrange_cutout[0] - yref
	SXADDPAR, Header, 'CDELT2', 1.0

	; delete WCS keywords set by the instrument
	SXDELPAR, Header, 'CD1_1'
	SXDELPAR, Header, 'CD1_2'
	SXDELPAR, Header, 'CD2_1'
	SXDELPAR, Header, 'CD2_2'

  ; write this out
  writefits, output_filename, this_cutout, header


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
    this_cutout_sigma = im_sigma[ * , slit_structure.yrange_cutout[0] : slit_structure.yrange_cutout[1] ]

    ; set to NaN the top and bottom edges of the cutout (to avoid extrapolations)
    this_cutout_sigma[*, 0] = !values.d_nan
    this_cutout_sigma[*,-1] = !values.d_nan

    ; write this out
    writefits, output_filename, this_cutout_sigma, /append

  endif


END


;******************************************************************


PRO flame_cutout_oneslit, fuel, i_slit

  print,'Working on slit ', fuel.slits[i_slit].name, ' - ', fuel.slits[i_slit].number

  ; create directory (and delete previous one if present)
  slitdir = fuel.util.intermediate_dir + $
    'slit' + string(fuel.slits[i_slit].number,format='(I02)') + '/'
  file_delete, slitdir, /allow_nonexistent, /recursive
  file_mkdir, slitdir

  ; file names for the cutouts
  output_filenames = strarr(fuel.util.science.n_frames)
  for i_frame=0, fuel.util.science.n_frames-1 do begin
    naked_filename = file_basename(fuel.util.science.corr_files[i_frame])
    output_filenames[i_frame] = flame_util_replace_string( slitdir + naked_filename, '.fits', '_slit' + string(fuel.slits[i_slit].number,format='(I02)') + '.fits'  )
  endfor

  ; read the x-axis and calculate slit edges
  im = readfits(fuel.util.science.corr_files[0])
  x_axis = dindgen( (size(im))[1] )
  bottom_y = poly(x_axis, fuel.slits[i_slit].bottom_poly)
  top_y = bottom_y + fuel.slits[i_slit].height

  ; set the vertical range for the cutout
  ymin = floor(min(bottom_y)) - 1
  ymax = floor(max(top_y)) + 2
  fuel.slits[i_slit].yrange_cutout = [ymin, ymax]

  ; get the vertical position of the reference star at each frame
  yref = fuel.diagnostics.position

  ; if the reference star is not measured, then set to zero
  if where(~finite(yref), /NULL) NE !NULL then yref[*] = 0.0

  ; extract slit
  print,'*** Cutting out slit ', fuel.slits[i_slit].name

  ; loop through science frames
  for i_frame=0, fuel.util.science.n_frames-1 do $
    flame_cutout_extract, fuel, fuel.slits[i_slit], fuel.util.science.corr_files[i_frame], output_filenames[i_frame], yref[i_frame]

  ; add the cutout filenames to the slits structure
  fuel.slits[i_slit].cutouts.filename = output_filenames

  ; cutout the corresponding arc
  if fuel.util.arc.n_frames gt 0 then begin

    print,'*** Cutting out slit from master arc frame'

    output_filename = slitdir + 'arc_slit' + string(fuel.slits[i_slit].number,format='(I02)') + '.fits'
    flame_cutout_extract, fuel, fuel.slits[i_slit], fuel.util.arc.master_file, output_filename, 0.0
    fuel.slits[i_slit].arc_cutout.filename = output_filename

  endif


END


;******************************************************************


PRO flame_cutout_slits, fuel

	flame_util_module_start, fuel, 'flame_cutout_slits'


  for i_slit=0, n_elements(fuel.slits)-1 do begin

    if fuel.slits[i_slit].skip then continue

    flame_cutout_oneslit, fuel, i_slit

  endfor


  flame_util_module_end, fuel

END
