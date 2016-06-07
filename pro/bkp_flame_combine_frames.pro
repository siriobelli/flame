;
;
; TO-DO:
;
; - frames are combined by shifting A-B-A by an integer number of pixels - can you do better?
; - dithering length is calculated from first A frame and first B frame - we need a better definition
; - handle better the case in which A or B are actually sky
;



PRO flame_combine_frames_slit, slit=slit, filenames=filenames, diagnostics=diagnostics, output_dir=output_dir

	print, 'Combining these files: '
	print, filenames

	; number of frames 
	N_frames = n_elements(filenames)

	; read first frame
	im_0 = readfits(filenames[0], header)

	; make big cube containing all frames
	im_cube = dblarr( N_frames, (size(im_0))[1], (size(im_0))[2] )

	; read all frames 
	for i_frame=0,N_frames-1 do im_cube[i_frame,*,*] = readfits(filenames[i_frame])

	; prefix for output file names
	filename_prefix = output_dir + 'slit' + string(slit.number, format='(I02)') + '-' + slit.name

	; select all A frames
	w_A = where(diagnostics.offset_pos eq 'A', /null)
	
	; select all B frames
	w_B = where(diagnostics.offset_pos eq 'B', /null)

	; if there are no A frames, select the non-B as A, and viceversa
	w_all = indgen(N_frames)
	if w_A eq !NULL then w_A = cgsetdifference(w_all, w_B)
	if w_B eq !NULL then w_B = cgsetdifference(w_all, w_A)

	; write raw sum of As and Bs 
	A_sum = mean(im_cube[w_A,*,*], dimension=1, /nan)
	B_sum = mean(im_cube[w_B,*,*], dimension=1, /nan)
	writefits, filename_prefix + '_A.fits', A_sum, header
	writefits, filename_prefix + '_B.fits', B_sum, header


	; first, let's sigma-clip so to eliminate cosmic rays
	; *******************************************************************

	; decide the amount of clipping
	threshold = 5.0		; 5-sigma clipping

	; calculate median at each pixel of the image, for the A and the B frames
	median_im_A = median(im_cube[w_A,*,*], dimension=1)
	median_im_B = median(im_cube[w_B,*,*], dimension=1)

	; calculate standard deviation at each pixel of the image, for the A and B frames
	sigma_im_A = stddev(im_cube[w_A,*,*], dimension=1, /nan)
	sigma_im_B = stddev(im_cube[w_B,*,*], dimension=1, /nan)

	; make a cube with the corresponding median at each (pixel,frame) position
	median_cube = im_cube
	median_cube[*]=0.0
	for i_A=0,n_elements(w_A)-1 do median_cube[w_A[i_A], *, *] = median_im_A
	for i_B=0,n_elements(w_B)-1 do median_cube[w_B[i_B], *, *] = median_im_B

	; make a cube with the corresponding sigma at each (pixel,frame) position
	sigma_cube = im_cube
	sigma_cube[*]=0.0
	for i_A=0,n_elements(w_A)-1 do sigma_cube[w_A[i_A], *, *] = sigma_im_A
	for i_B=0,n_elements(w_B)-1 do sigma_cube[w_B[i_B], *, *] = sigma_im_B

	; for each pixel in the cube, calculate the deviation from the median in units of sigma
	deviation = (im_cube - median_cube)/sigma_cube 

	; and now make a mask that indicates those pixels that are more than "threshold" times sigma from the median
	mask_cube = im_cube
	mask_cube[*] = 0.0
	mask_cube[where( abs(deviation) GT threshold, /null )] = 1.0

	; turn masked pixels into NaNs
	im_cube[where(mask_cube, /null)] = !values.d_nan

	; write the mask image
	writefits, filename_prefix + '_Npixels_rejected.fits', total(mask_cube, 1), header


	; now that we cleaned the cube, let's combine all frames 
	; *******************************************************************

	; these are clean sums
	A_cleansum = mean(im_cube[w_A,*,*], dimension=1, /nan)
	B_cleansum = mean(im_cube[w_B,*,*], dimension=1, /nan)

	; write them out 
	writefits, filename_prefix + '_Aclean.fits', A_cleansum, header
	writefits, filename_prefix + '_Bclean.fits', B_cleansum, header
	
	; A-B
	AB_cleansum = A_cleansum - B_cleansum
	writefits, filename_prefix + '_A-B.fits', AB_cleansum, header
	print, 'I wrote ', filename_prefix + '_A-B.fits'


	; let's now make an output with the negative-positive-negative traces
	; *******************************************************************

	; invert the master A-B
	BA_cleansum = -AB_cleansum

	; find the dithering length, which for now we take from the first A and the first B frames
	; (A is always the one on top)
	dithering_length = round(diagnostics[w_A[0]].position - diagnostics[w_B[0]].position)

	; zero padding on top
	padding = dblarr( (size(im_0))[1], dithering_length )
	padding[*] = 0.0
	AB_cleansum_padded = [ [AB_cleansum], [padding] ]
	BA_cleansum_padded = [ [BA_cleansum], [padding] ]

	BA_cleansum_shifted = shift(BA_cleansum_padded, 0, dithering_length)

	AB_combined = mean( [ [[AB_cleansum_padded]], [[BA_cleansum_shifted]] ], dimension=3, /nan )

	; last detail: the beginning and the end of the spectrum, outside the observed range, are NaN.
	; but the zero padding added zeroes there. Let's delete them

	; first find the x coordinates of the NaN regions, defined as columns where >90% of pixels (pre-padding) are NaNs
	number_NaNs = total( ~finite(AB_cleansum), 2)
	w_nans = where(number_NaNs GT 0.9 * (size(AB_cleansum))[2], /null )

	; now let's set all those columns equal to NaN
	if w_nans NE !NULL then $
		for i=0, n_elements(w_nans)-1 do AB_combined[w_nans[i],*] = !values.d_NaN

	; output final result
	writefits, filename_prefix + '_ABcombined.fits', AB_combined, header

	print, 'I wrote ', filename_prefix + '_ABcombined.fits'


END



; ---------------------------------------------------------------------------------------------------------------------------



PRO flame_combine_frames, fuel=fuel
 

	print, ' '
	print, 'Combine frames'
	print, '****************'

	; extract the slits 
	slits = *fuel.slits

 	; loop through all slits
	for i_slit=0, n_elements(slits)-1 do begin

		this_slit = slits[i_slit]
		print, 'Combining slit ', this_slit.number, ' - ', this_slit.name

		; combine the rectified sky-subtracted frames
		flame_combine_frames_slit, slit=this_slit, filenames=flame_util_replace_string(*this_slit.filenames, '.fits', '_skysub_rectified.fits'), diagnostics=*fuel.diagnostics, output_dir=fuel.output_dir

	endfor


END



