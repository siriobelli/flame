
PRO flame_spatialcal_one, fuel=fuel, slit=slit, cutout=cutout, $
		diagnostics=diagnostics, this_diagnostics=this_diagnostics

;
; The result of this routine is the set of coefficients gamma_coeff, saved in fuel.slits,
; that describe the transformation from observed (x,y) coordinates to the
; rectified vertical coordinate gamma
;


	; find the gamma coefficients -------------------------------------------------------


	; first, find the vertical offset to the reference star
	w_this_offset = where(fuel.diagnostics.offset_pos eq this_diagnostics.offset_pos, /null)
	ref_diagnostics = fuel.diagnostics[w_this_offset[0]]
	vertical_offset = this_diagnostics.position - ref_diagnostics.position

	; find the average x coordinate of the reference star trace
	xref = mean(fuel.settings.star_x_range)

	; set the size of the matrix: linear in y, and same degree as slit edge in x
	N_y = 2
	N_x = n_elements(slit.bottom_poly)

	; make the matrix of coefficients
	gamma_coeff = dblarr(N_y, N_x)

	; by definition, dgamma/dy = 1
	gamma_coeff[1,0] = 1.0

	; use the definition of the bottom edge of the slit to build the gamma matrix
	gamma_coeff[0,*] = -slit.bottom_poly

	; but the zero point is such that gamma at xref is the vertical distance to the star trace
	gamma_coeff[0,0] += slit.yrange_cutout[0] - this_diagnostics.position + poly(xref, slit.bottom_poly)


	; save and output -------------------------------------------------------

	; save into slit structure
	*cutout.gamma_coeff = gamma_coeff

END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************




PRO flame_spatialcal, fuel

	flame_util_module_start, fuel, 'flame_spatialcal'


	; loop through all slits
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		; this slit
		this_slit = fuel.slits[i_slit]

		if this_slit.skip then continue

	  print, 'Spatial calibration for slit ', strtrim(this_slit.number,2), ' - ', this_slit.name
		print, ' '

		; handle errors by ignoring that slit
		if fuel.settings.debugging eq 0 then begin
			catch, error_status
			if error_status ne 0 then begin
				print, ''
		    print, '**************************'
		    print, '***       WARNING      ***'
		    print, '**************************'
		    print, 'Error found. Skipping slit ' + strtrim(this_slit.number,2), ' - ', this_slit.name
				fuel.slits[i_slit].skip = 1
				catch, /cancel
				continue
			endif
		endif


		if fuel.util.arc.n_frames GT 0 then begin

				; calculate the polynomial transformation between observed and rectified frame, for the arcs
				flame_spatialcal_one, fuel=fuel, slit=this_slit, cutout=this_slit.arc_cutout, $
					diagnostics=fuel.diagnostics, this_diagnostics=(fuel.diagnostics)[0] 	; assume the dithering of the first frame

		endif


		for i_frame=0, n_elements(fuel.slits[i_slit].cutouts)-1 do begin

				; calculate the polynomial transformation between observed and rectified frame
				flame_spatialcal_one, fuel=fuel, slit=this_slit, cutout=this_slit.cutouts[i_frame], $
					diagnostics=fuel.diagnostics, this_diagnostics=(fuel.diagnostics)[i_frame]

		endfor

	endfor


  flame_util_module_end, fuel

END
