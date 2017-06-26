
PRO flame_checkdata_seeing, fuel

	; check if the reference star has been specified
	if fuel.input.star_y_A eq 0.0 then return

	; x coordinate where the star trace is certainly visible
	star_x = mean(fuel.util.star_x_range)

	; identify the slit with the reference star
	i_ref = -1
	for i_slit=0, n_elements(fuel.slits)-1 do $
		if fuel.input.star_y_A GE poly(star_x, fuel.slits[i_slit].bottom_poly) $
			and fuel.input.star_y_A LE poly(star_x, fuel.slits[i_slit].bottom_poly) + $
			fuel.slits[i_slit].height then $
				i_ref = i_slit

	; if there is no slit with the reference star, then exit
	if i_ref eq -1 then begin
		print, 'Did not find the slit with the reference star.'
		return
	endif

	print, 'Reference star is in slit number ', strtrim(fuel.slits[i_ref].number, 2)
	



END


;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_checkdata, fuel

	flame_util_module_start, fuel, 'flame_checkdata'


	; calculate seeing from the reference star
	flame_checkdata_seeing, fuel


  flame_util_module_end, fuel


  print, '-------------------------------------'
	print, '-------------------------------------'
	print, '-------------------------------------'

	; print total execution time
	print, ' '
	print, 'The data reduction took a total of ', $
		cgnumber_formatter((systime(/seconds) - fuel.util.start_time)/60.0, decimals=2), ' minutes.'
		print, ' '


END
