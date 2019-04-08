

;*******************************************************************************
;*******************************************************************************
;*******************************************************************************



PRO flame_extract_slit, fuel, slit, skysub=skysub


	; set up directories
	; ----------------------------------------------------------------------------

	; easy way to deal with names of directories
	if keyword_set(skysub) then skysub_string = '_skysub' else skysub_string = ''

	; directory from which to extract the spectra
	spec2d_dir = fuel.util.output_dir + 'spec2d' + skysub_string + path_sep()

	; directory where to save the 1D spectra
	extraction_dir = fuel.util.output_dir + 'spec1d' + skysub_string + path_sep()

	; if needed, create extraction directory in the output directory
  if ~file_test(extraction_dir) then file_mkdir, extraction_dir


	; extract spectrum
	; ----------------------------------------------------------------------------

	; FITS file to be used for the extraction
	filename_spec2d = spec2d_dir + slit.output_combined_file

	; filename for the output FITS file
	if fuel.settings.extract_optimal then extraction_type = 'optimal.' $
		else extraction_type = 'boxcar.'
	filename_output = extraction_dir + 'spec1d.' + extraction_type + string(slit.number, format='(I03)') + '.' + slit.output_combined_file

	; determine what type of extraction
	if fuel.settings.extract_optimal then begin
		if fuel.settings.extract_gaussian_profile then optimal_gaussian=1 else $
			optimal_profile=1
	endif else boxcar_auto=1

	; perform extraction
	flame_util_extract_spectrum, filename_spec2d, $
		boxcar_auto=boxcar_auto, optimal_profile=optimal_profile, optimal_gaussian=optimal_gaussian, $
		output_filename=filename_output


END



;*******************************************************************************
;*******************************************************************************
;*******************************************************************************


PRO flame_extract, fuel

	flame_util_module_start, fuel, 'flame_extract'

	; keep track of the spectra already extracted (useful for the "paired slits" nodding)
	slits_extracted = []

	; extract 1D spectrum for each slit
	for i_slit=0, n_elements(fuel.slits)-1 do begin

		if fuel.slits[i_slit].skip then continue

		; handle errors by ignoring that slit
		if fuel.settings.stop_on_error eq 0 then begin
			catch, error_status
			if error_status ne 0 then begin
				print, ''
		    print, '**************************'
		    print, '***       WARNING      ***'
		    print, '**************************'
		    print, 'Error found. Skipping slit ' + strtrim(fuel.slits[i_slit].number,2), ' - ', fuel.slits[i_slit].name
				fuel.slits[i_slit].skip = 1
				catch, /cancel
				continue
			endif
		endif

		; check if this slit has already been combined with another, and extracted
		if total(slits_extracted eq fuel.slits[i_slit].output_combined_file) GT 0 then continue

		; extract the 1D spectrum for this slit
		flame_extract_slit, fuel, fuel.slits[i_slit]
		if fuel.settings.skysub then $
			flame_extract_slit, fuel, fuel.slits[i_slit], /skysub

		; add this slit to the list of the ones that have been extracted
		slits_extracted = [slits_extracted, fuel.slits[i_slit].output_combined_file]

	endfor

	; save fuel structure to output directory
	filename = fuel.util.output_dir + 'fuel.sav'
	save, fuel, filename=filename
	print, 'fuel structure saved to ' + filename


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
