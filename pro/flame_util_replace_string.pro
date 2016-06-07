FUNCTION flame_util_replace_string, input, pattern, replacement

	; take care of recursive call for arrays
	if n_elements(input) GT 1 then begin
		output = strarr(n_elements(input))
		for i=0,n_elements(output)-1 do output[i] = flame_util_replace_string(input[i], pattern, replacement)
		return, output
	endif

	; find the pattern in the input string (find only the last occurrence)
	pattern_start = strpos(input, pattern, /reverse_search )
	pattern_length = strlen(pattern)

	; check that the pattern exist at least once
	if pattern_start eq -1 then begin
		;print, 'flame_util_replace_string: did not find the pattern in the input string'
		return, input 
	endif

	; identify the part before the pattern
	if pattern_start eq 0 then output_initial = '' $
		else output_initial = strmid(input, 0, pattern_start)

	; identify the part after the pattern
	if pattern_start + pattern_length eq strlen(input) then output_final = '' $
		else output_final = strmid(input, pattern_start + pattern_length, strlen(input))

	; return the sum of the first part, the replacement, and the final part
	return, output_initial + replacement + output_final


END