;
; adapt the line list from Rousselot et al. to the LUCAS needs
;

; read in line list
readcol, '~/idl/luci/LUCAS2/lib/line_list.dat', line_lambda, line_flux

;*********************
;        J band
;*********************

; select wavelength range and intensity
w_ok = where( line_lambda GT 1.0d4 and line_lambda LT 1.4d4 AND line_flux GT 40.0)

; extract lines of interests (in angstrom)
x = line_lambda[w_ok]

; round to the nearest 3 angstroms
x_approx = fix(x+0.5) / 3

; if there are doublets to within 1 angstrom, take only one of the two lines
w_uniq = uniq(x_approx)

; write line list as text file
writecol, 'line_list_J.dat', x[w_uniq]*1d-4



;*********************
;        K band
;*********************

; select wavelength range and intensity
w_ok = where( line_lambda GT 1.9d4 and line_lambda LT 2.5d4 AND line_flux GT 40.0)

; extract lines of interests (in angstrom)
x = line_lambda[w_ok]

; round to the nearest 3 angstroms
x_approx = fix(x+0.5) / 3

; if there are doublets to within 1 angstrom, take only one of the two lines
w_uniq = uniq(x_approx)

; write line list as text file
writecol, 'line_list_K.dat', x[w_uniq]*1d-4
