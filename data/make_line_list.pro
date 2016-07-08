;
; adapt the line list from Rousselot et al. to the flame needs
;
;************************************************************
;
; take successive windows of 0.10 microns and select the locally brightest lines
; (flux > 2% of flux of the brightest line), then merge narrow doublets that are within
; a resolveing power of 10000, and finally exclude the other doublets that are within a
; resolving power of 5000
;
;************************************************************



; read in line list
readcol, 'line_list.dat', line_lambda, line_flux

; transform to microns
line_lambda *= 1d-4

; load sky spectrum model
readcol, 'sky_emission_model.dat', sky_lambda, sky_flux


; set a small lambda range with which we can work
l0 = 2.4
l1 = 2.5
l_range = [l0 - 0.05*(l1-l0) , l1 + 0.05*(l1-l0) ]
wline = where(line_lambda GE l0 and line_lambda LT l1, /null)
wsky =  where(sky_lambda GE l0 and sky_lambda LT l1, /null)

; normalize line fluxes
norm_line = max(line_flux[wline])
norm_sky = max(sky_flux[wsky])

; discard lines that are less than p*100% of the brightest one in flux
p = 0.02
wline = where(line_lambda GE l0 and line_lambda LT l1 and line_flux GT p*norm_line, /null)

; show spectrum and line fluxes
cgplot, sky_lambda[wsky], sky_flux[wsky]/norm_sky, xra=l_range, /xsty
cgplot, line_lambda[wline], line_flux[wline]/norm_line, /overplot, psym=16, color='red'
for i=0,n_elements(wline)-1 do cgplot, line_lambda[wline[i]] + [0.,0.], [0,1], /overplot, color='red'



; create new list of lines 
goodline_lambda = line_lambda[wline]
goodline_flux = line_flux[wline]

; print lines 
for i=0,n_elements(goodline_lambda)-1 do print, i, goodline_lambda[i], goodline_flux[i]

; identify doublets that are so tight (less than R=10000) that we consider them single lines 
separation = (goodline_lambda - shift(goodline_lambda,-1)) / goodline_lambda
w_doub = where( abs(separation) LT 1d-4, /null)

; replace doublets with their average values
goodline_lambda[w_doub] = 0.5*(goodline_lambda[w_doub] + goodline_lambda[w_doub+1])
goodline_flux[w_doub] = 0.5*(goodline_flux[w_doub] + goodline_flux[w_doub+1])

; get rid of the second line in each doublet
goodline_lambda[w_doub+1] = !values.d_nan

; keep the good values
w_ok = where(finite(goodline_lambda), /null)
goodline_lambda = goodline_lambda[w_ok]
goodline_flux = goodline_flux[w_ok]

; show spectrum and line fluxes
cgplot, sky_lambda[wsky], sky_flux[wsky]/norm_sky, xra=l_range, /xsty
cgplot, goodline_lambda, goodline_flux/norm_line, /overplot, psym=16, color='blue'
for i=0,n_elements(goodline_lambda)-1 do cgplot, goodline_lambda[i] + [0.,0.], [0,1], /overplot, color='blue'



; print good lines 
for i=0,n_elements(goodline_lambda)-1 do print, i, goodline_lambda[i], goodline_flux[i]

; now identify the 'bad' doublets, that could be semi-resolved in the data (R~5000)
separation = (goodline_lambda - shift(goodline_lambda,-1)) / goodline_lambda
w_doub = where( abs(separation) LT 1.0/5d3, /null)

; show them in the plot
for i=0,n_elements(w_doub)-1 do cgplot, goodline_lambda[w_doub[i]] + [0.,0.], [0,1], /overplot, color='green'
for i=0,n_elements(w_doub)-1 do cgplot, goodline_lambda[w_doub[i]+1] + [0.,0.], [0,1], /overplot, color='green'

; eliminate them
goodline_lambda[w_doub] = !values.d_nan
goodline_lambda[w_doub+1] = !values.d_nan

; final selection
w_ok = where(finite(goodline_lambda), /null)
goodline_lambda = goodline_lambda[w_ok]
goodline_flux = goodline_flux[w_ok]

; print and show final selection
for i=0,n_elements(goodline_lambda)-1 do cgplot, goodline_lambda[i] + [0.,0.], [0,1], /overplot, color='red'
for i=0,n_elements(goodline_lambda)-1 do cgtext, goodline_lambda[i], 0.7, strtrim(i,2), color='red', alignment=-0.2, charsize=1.0
for i=0,n_elements(goodline_lambda)-1 do print, i, goodline_lambda[i], goodline_flux[i]

; print the final output that will go in the line list (after commenting out potential bad lines by hand)
for i=0,n_elements(goodline_lambda)-1 do print, goodline_lambda[i]


