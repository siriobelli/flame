FUNCTION flame_util_transform_coord, x, y, coeff
;
; given the coefficients Lambda_ij or Gamma_ij, transform the observed coordinates (x, y) into one of the rectified coordinates (lambda, gamma)
;

  ; check that the input arrays have the same length
  if n_elements(x) ne n_elements(y) then message, 'x and y must have the same number of elements!'

	; order of polynomials
	Nordx = (size(coeff))[2]
  Nordy = (size(coeff))[1]

  ; make empty array with all zeros
  new_coordinate = double(x)*0.0

  ; calculate the new, rectified, coordinate
	for i=0,Nordx-1 do for j=0,Nordy-1 do new_coordinate += double(coeff[j,i]) * double(x)^i * double(y)^j

  return, new_coordinate

END
