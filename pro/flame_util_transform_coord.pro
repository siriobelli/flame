FUNCTION flame_util_transform_coord, x, y, coeff
;
; given the coefficients Lambda_ij or Gamma_ij, transform the observed coordinates (x, y) into one of the rectified coordinates (lambda, gamma)
;

  ; check that the input arrays have the same length
  if n_elements(x) ne n_elements(y) then message, 'x and y must have the same number of elements!'

	; dimensionality of coeff
  N_dimensions = (size(coeff))[0]

  ; make sure that coeff is a vector or a matrix
  if N_dimensions NE 1 AND N_dimensions NE 2 then $
    message, 'the coefficient matrix must be a 1D or 2D array!'

  ; if coeff is a 1D vector
  if N_dimensions EQ 1 then begin
    Nordx = 1
    Nordy = (size(coeff))[1]
  endif

  ; if coeff is a 2D matrix
  if N_dimensions EQ 2 then begin
  	Nordx = (size(coeff))[2]
    Nordy = (size(coeff))[1]
  endif

  ; make empty array with all zeros
  new_coordinate = double(x)*0.0

  ; calculate the new, rectified, coordinate
	for i=0,Nordx-1 do for j=0,Nordy-1 do new_coordinate += double(coeff[j,i]) * double(x)^i * double(y)^j

  return, new_coordinate

END
