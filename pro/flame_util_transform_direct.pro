PRO flame_util_transform_direct, rectification, x=x, y=y, lambda=lambda, gamma=gamma
;
; given the rectification structure, transform the observed coordinates (x, y) into the rectified coordinates (lambda, gamma)
; x, y: input arrays
; lambda, gamma: output arrays
;

  ; check that the input arrays have the same length
  if n_elements(x) ne n_elements(y) then message, 'x and y must have the same number of elements!'


  ; calculate lambda -----------------------------------------------------------

	; order of polynomials
	Nordx = (size(rectification.Klambda))[2]
  Nordy = (size(rectification.Klambda))[1]

  ; make empty array with all zeros
  lambdax = double(x)*0.0

  ; calculate normalized lambda
	for i=0,Nordx-1 do for j=0,Nordy-1 do lambdax += double(rectification.Klambda[j,i]) * double(x)^i * double(y)^j

  ; transform into lambda
  lambda = rectification.lambda_min + rectification.lambda_delta * lambdax


  ; calculate gamma  -----------------------------------------------------------

  ; order of polynomials
	Nordx = (size(rectification.Kgamma))[2]
  Nordy = (size(rectification.Kgamma))[1]

  ; make empty array with all zeros
  gamma = double(x)*0.0

  ; calculate gamma
	for i=0,Nordx-1 do for j=0,Nordy-1 do gamma += double(rectification.Kgamma[j,i]) * double(x)^i * double(y)^j


END
