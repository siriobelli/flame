PRO flame_util_transform_direct, rectification, x=x, y=y, lambda=lambda, gamma=gamma
;
; given the rectification structure, transform the observed coordinates (x, y) into the rectified coordinates (lambda, gamma)
; x, y: input arrays
; lambda, gamma: output arrays
;

  ; check that the input arrays have the same length
  if n_elements(x) ne n_elements(y) then message, 'x and y must have the same number of elements!'

  ; make the arrays that will have the rectified coordinates
  lambda = double(x)*0.0
  gamma = double(x)*0.0

	; order of polynomial
	Nord = (size(rectification.Klambda))[1]

  ; exponents for the polynomial transformation
	xexp  = findgen(Nord)
	yexp  = findgen(Nord)

  ; calculate lambda
	for i=0, n_elements(x)-1 do lambda[i] = $
		rectification.lambda_min + rectification.lambda_delta * $
    total(((y[i])^xexp # (x[i])^yexp ) * rectification.Klambda)

  ; calculate gamma
	for i=0, n_elements(x)-1 do gamma[i] = $
    total(((y[i])^xexp # (x[i])^yexp ) * rectification.Kgamma)


END
