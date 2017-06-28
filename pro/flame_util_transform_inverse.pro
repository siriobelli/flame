PRO flame_util_transform_inverse, rectification, lambda=lambda, gamma=gamma, x=x, y=y
;
; given the rectification structure, transform the rectified coordinates (lambda, gamma) into the observed coordinates (x, y)
; lambda, gamma: input arrays
; x, y: output arrays
;

  ; check that the input arrays have the same length
  if n_elements(lambda) ne n_elements(gamma) then message, 'lambda and gamma must have the same number of elements!'

  ; make the arrays that will have the observed coordinates
  x = dblarr(n_elements(lambda))
  y = dblarr(n_elements(lambda))

	; order of polynomial
	Nord = (size(rectification.Kx))[1]

  ; exponents for the polynomial transformation
	lexp  = findgen(Nord)
	gexp  = findgen(Nord)

  ; instead of using lambda, let's use the "normalized" version
  lambda_norm = ( lambda - rectification.lambda_min ) / rectification.lambda_delta

  ; calculate x
	for i=0, n_elements(lambda)-1 do x[i] = $
    total(((gamma[i])^lexp # (lambda_norm[i])^gexp ) * rectification.Kx)

  ; calculate y
	for i=0, n_elements(lambda)-1 do y[i] = $
    total(((gamma[i])^lexp # (lambda_norm[i])^gexp ) * rectification.Ky)

END
