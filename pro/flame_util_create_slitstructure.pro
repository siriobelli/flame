FUNCTION flame_util_create_slitstructure, number=number, name=name, PA=PA, $
    approx_bottom=approx_bottom, approx_top=approx_top, approx_target=approx_target, $
    width_arcsec=width_arcsec, approx_R=approx_R, $
    range_lambda0=range_lambda0, range_delta_lambda=range_delta_lambda
;
; This function creates one element of the slit structure, and is called
; during the initialization of an instrument.
;


  ; create 'empty' slit structure
  slit = { $
      number:0, $
      name:'NN', $
      skip:0, $
      PA:!values.d_nan, $
      approx_bottom:0, $
      approx_top:0, $
      approx_target:0, $
      width_arcsec:!values.d_nan, $
      approx_R:0.0, $
      range_lambda0:[0.0, 0.0], $
      range_delta_lambda:[0.0, 0.0] }


    ; fill the structure with the values passed as arguments
    if n_elements(number) GT 0 then $
      slit.number = number

    if n_elements(name) GT 0 then $
      slit.name = name

    if n_elements(PA) GT 0 then $
      slit.PA = PA

    if n_elements(approx_bottom) GT 0 then $
      slit.approx_bottom = approx_bottom

    if n_elements(approx_top) GT 0 then $
      slit.approx_top = approx_top

    if n_elements(approx_target) GT 0 then $
      slit.approx_target = approx_target

    if n_elements(width_arcsec) GT 0 then $
      slit.width_arcsec = width_arcsec

    if n_elements(approx_R) GT 0 then $
      slit.approx_R = approx_R

    if n_elements(range_lambda0) GT 0 then $
      slit.range_lambda0 = range_lambda0

    if n_elements(range_delta_lambda) GT 0 then $
      slit.range_delta_lambda = range_delta_lambda


    ; return the slit structure
    return, slit


END
