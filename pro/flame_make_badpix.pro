;
; original version 14 Jan 2015, Eva Wuyts
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO flame_make_badpix, fuel=fuel


; identify the darks files
  readcol, fuel.darks_filelist, darklist, format='A'
  
  ndarks = n_elements(darklist)
  print, 'Number of darks: ',ndarks
  foo = readfits(darklist[0],/silent)
  s = size(foo, /dim)
  dark = fltarr(ndarks,s(0),s(1))
  for j=0, ndarks-1 do dark(j,*,*) = readfits(darklist[j],/silent)

  darktot = fltarr(s(0),s(1))
  for j=0, s(0)-1 do begin
     for k=0, s(1)-1 do begin
        darktot(j,k) = median(dark(*,j,k))
     endfor
  endfor
  
  writefits, fuel.intermediate_dir + 'darktot.fits', darktot
 
  ;
  ; need to get a better way to identify bad pixels!
  ; 
  ; find badpixels
  badpix = fltarr(s(0),s(1))
  bad = where(darktot lt -50 or darktot gt 50., /null)
  badpix(bad)=1
  
  writefits, fuel.intermediate_dir + fuel.badpix_filename, badpix
  
END
