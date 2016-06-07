PRO flame_test, fuel=fuel

; let's stack all FITS files in order to get good SNR on the OH lines

stop

i_slit=0


this_slit = (*fuel.slits)[i_slit]

filenames = *this_slit.filenames

im0 = readfits(filenames[0])
im = im0

for i_frame=1, n_elements(filenames)-1 do im += readfits(filenames[i_frame])


END