
PRO flame_util_extract_spectrum, filename, yrange=yrange, output_filename=output_filename
;
; Boxcar extraction of a 1D spectrum from a 2D spectrum output by Flame (i.e.,
; FITS file with six extensions: flux, noise, sigma, sky, exptime, weight)
; The vertical pixel range to be considered for the extraction must be specified with
; the yrange keyword (starting from zero; extrema are included in the extraction)
; The output is a FITS file with a structure containing the fields lambda, flux,
; ivar (inverse variance), ivar_empirical, and sky, and is compatible with SpecPro.
;


  ; test input
  ; ----------------------------------------------------------------------------

  ; number of input files
  if n_elements(filename) NE 1 then message, 'Input must be one and only one filename'

  ; set the output filename
  if ~keyword_set(output_filename) then begin
    file_extension = (strsplit(filename, '.', /extract))[-1]
    output_filename = flame_util_replace_string(filename, '.' + file_extension, '_spec1d.fits')
  endif

  if n_elements(yrange) NE 2 then message, 'yrange = [y0, y1] must be specified'


  ; read in 2D spectrum
  ; ----------------------------------------------------------------------------

  ; read in 2d spectrum
  spec2d = mrdfits(filename, 0, header, /silent)

  ; read in 2d theoretical error spectrum
  err2d = mrdfits(filename, 1, /silent)

	; read in 2d empirical error spectrum
	sig2d = mrdfits(filename, 2, /silent)

	; read in 2d sky
  sky2d = mrdfits(filename, 3, /silent)

  ; create ivar images
  ivar2d_th = 1d/err2d^2
	ivar2d_emp = 1d/sig2d^2

  ; read in wavelength axis
	lambda_1d = sxpar(header,'CRVAL1') + (findgen(sxpar(header,'NAXIS1')) - sxpar(header,'CRPIX1') + 1d) * sxpar(header,'CDELT1')

  ; spatial axis
  y_1d = findgen(sxpar(header,'NAXIS2'))

	; obtain the median profile
  profile = median(spec2d, dimension=1)


	; boxcar extraction
  ; ----------------------------------------------------------------------------

	; define boxcar aperture (include also pixels at the edge)
	w_boxcar = where( y_1d GE min(yrange) and y_1d LE max(yrange), /null )

  ; calculate boxcar extraction
  spec1d_boxcar = total(spec2d[*,w_boxcar], 2, /nan)
  ivar1d_boxcar = 1. / total(1./ivar2d_th[*,w_boxcar], 2, /nan)
	ivar1d_boxcar_emp = 1. / total(1./ivar2d_emp[*,w_boxcar], 2, /nan)

	; boxcar extraction for the sky
	sky2d_boxcar = sky2d[*,w_boxcar]
  sky1d_boxcar = total(sky2d_boxcar, 2, /nan)

  ; make nice output structure
	output_structure = { $
		lambda: lambda_1d, $
		flux: spec1d_boxcar, $
		ivar: ivar1d_boxcar, $
		ivar_empirical: ivar1d_boxcar_emp, $
		sky: sky1d_boxcar }


	; write FITS file
	; ----------------------------------------------------------------------------

	; convert wavelength to angstrom in order to be compatible with SpecPro
	output_structure.lambda *= 1d4

	; make new FITS header, with units
	header_output = header
 	sxdelpar, header_output, 'naxis'
	sxdelpar, header_output, 'naxis1'
	sxdelpar, header_output, 'naxis2'
	sxdelpar, header_output, 'bitpix'
	sxaddpar, header_output, 'TUNIT1', 'Angstrom'
	sxaddpar, header_output, 'TUNIT2', 'electron / (pix s)'

	; write structure to FITS file
	mwrfits, output_structure, output_filename, header_output, /create

	print, ''
  print, '1D spectrum extracted: ' + output_filename
	print, ''


END
