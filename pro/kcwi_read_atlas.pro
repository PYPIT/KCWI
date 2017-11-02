;
; Copyright (c) 2014, California Institute of Technology. All rights
;	reserved.
;+
; NAME:
;	KCWI_READ_ATLAS
;
; PURPOSE:
;	Read the atlas spectrum and convolve to nominal KCWI resolution
;
; CATEGORY:
;	Data reduction for the Keck Cosmic Web Imager (KCWI).
;
; CALLING SEQUENCE:
;	KCWI_READ_ATLAS, Kgeom, Ppar, Refspec, Refwave, Refdisp
;
; INPUTS:
;	Kgeom	- KCWI_GEOM struct from KCWI_TRACE_CBARS and KCWI_EXTRACT_ARCS
;	Ppar	- KCWI_PPAR pipeline parameter struct
;
; OUTPUTS:
;	Refspec	- Atlas reference spectrum
;	Refwave	- Atlas reference spectrum wavelengths
;	Refdisp	- Atlas reference spectrum dispersion in Ang/px
;
; INPUT KEYWORDS:
;
; SIDE EFFECTS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;	Written by:	Don Neill
;	2014-SEP-18	Initial Revision
;-
;
pro kcwi_read_atlas, kgeom, ppar, refspec, refwave, refdisp

pre = 'KCWI_READ_ATLAS'
q=''
;
; init
refspec = -1.
refwave = -1.
refdisp = -1.
;
; check inputs
if kcwi_verify_geom(kgeom,/init) ne 0 then return
if kcwi_verify_ppar(ppar,/init) ne 0 then return
;
; canonical resolution?
sigma = kgeom.atsig
;
; check if file is available
if not file_test(kgeom.refspec,/read,/regular) then begin
	kcwi_print_info,ppar,pre,'Atlas spectrum file not found',kgeom.refspec,$
		format='(a,a)',/error
	return
endif
;
; report the read
kcwi_print_info,ppar,pre,'Reading atlas spectrum in',kgeom.refspec, $
	format='(a,1x,a)'
;
; load the reference atlas spectrum.
rdfits1dspec,kgeom.refspec,refwave,atlas, $
	wavezero=refw0, deltawave=refdisp, refpix=refpix
refspec = atlas>0  
;
; we want to degrade this spectrum to the instrument resolution
xx = findgen(99)-50.0d
gaus = gaussian(xx,[1.0,0.0,sigma])
gaus /= total(gaus)
refspec = convolve(refspec,gaus)
;
return
end		; kcwi_read_atlas
