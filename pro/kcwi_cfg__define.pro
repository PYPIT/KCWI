;
; Copyright (c) 2013, California Institute of Technology. All rights
;	reserved.
;+
; NAME:
;	KCWI_CFG__DEFINE
;
; PURPOSE:
;	This procedure defines the structure for all the observation
;	information for KCWI raw data read from the raw headers.
;
; CATEGORY:
;	Data reduction for the Keck Cosmic Web Imager (KCWI).
;
; CALLING SEQUENCE:
;	kcfg = {KCWI_CFG}
;
; INPUTS:
;	None
;
; OUTPUTS:
;	None
;
; SIDE EFFECTS:
;	defines the KCWI_CFG data structure
;
; PROCEDURE:
;	Provide the automatic structure definition for the KCWI raw header
;	observation structure KCWI_CFG.
;
; EXAMPLE:
;	Instantiate and initialize a single cfg struct for KCWI:
;
;	kcfg = {kcwi_cfg}
;	kcfg = struct_init(kcfg)
;
; NOTES:
;	Keep structure name on a line by itself (see struct_init.pro).
;	Keep tag names 15 chars or less.
;	All tags above the DRIVED PROPERTIES should exactly match header
;		keywords (except case).
;
; MODIFICATION HISTORY:
;	Written by:	Don Neill (neill@caltech.edu)
;	2013-MAY-03	Initial version
;-
pro kcwi_cfg__define
;
tmp = { kcwi_cfg, $
;
; Observation properties
	observer:'', $		; observer
	telescop:'', $		; telescope
	instrume:'', $		; instrument
	object:'', $		; object name
	targname:'', $		; target name
	datepclr:'', $		; UT start of observation (YYYY-MM-DDTHH:MM:SS)
	daterend:'', $		; UT end of observation (YYYY-MM-DDTHH:MM:SS)
	el:-99.d0, $		; Telescope elevation (deg)
	az:-99.d0, $		; Telescope azimuth (deg)
	parang:-99.d0, $	; Parallactic angle (deg)
	epoch:-99.d0, $		; Coordinate epoch
	rotposn:-999.d0, $	; Rotator user angle (deg)
	rotrefan:-999.d0, $	; Rotator reference angle (deg)
	rotmode:'', $		; Rotator mode
	xposure:-9.0, $		; shutter open duration (s)
	telapse:-9.0, $		; dark current duration (s)
	airmass:-1.0, $		; Airmass
	ofname:'', $		; Original file name

;
; Configuration properties
	hatpos:'', $		; hatch position: Open or Closed
	flimagin:'', $		; dome lamp imaging mode: on or off
	flspectr:'', $		; dome lamp spectra mode: on or off
	caltype:'', $		; calibration type: bias, dark, arc, etc.
	frameno:0l, $		; image number
	skyobs:0, $		; sky observation? 0 - object, 1 - sky
	shuffmod:0, $		; is this a Nod & Shuffle observation?
	bgratnam:'', $		; Blue grating id
	bgratnum:0, $		; Blue graing number
	bgrangle:0., $		; Blue grating angle (degrees)
	bgrenc:0l, $		; Blue grating rotator encoder steps
	bfiltnam:'', $		; Blue filter id
	bfiltnum:0, $		; Blue filter number
	bartang:0., $		; Blue articulation angle
	bartenc:0l, $		; Blue cmaera articulation encoder steps
	bcwave:0., $		; Blue central wavelength (Ang)
	bpwave:0., $		; Blue peak wavelength (Ang)
	bfocpos:0l, $		; Blue focus stage encoder steps
	bfocus:0., $		; Blue focus (mm)
	bnasnam:'', $		; Blue mask position name
	bnaspos:0, $		; Blue mask position
	shufrows:0, $		; Number of CCD rows shuffled
	calpnam:'', $		; Cal position name ('Sky', 'Polar', 'Lens')
	callang:0., $		; Cal polarizer angle
	ifunum:0, $		; Slicer number (0-5, -1=unknown)
	ifunam:'', $		; Slicer name ("Small", etc.)
	cwave:0., $		; central wavelength (Ang)
	gratanom:0., $		; grating angle anomoly (degrees)
	wave0:0., $		; blue end  of wavelength range (Ang)
	wave1:0., $		; red end of wavelength range (Ang)
	dwav:0., $		; average dispersion (Ang/pix)
;
; Temperatures
	tmpa1:0., $		; Blue CCD temperature (K)
	tmpa8:0., $		; Optical bench temperature (K)
;
; Cal Lamp properties
	lmp0nam:'', $		; Lamp 0 name
	lmp0stat:0, $		; Lamp 0 status (0 - off, 1 - on)
	lmp0shst:0, $		; Lamp 0 shutter (0 - closed, 1 - open)
	lmp1nam:'', $		; Lamp 0 name
	lmp1stat:0, $		; Lamp 0 status (0 - off, 1 - on)
	lmp1shst:0, $		; Lamp 0 shutter (0 - closed, 1 - open)
	lmp3nam:'', $		; Lamp 0 name
	lmp3stat:0, $		; Lamp 0 status (0 - off, 1 - on)
;
; CCD properties
	biasrn1:3., $		; bias read noise for amp 1 in electrons/pixel
	biasrn2:3., $		; bias read noise for amp 2 in electrons/pixel
	biasrn3:3., $		; bias read noise for amp 3 in electrons/pixel
	biasrn4:3., $		; bias read noise for amp 4 in electrons/pixel
	gain1:0.145, $		; gain for amp 1 in electrons/DN
	gain2:0.145, $		; gain for amp 2 in electrons/DN
	gain3:0.145, $		; gain for amp 3 in electrons/DN
	gain4:0.145, $		; gain for amp 4 in electrons/DN
	gainmul:10, $		; gain multiplier (1, 2, 5, 10)
	nampsxy:'1 1', $	; number of amplifiers in x and y
	ampmode:'LL', $		; amplifier mode: LL,LR,UL,UR,DUP,DLO,QUAD
	nvidinp:0, $		; number of amps (video inputs)
	ccdmode:0, $		; ccd mode: 0 - slow, 1 - fast
;
; Image geometry properties
	naxis:0, $		; number of data axes
	naxis1:0, $		; length of data axis 1
	naxis2:0, $		; length of data axis 2
;
; Wavelength values for reduced data
	crval3:0., $		; wavelength zeropoint
	crpix3:0., $		; wavelength reference pixel
	cdelt3:0., $		; wavelength dispersion Ang/px
;
; DERIVED PROPERTIES (above should match real FITS header keywords)
	ra:-99.d0, $		; RA
	dec:-99.d0, $		; Dec
	juliandate:0.d0, $	; Julian date of observation
	date:'', $		; UT start of observation (YYYY-MM-DDTHH:MM:SS)
	binning:0, $		; is image binned?
	xbinsize:1, $		; binning in x
	ybinsize:1, $		; binning in y
	exptime:-9., $		; exposure time (either xposure or telapse)
;
; Nod-and-shuffle rows (in trimmed image)
	nsskyr0:0, $		; Sky region row 0 (bottom, pix)
	nsskyr1:0, $		; Sky region row 1 (top, pix)
	nsobjr0:0, $		; Object region row 0 (bottom, pix)
	nsobjr1:0, $		; Object region row 1 (top, pix)
;
; Input file properties
	obsfname:'', $		; input observation FITS file name (sans dir)
	obsdir:'', $		; input directory for observation
	obstype:'', $		; observation type: 'zero', 'cal', 'obj', 'std'
	imgnum:0L, $		; image number
;
; Instrument state
	imgtype:'', $		; observation type: bias, dark, arc, etc.
	gratid:'', $		; Grating id name
	gratnum:-1, $		; Gratind id number
	grangle:0., $		; Grating angle (degrees)
	grenc:0l, $		; Grating rotator encoder steps
	filter:'', $		; Filter id name
	filtnum:-1, $		; filter id number
	camang:0., $		; Camera articulation angle (degrees)
	campos:0l, $		; Camera articulation encoder steps
	focpos:0l, $		; Camera focus encoder steps
	focus:0., $		; Camera focus (mm)
	nasmask:0, $		; is the Nod & Shuffle mask deployed?
	lampname:'', $		; illumination source
;
; Master group properties
	groupnum:-1l, $		; group image number
	nimages:0, $		; number of images in group
	grouplist:'', $		; range list of images in group
	groupfile:'', $		; filename of grouped image
	grouppar:'', $		; parameter file for grouped image
;
; Status
	initialized:0, $	; 0 - no, 1 - yes
;
; timestamp
	timestamp:0.d0 $	; timestamp for struct (unix seconds)
	}
end
