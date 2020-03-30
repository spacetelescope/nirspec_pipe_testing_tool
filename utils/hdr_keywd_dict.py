import collections

'''
This script contains an ordered dictionary of the key words and their corresponding expected format and/or
allowed values.

* The information for the keywords was obtained from emails in the pipeline group (i.e. basecamp) and the information
contained in the online dictionary located at:
https://iwjwdmsdauiwebv.stsci.edu/portal/Mashup/Clients/jwkeywords/
and the pipeline GitHub page:
https://github.com/STScI-JWST/jwst/blob/master/jwst/datamodels/schemas/wcsinfo.schema.yaml
and the Confluence page:
https://confluence.stsci.edu/display/JWSTPWG/Build+7.1+Updates+for+Testing+and+Reports

'''

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed


# define dictionary
keywd_dict = collections.OrderedDict()

# Standard parameters
keywd_dict['SIMPLE']  = ['T', 'F'] # Written by IDL
keywd_dict['BITPIX']  = [int] # Bits per data value, e.g. 8
keywd_dict['NAXIS']   = [int] # Number of data array dimensions,  e.g. 0
keywd_dict['EXTEND']  = ['T', 'F'] # File may contain standard extensions
keywd_dict['NEXTEND'] = [int] # Number of standard extensions, e.g. 3

keywd_dict['DATAMODL']= ['Level1bModel'] # Type of data model
keywd_dict['TELESCOP']= ['JWST'] # Telescope used to acquire data
keywd_dict['INSTRUME']= ['NIRSPEC'] # Identifier for instrument used to acquire data
keywd_dict['RADESYS'] = [str] # Coordinate reference frame for RA and Dec, e.g. 'ICRS'

keywd_dict['DATE']    = [str] # UTC date file created, e.g. 2016-07-12T12:43:27
keywd_dict['ORIGIN']  = [str] # Institution responsible for creating FITS file, e.g. STScI
keywd_dict['FILENAME']= [str] # Name of file, e.g. jwtest1001001_01101_00001_NRS1_uncal.fits
keywd_dict['FILETYPE']= [str] # Type of data found in file, e.g. UNCALIBRATED
keywd_dict['SDP_VER']= [str] # Data processing software version number, e.g. 0.1.1

# Observation identifiers
keywd_dict['DATE-OBS']= [str] # UTC date at start of exposure, e.g. 2013-01-19
keywd_dict['TIME-OBS']= [str] # UTC time at start of exposure, e.g. 18:23:34.230
keywd_dict['OBS_ID']  = [str] # full programmatic observation identifier, e.g. SLIT-COMBO-001
keywd_dict['VISIT_ID']= [str] # visit identifier
keywd_dict['PROGRAM'] = [str] # program number, e.g. test1
keywd_dict['OBSERVTN']= [str] # observation number, e.g. 001
keywd_dict['VISIT']   = [str] # visit number, e.g. 001
keywd_dict['VISITGRP']= [str] # visit group identifier, e.g. 01
keywd_dict['SEQ_ID']  = [str] # parallel sequence identifier, e.g. 1
keywd_dict['ACT_ID']  = [str] # activity identifier, e.g. 01
keywd_dict['EXPOSURE']= [str] # exposure request number, e.g. 00001

# Exposure parameters
keywd_dict['DETECTOR']= [str] # name of detector used to acquire data, e.g. 'NRS1'
keywd_dict['NINTS']   = [int] # number of integrations within exposure, e.g. 1
keywd_dict['NGROUPS'] = [int] # number of groups within integration, e.g. 20
keywd_dict['ZEROFRAM']= ['T', 'F'] # boolean represented with string, T if a zero frame was read separately
keywd_dict['READPATT']= [str] # readout pattern, e.g. 'NRSRAPID'
keywd_dict['DATAPROB']= ['T', 'F'] # boolean represented with string, T if science telemetry indicated any problems

# Program information
keywd_dict['TITLE']   = [str] # proposal title
keywd_dict['PI_NAME'] = [str] # Name of principal investigator, e.g. UNKNOWN
keywd_dict['CATEGORY']= [str] # program category, e.g. N/A
keywd_dict['SUBCAT']  = [str] # program sub-category, e.g. N/A
keywd_dict['SCICAT']  = [str] # science category assigned during TAC process, e.g. N/A
keywd_dict['CONT_ID'] = [int] # continuation of the specified previous program, e.g. 1

# Observation information
keywd_dict['TEMPLATE']= [str] # proposal instruction template used, e.g. N/A
keywd_dict['OBSLABEL']= [str] # proposer label for observation, e.g. #TODO

# Visit information
keywd_dict['VISITYPE']= ['PRIME', 'PARALLEL', 'GENERIC'] # type of visit (prime or parallel)
keywd_dict['VSTSTART']= [str] # UTC visit start time, e.g. 2013-01-19T18:27:22
keywd_dict['NEXPOSUR']= [int] # total number of exposures in visit, e.g. 1
keywd_dict['INTARGET']= ['T', 'F'] # boolean represented with string, T if at least one exposure in visit is internal
keywd_dict['TARGOOPP']= ['T', 'F'] # boolean represented with string, visit scheduled as target of opportunity

# Exposure information
keywd_dict['PNTG_SEQ']= [int] # pointing sequence number, e.g. 1
keywd_dict['EXPCOUNT']= [int] # count of the exposures in visit, e.g. 1
keywd_dict['EXP_TYPE']= ['NRS_WATA', 'NRS_MSATA', 'NRS_TACONFIRM', 'NRS_CONFIRM', ' NRS_IMAGE', ' NRS_FOCUS',
                         'NRS_MIMF', 'NRS_FIXEDSLIT', 'NRS_AUTOWAVE',
                         'NRS_IFU', 'NRS_MSASPEC', ' NRS_AUTOFLAT', ' NRS_DARK',
                         'NRS_LAMP', 'NRS_BOTA', 'NRS_BRIGHTOBJ'] # type of data in exposure

# Target information
keywd_dict['TARGPROP']= [str] # proposer's name for the target
keywd_dict['TARGNAME']= [str] # standard astronomical catalog name for target, e.g. 'NGC 104 '
keywd_dict['TARGTYPE']= ['FIXED', 'MOVING', 'GENERIC'] # fixed target, moving target, or generic target
keywd_dict['TARG_RA'] = [float] # target RA computed at time of exposure, e.g. 0.0
keywd_dict['TARGURA'] = [float] # target RA uncertainty, e.g. 0.0
keywd_dict['TARG_DEC']= [float] # target DEC computed at time of exposure, e.g. 0.0
keywd_dict['TARRUDEC']= [float] # target Dec uncertainty, e.g. 0.0
keywd_dict['PROP_RA'] = [float] # proposer specified RA for the target, e.g. 0.0
keywd_dict['PROP_DEC']= [float] # proposer specified Dec for the target, e.g. 0.0
keywd_dict['MU_EPOCH']= [str] # epoch of proper motion values for RA and Dec, proposer specified, e.g. 2000.0

# Exposure times
keywd_dict['EXPSTART']= [float] # UTC exposure start time (MJD), e.g. 56311.76636840278
keywd_dict['EXPMID']  = [float] # UTC exposure mid time (MJD), e.g. 56311.76636840278
keywd_dict['EXPEND']  = [float] # UTC exposure end time (MJD), e.g. 56311.76763953704

# Exposure time parameters
keywd_dict['NSAMPLES']= [int] # number of A/D samples per pixel, e.g. 1
keywd_dict['NFRAMES'] = [int] # number of frames in group, e.g. 1
keywd_dict['GROUPGAP']= [int] # number of frames dropped between groups, e.g. 10
keywd_dict['TSAMPLE'] = [int] # delta time between samples in microseconds, e.g. 10
keywd_dict['NRESETS']  = [int] # number of resets between integrations, e.g. 1
keywd_dict['NRSTSTRT']= [int] # number of extra resets at start of exposure, e.g. 1
keywd_dict['TFRAME']  = [float] # [seconds] time between frames, e.g. 5.49132
keywd_dict['TGROUP']  = [float] # [seconds] time between groups, e.g 5.49132
keywd_dict['EFFINTTM']= [float] # [seconds] effective integration time, e.g. 104.33508
keywd_dict['EFFEXPTM']= [float] # [seconds] effective exposure time, e.g. 104.33508
keywd_dict['DURATION']= [float] # [seconds] total duration of exposuree.g. -1.0

# Subarray parameters
keywd_dict['SUBARRAY']= ['1024X16', '128X128', '128X2048', '2048X128', '2048X64', '32X32', '64X2048', '8X8',
                         'ALLSLITS', 'BRIGHTSKY', 'FULL', 'GENERIC', 'MASK1065', 'MASK1140', 'MASK1550',
                         'MASKLYOT', 'S1600A1', 'S200A1', 'S200A2', 'S200B1', 'S400A1', 'SLITLESSPRISM',
                         'STRIPE', 'SUB1024A', 'SUB1024B', 'SUB128', 'SUB16', 'SUB160', 'SUB160P',
                         'SUB1A', 'SUB1B', 'SUB2048', 'SUB256', 'SUB32', 'SUB320', 'SUB400P', 'SUB512',
                         'SUB64', 'SUB640', 'SUB64P', 'SUB80', 'SUB96', 'SUBGRISM128', 'SUBGRISM256',
                         'SUBGRISM64', 'SUBSTRIP256', 'SUBSTRIP96', 'SUBPRISM', 'WFSS128C', 'WFSS128R',
                         'WFSS64C', 'WFSS64R', 'N/A'] # name of subarray used
keywd_dict['SUBSTRT1']= [int] # starting pixel number in the SIAS x direction, e.g. 1
keywd_dict['SUBSIZE1']= [int] # number of pixels in the SIAS x direction, e.g. 2048
keywd_dict['SUBSTRT2']= [int] # starting pixel number in the SIAS y direction, e.g. 897
keywd_dict['SUBSIZE2']= [int] # number of pixels in the SIAS y direction, e.g. 256
keywd_dict['FASTAXIS']= [1, 2, -1, -2] # Direction of fast readout (+/-2 X direction, +/-
keywd_dict['SLOWAXIS']= [1, 2, -1, -2] # Direction of slow readout (+/-1 X direction, +/-

# NIRSpec configuration (NIRSpec only)
keywd_dict['FILTER']  = ['CLEAR', 'F070LP', 'F100LP', 'F110W', 'F140X', 'F170LP', 'F290LP',
                         'P750L', 'NULL', 'OPAQUE'] # name of the filter element used
keywd_dict['GRATING'] = ['G140M', 'G235M', 'G395M', 'G140H', 'G235H', 'G395H', 'PRISM',
                         'MIRROR', 'NULL', 'N/A', 'ANY'] # name of the grating element used
keywd_dict['GWA_XTIL']= [float] # grating y tilt, e.g. 0.3622055649757385
keywd_dict['GWA_YTIL']= [float] # grating y tilt, e.g. 0.1346436440944672
keywd_dict['GWA_TILT']= [float] # GWA temperature ave [K], e.g. 4.028447479156018e+01
keywd_dict['FXD_SLIT']= ['NONE', 'S200A1', 'S200A2', 'S200B1', 'S400A1', 'S1600A1', 'NULL'] # name of fixed slit aperture used
keywd_dict['MSASTATE']= ['CONFIGURED', 'LAUNCHLOCK_ALLCLOSED', 'PRIMARYPARK_ALLOPEN', 'PRIMARYPARK_ALLCLOSED',
                         'PRIMARYPARK_CONFIGURED'] # state of MSA: all open, all closed, configured
keywd_dict['FOCUSPOS']= [int] # [mm] focus position for NIRSpec, e.g. 1

# keywords added for build 7.1 from https://confluence.stsci.edu/display/JWSTPWG/Build+7.1+Updates+for+Testing+and+Reports
keywd_dict['NRS_NORM']= [int]  # Number of normal pixels in IRS2 readout, used also in IFU data
keywd_dict['NRS_REF']= [int]  # Number of reference pixels in IRS2 readout, used also in IFU data
keywd_dict['FRMDIVSR']= [int]  # Divisor applied to frame-averaged groups
keywd_dict['XOFFSET']= [float] # x offset from pattern starting position
keywd_dict['YOFFSET']= [float] # y offset from pattern starting position

# added for build 7.3
keywd_dict['PATT_NUM']= [int] # position number within primary pattern
keywd_dict['PATTSIZE']= [str, 'NoneType'] # [arcsec] primary dither pattern size: SMALL, MEDIUM, LARGE, None

# NIRSpec MSA supporting files (NIRSpec MSA only)
keywd_dict['MSAMETFL']= [str] # MSA configuration file name, e.g. 'N/A'
keywd_dict['MSAMETID']= [int] # MSA meta data ID for the exposure, e.g. 1

# lamp configuration
keywd_dict['LAMP']    = [str] # internal lamp state, e.g. 'ARGON'
keywd_dict['LAMPMODE']= [str] # NIRSpec internal lamp exposures, possible values: BRIGHTOBJ, FIXEDSLIT,
                              # GRATING-ONLY, IFU, MSASPEC, NULL

# Guide star information
keywd_dict['GS_ORDER']= [str, int] # index of guide star, e.g. 'N/A'
keywd_dict['GSSTRTTM']= [str, float] # UTC start time of guide star acquisition, e.g. 'N/A'
keywd_dict['GSENDTIM']= [str, float] # UTC end time of guide star acquisition, e.g. 'N/A'
keywd_dict['GDSTARID']= [str] # guide star identifier, e.g. 'N/A'
keywd_dict['GS_RA']   = [str, float] # guide star right ascension, e.g. 'N/A'
keywd_dict['GS_DEC']  = [str, float] # guide star declination, e.g. 'N/A'
keywd_dict['GSURA']   = [str, float] # guide star right ascension uncertainty, e.g. 'N/A'
keywd_dict['GSUDEC']  = [str, float] # guide star declination uncertainty, e.g. 'N/A'
keywd_dict['GS_MAG']  = [str, float] # guide star magnitude in FGS detector, e.g. 'N/A'
keywd_dict['GS_UMAG']  = [str, float] # guide star magnitude uncertainty, e.g. 'N/A'
keywd_dict['PCS_MODE']= [str, float] # Pointing Control System mode, e.g. 'N/A'

# JWST ephemeris information
keywd_dict['REFFRAME']= [str] # ephemeris coordinate system, e.g. 'N/A'
keywd_dict['EPH_TIME']= [float] # [sec] UTC time from ephemeris start time, e.g. 0.0
keywd_dict['JWST_X']  = [float] # [km] X spatial coordinate of JWST, e.g. 0.0
keywd_dict['JWST_Y']  = [float] # [km] Y spatial coordinate of JWST, e.g. 0.0
keywd_dict['JWST_Z']  = [float] # [km] Z spatial coordinate of JWST, e.g. 0.0
keywd_dict['JWST_DX'] = [float] # [km/sec] X component of JWST velocity vector, e.g. 0.0
keywd_dict['JWST_DY'] = [float] # [km/sec] Y component of JWST velocity vector, e.g. 0.0
keywd_dict['JWST_DZ'] = [float] # [km/sec] Z component of JWST velocity vector, e.g. 0.0

# Spacecraft pointing information
keywd_dict['PA_V3']   = [str, float] # [deg] position angle of V3-axis of JWST, e.g. 'N/A'
keywd_dict['RA_V1']   = [str, float] # [deg] RA of telescope V1 axis, e.g. 'N/A'
keywd_dict['DEC_V1']  = [str, float] # [deg] Dec of telescope V1 axis, e.g. 'N/A'

# Aperture pointing information
keywd_dict['APERNAME']= [str] # mnemonic for PDB science aperture used, e.g. #TODO
keywd_dict['PA_APER'] = [float] # [deg] position angle of aperture used, e.g. -999.0

# WCS parameters, these will be added to the science extension
keywd_dict['wcsinfo'] = {
                            'WCSAXES' : [int], # number of World Coordinate System axes, e.g. 3
                            'CRPIX1' : [int, float], # x-coordinate of the reference pixel, e.g. 1024.0
                            'CRPIX2' : [int, float], # y-coordinate of the reference pixel, e.g. 128.0
                            'CRPIX3' : [int, float], # z-coordinate of the reference pixel
                            'CRVAL1' : [float], # RA at the reference pixel (degrees), 5.3196
                            'CRVAL2' : [float], # Dec at the reference pixel (degrees), e.g. -72.98605000000001
                            'CRVAL3' : [float], # Wavelength at the reference pixel (microns), e.g. 2.5
                            'CTYPE1' : ['RA---TAN'], # first axis coordinate type
                            'CTYPE2' : ['DEC--TAN'], # second axis coordinate type
                            'CTYPE3' : [str], # third axis coordinate type, e.g. WAVE
                            'CUNIT1' : ['deg'], # units for first axis
                            'CUNIT2' : ['deg'], # units for seconds axis
                            'CUNIT3' : ['um'], # units for third axis, e.g. 'micron'
                            'CDELT1' : [float], # increment per pixel, axis 1, e.g. 0.12
                            'CDELT2' : [float], # increment per pixel, axis 2, e.g. 0.12
                            'CDELT3' : [float], # increment per pixel, axis 3, e.g. 0.000672
                            'PC1_1' : [float], # linear transformation matrix element, e.g. -1.0
                            'PC1_2' : [float], # linear transformation matrix element, e.g. 0.0
                            'PC1_3' : [float], # linear transformation matrix element, e.g. 0.0
                            'PC2_1' : [float], # linear transformation matrix element, e.g. 0.0
                            'PC2_2' : [float], # linear transformation matrix element, e.g. 1.0
                            'PC3_1' : [float], # linear transformation matrix element, e.g. 1.0
                            'PC3_2' : [float], # linear transformation matrix element, e.g. 0.0
                            'PC3_3' : [float], # linear transformation matrix element, e.g. 0.0
                            'S_REGION' : [str], # spatial extent of the observation, e.g. 'N/A'
                            'WAVSTART' : [float], # lower bound of the default wavelength range
                            'WAVEND' : [float], # upper bound of the default wavelength range
                            'SPORDER' : [int], # default spectral order
                            'V2_REF' : [float], # location of the aperture reference point in V2 (arcsec): 100-400 arcsec
                            'V3_REF' : [float], # location of the aperture reference point in V3 (arcsec): -100 to -400 arcsec
                            'VPARITY' : [int], # Relative sense of rotation between Ideal xy and V2V3
                            'V3I_YANG' : [float], # Angle from V3 axis to Ideal y axis (deg)
                            'RA_REF' : [float], # RA at the reference point (deg): 0 < RA < 360
                            'DEC_REF' : [float], # Dec at the reference point (deg): -90 < Dec < +90
                            'ROLL_REF' : [float] # Roll angle at the reference point (deg), e.g. 5.3196
}


# Velocity aberration correction
keywd_dict['DVA_RA']  = [float] # velocity aberration correction RA offset, e.g. 0.0
keywd_dict['DVA_DEC'] = [float] # velocity aberration correction Dec offset, e.g. 0.0
keywd_dict['VA_SCALE']=  [float] # velocity aberration scale factor, e.g. 1.0

#  Time related keywords
keywd_dict['BARTDELT']= [float] # calculated Barycentric time correction, e.g. 0.0
keywd_dict['BSTRTIME']= [float] # Solar System Barycentric exposure start time, e.g. 0.0
keywd_dict['BENDTIME']= [float] # Solar System Barycentric exposure end time, e.g. 0.0
keywd_dict['BMIDTIME']= [float] # Solar System Barycentric exposure mid time, e.g. 0.0
keywd_dict['HELIDELT']= [float] # calculated Heliocentric time correction, e.g. 0.0
keywd_dict['HSTRTIME']= [float] # Heliocentric exposure start time in MJD, e.g. 0.0
keywd_dict['HENDTIME']= [float] # Heliocentric exposure end time in MJD, e.g. 0.0
keywd_dict['HMIDTIME']= [float] # Heliocentric exposure mid time in MJD, e.g. 0.0
