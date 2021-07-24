import argparse
import collections
import sys
import pysiaf
import numpy as np
from astropy.io import fits
from astropy.time import Time
from datetime import timedelta

# import subarray dictionary
from .dict_info import subarray_dict as subdict

'''
This script contains mapping of pipeline keywords to the corresponding header keyword of the simulations file. The
script modifies the uncommented keywords in the dictionary.

This script is to be used AFTER the dummy header keyword values have been added, and it ONLY uses the keywords
that require a change in value, the others are commented out (but remain here for completion).

Keyword value options used in the dictionary:
    - 'specific_string' = this is a string that will not change from simulation to simulation
    - 'set_to_given_string' = will look for the given value when script is run
    - 'N/A' = this is a non applicable (not too important) keyword for the simulations
    - 'header_ext:keyword' = get the value from the header extension of the file with the keyword after the :
    - 'primary_ext:keyword' = this keyword is spelled differently in the IPS primary header, code will look for the string after the :
    - 'calculation' = these are numbers that need to be calculated

'''

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Jan 2020 - Version 1.0: initial version completed

# Define dictionary
stsci2ips_dict = collections.OrderedDict()

# Standard parameters
stsci2ips_dict['SIMPLE'] = 'primary_ext:SIMPLE'  # conforms to FITS standards
stsci2ips_dict['BITPIX'] = 'primary_ext:BITPIX'  # array data type
stsci2ips_dict['NAXIS'] = 'primary_ext:NAXIS'  # sumber of data array dimensions
stsci2ips_dict['EXTEND'] = 'primary_ext:EXTEND'  # File may contain standard extensions

# Level 3 Schema Metadata
stsci2ips_dict['DATE'] = 'primary_ext:DATE'  # UTC date file was created
stsci2ips_dict['ORIGIN'] = 'specific_string:NIRSpec Instrument Performance Simulator'  # Institution responsible for
# creating FITS file
stsci2ips_dict['FILENAME'] = 'specific_string:N/A'  # Name of uncal file
stsci2ips_dict['FILETYPE'] = 'specific_string:UNCALIBRATED'  # Type of data found in file
stsci2ips_dict['SDP_VER'] = 'specific_string:B7.4'  # Data processing software version number
stsci2ips_dict['SDP_VCS'] = 'specific_string:N/A'  # Calibration software version control sys number
# stsci2ips_dict['DATAMODL']= 'ImageModel' # Type of data model
# stsci2ips_dict['TELESCOP']= 'JWST' # Telescope used to acquire data

# Program information
stsci2ips_dict['TITLE'] = 'set_to_given_string'  # proposal title
# stsci2ips_dict['PI_NAME'] = 'UNKNOWN' # Name of principal investigator, e.g. UNKNOWN
# stsci2ips_dict['CATEGORY']= 'N/A' # program category, e.g. N/A
# stsci2ips_dict['SUBCAT']  = 'N/A' # program sub-category, e.g. N/A
# stsci2ips_dict['SCICAT']  = 'N/A' # science category assigned during TAC process, e.g. N/A
# stsci2ips_dict['CONT_ID'] = 1 # int, continuation of the specified previous program, e.g. N/A

# Observation identifiers
stsci2ips_dict['DATE-OBS'] = 'header_ext:DATE-OBS'  # UTC date at start of exposure, e.g. 2013-01-19
stsci2ips_dict['TIME-OBS'] = 'header_ext:TIME-OBS'  # UTC time at start of exposure, e.g. 18:23:34.230
stsci2ips_dict['OBS_ID'] = 'primary_ext:OBSID'  # full programmatic observation identifier
# stsci2ips_dict['VISIT_ID']= '' # visit identifier
# stsci2ips_dict['PROGRAM'] = '12345' # program number, e.g. '12345'
# stsci2ips_dict['OBSERVTN']= '001' # observation number, e.g. 001
# stsci2ips_dict['VISIT']   = '001' # visit number, e.g. 001
# stsci2ips_dict['VISITGRP']= '01' # visit group identifier, e.g. 01
# stsci2ips_dict['SEQ_ID']  = '1' # parallel sequence identifier, e.g. 1
# stsci2ips_dict['ACT_ID']  = '01' # activity identifier, e.g. 01
# stsci2ips_dict['EXPOSURE']= '00001' # exposure request number, e.g. 00001
# stsci2ips_dict['TEMPLATE']= 'N/A' # proposal instruction template used, e.g. N/A
# stsci2ips_dict['OBSLABEL']= '#TODO' # proposer label for observation, e.g. #TODO

# Visit information
# stsci2ips_dict['VISITYPE']= 'GENERIC' # type of visit (prime or parallel)
stsci2ips_dict['VSTSTART'] = 'primary_ext:DATE'  # UTC visit start time, e.g. 2013-01-19T18:27:22
# stsci2ips_dict['NEXPOSUR']= 1 # total number of exposures in visit, e.g. 1
# stsci2ips_dict['INTARGET']= True # boolean represented with string, T if at least one exposure in visit is internal
# stsci2ips_dict['TARGOOPP']= False # boolean represented with string, visit scheduled as target of opportunity

# Target information
# stsci2ips_dict['TARGPROP']= 'UNKNOWN' # proposer's name for the target
stsci2ips_dict['TARGNAME'] = 'set_to_given_string'  # standard astronomical catalog name for target, e.g. 'NGC 104 '
# stsci2ips_dict['TARGTYPE']= 'FIXED' # fixed target, moving target, or generic target, options are:
# 'FIXED', 'MOVING', 'GENERIC'
stsci2ips_dict['TARG_RA'] = 'primary_ext:RA_REF'  # target RA computed at time of exposure, e.g. 0.0
# stsci2ips_dict['TARGURA'] = 0.0 # target RA uncertainty, e.g. 0.0
stsci2ips_dict['TARG_DEC'] = 'primary_ext:DEC_REF'  # target DEC computed at time of exposure, e.g. 0.0
# stsci2ips_dict['TARRGURA']= 0.0 # target Dec uncertainty, e.g. 0.0
# stsci2ips_dict['MU_EPOCH']= '2000-01-01T00:00:00.000' # epoch of proper motion values for RA and Dec, proposer
# specified, e.g. 2000.0
# stsci2ips_dict['PROP_RA'] = 0.0 # proposer specified RA for the target, e.g. 0.0
# stsci2ips_dict['PROP_DEC']= 0.0 # proposer specified Dec for the target, e.g. 0.0

# Instrument configuration information
# stsci2ips_dict['INSTRUME']= 'NIRSPEC' # Identifier for instrument used to acquire data
stsci2ips_dict['DETECTOR'] = 'primary_ext:DET'  # name of detector used to acquire data
stsci2ips_dict['FILTER'] = 'primary_ext:FWA_POS'  # name of the grating element used, options are:
#                                'CLEAR', 'F070LP', 'F100LP', 'F110W', 'F140X', 'F170LP', 'F290LP', 'P750L', 'NULL'
stsci2ips_dict['GRATING'] = 'primary_ext:GWA_POS'  # name of grating used, options are:
#                                   'G140M', 'G235M', 'G395M', 'G140H', 'G235H', 'G395H', 'PRISM',
#                                   'MIRROR', 'NULL', 'N/A', 'ANY'
# stsci2ips_dict['FXD_SLIT']= 'NONE' # name of fixed slit aperture used, options are: 'NONE', 'S200A1',
# 'S200A2', 'S200B1', 'S400A1', 'S1600A1', 'NULL'
# stsci2ips_dict['FOCUSPOS']= 1 # int, [mm] focus position for NIRSpec, e.g. 0
# stsci2ips_dict['MSASTATE']= 'PRIMARYPARK_CONFIGURED' # state of MSA, options are:
#                                       'CONFIGURED', 'LAUNCHLOCK_ALLCLOSED', 'PRIMARYPARK_ALLOPEN',
#                                       'PRIMARYPARK_ALLCLOSED', 'PRIMARYPARK_CONFIGURED'
# stsci2ips_dict['MSAMETFL']= 'N/A' # MSA configuration file name, e.g. 'blah.fits'
# stsci2ips_dict['MSAMETID']= 1 # int, MSA meta data ID for the exposure, e.g. 1
stsci2ips_dict['LAMP'] = 'primary_ext:CAA_LAMP'  # internal lamp state, e.g. 'ARGON', 'NONE' the combinations are:
# F100LP LINE1  FLAT1
# F170LP LINE2  FLAT2
# F290LP LINE3  FLAT3
# CLEAR  LINE4  FLAT5
# F070LP FLAT4
# stsci2ips_dict['LAMPMODE']= 'NULL' # NIRSpec internal lamp exposures, possible values: BRIGHTOBJ, FIXEDSLIT,
# GRATING-ONLY, IFU, MSASPEC, NULL
stsci2ips_dict['GWA_XTIL'] = 'header_ext:GWA_XTIL'  # grating y tilt, e.g. 0.3622055649757385
stsci2ips_dict['GWA_YTIL'] = 'header_ext:GWA_YTIL'  # grating y tilt, e.g. 0.1346436440944672
# stsci2ips_dict['GWA_TILT']= 4.028447479156018e+01 # GWA temperature ave [K], e.g. 4.028447479156018e+01

# Exposure parameters
# stsci2ips_dict['PNTG_SEQ']= 1 # pointing sequence number, e.g. 1
# stsci2ips_dict['EXPCOUNT']= 1 # count of the exposures in visit, e.g. 1
# stsci2ips_dict['EXP_TYPE']= 'NRS_MSASPEC' # type of data in exposure, options are:
#                         'NRS_TASLIT', 'NRS_TACQ', 'NRS_TACONFIRM', 'NRS_CONFIRM', 'NRS_FIXEDSLIT', 'NRS_AUTOWAVE',
#                         'NRS_IFU', 'NRS_MSASPEC', ' NRS_AUTOFLAT', ' NRS_IMAGE', ' NRS_FOCUS', ' NRS_DARK',
#                         'NRS_LAMP', 'NRS_BOTA', 'NRS_BRIGHTOBJ'
stsci2ips_dict['EXPSTART'] = 'calculation'  # UTC exposure start time (MJD), e.g. 56311.76636840278
stsci2ips_dict['EXPMID'] = 'calculation'  # UTC exposure mid time (MJD), e.g. 56311.76636840278
stsci2ips_dict['EXPEND'] = 'calculation'  # UTC exposure end time (MJD), e.g. 56311.76763953704
stsci2ips_dict['READPATT'] = 'primary_ext:READOUT'  # readout pattern
# stsci2ips_dict['NINTS']   = 1 # number of integrations within exposure, e.g. 1
stsci2ips_dict['NGROUPS'] = 'header_ext:NGROUP'  # number of groups within integration, e.g. 20
stsci2ips_dict['NFRAMES'] = 'header_ext:NFRAME'  # number of frames in group, e.g. 1
# stsci2ips_dict['FRMDIVSR']= 1  # integer, Divisor applied to frame-averaged groups
# stsci2ips_dict['GROUPGAP']= 0 # number of frames dropped between groups, e.g. 10
# stsci2ips_dict['NSAMPLES']= 1 # number of A/D samples per pixel, e.g. 1
stsci2ips_dict['TSAMPLE'] = 'calculation'  # delta time between samples in microseconds, e.g. 10
# stsci2ips_dict['TFRAME']  = 10.73676 # [seconds] time between frames, e.g. 5.49132
stsci2ips_dict['TGROUP'] = 'calculation'  # [seconds] time between groups, e.g 5.49132
stsci2ips_dict['EFFINTTM'] = 'calculation'  # [seconds] effective integration time, e.g. 104.33508
stsci2ips_dict['EFFEXPTM'] = 'calculation'  # [seconds] effective exposure time, e.g. 104.33508
stsci2ips_dict['DURATION'] = 'calculation'  # [seconds] total duration of exposure.g. -1.0
# stsci2ips_dict['NRSTSTRT']= 1 # number of extra resets at start of exposure, e.g. 1
# stsci2ips_dict['NRESETS']  = 1 # number of resets between integrations, e.g. 1
# stsci2ips_dict['ZEROFRAM']= False # boolean represented with string, T if a zero frame was read separately
# stsci2ips_dict['DATAPROB']= False # boolean represented with string, T if science telemetry indicated any problems
# stsci2ips_dict['NRS_NORM']= 16  # int, Number of normal pixels in IRS2 readout, used also in IFU data
# stsci2ips_dict['NRS_REF']= 4  # int, Number of reference pixels in IRS2 readout, used also in IFU data

# Subarray parameters
stsci2ips_dict['SUBARRAY'] = 'primary_ext:SUBARRAY'  # name of subarray used, options are:
#                          '1024X16', '128X128', '128X2048', '2048X128', '2048X64', '32X32', '64X2048', '8X8',
#                          'ALLSLITS', 'BRIGHTSKY', 'FULL', 'GENERIC', 'MASK1065', 'MASK1140', 'MASK1550',
#                          'MASKLYOT', 'S1600A1', 'S200A1', 'S200A2', 'S200B1', 'S400A1', 'SLITLESSPRISM',
#                          'STRIPE', 'SUB1024A', 'SUB1024B', 'SUB128', 'SUB16', 'SUB160', 'SUB160P',
#                          'SUB1A', 'SUB1B', 'SUB2048', 'SUB256', 'SUB32', 'SUB320', 'SUB400P', 'SUB512',
#                          'SUB64', 'SUB640', 'SUB64P', 'SUB80', 'SUB96', 'SUBGRISM128', 'SUBGRISM256',
#                          'SUBGRISM64', 'SUBSTRIP256', 'SUBSTRIP96', 'SUBPRISM', 'WFSS128C', 'WFSS128R',
#                          'WFSS64C', 'WFSS64R', 'N/A']
# stsci2ips_dict['SUBSTRT1']= 1 # starting pixel number in the SIAS x direction, e.g. 1
# stsci2ips_dict['SUBSIZE1']= 2048 # number of pixels in the SIAS x direction, e.g. 2048
# stsci2ips_dict['SUBSTRT2']= 1 # starting pixel number in the SIAS y direction, e.g. 897
# stsci2ips_dict['SUBSIZE2']= 2048 # number of pixels in the SIAS y direction, e.g. 256
# stsci2ips_dict['FASTAXIS']= 2 # Direction of fast readout, options are: 1, 2, -1, -2
# stsci2ips_dict['SLOWAXIS']= 1 # Direction of slow readout, options are: 1, 2, -1, -2

# Dither information
# stsci2ips_dict['XOFFSET']= 0.0  # float, x offset from pattern starting position
# stsci2ips_dict['YOFFSET']= 0.0  # float, y offset from pattern starting position

# added for build 7.3
# stsci2ips_dict['PATT_NUM']= 1 # position number within primary pattern
# stsci2ips_dict['PATTSIZE']= "SMALL"  # [arcsec] primary dither pattern size: SMALL, MEDIUM, LARGE, None

# JWST ephemeris information
# stsci2ips_dict['REFFRAME']= 'N/A' # ephemeris coordinate system, e.g. 'N/A'
# stsci2ips_dict['EPH_TIME']= 0.0 # [sec] UTC time from ephemeris start time, e.g. 0.0
# stsci2ips_dict['JWST_X']  = 0.0 # [km] X spatial coordinate of JWST, e.g. 0.0
# stsci2ips_dict['JWST_Y']  = 0.0 # [km] Y spatial coordinate of JWST, e.g. 0.0
# stsci2ips_dict['JWST_Z']  = 0.0 # [km] Z spatial coordinate of JWST, e.g. 0.0
# stsci2ips_dict['JWST_DX'] = 0.0 # [km/sec] X component of JWST velocity vector, e.g. 0.0
# stsci2ips_dict['JWST_DY'] = 0.0 # [km/sec] Y component of JWST velocity vector, e.g. 0.0
# stsci2ips_dict['JWST_DZ'] = 0.0 # [km/sec] Z component of JWST velocity vector, e.g. 0.0

# Aperture pointing information
stsci2ips_dict['APERNAME'] = 'primary_ext:APERTURE'  # mnemonic for PDB science aperture used, e.g. #TODO

# Velocity aberration correction
# stsci2ips_dict['DVA_RA']  = 0.0 # velocity aberration correction RA offset, e.g. 0.0
# stsci2ips_dict['DVA_DEC'] = 0.0 # velocity aberration correction Dec offset, e.g. 0.0

#  Time information
# stsci2ips_dict['BARTDELT']= 0.0 # calculated Barycentric time correction, e.g. 0.0
# stsci2ips_dict['BSTRTIME']= 0.0 # Solar System Barycentric exposure start time, e.g. 0.0
# stsci2ips_dict['BENDTIME']= 0.0 # Solar System Barycentric exposure end time, e.g. 0.0
# stsci2ips_dict['BMIDTIME']= 0.0 # Solar System Barycentric exposure mid time, e.g. 0.0
# stsci2ips_dict['HELIDELT']= 0.0 # calculated Heliocentric time correction, e.g. 0.0
# stsci2ips_dict['HSTRTIME']= 0.0 # Heliocentric exposure start time in MJD, e.g. 0.0
# stsci2ips_dict['HENDTIME']= 0.0 # Heliocentric exposure end time in MJD, e.g. 0.0

# Guide star information
# stsci2ips_dict['GS_ORDER']= 1 # index of guide star, e.g. 'N/A'
# stsci2ips_dict['GSSTRTTM']= 'N/A' # UTC start time of guide star acquisition, e.g. 'N/A'
# stsci2ips_dict['GSENDTIM']= 'N/A' # UTC end time of guide star acquisition, e.g. 'N/A'
# stsci2ips_dict['GDSTARID']= 'N/A' # guide star identifier, e.g. 'N/A'
# stsci2ips_dict['GS_RA']   = 0.0 # guide star right ascension, e.g. 'N/A'
# stsci2ips_dict['GS_DEC']  = 0.0 # guide star declination, e.g. 'N/A'
# stsci2ips_dict['GS_MAG']  = 0.0 # guide star magnitude in FGS detector, e.g. 'N/A'
# stsci2ips_dict['GS_UMAG']  = 0.0 # guide star magnitude uncertainty, e.g. 'N/A'
# stsci2ips_dict['GSCENTX']   = 0.0 # guide star centroid x position in FGS ideal fra
# stsci2ips_dict['GSCENTY']   = 0.0 # guide star centroid y position in FGS ideal fra
# stsci2ips_dict['JITTERMS']   = 0.0 # RMS jitter over the exposure (arcsec)

# Reference file information

# CRDS parameters
stsci2ips_dict['CRDS_VER'] = 'specific_string:7.4.0'  # Version of CRDS file selection software used
stsci2ips_dict['CRDS_CTX'] = 'specific_string:jwst_0570.pmap'  # CRDS context (.pmap) used to select ref files

# in this section of the header the pipeline will store all reference files used through cal_detector1
# should information exist, keywords would look like the following commented section:
"""
 Dark reference file information

R_DARK  = 'crds://jwst_nirspec_dark_0086.fits' / Dark reference file name

        Gain reference file information

R_GAIN  = 'crds://jwst_nirspec_gain_0019.fits' / Gain reference file name

        Linearity reference file information

R_LINEAR= 'crds://jwst_nirspec_linearity_0018.fits' / Linearity reference file n

        Mask reference file information

R_MASK  = 'crds://jwst_nirspec_mask_0010.fits' / Mask reference file name

        Read noise reference file information

R_READNO= 'crds://jwst_nirspec_readnoise_0018.fits' / Read noise reference file

        Reference pixels reference file information

R_REFPIX= 'crds://jwst_nirspec_refpix_0004.fits' / Reference pixels reference fi

        Saturation reference file information

R_SATURA= 'crds://jwst_nirspec_saturation_0020.fits' / Saturation reference file

        Superbias reference file information

R_SUPERB= 'crds://jwst_nirspec_superbias_0113.fits' / Superbias reference file n

"""

# Calibration step information
# in this section, cal_detector1 should have added the following keywords:
# stsci2ips_dict['S_DARK']   = 'COMPLETE' # Dark Subtraction
# stsci2ips_dict['S_DQINIT'] = 'COMPLETE' # Data Quality Initialization
# stsci2ips_dict['S_GANSCL'] = 'SKIPPED ' # Gain Scale Correction
# stsci2ips_dict['S_JUMP']   = 'COMPLETE' # Jump Detection
# stsci2ips_dict['S_LINEAR'] = 'COMPLETE' # Linearity Correction
# stsci2ips_dict['S_RAMP']   = 'COMPLETE' # Ramp Fitting
# stsci2ips_dict['S_REFPIX'] = 'COMPLETE' # Reference Pixel Correction
# stsci2ips_dict['S_SATURA'] = 'COMPLETE' # Saturation Checking
# stsci2ips_dict['S_SUPERB'] = 'COMPLETE' # Superbias Subtraction
# stsci2ips_dict['NEXTEND'] = 3 # Number of standard extensions, e.g. 3
# stsci2ips_dict['RADESYS'] = 'ICRS' # Coordinate reference frame for RA and Dec
# stsci2ips_dict['DPSW_VER'] = '0.8.0' # Data processing software version number
# stsci2ips_dict['WFSVISIT'] = 'NO' # wavefront sensing and control visit indicator
# stsci2ips_dict['EXTARGET'] = 'T' # T if at least one exposure in visit is external
# stsci2ips_dict['TARRUDEC'] = 0.0 # target Dec uncertainty
# stsci2ips_dict['PROPEPOC'] = 2000.0 # proposer specified epoch for RA and Dec
# stsci2ips_dict['NRESET']   = 1 # number of resets between integrations
# stsci2ips_dict['CHRGTIME'] = -1.0 # [seconds] charge accumulation time
# stsci2ips_dict['NXLIGHT']  = '#TODO' # number of light sensitive x values (columns)
# stsci2ips_dict['MSACONFG'] = 'N/A' # MSA configuration file name
# stsci2ips_dict['GSURA']    = 0.0 # guide star right ascension uncertainty
# stsci2ips_dict['GSUDEC']   = 0.0 # guide star declination uncertainty
# stsci2ips_dict['GSUMAG']   = 0.0 # guide star magnitude uncertainty
# stsci2ips_dict['COORDSYS'] = 'N/A' # ephemeris coordinate system
# stsci2ips_dict['PA_V3'] = 'primary_ext:ROLL_REF'  # [deg] position angle of V3-axis of JWST, e.g. 'N/A'
# stsci2ips_dict['RA_V1'] = 0.0  # [deg] RA of telescope V1 axis, e.g. 'N/A'
# stsci2ips_dict['DEC_V1'] = 0.0  # [deg] Dec of telescope V1 axis, e.g. 'N/A'
# stsci2ips_dict['PA_APER'] = -999.0 # [deg] position angle of aperture used, e.g. -999.0
# stsci2ips_dict['VA_SCALE']=  1.0 # velocity aberration scale factor, e.g. 1.0

# WCS parameters in science extension
stsci2ips_dict['wcsinfo'] = {
                            # 'WCSAXES': 3,  # number of World Coordinate System axes, e.g. 3
                            # 'CRPIX1': 1024,  # x-coordinate of the reference pixel, e.g. 1024.0
                            # 'CRPIX2': 1024,  # y-coordinate of the reference pixel, e.g. 128.0
                            # 'CRPIX3': 1024,  # z-coordinate of the reference pixel
                            # 'CRVAL1': 5.3196,  # RA at the reference pixel (degrees), e.g. 5.3196
                            # 'CRVAL2': -72.98605000000001,  # Dec at reference pixel (degrees), e.g. -72.98605000000001
                            # 'CRVAL3': 2.5,  # Wavelength at the reference pixel (microns), e.g. 2.5
                            # 'CTYPE1': 'RA---TAN',  # first axis coordinate type
                            # 'CTYPE2': 'DEC--TAN',  # second axis coordinate type
                            # 'CTYPE3': 'WAVE',  # third axis coordinate type, e.g. WAVE
                            # 'CUNIT1': 'deg',  # units for first axis
                            # 'CUNIT2': 'deg',  # units for seconds axis
                            # 'CUNIT3': 'um',  # units for third axis, e.g. 'micron'
                            # 'CDELT1': 8.69045277777777E-06,  # increment per pixel, axis 1, e.g. 0.0000012
                            # 'CDELT2': 8.73766666666665E-06,  # increment per pixel, axis 2, e.g. 0.0000012
                            # 'CDELT3': 0.00000672,  # increment per pixel, axis 3, e.g. 0.000672
                            # 'PC1_1': -1.0,  # linear transformation matrix element, e.g. 1.0
                            # 'PC1_2': 0.0,  # linear transformation matrix element, e.g. 0.0
                            # 'PC1_3': 0.0,  # linear transformation matrix element, e.g. 0.0
                            # 'PC2_1': 0.0,  # linear transformation matrix element, e.g. 0.0
                            # 'PC2_2': 1.0,  # linear transformation matrix element, e.g. 1.0
                            # 'PC3_1': 1.0,  # linear transformation matrix element, e.g. 1.0
                            # 'PC3_2': 1.0,  # linear transformation matrix element, e.g. 0.0
                            # 'PC3_3': 0.0,  # linear transformation matrix element, e.g. 0.0
                            # 'S_REGION': 'N/A',  # spatial extent of the observation, e.g. 'N/A'
                            # 'WAVSTART': 1.0,  # lower bound of the default wavelength range
                            # 'WAVEND': 2.0,  # upper bound of the default wavelength range
                            # 'SPORDER': 1,  # default spectral order
                            # 'V2_REF': 'primary_ext:V2_REF',  # aperture reference V2 point (arcsec): 100 to -400
                            # 'V3_REF': 'primary_ext:V3_REF',  # aperture reference V3 point (arcsec): -100 to -400
                            # 'VPARITY': 'pysiaf',  # Relative sense of rotation between Ideal xy and V2V3
                            'V3I_YANG': 'pysiaf',  # Angle from V3 axis to Ideal y axis (deg)
                            'RA_REF': 'primary_ext:RA_REF',  # RA at the reference point (deg): 0 < RA < 360
                            'DEC_REF': 'primary_ext:DEC_REF',  # Dec at the reference point (deg): -90 < Dec < +90
                            'ROLL_REF': 'primary_ext:ROLL_REF'  # Roll angle at the reference point (deg), e.g. 5.3196
}


# Functions

def change_keyword2ips_value(ips_keywd_dict, st_pipe_ready_dict, ips_keywd, st_pipe_ready_keywd, st_pipe_ready_file,
                             verbose=False):
    """
    Function changes the keyword value to match IPS value.
    :param ips_keywd_dict: dictionary, IPS header
    :param st_pipe_ready_dict: dictionary, STScI pipeline-ready file header
    :param st_pipe_ready_keywd: string, keyword with STScI spelling
    :param ips_keywd: string, keyword with IPS spelling
    :param st_pipe_ready_file: string, path and name of the fits file that is STScI pipeline-ready
    :param verbose: boolean
    :return: nothing
    """
    for key, val in st_pipe_ready_dict.items():
        if key == st_pipe_ready_keywd:
            if ips_keywd in ips_keywd_dict:
                fits.setval(st_pipe_ready_file, st_pipe_ready_keywd, value=ips_keywd_dict[ips_keywd])
                if verbose:
                    print('Modified keyword: ', key, '   old_value=', val, '   new_value=', ips_keywd_dict[ips_keywd])
            else:
                if verbose:
                    print('IPS keyword ', ips_keywd, ', corresponding to ', st_pipe_ready_keywd,
                          'in STScI pipeline, not found in header of IPS file.')
            break


def find_subarray_size(subsize1, subsize2, detector, grating):
    stpipe_key, ssz1, ssz2, sst1, sst2 = None, None, None, None, None
    for sa in subdict.subarray_dict:
        ssz1 = subdict.subarray_dict[sa]["subsize1"]
        ssz2 = subdict.subarray_dict[sa]["subsize2"]
        if subsize1 == ssz1:
            if subsize2 == ssz2:
                stpipe_key = sa
                sst1 = subdict.subarray_dict[sa]["substrt1"]
                sst2_dict = subdict.subarray_dict[sa]["substrt2"]
                for grat, sst2_tuple in sst2_dict.items():
                    if grat.lower() == grating.lower():
                        if "1" in detector:
                            sst2 = sst2_tuple[0]
                        elif "2" in detector:
                            sst2 = sst2_tuple[1]
                        break
    return stpipe_key, ssz1, ssz2, sst1, sst2


def set_subarray_and_size_keywds(ips_keywd_dict, st_pipe_ready_dict, st_pipe_ready_file, verbose=False):
    """
    Set the subarray keyword to one of the values that the STScI pipeline excepts, with the corresponding size keywords.
    :param ips_keywd_dict: dictionary, IPS header
    :param st_pipe_ready_dict: dictionary, STScI pipeline-ready file header
    :param st_pipe_ready_file: string, path and name of the fits file that is STScI pipeline-ready
    :return: nothing
    """
    # get the information from the pipeline-ready header
    exp_type_keywd_value = st_pipe_ready_dict['EXP_TYPE']
    grating = st_pipe_ready_dict['GRATING']
    detector = st_pipe_ready_dict['DETECTOR']

    # set the subarray keyword to a value the STScI pipeline can process
    if 'MSA' in exp_type_keywd_value or 'IFU' in exp_type_keywd_value:
        if verbose:
            print('Modified keyword: SUBARRAY', '   old_value=', st_pipe_ready_dict['SUBARRAY'],
                  '   new_value= N/A')
        st_pipe_ready_dict['SUBARRAY'] = 'N/A'
    else:
        if ips_keywd_dict['SUBARRAY'] == 'F':
            if verbose:
                print('Modified keyword: SUBARRAY', '   old_value=', st_pipe_ready_dict['SUBARRAY'],
                      '   new_value= FULL')
            st_pipe_ready_dict['SUBARRAY'] = 'FULL'
        else:
            # match the subarray used to the sizes keywords
            data = fits.getdata(st_pipe_ready_file, 'SCI')
            shape = np.shape(data)
            subsize2, subsize1 = shape[-2], shape[-1]
            stpipe_key, ssz1, ssz2, sst1, sst2 = find_subarray_size(subsize1, subsize2, detector, grating)
            if stpipe_key is None:
                # try with inverted subsize 1 and 2
                subsize1, subsize2 = shape[-2], shape[-1]
                stpipe_key, ssz1, ssz2, sst1, sst2 = find_subarray_size(subsize1, subsize2, detector, grating)

            if stpipe_key is not None:
                fits.setval(st_pipe_ready_file, 'SUBARRAY', value=stpipe_key)
                if verbose:
                    print('Modified keyword: SUBARRAY', '   old_value=', st_pipe_ready_dict['SUBARRAY'],
                          '   new_value=', stpipe_key)
                    print('Modified keyword: SUBSTRT1', '   old_value=', st_pipe_ready_dict['SUBSTRT1'],
                          '   new_value=', sst1)
                    print('Modified keyword: SUBSTRT2', '   old_value=', st_pipe_ready_dict['SUBSTRT2'],
                          '   new_value=', sst2)
                    print('Modified keyword: SUBSIZE1', '   old_value=', st_pipe_ready_dict['SUBSIZE1'],
                          '   new_value=', ssz1)
                    print('Modified keyword: SUBSIZE2', '   old_value=', st_pipe_ready_dict['SUBSIZE2'],
                          '   new_value=', ssz2)
                fits.setval(st_pipe_ready_file, 'SUBSTRT1', value=sst1)
                fits.setval(st_pipe_ready_file, 'SUBSTRT2', value=sst2)
                fits.setval(st_pipe_ready_file, 'SUBSIZE1', value=ssz1)
                fits.setval(st_pipe_ready_file, 'SUBSIZE2', value=ssz2)
            else:
                print('WARNING! Unable to determine correct value for the following keywords: '
                      'SUBARRAY, SUBSTRT1. SUBSTRT2, SUBSIZE1, and SUBSIZE2')


def get_pysiaf_aperture(st_pipe_ready_dict, verbose):
    subarray = st_pipe_ready_dict['SUBARRAY']
    detector = st_pipe_ready_dict['DETECTOR']
    if 'ifu' in st_pipe_ready_dict['EXP_TYPE'].lower():
        aperture_name = detector + '_FULL_IFU'
        return aperture_name
    if 'full' in subarray.lower():
        aperture_name = detector + '_FULL'
    elif '200' in subarray or '400' in subarray:
        aperture_name = 'NRS_' + subarray + '_SLIT'
    else:
        aperture_name = 'NRS_S1600A1_SLIT'
    if verbose:
        print('PySIAF aperture name: ', aperture_name)
    return aperture_name


def match_IPS_keywords(stsci_pipe_ready_file, ips_file, additional_args_dict=None, verbose=False):
    """
    This function performs the change of keyword values for the STScI pipeline-ready file to match the values in the
    IPS file.
    :param stsci_pipe_ready_file: string, path and name of the STScI pipeline-ready file
    :param ips_file: string, path and name of the IPS file
    :param additional_args_dict: dictionary, keywords and corresponding values indicated in the command line
    :param verbose: boolean
    :return: nothing; the input fits file with the modified keyword values
    """
    # get the headers from the IPS file
    primary_ext_ips_keywd_dict = fits.getheader(ips_file, 0)
    header_ext_ips_keywd_dict = fits.getheader(ips_file, extname='header')

    # get the header from the STScI pipeline-ready file
    st_pipe_ready_dict = fits.getheader(stsci_pipe_ready_file, 0)

    # iterate over the map of STScI to IPS keywords dictionary (defined in this script)
    for key2modify, val2modify in stsci2ips_dict.items():
        # check if this is for the science extension in the ST ready file
        if key2modify == 'wcsinfo':
            # all these keywords do not need to be calculated
            for key, val in val2modify.items():
                if 'pysiaf' in val:
                    # get this value from pysiaf
                    NIRSpec_SIAF = pysiaf.Siaf('NIRSpec')
                    aperture_name = get_pysiaf_aperture(st_pipe_ready_dict, verbose=verbose)
                    refpoint = NIRSpec_SIAF[aperture_name].reference_point('tel')
                    V2_REF, V3_REF = refpoint[0], refpoint[1]
                    V3IdlYAngle = NIRSpec_SIAF[aperture_name].V3IdlYAngle
                    VIdlParity = NIRSpec_SIAF[aperture_name].VIdlParity
                    fits.setval(stsci_pipe_ready_file, 'V3I_YANG', value=V3IdlYAngle, extname='SCI')
                    fits.setval(stsci_pipe_ready_file, 'VPARITY', value=VIdlParity, extname='SCI')
                    # check if these are in the header
                    if 'V2_REF' in primary_ext_ips_keywd_dict:
                        V2_REF = primary_ext_ips_keywd_dict['V2_REF']
                    if 'V3_REF' in primary_ext_ips_keywd_dict:
                        V3_REF = primary_ext_ips_keywd_dict['V3_REF']
                    fits.setval(stsci_pipe_ready_file, 'V2_REF', value=V2_REF, extname='SCI')
                    fits.setval(stsci_pipe_ready_file, 'V3_REF', value=V3_REF, extname='SCI')
                    if verbose:
                        print('Modified keyword in science extension: V2_REF', '   new_value=', V2_REF)
                        print('Modified keyword in science extension: V3_REF', '   new_value=', V3_REF)
                        print('Modified keyword in science extension: V3I_YANG', '   new_value=', V3IdlYAngle)
                        print('Modified keyword in science extension: VPARITY', '   new_value=', VIdlParity)
                else:
                    ips_key = val.split(':')[-1]
                    fits.setval(stsci_pipe_ready_file, key, value=primary_ext_ips_keywd_dict[ips_key], extname='SCI')
                    if verbose:
                        print('Modified keyword in science extension: ', key,
                              '   new_value=', primary_ext_ips_keywd_dict[ips_key])
        if val2modify == 'N/A':
            # for a non applicable (not too important) keyword for the simulations, set value to N/A
            st_pipe_ready_dict[key2modify] = 'N/A'

        if 'primary_ext' in val2modify:
            # look for the same keyword in IPS file in the primary extension header
            ips_key = val2modify.split(':')[-1]
            if 'SUBARRAY' in key2modify:
                set_subarray_and_size_keywds(primary_ext_ips_keywd_dict, st_pipe_ready_dict, stsci_pipe_ready_file,
                                             verbose=verbose)
            else:
                change_keyword2ips_value(primary_ext_ips_keywd_dict, st_pipe_ready_dict, ips_key, key2modify,
                                         stsci_pipe_ready_file, verbose=verbose)

        if 'header_ext' in val2modify:
            # look for the same keyword in IPS file in the header extension
            ips_key = val2modify.split(':')[-1]
            change_keyword2ips_value(header_ext_ips_keywd_dict, st_pipe_ready_dict, ips_key, key2modify,
                                     stsci_pipe_ready_file, verbose=verbose)

        if 'set_to_given_string' in val2modify:
            # change the keyword value to that given in the command line - this is optional
            if additional_args_dict is None:
                continue
            else:
                if key2modify in additional_args_dict:
                    if verbose:
                        print('Modified keyword: ', key2modify, '   old_value=', st_pipe_ready_dict[key2modify],
                              '   new_value=', additional_args_dict[key2modify])
                    fits.setval(stsci_pipe_ready_file, key2modify, value=additional_args_dict[key2modify])
                else:
                    if verbose:
                        print('Value for keyword=', key2modify, ' not provided with line command.')
                    continue

        if 'calculation' in val2modify:
            if verbose:
                print('Value for keyword ', key2modify, ' will be calculated...')

            if key2modify == 'EXPSTART' or key2modify == 'EXPEND' or key2modify == 'EXPMID':
                continue

            if key2modify == 'DURATION':
                # calculate also EXPSTART, EXPEND and EXPMID

                # this the calculation follows the JWST keyword dictionary calculation:
                # https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/
                # duration = TFRAME*((NGROUPS*NFRAMES+(NGROUPS-1)*GROUPGAP+DROPFRAMES1*NINTS)) where DROPFRAMES1 is
                # a lookup in the PRD DataModes table.
                # However, NIRSpec only drops frames in TA, hence the calculation simplifies to:
                # duration = TFRAME*((NGROUPS*NFRAMES))
                new_val = st_pipe_ready_dict['TFRAME'] * st_pipe_ready_dict['NGROUPS'] * st_pipe_ready_dict['NFRAMES']
                # change value in dictionary for use within the script
                st_pipe_ready_dict['DURATION'] = new_val

                # calculate EXPSTART
                # this the calculation follows the JWST keyword dictionary calculation:
                # https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/
                # expstart = input('DATE-OBS') + 'T' + input('TIME-OBS') / UTC exposure start time (MJD)
                dateobs_string = st_pipe_ready_dict['DATE-OBS'].replace('T', ' ')
                timeobs_string = st_pipe_ready_dict['TIME-OBS']
                expstart_string = dateobs_string + 'T' + timeobs_string
                expstart = Time(expstart_string, format='fits', scale='utc')
                expstart_mjd = expstart.mjd
                # change value in dictionary for use within the script
                st_pipe_ready_dict['EXPSTART'] = expstart_mjd

                # calculate EXPEND
                duration = st_pipe_ready_dict['DURATION']
                expend = expstart + timedelta(seconds=duration)
                expend_mjd = expend.mjd
                # change value in dictionary for use within the script
                st_pipe_ready_dict['EXPEND'] = expend_mjd

                # calculate EXPMID
                expmid = expstart + timedelta(seconds=duration/2)
                expmid_mjd = expmid.mjd
                # change value in dictionary for use within the script
                st_pipe_ready_dict['EXPMID'] = expmid_mjd

            if key2modify == 'TSAMPLE':
                # this the calculation follows the JWST keyword dictionary calculation:
                # https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/
                # tsample = readout pattern lookup
                # table taken from:
                # https://jwst-docs.stsci.edu/near-infrared-spectrograph/nirspec-instrumentation/nirspec-detectors/nirspec-detector-readout-modes-and-patterns
                readout_patterns = {'NRSRAPID': 10.737,  # frames=1
                                    'NRSRAPIDD1': 21.474,  # frames=1
                                    'NRSRAPIDD2': 32.210,  # frames=1
                                    'NRSRAPIDD6': 75.159,  # frames=1
                                    'NRS': 42.947,  # frames=4
                                    'NRSIRS2RAPID': 14.589,  # frames=1
                                    'NRSIRS2': 72.944}  # frames=5
                new_val = readout_patterns[st_pipe_ready_dict['READPATT']]
                # change value in dictionary for use within the script
                st_pipe_ready_dict['TSAMPLE'] = new_val

            if key2modify == 'TGROUP':
                # this the calculation follows the JWST keyword dictionary calculation:
                # https://mast.stsci.edu/portal/Mashup/Clients/jwkeywords/
                # tgroup = (GROUPGAP+NFRAMES)*TFRAME
                # However, NIRSpec only drops frames in TA, hence the calculation simplifies to:
                # tgroup = NFRAMES*TFRAME
                # reference for GROUPGAP: http://www.stsci.edu/~tumlinso/nirspec_ocd_v6_DRAFT.pdf
                new_val = st_pipe_ready_dict['NFRAMES'] * st_pipe_ready_dict['TFRAME']
                # change value in dictionary for use within the script
                st_pipe_ready_dict['TGROUP'] = new_val

            if verbose:
                print('Modified keyword: ', key2modify, '   old_value=', st_pipe_ready_dict[key2modify],
                      '   new_value=', new_val)
                print('    * WARNING: This calculation needs to be verified ')
            fits.setval(stsci_pipe_ready_file, key2modify, value=new_val)

        elif 'specific_string' in val2modify:
            # now set all the other keyword whose string value will not change from simulation from simulation, this
            # is the case for the 'specific_string' in the map of STScI to IPS keywords dictionary
            new_val = val2modify.split(':')[-1]
            if verbose:
                print('Modified keyword: ', key2modify, '   old_value=', st_pipe_ready_dict[key2modify],
                      '   new_value=', new_val)
            fits.setval(stsci_pipe_ready_file, key2modify, value=val2modify)


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("stsci_pipe_ready_file",
                        action='store',
                        default=None,
                        help='Path and name of fits file that has the STScI pipeline-ready header, '
                             'i.e. blah_modified.fits')
    parser.add_argument("ips_file",
                        action='store',
                        default=None,
                        help='Path and name of the simulations crm file, i.e. blah.cts.fits')
    parser.add_argument("-p",
                        dest="proposal_title",
                        action='store',
                        default=None,
                        help='Add the proposal title to the header keyword, i.e. -p=some_title')
    parser.add_argument("-t",
                        dest="target_name",
                        action='store',
                        default=None,
                        help='Add the target name to the header keyword, i.e. -t=some_target')
    parser.add_argument("-n",
                        dest="new_file",
                        action='store_true',
                        default=True,
                        help='Use -n if wanting to create a new file with updated header. Default is to '
                             'update header without creating a new file')
    parser.add_argument("-v",
                        dest="verbose",
                        action='store_true',
                        default=True,
                        help='Use -v to print keywords values and script messages.')
    args = parser.parse_args()

    # Set the variables
    stsci_pipe_ready_file = args.stsci_pipe_ready_file
    ips_file = args.ips_file
    proposal_title = args.proposal_title
    target_name = args.target_name
    new_file = args.new_file
    verbose = args.verbose

    # create dictionary of command-line arguments
    additional_args_dict = {'TITLE': proposal_title,
                            'TARGNAME': target_name,
                            'new_file': new_file
                            }

    # modify the keyword values to match IPS information
    match_IPS_keywords(stsci_pipe_ready_file, ips_file, additional_args_dict=additional_args_dict,
                       verbose=verbose)

    print('\n * Script  level2b_hdr_keywd_dict_map2sim.py  finished * \n')


if __name__ == '__main__':
    sys.exit(main())
