import argparse
import pysiaf
import os
import re
import subprocess
import sys
import numpy as np
from datetime import datetime
from astropy.io import fits
from collections import OrderedDict
from glob import glob

# import the header keyword dictionaries
from .dict_info import hdr_keywd_dict as hkwd
from .dict_info import hdr_keywd_dict_sample as shkvd
# import subarray dictionary
from .dict_info import subarray_dict as subdict

'''
This script checks that the fits files to be used as input for the pipeline build 7.1, have the expected keywords in
the primary and science headers.

Example usage:
    The code works from the terminal.
    To create a NEW FS fits file with the updated header type:
        $ nptt_hdr_keywd_check.py blah.fits

    To simply update the header of the existing fits file type:
        $ nptt_hdr_keywd_check.py blah.fits -u -m=FS

where the -m flag is the mode used, i.e. FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf.

The code is NOT case-sensitive.

If a mode is not provided, the code will look for a mode_used variable in the pytests configuration file, and it
will crash if this config file does not exist.

'''

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.4"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Apr 2019 - Version 1.1: added dictionary to choose right GWA_XTIL keyword value according to GRATING
# May 2019 - Version 1.2: added logic for dark processing
# Jul 2020 - Version 1.3: changed default value of SUBARRAY according to CRDS rules
# Feb 2021 - Version 1.4: implemented adding the MSA metafile name to the header


# Paths
wit4_path = os.environ.get('WIT4_PATH')
if wit4_path is None:
    print("(hdr_keywd_check): The environment variable WIT4_PATH is not defined. To set it, follow the "
          "instructions at: \n"
          "                   https://github.com/spacetelescope/nirspec_pipe_testing_tool")
    exit()
nirspec_cdp3 = "nirspec/CDP3/03_Instrument_model/3.1_Files/NIRS_FM2_05_CV3_FIT1/Description"
path_to_tilt_files = os.path.join(wit4_path, nirspec_cdp3)


# General functions

def get_modeused_PTT_cfg_file(detector):
    # get script directory and config name
    ptt_cfg_file = "PTT_config.cfg"
    get_mode = False
    mode_used, msa_metafile = 'N/A', 'N/A'
    if os.path.isfile(ptt_cfg_file):
        get_mode = True
    else:
        ptt_cfg_file = ptt_cfg_file.replace(".cfg", "_"+detector+".cfg")
    if os.path.isfile(ptt_cfg_file):
        get_mode = True
    else:
        print("Unable to obtain instrument mode from the PTT configuration file.")
    if get_mode:
        with open(ptt_cfg_file, "r") as cfg:
            for i, line in enumerate(cfg.readlines()):
                if "mode_used" in line:
                    mode_used = line.split()[2]
                if "msa_conf_name" in line:
                    msa_metafile = line.split()[2]
                    msa_metafile = os.path.basename(msa_metafile)
    return mode_used, msa_metafile


def read_hdrfits(fits_file_name):
    """
    This function reads the header fits file and returns a dictionary of the keywords with
    corresponding values. Keywords will be stored in the order they are read.
    Args:
        fits_file_name: full path with name of the fits file

    Returns:
        A dictionary of keywords with corresponding values
    """
    #  Read the fits file
    hdulist = fits.open(fits_file_name)
    # print on screen what extensions are in the file
    # print ('\n FILE INFORMATION: \n')
    hdulist.info()
    # get and print header
    hdr = hdulist[0].header
    # close the fits file
    hdulist.close()
    return hdr  # this is really what we should be using...


def create_addedkeywds_file(fits_file):
    """
    This function create text file to log added keywords.
    Args:
        fits_file: string, name of the fits file to be keyword checked

    Returns:
        addedkeywds_file_name: string, the file name where all added keywords were saved
    """
    addedkeywds_file_name = fits_file.replace(".fits", "_addedkeywds.txt")
    print('Name of text file containing all keywords added:  ', addedkeywds_file_name)
    tf = open(addedkeywds_file_name, 'w')
    tf.write('### The following keywords were added or have the wrong format: \n')
    tf.write('# {:<12} {:<10} {:<25} \n'.format('Keyword', 'Extension', 'Comments'))
    tf.close()
    return addedkeywds_file_name


def check_addedkeywds_file(addedkeywds_file_name):
    """
    If no new keywords were added, the text file is set to empty, and this function will erase it.
    Args:
        addedkeywds_file_name: string, name of text file to check

    Returns:
        nothing

    """
    with open(addedkeywds_file_name, 'r') as wf:
        lines = wf.readlines()
        if len(lines) == 1:
            os.system('rm ' + addedkeywds_file_name)


# Functions to check specific keyword values

def check_value_type(key, val, hkwd_val, ext='primary'):
    """
    Check if the keyword has the right type of value.
    Args:
        key: string, keyword
        val: string, integer, or float, value of the keyword given
        hkwd_val: list of the allowed keyword values
        ext: string, extension number (default value is 'primary')

    Returns:
        warning: string or None, if string this is a sentence saying that either the keyword was empty or that its
                 value does not correspond to one of the expected values.
    """

    # Check if type of value matches expected
    valtype = type(val)

    # Store type that value should have, checking for option lists
    dict_type = type(hkwd_val[0]) if len(hkwd_val) > 1 else hkwd_val[0]

    # If we don't want a string, check for numerics
    if dict_type != str:
        # regex patterns
        float_pattern = re.compile(r"\b\s*[+-]?(?=\d*[.eE])([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?(?=\s|$)")
        int_pattern = re.compile(r"\b\s*[-+]?\d+(?=\s|$)")

        if val is None or valtype == '':
            warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'This keyword has an empty value')
            print(warning)
            return None, warning, dict_type
        if valtype == str:
            if not float_pattern.fullmatch(val) is None:
                val = float(val)
                valtype = float
            elif not int_pattern.fullmatch(val) is None:
                val = int(val)
                valtype = int

    if valtype == bool:
        return val, None, valtype

    if (valtype in hkwd_val) or (val in hkwd_val) or (valtype == type(dict_type)):
        print('{:<15} {:<9} {:<28} {:<16}'.format(key, ext, 'Allowed value type', val))
        warning = None
    else:
        warning = '{:<15} {:<9} {:<28}'.format(key, ext, 'Incorrect value type. Expected e.g. ' + repr(hkwd_val[0]) +
                                               ', got: ' + repr(val))
        print(warning)

    return val, warning, dict_type


def check3numbers(key, val, ext='primary'):
    """
    Check if this keyword value has a format like: 0.1.1
    Args:
        key: keyword
        val: value of that keyword
        ext: string, extension number (default value is 'primary')

    Returns:
        A warning (if the format is not what is expected) or nothing if it is.

    """
    warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect value format. Expected: 0.0.0, got: ' + str(val))
    r = re.compile('\d.\d.\d')  # check that the string has a structure like 0.1.1
    if r.match(val) is None:
        print(warning)
        return warning
    else:
        print('{:<15} {:<9} {:<25}'.format(key, ext, 'Matches expected format'))


def check_len(key, val, val_len=2, ext='primary'):
    """
    Check if the length of the keyword value has a a given length, default value for the length check is 2.
    Args:
        key: keyword
        val: keyword value
        val_len: length to be checked against
        ext: string, extension number (default value is 'primary')

    Returns:
        A warning (if the length was not what was expected) or nothing.

    """
    if isinstance(val, int):
        string_length = len(str(val))
    else:
        string_length = len(val)
    warning = '{:<15} {:<9} {:<25}'.format(key, ext,
                                           'Incorrect length of value. Expected: ' + repr(val_len) + ', got ' + repr(
                                               string_length))
    if not string_length == val_len:
        print(warning)
        return warning
    # else:
    #    print('{:<15} {:<9} {:<25}'.format(key, ext, 'Correct format'))


def check_datetimeformat(key, val, check_time, check_date, check_datetime, ext='primary'):
    """
    Check if the date and/or time has the expected format.
    Args:
        key: keyword
        val: keyword value
        check_time: boolean, if true check against this format hr:min:sec
        check_date: boolean, if true check against this format year:month:day
        check_datetime: boolean, if true check against this format year-month-dayThr:min:sec
        ext: string, extension number (default value is 'primary')

    Returns:
        A warning that the format is not correct, or nothing.
    """
    warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect value format')
    if '.' in val:
        vlist = val.split(':')
        v = float(vlist[-1])
        v = str(int(np.round(v, decimals=0)))
        val = vlist[0] + ':' + vlist[1] + ':' + v
    if check_time:
        val = datetime.strptime(val, '%H:%M:%S')
    if check_date:
        val = datetime.strptime(val, '%Y-%m-%d')
    if check_datetime:
        val = datetime.strptime(val, '%Y-%m-%dT%H:%M:%S')
    if not isinstance(val, datetime):
        print(warning)
        return warning
    # else:
    #    print('{:<15} {:<9} {:<28} {:<16}'.format(key, ext, 'Correct value format', val))


def get_gwa_Xtil_val(grating, path_to_tilt_files):
    """
    This function gets the right GWA_XTIL value according to the grating given. The reference file has a different
    reference frame, where X and Y are inverted with respect to the pipeline.
    Args:
        grating: string, value from GRATING keyword

    Returns:
        gwaxtil: float, corresponding GWA_XTIL value to the grating
    """
    dispersion_files_list = glob(os.path.join(path_to_tilt_files, "disperser_*_TiltY.gtp"))
    gwa_xtil_found = False
    for dispersion_file in dispersion_files_list:
        if grating in dispersion_file:
            print("Using this file for setting the GWA_XTIL keyword: \n", dispersion_file)
            with open(dispersion_file, "r") as df:
                for line in df.readlines():
                    if "*Zeroreadings 1" in line:
                        gwa_xtil_found = True
                        continue
                    if gwa_xtil_found:
                        line = line.replace("\n", "")
                        gwa_xtil = float(line)
                        break
    return gwa_xtil


def get_gwa_Ytil_val(grating, path_to_tilt_files):
    """
    This function gets the right GWA_YTIL value according to the grating given. The reference file has a different
    reference frame, where X and Y are inverted with respect to the pipeline.
    Args:
        grating: string, value from GRATING keyword

    Returns:
        gwaxtil: float, corresponding GWA_YTIL value to the grating
    """
    dispersion_files_list = glob(os.path.join(path_to_tilt_files, "disperser_*_TiltX.gtp"))
    gwa_ytil_found = False
    for dispersion_file in dispersion_files_list:
        if grating in dispersion_file:
            print("Using this file for setting the GWA_YTIL keyword: \n", dispersion_file)
            with open(dispersion_file, "r") as df:
                for line in df.readlines():
                    if "*Zeroreadings 1" in line:
                        gwa_ytil_found = True
                        continue
                    if gwa_ytil_found:
                        line = line.replace("\n", "")
                        gwa_ytil = float(line)
                        break
    return gwa_ytil


def get_pysiaf_aperture(mode, subarray, detector):
    if 'ifu' in mode.lower():
        aperture_name = detector + '_FULL_IFU'
        return aperture_name
    if 'full' in subarray.lower():
        aperture_name = detector + '_FULL'
    elif '200' in subarray or '400' in subarray:
        aperture_name = 'NRS_' + subarray + '_SLIT'
    else:
        aperture_name = 'NRS_S1600A1_SLIT'
    print('PySIAF aperture name: ', aperture_name)
    return aperture_name


def get_pipe_subarray_name(subarray):
    print('changing subarray from value: ', subarray)
    pipe_subarr_val = 'NOT_FOUND'
    if 'FULL' in subarray:
        pipe_subarr_val = 'FULL'
    elif '200A1' in subarray:
        pipe_subarr_val = 'S200A1'
    elif '200A2' in subarray:
        pipe_subarr_val = 'S200A2'
    elif '200B1' in subarray:
        pipe_subarr_val = 'S200B1'
    elif '400A1' in subarray:
        pipe_subarr_val = 'S400A1'
    elif '1600' in subarray:
        pipe_subarr_val = 'S1600A1'
    print('  to value: ', pipe_subarr_val)
    return pipe_subarr_val


# keyword and format check

def check_keywds(file_keywd_dict, warnings_file_name, warnings_list, missing_keywds, mode_used):
    """
    This function will check keywords against those in hdr_keywod_dict.py
    Args:
        file_keywd_dict: dictionary of the original keyword header
        warnings_file_name: string, name of the file to reccord added keywords
        warnings_list: list of the warnings to be written into the file
        missing_keywds: list of the keywords not in the original the header
        mode_used: str or None, observation mode used FS, MOS, or IFU (if None then a configuration file
                    is expected to exist and contain a variable named mode_used)

    Returns:
        specific_keys_dict: dictionary with specific keys and values that need to be changed
    """
    # get the name and path of the input fits file
    ff = warnings_file_name.replace('_addedkeywds.txt', '.fits')

    # get the detector and grating from the input file
    detector = fits.getval(ff, "DETECTOR", 0)
    grating = fits.getval(ff, "GRATING", 0)

    # get the instrument mode used either from the input of the script or from the configuration file
    mode, msa_metafile = get_modeused_PTT_cfg_file(detector)
    if mode_used is None:
        mode_used = mode

    # initialize the dictionary to hold keyword values that will be changed and written to the file
    specific_keys_dict = {}

    # loop through the keywords and values of the PTT dictionary and add keywords that are not in input file
    for hkwd_key, hkwd_val in hkwd.keywd_dict.items():
        # start by making the warning for each keyword None and assigning key and val
        key = hkwd_key

        # Check if keyword is in the file, expect for wcs info (these go in the science extension)
        if key != 'wcsinfo':
            ext = 'primary'
            if key not in file_keywd_dict:
                # make sure the MSA metafile is pointing to the right place
                if key == 'MSAMETFL':
                    if ('mos' in mode_used.lower()) or ("msa" in mode_used.lower()):
                        val = msa_metafile
                        specific_keys_dict[key] = val
                        print('     Setting value of ', key, ' to ', val)
                missing_keywds.append(key)
                warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'New keyword added to header')
                warnings_list.append(warning)
            else:
                val = file_keywd_dict[hkwd_key]

                # Check simple standard keyword values
                if val in hkwd_val:
                    print('{:<15} {:<9} {:<28} {:<16}'.format(key, ext, 'Has correct format', val))
                    warning = None
                else:
                    new_val, warning, dict_type = check_value_type(key, val, hkwd_val)
                    if warning is not None and "Incorrect value type" in warning:
                        if dict_type == int:
                            new_val = int(float(new_val))
                        else:
                            new_val = dict_type(new_val)
                        print('     Setting value of ', key, ' to type ', dict_type, ' and value ', val)
                        warning = None
                    if new_val != val:
                        specific_keys_dict[key] = new_val
                        missing_keywds.append(key)

                # Check for specific keywords
                if key == 'DPSW_VER':
                    warning = check3numbers(key, val)
                elif (key == 'VISITGRP') or (key == 'ACT_ID'):
                    warning = check_len(key, val, val_len=2)
                elif (key == 'OBSERVTN') or (key == 'VISIT'):
                    warning = check_len(key, val, val_len=3)
                elif key == 'EXPOSURE':
                    warning = check_len(key, val, val_len=5)
                elif (key == 'DATE') or (key == 'VSTSTART'):
                    warning = check_datetimeformat(key, val, check_date=False, check_datetime=True,
                                                   check_time=False)
                elif key == 'DATE-OBS':
                    warning = check_datetimeformat(key, val, check_date=True, check_datetime=False,
                                                   check_time=False)
                    """ THIS ISSUE SEEMS TO HAVE BEEN SOLVED FOR NOW, BUT REMAINS HERE IN CASE IT RESURFACES
                    # check that for IFU data date is later than 2016-5-11, this is so that the latest IFU_post
                    # file is used, otherwise assign_wcs will crash 
                    if mode_used.lower() == "ifu":
                        yr = val.split("-")[0]
                        mo = val.split("-")[1]
                        if int(yr) < 2016:
                            val = "2016-05-12"
                            specific_keys_dict[key] = val
                            missing_keywds.append(key)
                        elif int(mo) < 5:
                            val = "2016-05-12"
                            specific_keys_dict[key] = val
                            missing_keywds.append(key)
                    """
                elif key == 'TIME-OBS':
                    warning = check_datetimeformat(key, val, check_date=False, check_datetime=False,
                                                   check_time=True)

                # specific check for VISITYPE, set to GENERIC
                if key == 'VISITYPE':
                    if val not in hkwd_val:
                        # for now always set this keyword to generic
                        print("Replacing ", key, fits.getval(ff, "VISITYPE", 0), "for GENERIC")
                        # fits.setval(ff, key, 0, value='GENERIC')
                        specific_keys_dict[key] = 'GENERIC'
                        missing_keywds.append(key)

                # specific check for SUBARRAY
                if key == 'SUBARRAY':
                    if 'ifu' in mode_used.lower() or "mos" in mode_used.lower():
                        if val not in hkwd_val:
                            print("Replacing ", key, fits.getval(ff, "SUBARRAY", 0), "for N/A")
                            specific_keys_dict[key] = 'N/A'
                            missing_keywds.append(key)
                    else:  # set SUBARRAY for anything else other than IFU or MOS
                        # set the subarray according to size
                        substrt1 = fits.getval(ff, "substrt1", 0)
                        substrt2 = fits.getval(ff, "substrt2", 0)
                        subsize1 = fits.getval(ff, "subsize1", 0)
                        subsize2 = fits.getval(ff, "subsize2", 0)
                        for subarrd_key, subarrd_vals_dir in subdict.subarray_dict.items():
                            sst1 = subarrd_vals_dir["substrt1"]
                            sst2_dict = subarrd_vals_dir["substrt2"]
                            ssz1 = subarrd_vals_dir["subsize1"]
                            ssz2 = subarrd_vals_dir["subsize2"]
                            # print ("sst1=", sst1, "  sst2_list=", sst2_list, "  ssz1=", ssz1, "  ssz2=", ssz2)
                            if substrt1 == sst1 and subsize1 == ssz1 and subsize2 == ssz2:
                                for grat, sst2_tuple in sst2_dict.items():
                                    if grat.lower() == grating.lower():
                                        if 'FULL' in subarrd_key:
                                            subarrd_key = 'FULL'
                                        elif '200A1' in subarrd_key:
                                            subarrd_key = 'S200A1'
                                        elif '200A2' in subarrd_key:
                                            subarrd_key = 'S200A2'
                                        elif '200B1' in subarrd_key:
                                            subarrd_key = 'S200B1'
                                        elif '400A1' in subarrd_key:
                                            subarrd_key = 'S400A1'
                                        elif '1600' in subarrd_key:
                                            if mode_used.lower() == "fs":
                                                subarrd_key = 'SUB2048'
                                            elif mode_used.lower() == "bots":
                                                # determine which 1600 subarray is it
                                                data = fits.getdata(ff, 1)
                                                subsize1, subsize2 = np.shape(data)
                                                for sa in subdict.subarray_dict:
                                                    ssz1 = subdict.subarray_dict[sa]["subsize1"]
                                                    ssz2 = subdict.subarray_dict[sa]["subsize2"]
                                                    if subsize1 == ssz1:
                                                        if subsize2 == ssz2:
                                                            subarrd_key = sa
                                                            break
                                        specific_keys_dict[key] = subarrd_key
                                        print("changing subarray keyword to ", subarrd_key)
                                        missing_keywds.append(key)
                                        # and make sure to change the primary slit keyword accordingly
                                        if mode_used.lower() == "fs":
                                            FXD_SLIT = subarrd_key
                                            if '2048' or '1024' or '512' or '32' in subarrd_key:
                                                FXD_SLIT = 'S1600A1'
                                            specific_keys_dict['FXD_SLIT'] = FXD_SLIT
                                            print("changing primary slit keyword to FXD_SLIT=", subarrd_key)
                                            missing_keywds.append('FXD_SLIT')
                                        # this part is simply to check that the subarray values are correct
                                        # but no values will be changed in the input file
                                        if "1" in detector:
                                            sst2 = sst2_tuple[0]
                                        elif "2" in detector:
                                            sst2 = sst2_tuple[1]
                                        print("Subarray values in input file: \n", )
                                        print("substrt1=", substrt1, " substrt2=", substrt2, " subsize1=", subsize1,
                                              " subsize2=", subsize2)
                                        print("Subarray values in PTT dictionary: \n", )
                                        print("substrt1=", sst1, " substrt2=", sst2, " subsize1=", ssz1, " subsize2=",
                                              ssz2)

                # check for right value for EXP_TYPE, default will be to add the sample value: NRS_MSASPEC
                if key == 'EXP_TYPE':
                    # make sure there are no white spaces
                    if " " in mode_used:
                        mode_used = mode_used.replace(" ", "")
                    print('   * MODE_USED  = ', mode_used)
                    if mode_used.lower() == "fs":
                        val = 'NRS_FIXEDSLIT'
                        # set the lamp mode to the correct value
                        specific_keys_dict['LAMPMODE'] = 'FIXEDSLIT'
                        missing_keywds.append('LAMPMODE')
                    if mode_used.lower() == "ifu":
                        val = 'NRS_IFU'
                        # set the lamp mode to the correct value
                        specific_keys_dict['LAMPMODE'] = 'IFU'
                        missing_keywds.append('LAMPMODE')
                    if mode_used.lower() == "mos" or mode_used.lower() == "msa":
                        val = 'NRS_MSASPEC'
                        # set the lamp mode to the correct value
                        specific_keys_dict['LAMPMODE'] = 'MSASPEC'
                        missing_keywds.append('LAMPMODE')
                    if mode_used.lower() == "bots":
                        val = 'NRS_BRIGHTOBJ'
                        # set the lamp mode to the correct value
                        specific_keys_dict['LAMPMODE'] = 'BRIGHTOBJ'
                        missing_keywds.append('LAMPMODE')
                    if mode_used.lower() == "dark":
                        val = 'NRS_DARK'
                    if mode_used.lower() == "image":
                        val = 'NRS_IMAGE'
                    if mode_used.lower() == "confirm":
                        val = 'NRS_CONFIRM'
                    if mode_used.lower() == "taconfirm":
                        val = 'NRS_TACONFIRM'
                    if mode_used.lower() == "wata":
                        val = 'NRS_WATA'
                    if mode_used.lower() == "msata":
                        val = 'NRS_MSATA'
                    if mode_used.lower() == "focus":
                        val = 'NRS_FOCUS'
                    if mode_used.lower() == "mimf":
                        val = 'NRS_MIMF'
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)
                    print('     Setting value of ', key, ' to ', val)

                if key == 'TSOVISIT':
                    if mode_used.lower() == 'bots':
                        val = True
                    else:
                        val = False
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)

                # make sure the MSASTATE keyword is set correctly
                if key == 'MSASTATE':
                    if (mode_used.lower() == 'fs') or (mode_used.lower() == 'ifu'):
                        val = 'PRIMARYPARK_ALLCLOSED'
                        specific_keys_dict[key] = val
                        missing_keywds.append(key)

                # make sure the MSA metafile is pointing to the right place
                if key == 'MSAMETFL':
                    if (mode_used.lower() == 'mos') or ("msa" in mode_used.lower()):
                        val = msa_metafile
                        specific_keys_dict[key] = val
                        missing_keywds.append(key)
                        print('     Setting value of ', key, ' to ', val)

                # only modify these keywords if present
                if key == 'GWA_XP_V':
                    val = float(val)
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)
                if key == 'GWA_YP_V':
                    val = float(val)
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)
                if key == 'GWA_PXAV':
                    val = float(val)
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)
                if key == 'GWA_PYAV':
                    val = float(val)
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)
                if key == 'PHOTMJSR':
                    val = float(val)
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)

            if warning is not None:
                missing_keywds.append(key)
                warnings_list.append(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning + '\n')
        else:
            for subkey, _ in hkwd_val.items():
                # add the WCS keywords to science extension
                missing_keywds.append(key)
                # now add the keyword to in the list to be added into the science extension
                warning = '{:<15} {:<9} {:<25}'.format(subkey, 'sci', 'New keyword added to header')
                warnings_list.append(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning + '\n')

    print("keywords to be modified: ", list(OrderedDict.fromkeys(missing_keywds)))
    return specific_keys_dict


def add_keywds(fits_file, only_update, missing_keywds, specific_keys_dict):
    """
    This function adds the missing keywords from the hdr_keywords_dictionary.py (hkwd) file and gives
    the fake values taken from the dictionary sample_hdr_keywd_vals_dict.py (shkvd).
    Args:
        only_update: If false a copy of the original fits file will be created with the
                     updated header.
        missing_keywds: list, missing keywords will be appended here
        specific_keys_dict: dictionary with specific keys and values that need to be changed
    """
    missing_keywds = list(OrderedDict.fromkeys(missing_keywds))
    # print ("specific_keys_dict = ", specific_keys_dict)
    # create name for updated fits file
    updated_fitsfile = fits_file
    if not only_update:
        updated_fitsfile = fits_file.replace('.fits', '_updatedHDR.fits')
        subprocess.run(["cp", fits_file, updated_fitsfile])
    # add missimg keywords
    wcs_keywds_from_main_hdr = {}
    print('Saving keyword values in file: ', updated_fitsfile)
    ext = 0
    for i, key in enumerate(missing_keywds):
        if key in specific_keys_dict:
            # print ("found it in the dict: ", key, specific_keys_dict[key])
            if specific_keys_dict[key] == 'remove':
                # keyword to be deleted
                try:
                    fits.delval(updated_fitsfile, key, ext)
                except:
                    KeyError
            else:
                # change the keyword to specified value
                fits.setval(updated_fitsfile, key, value=specific_keys_dict[key])
            continue
        # get the index of the keyword previous to the one you want to add
        prev_key_idx = list(shkvd.keywd_dict.keys()).index(key) - 1
        # add the keyword in the right place from the right dictionary
        new_value = shkvd.keywd_dict[key]
        after_key = list(shkvd.keywd_dict.keys())[prev_key_idx]
        if after_key == 'wcsinfo':
            after_key = list(shkvd.keywd_dict.keys())[prev_key_idx - 1]
        if key != 'wcsinfo':
            print("adding keyword: ", key, " in extension: primary    after: ", after_key,
                  "    with value: ", new_value)
            # choose right value for GWA_XTIL according to the grating
            if key == "GWA_XTIL":
                grating = fits.getval(fits_file, "GRATING", 0)
                new_value = get_gwa_Xtil_val(grating, path_to_tilt_files)
                print("Replacing value of keyword ", key, " corresponding to GRATING=", grating,
                      " and value ", new_value)
            # choose right value for GWA_YTIL according to the grating
            if key == "GWA_YTIL":
                grating = fits.getval(fits_file, "GRATING", 0)
                new_value = get_gwa_Ytil_val(grating, path_to_tilt_files)
                print("Replacing value of keyword ", key, " corresponding to GRATING=", grating,
                      " and value ", new_value)
            fits.setval(updated_fitsfile, key, 0, value=new_value, after=after_key)
        else:
            # go into the subdictionary for WCS keywords
            extname = 'sci'
            sci_hdr = fits.getheader(updated_fitsfile, extname)
            main_hdr = fits.getheader(updated_fitsfile, 0)
            for subkey, new_value in shkvd.keywd_dict["wcsinfo"].items():
                # first remove these keywords from the main header
                if subkey in main_hdr:
                    # val = fits.getval(updated_fitsfile, subkey, 0)
                    # wcs_keywds_from_main_hdr[subkey] = val
                    fits.delval(updated_fitsfile, subkey, 0)
                # following commented lines are not working
                # uncomment this line if wanting to use the original keyword value given by create_data
                # new_value = wcs_keywds_from_main_hdr[subkey]
                if subkey not in sci_hdr:
                    print("adding keyword: ", subkey, " in extension: ", extname, " with value: ", new_value)
                    fits.setval(updated_fitsfile, subkey, 1, value=new_value, after='EXTNAME')
            print("Science header has been updated.")


def check_hdr_keywds(fits_file, only_update, mode_used):
    """
    This is the function that does all the work in this script (i.e. uses all other functions) to update the header
    Args:
        fits_file: string, name of the input file to be checked
        only_update: boolean, if False a new file will be created; True will only update

    Returns:
        Nothing. The outputs are a text file with all the added keywords and the new/updated fits file.
    """

    # read the keywords and corresponding values from fits file directly
    file_keywd_dict = read_hdrfits(fits_file)

    # create text file to log warnings
    print('')
    addedkeywds_file_name = create_addedkeywds_file(fits_file)

    # check the keywords
    print('\n   Starting keyword check...')
    warnings_list, missing_keywds = [], []
    specific_keys_dict = check_keywds(file_keywd_dict, addedkeywds_file_name, warnings_list, missing_keywds, mode_used)

    # if warnings text file is empty erase it
    check_addedkeywds_file(addedkeywds_file_name)

    # create new file with updated header or simply update the input fits file
    print('\n   Adding keywords...')
    add_keywds(fits_file, only_update, missing_keywds, specific_keys_dict)


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    parser.add_argument("-u",
                        dest="only_update",
                        action='store_true',
                        default=False,
                        help='Use -u if NOT wanting to create a new file with updated header.')
    parser.add_argument("-m",
                        dest="mode_used",
                        action='store',
                        default=None,
                        help='Observation mode used: FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, '
                             'msata, focus, mimf - The code is not case-sensitive.')
    args = parser.parse_args()

    # Set the variables
    fits_file = args.fits_file
    only_update = args.only_update
    mode_used = args.mode_used

    # Perform the keyword check
    check_hdr_keywds(fits_file, only_update, mode_used)

    print('\n * Script  hdr_keywd_check.py  finished * \n')


if __name__ == '__main__':
    sys.exit(main())
