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

# import the sample header keyword dictionary of level 2b
from .dict_info import level2b_hdr_keywd_dict_sample as lev2bdict
# import subarray dictionary
from .dict_info import subarray_dict as subdict


'''
This script checks that the fits files to be used as input for the pipeline build 7.3, have the expected keywords in
the primary and science headers.

Example usage:
    The code works from the terminal.
    To create a NEW FS fits file with the updated header type:
        $ nptt_level2b_hdr_keywd_check blah.fits IFU

    To simply update the header of the existing fits file type:
        $ nptt_level2b_hdr_keywd_check blah.fits IFU -u

where the mode is either FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf.

The code is NOT case-sensitive.

If a mode is not provided, the code will look for a mode_used variable in the pytests configuration file, and it
will crash if this config file does not exist.

'''

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.6"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Apr 2019 - Version 1.1: added dictionary to choose right GWA_XTIL keyword value according to GRATING
# May 2019 - Version 1.2: added logic for dark processing
# Dec 2019 - Version 1.3: added logic for image processing
# Jul 2020 - Version 1.4: changed default value of SUBARRAY according to CRDS rules
# Aug 2020 - Version 1.5: fixed bug with set_exp_type_value function
# Feb 2021 - Version 1.6: implemented adding the MSA metafile name to the header


# Paths
wit4_path = os.environ.get('WIT4_PATH')
nirspec_cdp3 = "nirspec/CDP3/03_Instrument_model/3.1_Files/NIRS_FM2_05_CV3_FIT1/Description"
path_to_tilt_files = os.path.join(wit4_path, nirspec_cdp3)


# General functions

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
    print('(level2b_hdr_keword_check.read_hdrfits:) File contents')
    hdulist.info()
    # get and print header
    hdr = hdulist[0].header
    hdulist.close()
    return hdr


def create_addedkeywds_file(fits_file, mktxt=True):
    """
    This function create text file to log added keywords.
    Args:
        fits_file: string, name of the fits file to be keyword checked

    Returns:
        addedkeywds_file_name: string, the file name where all added keywords were saved
    """
    addedkeywds_file_name = fits_file.replace(".fits", "_addedkeywds.txt")
    if mktxt:
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
            os.system('rm '+addedkeywds_file_name)


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

    # Store type that value should have
    dict_type = type(hkwd_val[0])

    # Check if type of value correspond to what is given, else change it
    count = 0
    if val is not None:
        if isinstance(val, bool):
            if valtype == dict_type:
                warning = None
            else:
                warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'This keyword is boolean but not expected to be.')
            print(warning)
            val_and_valtype = [val, dict_type]
            return warning, val_and_valtype

        else:
            if isinstance(val, float):
                valtype = type(val)
                is_float = True
            else:   # maybe a float disguised as a string?
                is_float = False
                for v in val:
                    # check if value is a float
                    if v == '.':
                        count += 1

        if not is_float:
            # check if value has a negative sign
            neg_in_value = False
            if (count == 0) and ('-' in val):
                val_list = val.split("-")
                if len(val_list) == 2:
                    neg_in_value = True

            # check for letters in value
            no_letters_in_string = True
            for char in val:
                if char.isalpha():
                    no_letters_in_string = False
            if no_letters_in_string:
                if (count == 0) and (':' not in val) and neg_in_value:
                    if val == '':
                        warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'This keyword has an empty value')
                        print(warning)
                        val_and_valtype = [val, dict_type]
                        return warning, val_and_valtype

                    val = int(val)
                if (count == 1) and (':' not in val):
                    val = float(val)
                valtype = type(val)

    if (valtype in hkwd_val) or (val in hkwd_val) or (valtype == dict_type):
        print('{:<15} {:<9} {:<25}'.format(key, ext, 'Allowed value type'))
        warning = None
    else:
        warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect value type. Expected e.g. ' + repr(hkwd_val[0])
                                               + ', got: '+repr(val))
        print(warning)
        # if the gotten value contains letters then chenge it for the dictionary value, otherwise just return the type
        # of value that it should be changed to
        if re.search('[a-zA-Z]', str(val)):
            val = hkwd_val[0]

    val_and_valtype = [val, dict_type]
    return warning, val_and_valtype


def check3numbers(key, val, ext='primary', verbose=False):
    """
    Check if this keyword value has a format like: 0.1.1
    Args:
        key: keyword
        val: value of that keyword
        ext: string, extension number (default value is 'primary')

    Returns:
        A warning (if the format is not what is expected) or nothing if it is.

    """
    warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect value format. Expected: 0.0.0, got: '+str(val))
    r = re.compile('\d.\d.\d')  # check that the string has a structure like 0.1.1
    if r.match(val) is None:
        if verbose:
            print(warning)
        return warning
    else:
        if verbose:
            print('{:<15} {:<9} {:<25}'.format(key, ext, 'Matches expected format'))


def check_len(key, val, val_len=2, ext='primary', verbose=False):
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
    warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect length of value. Expected: ' + repr(val_len) +
                                           ', got ' + repr(string_length))
    if string_length == val_len:
        if verbose:
            print('{:<15} {:<9} {:<25}'.format(key, ext, 'Correct format'))
    else:
        if verbose:
            print(warning)
        return warning


def check_datetimeformat(key, val, check_time, check_date, check_datetime, ext='primary', verbose=False):
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
        val = vlist[0]+':'+vlist[1]+':'+v
    if check_time:
        val = datetime.strptime(val, '%H:%M:%S')
    if check_date:
        val = datetime.strptime(val, '%Y-%m-%d')
    if check_datetime:
        val = datetime.strptime(val, '%Y-%m-%dT%H:%M:%S')
    if isinstance(val, datetime):
        if verbose:
            print('{:<15} {:<9} {:<25}'.format(key, ext, 'Correct value format'))
    else:
        if verbose:
            print(warning)
        return warning


def get_gwa_Xtil_val(grating, path_to_tilt_files, verbose=False):
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
            if verbose:
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


def get_gwa_Ytil_val(grating, path_to_tilt_files, verbose=False):
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
            if verbose:
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


# keyword and format check

def set_exp_type_value(mode_used):
    """
    This function sets the appropriate value according to the mode.
    Args:
        mode_used: string
    Returns:
        val: string, expected pipeline value
    """
    # make sure there are no white spaces
    if " " in mode_used:
        mode_used = mode_used.replace(" ", "")
    print('   * MODE_USED  = ', mode_used)
    val = None
    if "fs" in mode_used.lower():
        val = 'NRS_FIXEDSLIT'
    if "ifu" in mode_used.lower():
        val = 'NRS_IFU'
    if "mos" in mode_used.lower() or mode_used.lower() == "msa":
        val = 'NRS_MSASPEC'
    if "bots" in mode_used.lower():
        val = 'NRS_BRIGHTOBJ'
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
    if val is None:
        print("\nWARNING: Cannot determine the EXP_TYPE. Try again with one of these modes: FS, MOS, IFU, BOTS, \n"
              "         dark, image, confirm, taconfirm, wata, msata, focus, mimf")
        print("         Exiting script level2b_hdr_keywd_check. \n")
        exit()
    print('     corresponds to EXP_TYPE =', val)
    return val


def determine_subarray(key, ff, detector, grating, specific_keys_dict, missing_keywds):
    data = fits.getdata(ff, 1)
    subsize1, subsize2 = np.shape(data)
    pipe_subarr_val = None
    for sa in subdict.subarray_dict:
        ssz1 = subdict.subarray_dict[sa]["subsize1"]
        ssz2 = subdict.subarray_dict[sa]["subsize2"]
        if subsize1 == ssz1:
            if subsize2 == ssz2:
                pipe_subarr_val = sa
                sst1 = subdict.subarray_dict[sa]["substrt1"]
                sst2_dict = subdict.subarray_dict[sa]["substrt2"]
                for grat, sst2_tuple in sst2_dict.items():
                    if grat.lower() == grating.lower():
                        if "1" in detector:
                            sst2 = sst2_tuple[0]
                        elif "2" in detector:
                            sst2 = sst2_tuple[1]
                        break

    if pipe_subarr_val is not None:
        specific_keys_dict[key] = pipe_subarr_val
        specific_keys_dict['SUBSTRT1'] = sst1
        specific_keys_dict['SUBSIZE1'] = ssz1
        specific_keys_dict['SUBSTRT2'] = sst2
        specific_keys_dict['SUBSIZE2'] = ssz2
        missing_keywds.append(key)
        missing_keywds.append('SUBSTRT1')
        missing_keywds.append('SUBSIZE1')
        missing_keywds.append('SUBSTRT2')
        missing_keywds.append('SUBSIZE2')
    return pipe_subarr_val, specific_keys_dict, missing_keywds


def get_pipe_subarray_name(subarray):
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
    return pipe_subarr_val


def set_pysiaf_keywords(input_fits_file, mode, FXD_SLIT, detector):
    # set up these keywords from SIAF
    NIRSpec_SIAF = pysiaf.Siaf('NIRSpec')
    if 'ifu' in mode.lower():
        aperture_name = detector + '_FULL_IFU'
        return aperture_name
    if 'full' in FXD_SLIT.lower():
        aperture_name = detector + '_FULL'
    elif '200' in FXD_SLIT or '400' in FXD_SLIT:
        aperture_name = 'NRS_' + FXD_SLIT + '_SLIT'
    else:
        aperture_name = 'NRS_S1600A1_SLIT'
    print('PySIAF aperture name: ', aperture_name)
    refpoint = NIRSpec_SIAF[aperture_name].reference_point('tel')
    V2_REF, V3_REF = refpoint[0], refpoint[1]
    V3IdlYAngle = NIRSpec_SIAF[aperture_name].V3IdlYAngle
    VIdlParity = NIRSpec_SIAF[aperture_name].VIdlParity
    fits.setval(input_fits_file, 'V2_REF', value=V2_REF, extname='SCI')
    fits.setval(input_fits_file, 'V3_REF', value=V3_REF, extname='SCI')
    fits.setval(input_fits_file, 'V3I_YANG', value=V3IdlYAngle, extname='SCI')
    fits.setval(input_fits_file, 'VPARITY', value=VIdlParity, extname='SCI')


def check_keywds(file_keywd_dict, warnings_file_name, warnings_list, missing_keywds, mode_used, detector=None,
                 subarray=None, msa_metafile=None, verbose=False):
    """
    This function will check keywords against those in hdr_keywod_dict.py
    Args:
        file_keywd_dict: dictionary of the original keyword header
        warnings_file_name: string, name of the file to reccord added keywords
        warnings_list: list of the warnings to be written into the file
        missing_keywds: list of the keywords not in the original the header
        mode_used: str or None, observation mode used FS, MOS, or IFU (if None then a configuration file
                    is expected to exist and contain a variable named mode_used)
        detector: string, expects NRS1, NRS2, or None (in this case it will be read from the header)
        subarray: None or string, name of the subarray to use
        msa_metafile: None or string, name of the MSA metafile
        verbose: boolean

    Returns:
        specific_keys_dict: dictionary with specific keys and values that need to be changed
    """

    # get the name and path of the input fits file and its header
    ff = warnings_file_name.replace('_addedkeywds.txt', '.fits')
    original_sci_header = fits.getheader(ff, 1)

    # initialize the dictionary to hold keyword values that will be changed and written to the file
    specific_keys_dict = {}

    # get the detector and grating from the input file
    try:
        if detector is None:
            detector = file_keywd_dict["DETECTOR"]
        grating = file_keywd_dict["GRATING"]
    except KeyError:
        if detector is None:
            detector = file_keywd_dict["DET"]
        grating = file_keywd_dict["GWA_POS"]
    try:
        dateobs = file_keywd_dict["DATE-OBS"]
        timeobs = file_keywd_dict["TIME-OBS"]
    except KeyError:
        dateobs = file_keywd_dict["DATE"].split("T")[0]
        timeobs = file_keywd_dict["DATE"].split("T")[1]

    # loop through the keywords and values of the PTT dictionary and add keywords that are not in input file
    for lev2bdict_key, lev2bdict_val in lev2bdict.keywd_dict.items():
        # start by making the warning for each keyword None and assigning key and val
        key = lev2bdict_key
        val = lev2bdict_val

        # Check if keyword is in the file, expect for wcs info (these go in the science extension)
        if key != 'wcsinfo':
            ext = 'primary'
            if key not in file_keywd_dict:
                # make sure the MSA metafile is pointing to the right place
                if ('mos' in mode_used.lower()) or ("msa" in mode_used.lower()):
                    if key == 'MSAMETFL':
                            val = msa_metafile
                            if verbose:
                                print('     Setting value of ', key, ' to ', val)
                if key == 'EXP_TYPE':
                    val = set_exp_type_value(mode_used)
                if key == 'TSOVISIT':
                    if mode_used.lower() == 'bots':
                        val = True
                    else:
                        val = False
                if key == 'FXD_SLIT':
                    val = 'S200A1'
                    if mode_used.lower() == 'bots':
                        val = 'S1600A1'
                if key == 'DATE-OBS':
                    val = dateobs
                if key == 'TIME-OBS':
                    val = timeobs
                if key == 'SUBARRAY':
                    subar_data = determine_subarray(key, ff, detector, grating, specific_keys_dict, missing_keywds)
                    pipe_subarr_val, subarr_specific_keys_dict, subarr_missing_keywds = subar_data
                    specific_keys_dict.update(subarr_specific_keys_dict)
                    for item in subarr_missing_keywds:
                        missing_keywds.append(item)

                if key == 'SUBSTRT1' or key == 'SUBSTRT2' or key == 'SUBSIZE1' or key == 'SUBSIZE2':
                    continue

                specific_keys_dict[key] = val
                missing_keywds.append(key)
                warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'New keyword added to header')
                warnings_list.append(warning)

            else:
                # check if the keyword exists in the science header
                orig_val = None
                try:
                    orig_val = original_sci_header[key]
                except:
                    KeyError

                if orig_val is not None:
                    val = orig_val
                else:
                    val = file_keywd_dict[lev2bdict_key]

                # Check simple standard keyword values
                if type(val) == type(lev2bdict_val):
                    if verbose:
                        print('{:<15} {:<9} {:<25}'.format(key, ext, 'Has correct format'))
                    warning = None
                else:
                    if not isinstance(lev2bdict_val, list):
                        lev2bdict_val = [lev2bdict_val]
                    warning, val_and_valtype = check_value_type(key, val, lev2bdict_val)
                    val, dict_type = val_and_valtype
                    if warning is not None and "Incorrect value type" in warning:
                        if dict_type == int:
                            val = int(float(val))
                        elif dict_type == float:
                            val = float(val)
                        elif dict_type == str:
                            val = str(val)
                        specific_keys_dict[key] = val
                        missing_keywds.append(key)
                        if verbose:
                            print('     Setting value of ', key, ' to type ', dict_type, ' and value ', val)
                        warning = None

                # Check for specific keywords
                if key == 'DPSW_VER':
                    warning = check3numbers(key, val, verbose)
                elif (key == 'VISITGRP') or (key == 'ACT_ID'):
                    warning = check_len(key, val, val_len=2, verbose=verbose)
                elif (key == 'OBSERVTN') or (key == 'VISIT'):
                    warning = check_len(key, val, val_len=3, verbose=verbose)
                elif key == 'EXPOSURE':
                    warning = check_len(key, val, val_len=5, verbose=verbose)
                elif (key == 'DATE') or (key == 'VSTSTART'):
                    warning = check_datetimeformat(key, val, check_date=False, check_datetime=True,
                                                   check_time=False, verbose=verbose)
                elif key == 'DATE-OBS':
                    warning = check_datetimeformat(key, val, check_date=True, check_datetime=False,
                                                   check_time=False, verbose=verbose)
                elif key == 'TIME-OBS':
                    warning = check_datetimeformat(key, val, check_date=False, check_datetime=False,
                                                   check_time=True, verbose=verbose)

                # specific check for VISITYPE, set to GENERIC
                if key == 'VISITYPE':
                    if val != lev2bdict_val:
                        # for now always set this keyword to generic
                        if verbose:
                            print("Replacing ", key, fits.getval(ff, "VISITYPE", 0), "for GENERIC")
                        specific_keys_dict[key] = 'GENERIC'
                        missing_keywds.append(key)

                # make sure the MSA metafile is pointing to the right place
                if key == 'MSAMETFL':
                    if ('mos' in mode_used.lower()) or ("msa" in mode_used.lower()):
                        val = msa_metafile
                        specific_keys_dict[key] = val
                        missing_keywds.append(key)
                        if verbose:
                            print('     Setting value of ', key, ' to ', val)
                # specific check for SUBARRAY
                if key == 'SUBARRAY':
                    if 'ifu' in mode_used.lower() or "mos" in mode_used.lower():
                        if val != lev2bdict_val:
                            if verbose:
                                print("Replacing ", key, fits.getval(ff, "SUBARRAY", 0), "for N/A")
                            specific_keys_dict[key] = 'N/A'
                            missing_keywds.append(key)
                    else:  # set SUBARRAY for anything else other than IFU or MOS
                        if subarray is not None:
                            # force the subarray keyword to be set to input
                            if isinstance(subarray, bool):
                                subarray_info = determine_subarray(key, ff, detector, grating, specific_keys_dict,
                                                                   missing_keywds)
                                pipe_subarr_val, specific_keys_dict, missing_keywds = subarray_info
                            else:
                                pipe_subarr_val = get_pipe_subarray_name(subarray)
                                if '1600' in pipe_subarr_val:
                                    if mode_used.lower() == "fs":
                                        pipe_subarr_val = 'SUB2048'
                                    elif mode_used.lower() == "bots":
                                        # determine which 1600 subarray is it
                                        data = fits.getdata(ff, 1)
                                        subsize1, subsize2 = np.shape(data)
                                        for sa in subdict.subarray_dict:
                                            ssz1 = subdict.subarray_dict[sa]["subsize1"]
                                            ssz2 = subdict.subarray_dict[sa]["subsize2"]
                                            if subsize1 == ssz1:
                                                if subsize2 == ssz2:
                                                    pipe_subarr_val = sa
                                                    break
                                specific_keys_dict[key] = pipe_subarr_val
                            subarray = pipe_subarr_val
                            if verbose:
                                print("changing subarray keyword to ", pipe_subarr_val)
                            missing_keywds.append(key)
                            # and make sure to change the primary slit keyword accordingly
                            if mode_used.lower() == "fs" or mode_used.lower() == "bots":
                                if not pipe_subarr_val:
                                    FXD_SLIT = 'S200A1'
                                if '2048' or '1024' or '512' or '32' in pipe_subarr_val:
                                    FXD_SLIT = 'S1600A1'
                                specific_keys_dict['FXD_SLIT'] = FXD_SLIT
                                if verbose:
                                    print("changing primary slit keyword to FXD_SLIT=", pipe_subarr_val)
                                missing_keywds.append('FXD_SLIT')
                                specific_keys_dict[key] = pipe_subarr_val
                            # set the subarray sizes and start keywords accordingly
                            if subarray in subdict.subarray_dict:
                                ssz1 = subdict.subarray_dict[subarray]["subsize1"]
                                ssz2 = subdict.subarray_dict[subarray]["subsize2"]
                                sst1 = subdict.subarray_dict[subarray]["substrt1"]
                                sst2_dict = subdict.subarray_dict[subarray]["substrt2"]
                                for grat, sst2_tuple in sst2_dict.items():
                                    if grat.lower() == grating.lower():
                                        if "1" in detector:
                                            sst2 = sst2_tuple[0]
                                        elif "2" in detector:
                                            sst2 = sst2_tuple[1]
                                        break
                                specific_keys_dict['SUBSTRT1'] = sst1
                                specific_keys_dict['SUBSIZE1'] = ssz1
                                specific_keys_dict['SUBSTRT2'] = sst2
                                specific_keys_dict['SUBSIZE2'] = ssz2
                                missing_keywds.append('SUBSTRT1')
                                missing_keywds.append('SUBSIZE1')
                                missing_keywds.append('SUBSTRT2')
                                missing_keywds.append('SUBSIZE2')
                                if verbose:
                                    print("Subarray size and start keywords now set to: \n", )
                                    print("   substrt1=", sst1, " substrt2=", sst2, " subsize1=",
                                          ssz1, " subsize2=", ssz2)
                        else:
                            # determine subarray from size
                            subarray_info = determine_subarray(key, ff, detector, grating, specific_keys_dict,
                                                               missing_keywds)
                            pipe_subarr_val, specific_keys_dict, missing_keywds = subarray_info

                # check for right value for EXP_TYPE, default will be to add the sample value: NRS_MSASPEC
                if key == 'EXP_TYPE':
                    val = set_exp_type_value(mode_used)
                    if mode_used.lower() == "ifu":
                        specific_keys_dict['DATAMODL'] = 'IFUImageModel'
                        missing_keywds.append('DATAMODL')
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)
                    if verbose:
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
                    tf.write(warning+'\n')
        else:
            # add the WCS keywords to science extension
            for subkey, _ in lev2bdict_val.items():
                if subkey == 'V2_REF':
                    if isinstance(file_keywd_dict["SUBARRAY"], bool):
                        subarray_info = determine_subarray(key, ff, detector, grating, specific_keys_dict,
                                                           missing_keywds)
                        pipe_subarr_val, _, _ = subarray_info
                    else:
                        pipe_subarr_val = get_pipe_subarray_name(file_keywd_dict["SUBARRAY"])
                    set_pysiaf_keywords(ff, mode_used, pipe_subarr_val, detector)
                elif subkey == 'V3_REF' or subkey == 'V3I_YANG' or subkey == 'VPARITY':
                    continue
                else:
                    missing_keywds.append(key)
                    specific_keys_dict[key] = val
                # now add the keyword to in the list to be added into the science extension
                warning = '{:<15} {:<9} {:<25}'.format(subkey, 'sci', 'New keyword added to header')
                warnings_list.append(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')

    if verbose:
        print("keywords to be modified: ", list(OrderedDict.fromkeys(missing_keywds)))
    return specific_keys_dict, missing_keywds


def add_keywds(fits_file, only_update, missing_keywds, specific_keys_dict, mode_used, verbose=False):
    """
    This function adds the missing keywords from the hdr_keywords_dictionary.py (hkwd) file and gives
    the fake values taken from the dictionary sample_hdr_keywd_vals_dict.py (shkvd).
    Args:
        only_update: If false a copy of the original fits file will be created with the
                     updated header.
        missing_keywds: list, missing keywords will be appended here
        specific_keys_dict: dictionary with specific keys and values that need to be changed
        mode_used: str or None, observation mode used FS, MOS, or IFU (if None then a configuration file
                    is expected to exist and contain a variable named mode_used)
        verbose: boolean
    """
    missing_keywds = list(OrderedDict.fromkeys(missing_keywds))
    # create name for updated fits file
    updated_fitsfile = fits_file
    if not only_update:
        updated_fitsfile = fits_file.replace('.fits', '_updatedHDR.fits')
        subprocess.run(["cp", fits_file, updated_fitsfile])
    # add missimg keywords
    if verbose:
        print('Saving keyword values in file: ', updated_fitsfile)
    ext = 0
    for i, key in enumerate(missing_keywds):
        if key != 'wcsinfo':
            if key in specific_keys_dict:
                if specific_keys_dict[key] == 'remove':
                    # keyword to be deleted
                    try:
                        fits.delval(updated_fitsfile, key, ext)
                    except:
                        KeyError
                else:
                    # change the keyword to specified value
                    if verbose:
                        print('setting key ', key, ' to ', specific_keys_dict[key])
                    fits.setval(updated_fitsfile, key, value=specific_keys_dict[key])
            continue

        # get the index of the keyword previous to the one you want to add
        prev_key_idx = list(lev2bdict.keywd_dict.keys()).index(key) - 1
        # add the keyword in the right place from the right dictionary
        new_value = lev2bdict.keywd_dict[key]
        after_key = list(lev2bdict.keywd_dict.keys())[prev_key_idx]
        if after_key == 'wcsinfo':
            after_key = list(lev2bdict.keywd_dict.keys())[prev_key_idx-1]
        if key != 'wcsinfo':
            # the DATAMODL keyword will only be modified if mode is IFU
            if key == 'DATAMODL':
                if 'ifu' in mode_used.lower():
                    new_value = 'IFUImageModel'
                else:
                    new_value = 'ImageModel'
            # for crm simulation files, get keywords from ESA file
            if 'NRSV' not in updated_fitsfile:
                continue
            if key == 'DETECTOR':
                new_value = fits.getval(updated_fitsfile, 'DET', 0)
                grating = new_value
            if key == 'GRATING':
                new_value = fits.getval(updated_fitsfile, 'GWA_POS', 0)
                grating = new_value
            if key == 'FILTER':
                new_value = fits.getval(updated_fitsfile, 'FWA_POS', 0)
            # choose right value for GWA_XTIL according to the grating
            if key == "GWA_XTIL":
                new_value = get_gwa_Xtil_val(grating, path_to_tilt_files, verbose=verbose)
                if verbose:
                    print("Replacing value of keyword ", key, " corresponding to GRATING=", grating, "and value ", new_value)
            # choose right value for GWA_YTIL according to the grating
            if key == "GWA_YTIL":
                new_value = get_gwa_Ytil_val(grating, path_to_tilt_files, verbose=verbose)
                if verbose:
                    print("Replacing value of keyword ", key, " corresponding to GRATING=", grating, "and value ", new_value)
            if key == "EXP_TYPE":
                new_value = set_exp_type_value(mode_used)
            fits.setval(updated_fitsfile, key, 0, value=new_value, after=after_key)
            if verbose:
                print("adding keyword: ", key, " in extension: PRIMARY    after: ", after_key, "   value: ", new_value)
        else:
            # go into the sub-dictionary for WCS keywords
            extname = 'sci'
            sci_hdr = fits.getheader(updated_fitsfile, extname)
            main_hdr = fits.getheader(updated_fitsfile, 0)
            for subkey, new_value in lev2bdict.keywd_dict["wcsinfo"].items():
                # first remove these keywords from the main header
                if subkey in main_hdr:
                    # val = fits.getval(updated_fitsfile, subkey, 0)
                    # wcs_keywds_from_main_hdr[subkey] = val
                    try:
                        fits.delval(updated_fitsfile, subkey, 0)
                    except:
                        KeyError
                # following commented lines are not working
                # uncomment this line if wanting to use the original keyword value given by create_data
                # new_value = wcs_keywds_from_main_hdr[subkey]
                if subkey not in sci_hdr:
                    if verbose:
                        print("adding keyword: ", subkey, " in extension: ", extname, " with value: ", new_value)
                    fits.setval(updated_fitsfile, subkey, 1, value=new_value, after='EXTNAME')
    # if the keyword is not in the sample dictionary then remove it
    for key in lev2bdict.keywd_dict:
        if key != 'wcsinfo':
            try:
                fits.delval(updated_fitsfile, key, 1)
            except:
                KeyError
    if verbose:
        print("Main and science headers have been updated.")
    return updated_fitsfile


def check_lev2b_hdr_keywd(fits_file, only_update, mode_used, detector=None, subarray=None, msa_metafile=None,
                          mktxt=True, verbose=True):
    """
    This is the function that does all the work in this script (i.e. uses all other functions) to update the header
    Args:
        fits_file: string, name of the input file to be checked
        only_update: boolean, if False a new file will be created; True will only update
        mode_used: string, observation mode used FS, MOS, or IFU (if None then a configuration file
                    is expected to exist and contain a variable named mode_used)
        detector: string, expects NRS1, NRS2, or None (in this case it will be read from the header)
        subarray: None or string, name of subarray to use
        msa_metafile: None or string, name of the MSA metafile
        mktxt: boolean, create a text file with all keywords added
        verbose: boolean

    Returns:
        updated_fitsfile: string, path and name of the outputs are a text file with all the added keywords and the
                          new/updated fits file.
    """

    # read the keywords and corresponding values from fits file directly
    file_keywd_dict = read_hdrfits(fits_file)

    # create text file to log warnings
    addedkeywds_file_name = create_addedkeywds_file(fits_file, mktxt=mktxt)

    # check the keywords
    print('\n(level2b_hdr_keywd_check.check_lev2b_hdr_keywd:) Starting keyword check...')
    warnings_list, missing_keywds = [], []
    specific_keys_dict, missing_keywds = check_keywds(file_keywd_dict, addedkeywds_file_name, warnings_list,
                                                      missing_keywds, mode_used, detector, subarray,
                                                      msa_metafile, verbose=verbose)

    # if warnings text file is empty erase it
    check_addedkeywds_file(addedkeywds_file_name)

    # create new file with updated header or simply update the input fits file
    print('\n(level2b_hdr_keywd_check.check_lev2b_hdr_keywd:) Adding and/or changing keywords...')
    updated_fitsfile = add_keywds(fits_file, only_update, missing_keywds, specific_keys_dict, mode_used,
                                  verbose=verbose)
    if not mktxt:
        os.system('rm ' + addedkeywds_file_name)

    return updated_fitsfile


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    parser.add_argument("mode_used",
                        action='store',
                        default=None,
                        help='Observation mode used: FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata,'
                             ' focus, mimf - The code is not case-sensitive.')
    parser.add_argument("-u",
                        dest="only_update",
                        action='store_true',
                        default=False,
                        help='Use -u if NOT wanting to create a new file with updated header.')
    parser.add_argument("-d",
                        dest="detector",
                        action='store',
                        default=None,
                        help='Use -d to provide the detector: -d=NRS1 or -d=NRS2.')
    parser.add_argument("-m",
                        dest="msa_metafile",
                        action='store',
                        default=None,
                        help='Use -m to provide the msa metafile name, e.g. -m=V9621500100101_msa.fits.')
    parser.add_argument("-v",
                        dest="verbose",
                        action='store_true',
                        default=False,
                        help='Use -v to print on-screen keywords and values.')
    args = parser.parse_args()

    # Set the variables
    fits_file = args.fits_file
    mode_used = args.mode_used
    only_update = args.only_update
    detector = args.detector
    msa_metafile = args.msa_metafile
    verbose = args.verbose

    # Perform the keyword check
    check_lev2b_hdr_keywd(fits_file, only_update, mode_used, detector=detector, msa_metafile=msa_metafile,
                          verbose=verbose)

    print('\n * Script  level2b_hdr_keywd_check.py  finished * \n')


if __name__ == '__main__':
    sys.exit(main())

