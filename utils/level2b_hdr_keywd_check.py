import argparse
import collections
import os
import re
import subprocess
import numpy as np
from datetime import datetime
from astropy.io import fits
from collections import OrderedDict

# import the sample header keyword dictionary of level 2b
import level2b_hdr_keywd_dict_sample as lev2bdict
# import subarray dictionary
import subarray_dict as subdict


'''
This script checks that the fits files to be used as input for the pipeline build 7.1, have the expected keywords in
the primary and science headers.

Example usage:
    The code works from the terminal.
    To create a NEW FS fits file with the updated header type:
        > python /path_to_this_script/level2b_hdr_keywd_check.py blah.fits IFU

    To simply update the header of the existing fits file type:
        > python /path_to_this_script/level2b_hdr_keywd_check.py blah.fits IFU -u

where the mode is either FS, MOS, IFU, BOTS. If a mode is not provided, the code will look for a mode_used variable
in the pytests configuration file, and it will crash if this config file does not exist.

'''

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.1"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Apr 2019 - Version 1.1: added dictionary to choose right GWA_XTIL keyword value according to GRATING


### General functions

def read_hdrfits(fits_file_name):
    '''
    This function reads the header fits file and returns a dictionary of the keywords with
    corresponding values. Keywords will be stored in the order they are read.
    Args:
        hdr_txt_file: full path with name of the header text file

    Returns:
        A dictionary of keywords with corresponding values
    '''
    #  Read the fits file
    hdulist = fits.open(fits_file_name)
    # print on screen what extensions are in the file
    #print ('\n FILE INFORMATION: \n')
    hdulist.info()
    # get and print header
    #print ('\n FILE HEADER: \n')
    hdr = hdulist[0].header
    sci_hdr = hdulist[1].header
    #print (repr(hdr))
    # close the fits file
    hdulist.close()
    # set the name of the text file and save the header
    text_file_name = fits_file_name.replace('.fits', '_header.txt')
    tf = open(text_file_name, 'w')
    tf.write(repr(hdr))
    tf.close()
    # read the text file
    keywd_dict = read_hdrtxt(text_file_name)
    # remove the text file
    os.system("rm "+text_file_name)
    return keywd_dict


def read_hdrtxt(hdr_txt_file):
    '''
    This function reads the header text file and returns a dictionary of the keywords with
    corresponding values. Keywords will be stored in the order they are read.
    Args:
        hdr_txt_file: full path with name of the header text file

    Returns:
        A dictionary of keywords with corresponding values
    '''
    keywd_dict = collections.OrderedDict()
    with open(hdr_txt_file, 'r') as htf:
        for line in htf.readlines():   # identify keywords by lines containing a =
            if '=' in line:
                line_list = line.split('=')
                keywd = line_list[0].split()[0]   # remove the white spaces from the keyword
                keywd_val = line_list[1].split()[0]  # remove the white spaces from the keyword value
                if "'" in keywd_val:
                    keywd_val = keywd_val.replace("'", "")   # remove the extra '
                keywd_dict[keywd] = keywd_val   # add dictionary entry
    return keywd_dict


def create_addedkeywds_file(fits_file):
    """
    This function create text file to log added keywords.
    Args:
        fits_file: string, name of the fits file to be keyword checked

    Returns:
        addedkeywds_file_name: string, the file name where all added keywords were saved
    """
    addedkeywds_file_name = fits_file.replace(".fits", "_addedkeywds.txt")
    print ('Name of text file containing all keywords added:  ', addedkeywds_file_name)
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


### Functions to check specific keyword values

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
    if val is not None:
        count = 0
        for v in val:
            # check if value is a float
            if v=='.':
                count += 1

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
            if (count == 0) and (':' not in val) and (neg_in_value):
                if val=='':
                    warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'This keyword has an empty value')
                    print (warning)
                    val_and_valtype = [val, dict_type]
                    return warning, val_and_valtype

                val = int(val)
            if (count==1) and (':' not in val):
                val = float(val)
            valtype = type(val)

    if (valtype in hkwd_val) or (val in hkwd_val) or (valtype == dict_type):
        print ('{:<15} {:<9} {:<25}'.format(key, ext, 'Allowed value type'))
        warning = None
    else:
        warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect value type. Expected e.g. '+repr(hkwd_val[0])+', got: '+repr(val))
        print (warning)
        # if the gotten value contains letters then chenge it for the dictionary value, otherwise just return the type
        # of value that it should be changed to
        if re.search('[a-zA-Z]', str(val)):
            val = hkwd_val[0]

    val_and_valtype = [val, dict_type]
    return warning, val_and_valtype


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
    warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect value format. Expected: 0.0.0, got: '+str(val))
    r = re.compile('\d.\d.\d') # check that the string has a structure like 0.1.1
    if r.match(val) is None:
        print (warning)
        return warning
    else:
        print ('{:<15} {:<9} {:<25}'.format(key, ext, 'Matches expected format'))


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
    warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect length of value. Expected: '+repr(val_len)+', got '+repr(string_length))
    if string_length == val_len:
        print ('{:<15} {:<9} {:<25}'.format(key, ext, 'Correct format'))
    else:
        print (warning)
        return warning


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
        val = vlist[0]+':'+vlist[1]+':'+v
    if check_time:
        val = datetime.strptime(val, '%H:%M:%S')
    if check_date:
        val = datetime.strptime(val, '%Y-%m-%d')
    if check_datetime:
        val = datetime.strptime(val, '%Y-%m-%dT%H:%M:%S')
    if isinstance(val, datetime):
        print ('{:<15} {:<9} {:<25}'.format(key, ext, 'Correct value format'))
    else:
        print (warning)
        return warning


def get_gwaxtil_val(grating):
    """
    This function gets the right GWA_XTIL value according to the grating given
    Args:
        grating: string, value from GRATING keyword

    Returns:
        gwaxtil: float, corresponding GWA_XTIL value to the grating
    """
    # dictionary containing default XTILT values
    # (taken from "Zeroreadings 1" parameter in "disperser_*_TiltY.gtp" reference files)
    default_xtil = {
            "G140H" : 0.361956834793,
            "G140M" : 0.331874251366,
            "G235H" : 0.353678375483,
            "G235M" : 0.321180671453,
            "G395H" : 0.323275774717,
            "G395M" : 0.28645208478,
            "MIRROR" : 0.347495883703,
            "PRISM" : 0.336739093065
            }
    gwaxtil = default_xtil[grating]
    return gwaxtil


### keyword and format check

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
    ff = warnings_file_name.replace('_addedkeywds.txt', '.fits')
    original_header = fits.getheader(ff, 1)
    specific_keys_dict = {}
    for lev2bdict_key, lev2bdict_val in lev2bdict.keywd_dict.items():
        # start by making the warning for each keyword None and assigning key and val
        key = lev2bdict_key

        # Check if keyword is in the file, expect for wcs info (these go in the science extension)
        if key != 'wcsinfo':
            ext = 'primary'
            if key not in file_keywd_dict:
                missing_keywds.append(key)
                warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'New keyword added to header')
                warnings_list.append(warning)
            else:
                # check if the keyword exists in the science header
                orig_val = None
                try:
                    orig_val = original_header[key]
                except:
                    KeyError

                if orig_val is not None:
                    val = orig_val
                else:
                    val = file_keywd_dict[lev2bdict_key]

                # Check simple standard keyword values
                if type(val) == type(lev2bdict_val):
                    print ('{:<15} {:<9} {:<25}'.format(key, ext, 'Has correct format'))
                    warning = None
                else:
                    if not isinstance(lev2bdict_val, list):
                        lev2bdict_val = [lev2bdict_val]
                    warning, val_and_valtype = check_value_type(key, val, lev2bdict_val)
                    val, dict_type = val_and_valtype
                    if warning is not None  and  "Incorrect value type" in warning:
                        if dict_type == int:
                            val = int(float(val))
                        elif dict_type == float:
                            val = float(val)
                        elif dict_type == str:
                            val = str(val)
                        specific_keys_dict[key] = val
                        missing_keywds.append(key)
                        print('     Setting value of ', key, ' to type ', dict_type, ' and value ', val)
                        warning = None

                # Check for specific keywords
                if key=='DPSW_VER':
                    warning = check3numbers(key, val)
                elif (key=='VISITGRP') or (key=='ACT_ID'):
                    warning = check_len(key, val, val_len=2)
                elif (key=='OBSERVTN') or (key=='VISIT'):
                    warning = check_len(key, val, val_len=3)
                elif (key=='EXPOSURE'):
                    warning = check_len(key, val, val_len=5)
                elif (key=='DATE') or (key=='VSTSTART'):
                    warning = check_datetimeformat(key, val, check_date=False, check_datetime=True,
                                                   check_time=False)
                elif key=='DATE-OBS':
                    warning = check_datetimeformat(key, val, check_date=True, check_datetime=False,
                                                   check_time=False)
                elif key=='TIME-OBS':
                    warning = check_datetimeformat(key, val, check_date=False, check_datetime=False,
                                                   check_time=True)

                # specific check for VISITYPE, set to GENERIC
                if key == 'VISITYPE':
                    if val != lev2bdict_val:
                        # for now always set this keyword to generic
                        print ("Replacing ", key, fits.getval(ff, "VISITYPE", 0), "for GENERIC")
                        #fits.setval(ff, key, 0, value='GENERIC')
                        specific_keys_dict[key] = 'GENERIC'
                        missing_keywds.append(key)

                # choose right value for GWA_XTIL according to the grating
                if key == "GWA_XTIL":
                    grating = fits.getval(ff, "GRATING", 0)
                    gwa_xtil = get_gwaxtil_val(grating)
                    print("Replacing value of keyword ", key, " corresponding to GRATING=", grating,
                          "and value ", gwa_xtil)
                    specific_keys_dict[key] = gwa_xtil
                    missing_keywds.append(key)

                # specific check for SUBARRAY
                if key == 'SUBARRAY':
                    if 'IFU' in mode_used  or  "MOS" in mode_used:
                        if val != lev2bdict_val:
                            print ("Replacing ", key, fits.getval(ff, "SUBARRAY", 0), "for GENERIC")
                            #fits.setval(ff, key, 0, value='GENERIC')
                            specific_keys_dict[key] = 'GENERIC'
                            missing_keywds.append(key)
                    elif mode_used == "FS"  or  mode_used == "BOTS":
                        # set the subarray according to size
                        substrt1 = fits.getval(ff, "SUBSTRT1", 0)
                        substrt2 = fits.getval(ff, "SUBSTRT2", 0)
                        subsize1 = fits.getval(ff, "SUBSIZE1", 0)
                        subsize2 = fits.getval(ff, "SUBSIZE2", 0)
                        for subarrd_key, subarrd_vals_dir in subdict.subarray_dict.items():
                            sst1 = subarrd_vals_dir["substrt1"]
                            sst2_list = subarrd_vals_dir["substrt2"]
                            ssz1 = subarrd_vals_dir["subsize1"]
                            ssz2 = subarrd_vals_dir["subsize2"]
                            if substrt1 == sst1  and  subsize1 == ssz1  and  subsize2 == ssz2:
                                for sst2 in sst2_list:
                                    if substrt2 == sst2:
                                        specific_keys_dict[key] = subarrd_key
                                        print ("changing subarray keyword to ", subarrd_key)
                                        missing_keywds.append(key)

                # check for right value for EXP_TYPE, default will be to add the sample value: NRS_MSASPEC
                if key == 'EXP_TYPE':
                    print('   * MODE_USED  = ', mode_used)
                    if 'FS' in mode_used:
                        val = 'NRS_FIXEDSLIT'
                    if 'IFU' in mode_used:
                        val = 'NRS_IFU'
                        missing_keywds.append('DATAMODL')
                    if 'MOS' in mode_used:
                        val = 'NRS_MSASPEC'
                    if mode_used == "BOTS":
                        val = 'NRS_BRIGHTOBJ'
                    specific_keys_dict[key] = val
                    missing_keywds.append(key)
                    print('     Setting value of ', key, ' to ', val)

                # make sure the MSASTATE keyword is set correctly
                if key == 'MSASTATE':
                    #mode_used = fits.getval(ff, 'MODEUSED', 0)
                    if (mode_used == 'FS') or (mode_used == 'IFU'):
                        val = 'PRIMARYPARK_ALLCLOSED'
                        specific_keys_dict[key] = val
                        missing_keywds.append(key)

            if warning is not None:
                missing_keywds.append(key)
                warnings_list.append(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')
        else:
            # add the WCS keywords to science extension
            missing_keywds.append(key)
            for subkey, _ in lev2bdict_val.items():
                # now add the keyword to in the list to be added into the science extension
                warning = '{:<15} {:<9} {:<25}'.format(subkey, 'sci', 'New keyword added to header')
                warnings_list.append(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')

    print("keywords to be modified: ", list(OrderedDict.fromkeys(missing_keywds)))
    return specific_keys_dict


def add_keywds(fits_file, only_update, missing_keywds, specific_keys_dict):
    '''
    This function adds the missing keywords from the hdr_keywords_dictionary.py (hkwd) file and gives
    the fake values taken from the dictionary sample_hdr_keywd_vals_dict.py (shkvd).
    Args:
        only_update: If false a copy of the original fits file will be created with the
                     updated header.
        missing_keywds: list, missing keywords will be appended here
        specific_keys_dict: dictionary with specific keys and values that need to be changed
    '''
    missing_keywds = list(OrderedDict.fromkeys(missing_keywds))
    #print ("specific_keys_dict = ", specific_keys_dict)
    # create name for updated fits file
    updated_fitsfile = fits_file
    if not only_update:
        updated_fitsfile = fits_file.replace('.fits', '_updatedHDR.fits')
        subprocess.run(["cp", fits_file, updated_fitsfile])
    # add missimg keywords
    wcs_keywds_from_main_hdr = {}
    print ('Saving keyword values in file: ', updated_fitsfile)
    ext = 0
    for i, key in enumerate(missing_keywds):
        if key in specific_keys_dict:
            #print ("found it in the dict: ", key, specific_keys_dict[key])
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
        prev_key_idx = list(lev2bdict.keywd_dict.keys()).index(key) - 1
        # add the keyword in the right place from the right dictionary
        new_value = lev2bdict.keywd_dict[key]
        after_key = list(lev2bdict.keywd_dict.keys())[prev_key_idx]
        if after_key == 'wcsinfo':
            after_key = list(lev2bdict.keywd_dict.keys())[prev_key_idx-1]
        if key != 'wcsinfo':
            print("adding keyword: ", key, " in extension: primary    after: ", after_key)
            # the DATAMODL keyword will only be modified if mode is IFU
            if key == 'DATAMODL':
                if 'IFU' in mode_used:
                    new_value = 'IFUImageModel'
                else:
                    new_value = 'ImageModel'
            if key == 'GRATING':
                new_value = fits.getval(updated_fitsfile, 'GWA_POS', 0)
            if key == 'FILTER':
                new_value = fits.getval(updated_fitsfile, 'FWA_POS', 0)
            fits.setval(updated_fitsfile, key, 0, value=new_value, after=after_key)
        else:
            # go into the subdictionary for WCS keywords
            extname = 'sci'
            sci_hdr = fits.getheader(updated_fitsfile, extname)
            main_hdr = fits.getheader(updated_fitsfile, 0)
            for subkey, new_value in lev2bdict.keywd_dict["wcsinfo"].items():
                # first remove these keywords from the main header
                if subkey in main_hdr:
                    #val = fits.getval(updated_fitsfile, subkey, 0)
                    #wcs_keywds_from_main_hdr[subkey] = val
                    try:
                        fits.delval(updated_fitsfile, subkey, 0)
                    except:
                        KeyError
                # following commented lines are not working
                # uncomment this line if wanting to use the original keyword value given by create_data
                #new_value = wcs_keywds_from_main_hdr[subkey]
                if subkey not in sci_hdr:
                    print("adding keyword: ", subkey, " in extension: ", extname, " with value: ", new_value)
                    fits.setval(updated_fitsfile, subkey, 1, value=new_value, after='EXTNAME')
    # if the keyword is not in the sample dictionary then remove it
    for key in lev2bdict.keywd_dict:
        if key != 'wcsinfo':
            try:
                fits.delval(updated_fitsfile, key, 1)
            except:
                KeyError
    print ("Main and science headers have been updated.")



def perform_check(fits_file, only_update, mode_used):
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



if __name__ == '__main__':

    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    parser.add_argument("mode_used",
                        #dest="mode_used",
                        action='store',
                        default=None,
                        help='Observation mode used: FS, MOS, or IFU.')
    parser.add_argument("-u",
                        dest="only_update",
                        action='store_true',
                        default=False,
                        help='Use -u if NOT wanting to create a new file with updated header.')
    args = parser.parse_args()

    # Set the variables
    fits_file = args.fits_file
    mode_used = args.mode_used
    only_update = args.only_update

    # Perform the keyword check
    perform_check(fits_file, only_update, mode_used)

    print ('\n * Script  level2b_hdr_keywd_check.py  finished * \n')


