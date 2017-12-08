import argparse
import collections
import os
import re
import numpy as np
from datetime import datetime
from astropy.io import fits

# import the header keyword dictionaries
import hdr_keywd_dict as hkwd
import hdr_keywd_dict_sample as shkvd


'''
This script checks that the fits files to be used as input for the pipeline build 7.1, have the expected keywords in
the primary and science headers.

Example usage:
    The code works from the terminal.
    To create a NEW FS fits file with the updated header type:
        > python /path_to_this_script/hdr_keywd_check.py blah.fits

    To simply update the header of the existing fits file type:
        > python /path_to_this_script/hdr_keywd_check.py blah.fits -u

'''

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
    # Check if type of value correspond to what is given, else change it
    if val is not None:
        count = 0
        for v in val:
            if v=='.':
                count += 1
        no_letters_in_string = True
        for char in val:
            if char.isalpha():
                no_letters_in_string = False
        if no_letters_in_string:
            if (count == 0) and (':' not in val) and ('-' not in val):
                if val=='':
                    warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'This keyword has an empty value')
                    print (warning)
                    return warning
                val = int(val)
            if (count==1) and (':' not in val):
                val = float(val)
            valtype = type(val)
    if (valtype in hkwd_val) or (val in hkwd_val):
        print ('{:<15} {:<9} {:<25}'.format(key, ext, 'Allowed value type'))
        warning = None
    else:
        warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'Incorrect value type. Expected: '+str(val)+', got: '+str(valtype))
        print (warning)
    return warning


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


### keyword and format check

def check_keywds(file_keywd_dict, warnings_file_name, warnings_list, missing_keywds):
    """
    This function will check keywords against those in hdr_keywod_dict.py
    Args:
        file_keywd_dict: dictionary of the original keyword header
        warnings_file_name: string, name of the file to reccord added keywords
        warnings_list: list of the warnings to be written into the file
        missing_keywds: list of the keywords not in the original the header

    Returns:
        nothing
    """
    for hkwd_key, hkwd_val in hkwd.keywd_dict.items():
        # start by making the warning for each keyword None and assigning key and val
        key = hkwd_key

        # Check if keyword is in the file, expect for wcs info (these go in the science extension)
        if key != 'wcsinfo':
            ext = 'primary'
            if key not in file_keywd_dict:
                missing_keywds.append(key)
                warning = '{:<15} {:<9} {:<25}'.format(key, ext, 'New keyword added to header')
                warnings_list.append(warning)
            else:
                val = file_keywd_dict[hkwd_key]

                # Check simple standard keyword values
                if val in hkwd_val:
                    print ('{:<15} {:<9} {:<25}'.format(key, ext, 'Has correct format'))
                    warning = None
                else:
                    warning = check_value_type(key, val, hkwd_val)

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

            if warning is not None:
                missing_keywds.append(key)
                warnings_list.append(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')
        else:
            # add the WCS keywords to science extension
            missing_keywds.append(key)
            for subkey, _ in hkwd_val.items():
                # now add the keyword to in the list to be added into the science extension
                warning = '{:<15} {:<9} {:<25}'.format(subkey, 'sci', 'New keyword added to header')
                warnings_list.append(warning)
                with open(warnings_file_name, "a") as tf:
                    tf.write(warning+'\n')
    print("missing_keywords = ", missing_keywds)


def add_keywds(fits_file, only_update, missing_keywds):
    '''
    This function adds the missing keywords from the hdr_keywords_dictionary.py (hkwd) file and gives
    the fake values taken from the dictionary sample_hdr_keywd_vals_dict.py (shkvd).
    Args:
        only_update: If false a copy of the original fits file will be created with the
                     updated header.
        missing_keywds: list, missing keywords will be appended here
    '''
    # create name for updated fits file
    updated_fitsfile = fits_file
    if not only_update:
        updated_fitsfile = fits_file.replace('.fits', '_updatedHDR.fits')
        os.system('cp '+fits_file+' '+updated_fitsfile)
    # add missimg keywords
    wcs_keywds_from_main_hdr = {}
    print ('Path of updated file: ', updated_fitsfile)
    for i, key in enumerate(missing_keywds):
        # get the index of the keyword previous to the one you want to add
        prev_key_idx = list(hkwd.keywd_dict.keys()).index(key) - 1
        # add the keyword in the right place from the right dictionary
        new_value = shkvd.keywd_dict[key]
        after_key = list(shkvd.keywd_dict.keys())[prev_key_idx]
        if key != 'wcsinfo':
            print("adding keyword: ", key, " in extension: primary")
            fits.setval(updated_fitsfile, key, value=new_value, after=after_key)
            #print ('\n New header: ')
            #hdulist = fits.open(updated_fitsfile)
            #hdr = hdulist[0].header
            #print (repr(hdr))
        else:
            # go into the subdictionary for WCS keywords
            extname = 'sci'
            main_hdr = fits.getheader(updated_fitsfile, 0)
            for subkey, new_value in shkvd.keywd_dict["wcsinfo"].items():
                # first remove these keywords from the main header
                if subkey in main_hdr:
                    # uncomment this section if wanting to save the original keyword value
                    #val = file_keywd_dict[hkwd_key]
                    #wcs_keywds_from_main_hdr[hkwd_key] = val
                    fits.delval(updated_fitsfile, subkey, 0)
                # uncomment this line if wanting to use the original keyword value given by create_data
                #new_value = wcs_keywds_from_main_hdr[subkey]
                print("adding keyword: ", subkey, " in extension: ", extname, " with value: ", new_value)
                fits.setval(updated_fitsfile, subkey, 1, value=new_value, after='EXTNAME')
            print ("Science header has been updated.")



def perform_check(fits_file, only_update):
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
    check_keywds(file_keywd_dict, addedkeywds_file_name, warnings_list, missing_keywds)

    # if warnings text file is empty erase it
    check_addedkeywds_file(addedkeywds_file_name)

    # create new file with updated header or simply update the input fits file
    print('\n   Adding keywords...')
    add_keywds(fits_file, only_update, missing_keywds)



if __name__ == '__main__':

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
    args = parser.parse_args()

    # Set the variables
    fits_file = args.fits_file
    only_update = args.only_update

    # Perform the keyword check
    perform_check(fits_file, only_update)

    print ('\n * Script  hdr_keywd_check.py  finished * \n')


