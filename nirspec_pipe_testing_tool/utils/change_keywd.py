import argparse
import os
import sys
from astropy.io import fits

from .dict_info import hdr_keywd_dict_sample as shkvd

"""
This script changes the value of a keyword to the specified value. If no extension is given, the keyword will be
added/modified in the main header.

Usage:
    In a terminal, the general format for the script is as follows:
    > python change_keywd.py file.fits keyword_to_be_changed value

    Specific examples:
    > python change_keywd.py fits_file.fits subarray full

    Remember that if you are changing the value from True to False, the pipeline only takes these as either T or F:
    > python change_keywd.py fits_file.fits TARGOOPP F

    To modify a keyword in extension 1 add -e=1:
    > python change_keywd.py fits_file.fits V3_REF -202.2 -e=1

"""


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed


def keywd_in_dict_or_new_keywd(sample_kyewd_dict, extension_number, value):
    if keyword in sample_kyewd_dict:
        # obtain the type of the sample value and convert the input value
        dict_value_type = type(sample_kyewd_dict[keyword])
        new_value = dict_value_type(value)
    else:
        new_value = new_keywd(extension_number, value)
    return new_value


def new_keywd(extension_number, value):
    add_keywd = input("This keyword is not in the sample keyword dictionary in extension "+repr(extension_number)+". Would you like to add it?  [y]  n   ")
    if "n" in add_keywd:
        print ("Exiting the code.")
        exit()
    else:
        val_type = input("Please type the value type of the keyword value (use one of the following: [str], float, int)   ")
        if (val_type == "float") or (val_type == "int"):
            new_value = val_type(value)
        else:
            new_value = value
    return new_value



def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    parser.add_argument("keyword",
                        action='store',
                        default=None,
                        help='Name of the keyword to be changed, it does not need to be all capital letters.')
    parser.add_argument("value",
                        action='store',
                        default=None,
                        help='New value for the given keyword.')
    parser.add_argument("-e",
                        dest="ext_number",
                        action='store',
                        default=None,
                        help='Save a text file with the header information of the given extension.')
    args = parser.parse_args()

    # Set the variables
    fits_file = args.fits_file
    keyword = args.keyword
    value = args.value
    ext_number = args.ext_number
    if ext_number is None:
        extension_number = 0
    else:
        extension_number = int(ext_number)


    # convert the keyword value input to the right type
    keyword = keyword.upper()
    if extension_number == 0:
        sample_kyewd_dict = shkvd.keywd_dict
        new_value = keywd_in_dict_or_new_keywd(sample_kyewd_dict, extension_number, value)

    elif extension_number == 1:
        sample_kyewd_dict = shkvd.keywd_dict["wcsinfo"]
        new_value = keywd_in_dict_or_new_keywd(sample_kyewd_dict, extension_number, value)

    else:
        new_value = new_keywd(extension_number, value)


    # now change the the value
    fits.setval(fits_file, keyword, extension_number, value=new_value)

    print ('\n * Script  change_keywd.py  finished * \n')


if __name__ == '__main__':
    sys.exit(main())
