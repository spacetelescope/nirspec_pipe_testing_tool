#! /usr/bin/env python

from astropy.io import fits
import argparse
import sys

"""
This script reads the header of a fits file and gives info of its extensions.

Example usage:
    The code works from the terminal.
    To simply see the header on-screen type:
        > python /path_to_this_script/read_hdr.py blah.fits
    To see and save the header into a text file type:
        > python /path_to_this_script/read_hdr.py blah.fits -s
    To see and save the header of a different extension from the main, use -e=extension_number:
        > python /path_to_this_script/read_hdr.py blah.fits -s -e=1

"""


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed

def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file_name",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    parser.add_argument("-s",
                        dest="save_txt",
                        action='store_true',
                        default=False,
                        help='Save a text file with the header information.')
    parser.add_argument("-e",
                        dest="ext_number",
                        action='store',
                        default=None,
                        help='Save a text file with the header information of the given extension.')
    args = parser.parse_args()

    # Set the variables
    fits_file_name = args.fits_file_name
    save_txt = args.save_txt
    ext_number = args.ext_number
    if ext_number is None:
        extension_number = 0
    else:
        extension_number = int(ext_number)

    #  Read the fits file
    hdulist = fits.open(fits_file_name)

    # print on screen what extensions are in the file
    print ('\n FILE INFORMATION: \n')
    hdulist.info()
    print ()

    # get and print header
    print ('\n FILE HEADER FOR EXTENSION: ', extension_number, '\n')
    hdr = hdulist[extension_number].header
    print (repr(hdr))

    # close the fits file
    hdulist.close()

    # save the info into a text file
    if save_txt:
        # set the name of the text file
        text_file_name = fits_file_name.replace('.fits', '_hdr_ext'+repr(extension_number)+'.txt')
        tf = open(text_file_name, 'w')
        tf.write(repr(hdr))
        tf.close()


if __name__ == '__main__':
    sys.exit(main())
