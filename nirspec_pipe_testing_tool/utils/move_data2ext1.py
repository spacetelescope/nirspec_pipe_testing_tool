import argparse
import numpy as np
import sys
import os
from astropy.io import fits

'''
This script moves the science data to extension 1.

Example usage:
    The code works from the terminal.
    To create a NEW fits file with the updated header type:

        $ nptt_move_data2ext1.py blah.fits NRS1

        NRS1 (491) is the detector for which the data is being extracted, NRS2 (492) can also be used. To add
        the reference pixels add at the end of the command -rfpx, e.g.

        $ nptt_move_data2ext1.py blah.fits NRS1 -rfpx

'''

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.2"


# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Apr 2019 - Version 1.1: implemented option to add reference pixels
# Feb 2021 - Version 1.2: implemented option for output directory

# General functions


def move_data(input_fits_file, detector='NRS1', add_ref_pix=False, output_dir=None):
    """
    This function does the move of the data to the science extension.
    Args:
        input_fits_file: string, name of the input file
        detector: string, expects either NRS1 or NRS2
        add_ref_pix: boolean, if True the code will add the reference pixels
        output_dir: string, path to place the output file - if None output will be in same dir as input
    Output:
        Nothing. The function creates a new fits file with the primary and science data extensions
    """

    # create the fits list to hold the calculated flat values for each slit
    original_hdulist = fits.open(input_fits_file)
    outfile = fits.HDUList()
    outfile.append(original_hdulist[0])

    # move the the data to extension corresponding to de detector indicated in the function arguments
    det = 1  # default for NRS1=491
    if len(original_hdulist) > 6:
        if '2' in detector:
            det = 4
    input_fits_file_data = fits.getdata(input_fits_file, det)
    if add_ref_pix:
        print("(move_data2ext1.move_data:) Reference pixels will be added.")
        input_fits_file_data = add_ref_pixels(input_fits_file_data)
    outfile_ext = fits.ImageHDU(input_fits_file_data, name="SCI")
    outfile.append(outfile_ext)

    # do the same for the error extension
    input_fits_file_data = fits.getdata(input_fits_file, det + 1)
    if add_ref_pix:
        input_fits_file_data = add_ref_pixels(input_fits_file_data)
    outfile_ext = fits.ImageHDU(input_fits_file_data, name="ERR")
    outfile.append(outfile_ext)

    # do the same for the data quality extension
    input_fits_file_data = fits.getdata(input_fits_file, det + 2)
    if add_ref_pix:
        input_fits_file_data = add_ref_pixels(input_fits_file_data)
    outfile_ext = fits.ImageHDU(input_fits_file_data, name="DQ")
    outfile.append(outfile_ext)

    # write the new output file
    outfile_name = input_fits_file.replace(".fits", "_" + detector + "_modified.fits")
    if output_dir is not None:
        outfile_name = os.path.join(output_dir, os.path.basename(outfile_name))
        if ".fits" in output_dir:
            outfile_name = output_dir
    outfile.writeto(outfile_name, overwrite=True)
    original_hdulist.close()
    return outfile_name


def add_ref_pixels(data):
    """
    This function adds the reference pixels to the given data array.
    Args:
        data: the science, DQ, or error data

    Returns:

    """
    # define the correct size
    data_out = np.zeros((2048, 2048))

    # add the reference pixels
    for i in range(2040):
        for j in range(2040):
            data_out[j + 4, i + 1] = data[j, i]

    return data_out


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    parser.add_argument("detector",
                        action='store',
                        default=None,
                        help='The data for the given detector will be moved to a science extension.')
    parser.add_argument("-rfpx",
                        dest="add_ref_pix",
                        action='store_true',
                        default=False,
                        help='Add the reference pixels.')
    args = parser.parse_args()

    # Set the variables
    input_fits_file = args.fits_file
    detector = args.detector
    add_ref_pix = args.add_ref_pix

    # Perform data move to the science extension
    move_data(input_fits_file, detector, add_ref_pix)

    print('\n * Script  move_data2ext1.py  finished * \n')


if __name__ == '__main__':
    sys.exit(main())
