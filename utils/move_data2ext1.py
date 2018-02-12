import argparse
from astropy.io import fits



'''
This script moves the science data to extension 1.

Example usage:
    The code works from the terminal.
    To create a NEW fits file with the updated header type:
        > python /path_to_this_script/hdr_keywd_check.py blah.fits


'''

### General functions

def move_data(input_fits_file):
    """
    This function does the move of the data to the science extension.
    Args:
        input_fits_file: string, name of the input file
    Output:
        Nothing. The function creates a new fits file with the primary and science data extensions
    """
    
    # create the fits list to hold the calculated flat values for each slit
    hdu0 = fits.PrimaryHDU()
    outfile = fits.HDUList()
    outfile.append(hdu0)
    
    # move the the data to extension 1
    input_fits_file_data = fits.getdata(input_fits_file)
    outfile_ext = fits.ImageHDU(input_fits_file_data, name="SCI")
    outfile.append(outfile_ext)

    # write the new output file
    outfile_name = input_fits_file.replace(".fits", "_modified.fits")
    outfile.writeto(outfile_name, overwrite=True)





if __name__ == '__main__':

    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    args = parser.parse_args()

    # Set the variables
    input_fits_file = args.fits_file

    # Perform data move to the science extension
    move_data(input_fits_file)

    print ('\n * Script  move_data2ext1.py  finished * \n')


