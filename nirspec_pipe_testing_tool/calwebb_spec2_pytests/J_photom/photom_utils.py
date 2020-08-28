from astropy.io import fits

from .. auxiliary_code.reffile_test import create_rfile_test

"""
This file contains the functions which will be used to test the photom step
of the JWST Calibration Pipeline.

"""


# HEADER
__author__ = "M. A. Pena-Guerrero & Gray Kanarek"
__version__ = "2.1"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check
# Jun 2020 - Version 2.1: implemented change of units test


# VERIFICATION FUNCTIONS

def s_photom_exists(output_hdul):
    """
    This function checks that the keyword S_PHOTOM was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_PHOTOM" in output_hdul
    return result


def units_logic(output_hdul):
    """
    Check that for point sources the units are flux density (MJy) and for extended sources units are surface
    brightness (MJy/sr).
    :param output_hdul: list
    :return:
        result: boolean, true if the keyword was indeed added
        msg: string, output message
    """
    # get the units keyword from the science extension of the photom step output file
    step_output_file = output_hdul[1]
    srctype = fits.getval(step_output_file, 'SRCTYPE', 'SCI')
    units = fits.getval(step_output_file, 'BUNIT', 'SCI')
    result = False

    if srctype == 'POINT':
        if units == 'MJy':
            msg = "Photom units keyword in SCI extension are flux density, as expected for POINT source."
            result = True
        else:
            msg = "Photom units keyword in SCI extension are NOT flux density, as expected for POINT source."

    if srctype == 'EXTENDED':
        if units == 'MJy/sr':
            msg = "Photom units keyword in SCI extension are surface brightness, as expected for EXTENDED source."
            result = True
        else:
            msg = "Photom units keyword in SCI extension are NOT surface brightness, as expected for EXTENDED source."

    return result, msg


# Reference file checks
pthlos_rfile_is_correct = create_rfile_test("R_PTHLOS", "pathloss")

