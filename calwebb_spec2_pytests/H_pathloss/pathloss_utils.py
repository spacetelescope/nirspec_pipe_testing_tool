from .. import core_utils
from .. auxiliary_code.reffile_test import create_rfile_test

"""
This file contains the functions which will be used to test the pathloss step
of the JWST Calibration Pipeline.

"""



# HEADER
__author__ = "M. A. Pena-Guerrero & Gray Kanarek"
__version__ = "2.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check


### VERIFICATION FUNCTIONS

def s_pthlos_exists(output_hdul):
    """
    This function checks that the keyword S_PTHLOS was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_PTHLOS" in output_hdul
    return result


def r_pthlos_exists(output_hdul):
    """
    This function checks that the keyword R_PTHLOS was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "R_PTHLOS" in output_hdul
    if result:
        print (" Reference file used for pathloss step: ", output_hdul["R_PTHLOS"])
    return result


pthlos_rfile_is_correct = create_rfile_test("R_PTHLOS", "pathloss")
