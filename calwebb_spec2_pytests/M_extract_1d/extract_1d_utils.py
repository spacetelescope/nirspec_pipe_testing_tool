from .. import core_utils
from .. auxiliary_code.reffile_test import create_rfile_test

"""
This file contains the functions which will be used to test the extract_1d step
of the JWST Calibration Pipeline.

"""



# HEADER
__author__ = "M. A. Pena-Guerrero & Gray Kanarek"
__version__ = "2.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check


### VERIFICATION FUNCTIONS

def s_extr1d_exists(output_hdul):
    """
    This function checks that the keyword S_EXTR1D was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_EXTR1D" in output_hdul
    return result


extract1d_rfile_is_correct = create_rfile_test("R_EXTR1D", "extract 1D step")

