from .. auxiliary_code.reffile_test import create_rfile_test

"""
This file contains the functions which will be used to test the flat_field step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""



# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed


### VERIFICATION FUNCTIONS

def s_flat_exists(output_hdul):
    """
    This function checks that the keyword S_FLAT was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_FLAT" in output_hdul
    return result


fflat_rfile_is_correct = create_rfile_test("R_FFLAT", "F-Flat field")
sflat_rfile_is_correct = create_rfile_test("R_SFLAT", "S-Flat field")
dflat_rfile_is_correct = create_rfile_test("R_DFLAT", "D-Flat field")

