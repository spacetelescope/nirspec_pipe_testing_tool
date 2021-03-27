from .. auxiliary_code.reffile_test import create_rfile_test

"""
This file contains the functions which will be used to test the extract_1d step of the JWST Calibration Pipeline.

"""

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Mar 2021 - Version 1.0: initial version completed


# VERIFICATION FUNCTIONS

def s_wavecor_exists(output_hdul):
    """
    This function checks that the keyword S_WAVCOR was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_WAVCOR" in output_hdul
    return result


wavecor_rfile_is_correct = create_rfile_test("R_WAVCOR", "wavecorr step")

