# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Jan 2019 - Version 1.1: Maria modified and added Gray's code for validation tests
# Apr 2019 - Version 1.2: implemented logging capability
# Jul 2020 - Version 1.3: moved the corners comparison to benchmark data to an independent scipt


"""
This file contains the functions which will be used to test the extract_2d step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


# VERIFICATION FUNCTIONS

def s_ext2d_exists(output_hdul):
    """
    This function checks that the keyword S_EXTR2D was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_EXTR2D" in output_hdul
    return result


