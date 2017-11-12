"""
This file contains the functions which will be used to test the flat_field step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


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


