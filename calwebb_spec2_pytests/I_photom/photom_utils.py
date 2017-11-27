"""
This file contains the functions which will be used to test the photom step
of the JWST Calibration Pipeline.

"""


### VERIFICATION FUNCTIONS

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



