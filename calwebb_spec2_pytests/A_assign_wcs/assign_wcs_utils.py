"""
This file contains the functions which will be used to test the wcs_assign step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


### VERIFICATION FUNCTIONS

def wavstart_exists(output_hdul):
    """
    This function checks that the keyword WAVSTART was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "WAVSTART" in output_hdul
    return result


def wavend_exists(output_hdul):
    """
    This function checks that the keyword WAVEND was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "WAVEND" in output_hdul
    return result


def sporder_exists(output_hdul):
    """
    This function checks that the keyword SPORDER was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "SPORDER" in output_hdul
    return result


def s_wcs_exists(output_hdul):
    """
    This function checks that the keyword S_WCS was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_WCS" in output_hdul
    return result



