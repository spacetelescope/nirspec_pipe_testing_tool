"""
This file contains the functions which will be used to test the source_type step
of the JWST Calibration Pipeline.

"""


### VERIFICATION FUNCTIONS

def s_srctype_exists(output_hdul):
    """
    This function checks that the keyword SRCTYPE was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "SRCTYPE" in output_hdul
    return result


