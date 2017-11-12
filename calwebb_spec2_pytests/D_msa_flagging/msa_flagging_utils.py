"""
This file contains the functions which will be used to test the msa_flagging step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


### VERIFICATION FUNCTIONS

def msa_failed_open_exists(output_hdul):
    """
    This function checks that the keyword MSA_FAILED_OPEN was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "MSA_FAILED_OPEN" in output_hdul
    return result

