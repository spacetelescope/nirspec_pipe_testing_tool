from .. import core_utils

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


def pthlos_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the photom matches the expected one.
    Args:
        output_hdul: main header of the photom step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_PHOTOM"
    if key in output_hdul:
        rfile = output_hdul[key]
    else:
        rfile = "".join((key, " not in header"))

    # get the list of reference files used so far in the file header
    ref_files_used_so_far = core_utils.get_reffile_used(output_hdul)

    # determine if the reference file matches expected one
    if key in ref_files_used_so_far:
        expected_rfile = ref_files_used_so_far[key]
        result = rfile in expected_rfile

    return result

