from .. import core_utils

"""
This file contains the functions which will be used to test the pathloss step
of the JWST Calibration Pipeline.

"""


### VERIFICATION FUNCTIONS

def s_pthlos_exists(output_hdul):
    """
    This function checks that the keyword S_PTHLOS was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_PTHLOS" in output_hdul
    return result


def r_pthlos_exists(output_hdul):
    """
    This function checks that the keyword R_PTHLOS was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "R_PTHLOS" in output_hdul
    if result:
        print (" Reference file used for pathloss step: ", output_hdul["R_PTHLOS"])
    return result


def pthlos_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the pathloss matches the expected one.
    Args:
        output_hdul: main header of the pathloss step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_PTHLOS"
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

