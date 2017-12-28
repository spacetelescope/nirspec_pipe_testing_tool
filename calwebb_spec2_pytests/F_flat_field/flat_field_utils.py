from .. import core_utils

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


def fflat_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the F-Flat field matches the expected one.
    Args:
        output_hdul: main header of the flat field step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_FFLAT"
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


def sflat_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the S-Flat field matches the expected one.
    Args:
        output_hdul: main header of the flat field step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_SFLAT"
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


def dflat_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the D-Flat field matches the expected one.
    Args:
        output_hdul: main header of the flat field step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_DFLAT"
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

