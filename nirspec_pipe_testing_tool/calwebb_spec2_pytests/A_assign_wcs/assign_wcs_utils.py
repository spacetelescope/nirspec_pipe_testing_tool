import subprocess
import configparser

from .. auxiliary_code.reffile_test import create_rfile_test

"""
This file contains the functions which will be used to test the wcs_assign step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


# HEADER
__author__ = "M. Pena-Guerrero & G. Kanarek"
__version__ = "2.2"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check
# Jul 2018 - Version 2.1: Maria removed the function to specifically create the
#                         full_run_map file, now this txt file
#                         has the same format and information as the True_steps_map.txt file
# Apr 2023 - Version 2.2: Code clean up


# VERIFICATION FUNCTIONS

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


# REFERENCE FILES check

# from running calwebb_spec1

rmask_rfile_is_correct = create_rfile_test("R_MASK", "R mask")
saturation_rfile_is_correct = create_rfile_test("R_SATURA", "saturation step")
superbias_rfile_is_correct = create_rfile_test("R_SUPERB", "superbias step")
linearity_rfile_is_correct = create_rfile_test("R_LINEAR", "linearity step")
dark_rfile_is_correct = create_rfile_test("R_DARK", "dark")
readnoise_rfile_is_correct = create_rfile_test("R_READNO", "read noise")
gain_rfile_is_correct = create_rfile_test("R_GAIN", "gain")

# specific to the WCS step

camera_rfile_is_correct = create_rfile_test("R_CAMERA", "camera model")
colimator_rfile_is_correct = create_rfile_test("R_COLLIM", "collimator model")
disperser_rfile_is_correct = create_rfile_test("R_DISPER", "disperser model")
fore_rfile_is_correct = create_rfile_test("R_FORE", 
                                           "transform through the FORE optics")
fpa_rfile_is_correct = create_rfile_test("R_FPA", "FPA plane")
ifufore_rfile_is_correct = create_rfile_test("R_IFUFOR", 
                                              "transform from the MSA plane to the plane of the IFU slicer")
ifupost_rfile_is_correct = create_rfile_test("R_IFUPOS", 
                                              "transform from the slicer plane to the MSA plane")
ifuslicer_rfile_is_correct = create_rfile_test("R_IFUSLI", "metrology of the IFU slicer")
msa_rfile_is_correct = create_rfile_test("R_MSA", "metrology of the MSA plane")
ote_rfile_is_correct = create_rfile_test("R_OTE", 
                                          "transform through the optical telescope element")
wavran_rfile_is_correct = create_rfile_test("R_WAVRAN", "typical wavelength ranges")

