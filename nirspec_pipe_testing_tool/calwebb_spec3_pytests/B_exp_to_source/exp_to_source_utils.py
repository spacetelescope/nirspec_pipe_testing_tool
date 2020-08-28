from nirspec_pipe_testing_tool.calwebb_spec2_pytests.auxiliary_code.reffile_test import create_rfile_test

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Jul 2020 - Version 1.0: initial version completed


# VERIFICATION FUNCTIONS

# S_OUTLIR, R_RESAMP, S_IFUCUB, S_EXTR1D

def masterbg_exists(output_hdul):
    """
    This function checks that the keyword MASTERBG was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "MASTERBG" in output_hdul
    return result


# REFERENCE FILES checks

#camera_rfile_is_correct = create_rfile_test("R_CAMERA", "camera model")

