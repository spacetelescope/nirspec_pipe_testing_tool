from .. import core_utils

"""
This file contains the functions which will be used to test the wcs_assign step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


def create_completed_steps_txtfile(True_steps_suffix_map, step_input_file):
    """
    This function creates the completed steps along with the corresponding suffix of the output file name into a text file.
    Args:
        True_steps_suffix_map: string, full path of where the text file will be written into
        step_input_file: string, name of the input file for the pipeline step

    Returns:
        Nothing. A text file will be created in the pytests directory where all steps will be added
    """
    # name of the text file to collect the step name and suffix
    print ("Map created at: ", True_steps_suffix_map)
    line0 = "# {:<20}".format("Input file: "+step_input_file)
    line1 = "# {:<17} {:<20} {:<20} {:<20}".format("Step", "Added suffix", "Step complition", "Time to run [s]")
    with open(True_steps_suffix_map, "w+") as tf:
        tf.write(line0+"\n")
        tf.write(line1+"\n")


def create_map_from_full_run(full_run_map, step_input_file):
    """
    This function creates the map of fits file names of the intermediary step outputs when running the pipeline in
    a single complete run using the calwebb_spec2.cfg file.
    Args:
        full_run_map: string, full path of where the text file will be written into
        step_input_file: string, name of the input file for the pipeline step

    Returns:
        Nothing. A text file will be created in the pytests directory where all intermediary product files are mapped.
    """

    # list of the steps that produce output fits files (in order of completion)
    pipe_steps = ["assign_wcs", "bkg_subtract", "imprint_subtract", "msa_flagging", "extract_2d", "flat_field",
                  "pathloss", "barshadow", "photom", "resample_spec", "cube_build", "extract_1d"]

    # name of the text file to collect the step name and suffix
    print ("Map created at: ", full_run_map)
    line0 = "# {:<20}".format("Input file: "+step_input_file)
    line1 = "# {:<17} {:<20}".format("Step", "File name")
    with open(full_run_map, "w+") as tf:
        tf.write(line0+"\n")
        tf.write(line1+"\n")
        for stp in pipe_steps:
            line2 = "{:<20} {:<20}".format(stp, "".join((stp, ".fits")))
            tf.write(line2+"\n")




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


# REFERENCE FILES check

# from running calwebb_spec1

def rmask_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the R mask matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_MASK"
    # get the reference file used
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


def saturation_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the saturation step matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_SATURA"
    # get the reference file used
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


def superbias_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the superbias step matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_SUPERB"
    # get the reference file used
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


def linearity_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the linearity step matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_LINEAR"
    # get the reference file used
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


def dark_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the dark matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_DARK"
    # get the reference file used
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


def readnoise_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the read noise matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_READNO"
    # get the reference file used
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


def gain_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the gain matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_GAIN"
    # get the reference file used
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



# specific to the WCS step

def camera_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the camera model matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used for the camera model
    key = "R_CAMERA"
    if key in output_hdul:
        rfile = output_hdul[key]
    else:
        rfile = "".join((key, " not in header"))

    # get the list of reference files used so far in the file header
    ref_files_used_so_far = core_utils.get_reffile_used(output_hdul)

    # determine if the reference file for camera matches expected value
    if key in ref_files_used_so_far:
        expected_rfile = ref_files_used_so_far[key]
        result = rfile in expected_rfile

    return result


def colimator_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the colimator model matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_COLLIM"
    # get the reference file used
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


def disperser_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the disperser model matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    key = "R_DISPER"
    # get the reference file used
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


def fore_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the transform through the FORE optics matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_FORE"
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


def fpa_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the transform in the FPA plane matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_FPA"
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


def ifufore_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the transform from the MSA plane to the plane of the IFU
     slicer matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_IFUFOR"
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


def ifupost_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the transform from slicer plane to the MSA plane
    matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_IFUPOS"
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


def ifuslicer_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the metrology of the IFU slicer matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_IFUSLI"
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


def msa_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the metrology of the MSA plane matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_MSA"
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


def ote_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the transform through the optical telescope element
    matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_OTE"
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


def wavran_rfile_is_correct(output_hdul):
    """
    This function determines if the reference file for the typical wavelength ranges matches the expected one.
    Args:
        output_hdul: main header of the WCS step output

    Returns:
        result: boolean, true if the reference file matches expected value
    """
    # get the reference file used
    key = "R_WAVRAN"
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



