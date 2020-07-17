import subprocess
import configparser

from .. auxiliary_code.reffile_test import create_rfile_test

"""
This file contains the functions which will be used to test the wcs_assign step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


# HEADER
__author__ = "M. A. Pena-Guerrero & Gray Kanarek"
__version__ = "2.1"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check
# Jul 2018 - Version 2.1: Maria removed the function to specifically create the full_run_map file, now this txt file
#                         has the same format and information as the True_steps_map.txt file


def create_completed_steps_txtfile(txt_suffix_map, step_input_file):
    """
    This function creates the completed steps along with the corresponding suffix of the output file name into a text file.
    Args:
        txt_suffix_map: string, full path of where the text file will be written into
        step_input_file: string, name of the input file for the pipeline step

    Returns:
        Nothing. A text file will be created in the pytests directory where all steps will be added
    """
    # name of the text file to collect the step name and suffix
    # print ("Map created at: ", txt_suffix_map)
    line0 = "# {:<20}".format("Input file: "+step_input_file)
    line1 = "# {:<17} {:<20} {:<20} {:<20}".format("Step", "Added suffix", "Step complition", "Time to run [s]")
    print(line1)
    with open(txt_suffix_map, "w+") as tf:
        tf.write(line0+"\n")
        tf.write(line1+"\n")


def print_time2file(txt_name, end_time, string2print):
    """
    This function is only used in the case of running the pipeline in full. It changes the total running time to the
    appropriate units and returns a string of the right spacing.
    Args:
        txt_name: string, name of the suffix map file
        end_time: string, time it took to run
        string2print: string, phrase to be printed in the file before the time

    Returns:
        nothing
    """
    if float(end_time) > 60.0:
        total_time_min = repr(round(float(end_time)/60.0, 1))   # in minutes
        total_time = total_time_min+'min'
        if float(total_time_min) > 60.0:
            total_time_hr = repr(round(float(end_time)/60.0, 1))   # in hours
            total_time = total_time_hr+'hr'
        end_time = end_time+'  ='+total_time
    line2write = "{:<20} {:<20} {:<20} {:<20}".format('', '', string2print, end_time)
    print(line2write)
    with open(txt_name, "a") as tf:
        tf.write(line2write+"\n")


def set_pipe_log(calwebb_spec2_cfg, detector):
    # copy the configuration file to create the pipeline log
    stpipelogcfg = calwebb_spec2_cfg.replace("calwebb_spec2.cfg", "stpipe-log.cfg")
    subprocess.run(["cp", stpipelogcfg, "."])
    # make sure that the handler has the correct name for creating the pipeline log
    config = configparser.ConfigParser()
    config.read([stpipelogcfg])
    handler = config.get("*", "handler")
    # correct the right pipeline log name if needed
    pipelog = "pipeline.log"
    if pipelog in handler:
        # re-write the .cfg file
        pipelog = "pipeline_" + detector + ".log"
        config = configparser.ConfigParser()
        config.add_section("*")
        config.set("*", "handler", "file:" + pipelog)
        config.set("*", "level", "INFO")
        config.write(open("stpipe-log.cfg", "w"))
        print("Pipeline log name set to:  ", pipelog)
    return pipelog


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

