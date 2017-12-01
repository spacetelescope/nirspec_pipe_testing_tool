
"""
py.test module for unit testing the assign_wcs step.
"""

import pytest
import os
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep

from .. import core_utils
from . import assign_wcs_utils


def create_completed_steps_txtfile(True_steps_suffix_map, step_input_file, step, outstep_file_suffix, step_completed):
    """
    This function adds the completed steps along with the corresponding suffix of the output file name into a text file.
    Args:
        True_steps_suffix_map: string, full path of where the text file will be written into
        step_input_file: string, name of the input file for the pipeline step
        step: string, pipeline step just ran
        outstep_file_suffix: string, suffix added right before .fits to the input file
        step_completed: boolean, True if the step was completed and False if it was skiped

    Returns:
        nothing
    """
    # name of the text file to collect the step name and suffix
    print ("Map created at: ", True_steps_suffix_map)
    line0 = "# {:<20}".format("Input file: "+step_input_file)
    line1 = "# {:<17} {:<20} {:<20}".format("Step", "Added suffix", "Step complition")
    line2write = "{:<20} {:<20} {:<20}".format(step, outstep_file_suffix, str(step_completed))
    with open(True_steps_suffix_map, "w+") as tf:
        tf.write(line0+"\n")
        tf.write(line1+"\n")
        tf.write(line2write+"\n")



# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(config):
    step = "assign_wcs"
    step_dict = dict(config.items("steps"))
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    True_steps_suffix_map = config.get("calwebb_spec2_input_file", "True_steps_suffix_map")
    pytests_directory = os.getcwd()
    True_steps_suffix_map = os.path.join(pytests_directory, True_steps_suffix_map)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, step_dict)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    stp = AssignWcsStep()
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    skip_runing_pipe_step = config.getboolean("tests_only", "_".join((step, "tests")))
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    if run_calwebb_spec2:
        print ("*** Will run calwebb_spec2... ")
        calwebb_spec2_cfg = config.get("run_calwebb_spec2_in_full", "calwebb_spec2_cfg")
        final_output_name = step_input_file.replace(".fits", "_calwebb_spec2.fits")
        #result_level2B = Spec2Pipeline.call(step_input_file, config_file=calwebb_spec2_cfg)
        #result_level2B.save(final_output_name)
        hdul = core_utils.read_hdrfits(final_output_name, info=True, show_hdr=True)
        return hdul
    else:
        if config.getboolean("steps", step):
            print ("*** Step "+step+" set to True")
            if os.path.isfile(step_input_file):
                if not skip_runing_pipe_step:
                    result = stp.call(step_input_file)
                    result.save(step_output_file)
                step_completed = True
                create_completed_steps_txtfile(txt_name, step_input_file, step,
                               outstep_file_suffix, step_completed)
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                return hdul
            else:
                print("Skipping step. Intput file "+step_input_file+" does not exit.")
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                pytest.skip("Skipping "+step+" because the input file does not exist.")
        else:
            create_completed_steps_txtfile(txt_name, step_input_file, step,
                                           outstep_file_suffix, step_completed)
            pytest.skip("Skipping "+step+". Step set to False in configuration file.")



# Unit tests

def test_wavstart_exists(output_hdul):
    assert assign_wcs_utils.wavstart_exists(output_hdul), "The keyword WAVSTART was not added to the header."

def test_sporder_exists(output_hdul):
    assert assign_wcs_utils.sporder_exists(output_hdul), "The keyword SPORDER was not added to the header."

def test_s_wcs_exists(output_hdul):
    assert assign_wcs_utils.s_wcs_exists(output_hdul), "The keyword S_WCS was not added to the header --> extract_2d step was not completed."

