from __future__ import print_function, division
"""
py.test module for unit testing the bkg_subtract step.
"""

import pytest
import os
from jwst.pipeline import Spec2Pipeline
from jwst.background.background_step import BackgroundStep

from .. import core_utils
from . import bkg_subtract_utils


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
    step = "bkg_subtract"
    step_dict = dict(config.items("steps"))
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    True_steps_suffix_map = config.get("calwebb_spec2_input_file", "True_steps_suffix_map")
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, step_dict)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    step = set_inandout_filenames[0]
    step_input_filename = set_inandout_filenames[1]
    output_file = set_inandout_filenames[2]
    outstep_file_suffix = set_inandout_filenames[4]
    True_steps_suffix_map = set_inandout_filenames[5]
    txt_name = os.path.join(working_directory, True_steps_suffix_map)
    step_input_file = os.path.join(working_directory, step_input_filename)
    step_output_file = os.path.join(working_directory, output_file)
    stp = BackgroundStep()
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    if run_calwebb_spec2:
        print ("*** Will run calwebb_spec2... ")
        calwebb_spec2_cfg = config.get("run_calwebb_spec2_in_full", "calwebb_spec2_cfg")
        final_output_name = step_input_file.replace(".fits", "_calwebb_spec2.fits")
        result_level2B = Spec2Pipeline.call(step_input_file, config_file=calwebb_spec2_cfg)
        result_level2B.save(final_output_name)
        hdul = core_utils.read_hdrfits(final_output_name, info=True, show_hdr=True)
        return hdul
    else:
        if config.getboolean("steps", step):
            print ("*** Step "+step+" set to True")
            if os.path.isfile(step_input_file):
                print(" The input file ", step_input_filename,"exists... will run step "+step)
                bkg_list = core_utils.getlist("additional_arguments", "bkg_list")
                existing_bgfiles = 0
                for bg_file in bkg_list:
                    if os.path.isfile(bg_file):
                        existing_bgfiles += 1
                if existing_bgfiles == 0:
                    print (" Need at least one background file to continue. Step will be skipped.")
                    create_completed_steps_txtfile(txt_name, step_input_filename, step,
                                                   outstep_file_suffix, step_completed)
                    pytest.skip("Skipping "+step+" because files listed on bkg_list in the configuration file do not exist.")
                else:
                    result = stp.call(step_input_file, bkg_list)
                    if result is not None:
                        result.save(step_output_file)
                        hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                        step_completed = True
                    else:
                        hdul = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
                    create_completed_steps_txtfile(txt_name, step_input_filename, step,
                                                   outstep_file_suffix, step_completed)
                    return hdul
            else:
                print (" The input file does not exist. Skipping step.")
                create_completed_steps_txtfile(txt_name, step_input_filename, step,
                                               outstep_file_suffix, step_completed)
                pytest.skip("Skipping "+step+" because the input file does not exist.")
        else:
            create_completed_steps_txtfile(txt_name, step_input_filename, step,
                                           outstep_file_suffix, step_completed)
            pytest.skip("Skipping "+step+". Step set to False in configuration file.")



# Unit tests

def test_s_bkdsub_exists(output_hdul):
    assert bkg_subtract_utils.s_bkdsub_exists(output_hdul)

