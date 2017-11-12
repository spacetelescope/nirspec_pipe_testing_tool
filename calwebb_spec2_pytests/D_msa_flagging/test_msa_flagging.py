from __future__ import print_function, division
"""
py.test module for unit testing the msa_flagging step.
"""

import pytest
import os
from jwst.msaflagopen.msaflagopen_step import MSAFlagOpenStep

from .. import core_utils
from . import msa_flagging_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "msa_flagging"
    step_dict = dict(config.items("steps"))
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
    True_steps_suffix_map = config.get("calwebb_spec2_input_file", "True_steps_suffix_map")
    txt_name = os.path.join(working_directory, True_steps_suffix_map)
    steps_list, suffix_list, completion_list = core_utils.read_True_steps_suffix_map(txt_name)
    step_input_filename = core_utils.get_correct_input_step_filename(initial_input_file, steps_list,
                                                                     suffix_list, completion_list)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, step_dict)
    in_file_suffix, out_file_suffix, _, _ = suffix_and_filenames
    step_output_filename = step_input_filename.replace(".fits", out_file_suffix+".fits")
    print ("step_input_filename = ", step_input_filename)
    print ("step_output_filename = ", step_output_filename)
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
    stp = MSAFlagOpenStep()
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    if not run_calwebb_spec2:
        if config.getboolean("steps", step):
            print ("*** Step "+step+" set to True")
            if os.path.isfile(step_input_file):
                result = stp.call(step_input_file)
                result.save(step_output_file)
                step_completed = True
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                return hdul
            else:
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                pytest.skip("Skipping "+step+" because the input file does not exist.")
        else:
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
            pytest.skip("Skipping "+step+". Step set to False in configuration file.")



# Unit tests

def test_msa_failed_open_exists(output_hdul):
    assert msa_flagging_utils.msa_failed_open_exists(output_hdul)
