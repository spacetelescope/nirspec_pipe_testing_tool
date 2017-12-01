
"""
py.test module for unit testing the extract_1d step.
"""

import pytest
import os

from jwst.extract_1d.extract_1d_step import Extract1dStep
from .. import core_utils
from . import extract_1d_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "extract_1d"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    skip_runing_pipe_step = config.getboolean("tests_only", "_".join((step, "tests")))
    stp = Extract1dStep()
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    if not run_calwebb_spec2:
        if config.getboolean("steps", step):
            print ("*** Step "+step+" set to True")
            if os.path.isfile(step_input_file):
                if not skip_runing_pipe_step:
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

def test_s_extr1d_exists(output_hdul):
    assert extract_1d_utils.s_extr1d_exists(output_hdul), "The keyword S_EXTR1D was not added to the header --> Extract 1D step was not completed."

