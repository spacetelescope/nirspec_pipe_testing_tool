
"""
py.test module for unit testing the barshadow step.
"""

import pytest
import os

from jwst.barshadow.barshadow_step import BarShadowStep
from .. import core_utils
from . import barshadow_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "barshadow"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    stp = BarShadowStep()
    # Only run step if data is MOS
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if core_utils.check_MOS_true(inhdu):
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
    else:
        pytest.skip("Skipping "+step+" because data is not MOS.")



# Unit tests

def test_s_barsha_exists(output_hdul):
    assert barshadow_utils.s_barsha_exists(output_hdul), "The keyword S_BARSHA was not added to the header --> Barshadow step was not completed."

