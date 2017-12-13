
"""
py.test module for unit testing the pathloss step.
"""

import pytest
import os
import time
import subprocess

from jwst.pathloss.pathloss_step import PathLossStep
from .. import core_utils
from . import pathloss_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "pathloss"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    stp = PathLossStep()
    skip_runing_pipe_step = config.getboolean("tests_only", "_".join((step, "tests")))
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'
    if run_calwebb_spec2:
        # read the assign wcs fits file
        local_step_output_file = core_utils.read_completion_to_full_run_map("full_run_map.txt", step)
        hdul = core_utils.read_hdrfits(local_step_output_file, info=False, show_hdr=False)
        # move the output file into the working directory
        working_directory = config.get("calwebb_spec2_input_file", "working_directory")
        step_output_file = os.path.join(working_directory, local_step_output_file)
        print ("Step product was saved as: ", step_output_file)
        subprocess.run(["mv", local_step_output_file, step_output_file])
        return hdul
    else:
        if config.getboolean("steps", step):
            print ("*** Step "+step+" set to True")
            if os.path.isfile(step_input_file):
                if not skip_runing_pipe_step:
                    # start the timer to compute the step running time
                    start_time = time.time()
                    result = stp.call(step_input_file)
                    result.save(step_output_file)
                    # end the timer to compute calwebb_spec2 running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    print(" * calwebb_spec2 took "+end_time+" seconds to finish.")
                step_completed = True
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                return hdul
            else:
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because the input file does not exist.")
        else:
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+". Step set to False in configuration file.")



# Unit tests

def test_s_pthlos_exists(output_hdul):
    assert pathloss_utils.s_pthlos_exists(output_hdul), "The keyword S_PTHLOS was not added to the header --> Pathloss step was not completed."

def test_r_pthlos_exists(output_hdul):
    assert pathloss_utils.r_pthlos_exists(output_hdul), "The keyword R_PTHLOS was not added to the header --> Not sure what reference file was used."
