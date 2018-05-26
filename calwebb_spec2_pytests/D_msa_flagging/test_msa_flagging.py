
"""
py.test module for unit testing the msa_flagging step.
"""

import pytest
import os
import time
import subprocess
from jwst.msaflagopen.msaflagopen_step import MSAFlagOpenStep

from .. import core_utils
from . import msa_flagging_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "msa_flagging"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    skip_runing_pipe_step = config.getboolean("tests_only", "_".join((step, "tests")))
    # Only run step if data is MOS
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    end_time = '0.0'
    print ("Is data FS or BOTS?", core_utils.check_FS_true(inhdu), core_utils.check_BOTS_true(inhdu))
    if not core_utils.check_FS_true(inhdu) and not core_utils.check_BOTS_true(inhdu):
        stp = MSAFlagOpenStep()
        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False
        if run_calwebb_spec2:
            # read the assign wcs fits file
            input_file = config.get("calwebb_spec2_input_file", "input_file")
            local_step_output_file = input_file.replace(".fits", "_msa_flagging.fits")
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
                        # get the right configuration files to run the step
                        #local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                        # start the timer to compute the step running time
                        start_time = time.time()
                        #if local_pipe_cfg_path == "pipe_source_tree_code":
                        result = stp.call(step_input_file)
                        #else:
                        #    result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/NOCONFIGFI.cfg')
                        # end the timer to compute the step running time
                        end_time = repr(time.time() - start_time)   # this is in seconds
                        print("Step "+step+" took "+end_time+" seconds to finish")
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
    else:
        pytest.skip("Skipping "+step+" because data is FS or BOTS.")



# Unit tests

def test_msa_failed_open_exists(output_hdul):
    assert msa_flagging_utils.msa_failed_open_exists(output_hdul), "The keyword S_MSAFLG was not added to the header --> msa_flagging step was not completed."
