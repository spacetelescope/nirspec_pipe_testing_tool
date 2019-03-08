
"""
py.test module for unit testing the barshadow step.
"""

import os
import time
import pytest
from astropy.io import fits
from jwst.barshadow.barshadow_step import BarShadowStep

from . import barshadow_utils
from .. import core_utils
from .. auxiliary_code import change_filter_opaque2science



# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.1"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: separated completion from other tests


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
    run_pipe_step = config.getboolean("run_pipe_steps", step)
    # determine which tests are to be run
    barshadow_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    #barshadow_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    #barshadow_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [barshadow_completion_tests]#, barshadow_reffile_tests, barshadow_validation_tests]

    end_time = '0.0'

    # Only run step if data is MOS
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    print("core_utils.check_MOS_true(inhdu)=", core_utils.check_MOS_true(inhdu))
    if core_utils.check_MOS_true(inhdu):
        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False

        # check if the filter is to be changed
        change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
        if change_filter_opaque:
            is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file, step=step)
            if is_filter_opaque:
                print ("With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Barshadow pytest now set to Skip.")
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

        if run_calwebb_spec2:
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            return hdul, step_output_file, run_pytests
        else:
            if os.path.isfile(step_input_file):
                if run_pipe_step:
                    print ("*** Step "+step+" set to True")
                    stp = BarShadowStep()

                    # check that previous pipeline steps were run up to this point
                    core_utils.check_completed_steps(step, step_input_file)

                    # get the right configuration files to run the step
                    #local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                    # start the timer to compute the step running time
                    start_time = time.time()
                    #if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                    #else:
                    #    result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/NOCONFIGFILE.cfg')
                    result.save(step_output_file)
                    # end the timer to compute calwebb_spec2 running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    print(" * calwebb_spec2 took "+end_time+" seconds to finish.")
                else:
                    print("Skipping running pipeline step ", step)
                    # add the running time for this step
                    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
                    # Get the detector used
                    det = fits.getval(step_input_file, "DETECTOR", 0)
                    end_time = core_utils.get_stp_run_time_from_screenfile(step, det, working_directory)
                step_completed = True
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                return hdul, step_output_file, run_pytests

            else:
                print (" The input file does not exist. Skipping step.")
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because the input file does not exist.")

    else:
        pytest.skip("Skipping "+step+" because data is not MOS.")



# Unit tests

def test_s_barsha_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[2] = barshadow_completion_tests, barshadow_reffile_tests, barshadow_validation_tests
    run_pytests = output_hdul[2][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running completion pytest...\n")
        assert barshadow_utils.s_barsha_exists(output_hdul[0]), "The keyword S_BARSHA was not added to the header --> Barshadow step was not completed."

