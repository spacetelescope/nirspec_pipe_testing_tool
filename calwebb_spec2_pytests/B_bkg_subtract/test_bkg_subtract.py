
"""
py.test module for unit testing the bkg_subtract step.
"""

import pytest
import os
import time
from astropy.io import fits
from jwst.background.background_step import BackgroundStep

from .. import core_utils
from . import bkg_subtract_utils



# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.1"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# March 2019 - Version 1.1: separated completion from numerical tests


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "bkg_subtract"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    # determine which steps are to be run, if not run in full
    run_pipe_step = config.getboolean("run_pipe_steps", step)
    # determine which tests are to be run
    bkg_subtract_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    #bkg_subtract_numerical_tests = config.getboolean("run_pytest", "_".join((step, "numerical", "tests")))
    #bkg_subtract_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [bkg_subtract_completion_tests]#, bkg_subtract_numerical_tests, bkg_subtract_validation_tests]

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # skip if BOTS data
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if core_utils.check_BOTS_true(inhdu):

        if run_calwebb_spec2:
            if os.path.isfile(step_output_file):
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            else:
                pytest.skip("Skipping "+step+" because the output file does not exist.")
            return hdul, step_output_file, step_input_file, run_pytests
        else:
            if os.path.isfile(step_input_file):
                print(" The input file ", step_input_file,"exists... will run step "+step)
                bkg_list = core_utils.getlist("additional_arguments", "bkg_list")
                existing_bgfiles = 0
                for bg_file in bkg_list:
                    if os.path.isfile(bg_file):
                        existing_bgfiles += 1
                if existing_bgfiles == 0:
                    print (" Need at least one background file to continue. Step will be skipped.")
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Skipping "+step+" because files listed on bkg_list in the configuration file do not exist.")
                else:
                    if run_pipe_step:
                        print ("*** Step "+step+" set to True")
                        stp = BackgroundStep()

                        # check that previous pipeline steps were run up to this point
                        core_utils.check_completed_steps(step, step_input_file)

                        # get the right configuration files to run the step
                        local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                        # start the timer to compute the step running time
                        start_time = time.time()
                        if local_pipe_cfg_path == "pipe_source_tree_code":
                            result = stp.call(step_input_file, bkg_list)
                        else:
                            result = stp.call(step_input_file, bkg_list, config_file=local_pipe_cfg_path+'/background.cfg')
                        if result is not None:
                            result.save(step_output_file)
                            # end the timer to compute the step running time
                            end_time = repr(time.time() - start_time)   # this is in seconds
                            print("Step "+step+" took "+end_time+" seconds to finish")
                            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                            step_completed = True
                        else:
                            hdul = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
                    else:
                        print("Skipping running pipeline step ", step)
                        # add the running time for this step
                        working_directory = config.get("calwebb_spec2_input_file", "working_directory")
                        # Get the detector used
                        det = fits.getval(step_input_file, "DETECTOR", 0)
                        end_time = core_utils.get_stp_run_time_from_screenfile(step, det, working_directory)
                        hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                        step_completed = True
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return hdul, step_output_file, step_input_file, run_pytests

            else:
                print (" The input file does not exist. Skipping step.")
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because the input file does not exist.")

    else:
        msg = "Skipping "+step+" because data is BOTS."
        print (msg)
        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
        pytest.skip(msg)


### FUNCTION FOR VALIDATION

# fixture to validate the background substract
#@pytest.fixture(scope="module")
#def check_if_subtract_is_zero(step_input_file):
    """
    This function uses a copy of the background input file and runs it through the step to test if the subtraction
    is performed correctly, i.e. if the result is zero.
    Args:
        step_input_file: string, name of the step input file

    Returns:
        result: float, the result from the subtraction step.
    """
#    bgfile_copy = copy.deepcopy(step_input_file)
#    stp = BackgroundStep()
#    result = stp.call(step_input_file, bgfile_copy)




# Unit tests

def test_s_bkdsub_exists(output_hdul):
    # want to run this pytest?
    # output_hdu[3] = bkg_subtract_completion_tests, bkg_subtract_numerical_tests, bkg_subtract_validation_tests
    run_pytests = output_hdul[3][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running completion pytest...\n")
        assert bkg_subtract_utils.s_bkdsub_exists(output_hdul[0]), "The keyword S_BKDSUB was not added to the header --> background step was not completed."

