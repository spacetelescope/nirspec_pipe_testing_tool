
"""
py.test module for unit testing the bkg_subtract step.
"""

import pytest
import os
import time
import copy
import subprocess
from glob import glob
from astropy.io import fits
from jwst.background.background_step import BackgroundStep

from nirspec_pipe_testing_tool import core_utils
from . import bkg_subtract_utils


# HEADER
__author__ = "M. Pena-Guerrero"
__version__ = "1.4"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: separated completion from numerical tests
# Apr 2019 - Version 1.2: implemented nptt_log capability
# Jun 2020 - Version 1.3: implemented numerical test
# Apr 2023 - Version 1.4: Cleaned-up code


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
def output_vars(set_inandout_filenames, config):
    # determine if the pipeline is to be run in full, per steps, or skipped
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    if run_calwebb_spec2 == "skip":
        print('\n * NPTT finished processing run_calwebb_spec2 is set to skip. \n')
        pytest.exit("Skipping pipeline run and tests for spec2, run_calwebb_spec2 is set to skip in NPTT_config file.")
    elif "T" in run_calwebb_spec2:
        run_calwebb_spec2 = True
    else:
        run_calwebb_spec2 = False

    # get the general info
    set_inandout_filenames_info = core_utils.read_info4output_vars(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, outstep_file_suffix = set_inandout_filenames_info

    # determine which steps are to be run, if not run in full
    run_pipe_step = config.getboolean("run_spec2_steps", step)
    # determine which tests are to be run
    bkg_subtract_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    bkg_subtract_numerical_tests = config.getboolean("run_pytest", "_".join((step, "numerical", "tests")))
    #bkg_subtract_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [bkg_subtract_completion_tests, bkg_subtract_numerical_tests]#, bkg_subtract_validation_tests]

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # Get the main header and the detector used
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(output_directory, initial_input_file)
    if os.path.isfile(initial_input_file):
        inhdr = fits.getheader(step_input_file)
        detector = inhdr["DETECTOR"]
    else:
        msg = "Skipping "+step+" because the initial input file given in NPTT_config.cfg does not exist."
        pytest.skip(msg)

    # Get the logfile instance for NPTT created in the run.py script
    nptt_log = os.path.join(output_directory, 'NPTT_calspec2_' + detector + '.log')
    nptt_log = core_utils.mk_nptt_log(nptt_log, reset=False)

    # skip if BOTS data, else perform step
    if not core_utils.check_BOTS_true(inhdr):
        if run_calwebb_spec2:
            if os.path.isfile(step_output_file):
                outhdr = fits.getheader(step_output_file)
            else:
                msg = "Skipping "+step+" because the output file does not exist."
                nptt_log.info(msg)
                pytest.skip(msg)
            return outhdr, step_output_file, step_input_file, run_pytests, nptt_log
        else:
            if run_pipe_step:
                if os.path.isfile(step_input_file):
                    # Create the pipeline step log
                    stp_pipelog = "calspec2_" + step + "_" + detector + ".log"
                    core_utils.mk_stpipe_log_cfg(output_dir, stp_pipelog)
                    print("Pipeline step screen output will be logged in file: ", stp_pipelog)

                    msg = " The input file "+step_input_file+" exists... will run step "+step
                    print(msg)
                    nptt_log.info(msg)
                    bkg_list = core_utils.getlist("additional_arguments", "bkg_list")
                    existing_bgfiles = 0
                    for bg_file in bkg_list:
                        if os.path.isfile(bg_file):
                            existing_bgfiles += 1
                    if existing_bgfiles == 0:
                        msg = " Need at least one background file to continue. Step will be skipped."
                        print(msg)
                        nptt_log.info(msg)
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                        pytest.skip("Skipping "+step+" because files listed on bkg_list in the configuration "
                                                     "file do not exist.")
                    else:
                        # continue since the file(s) exist
                        msg = "*** Step "+step+" set to True"
                        print(msg)
                        nptt_log.info(msg)
                        stp = BackgroundStep()

                        # check that previous pipeline steps were run up to this point
                        core_utils.check_completed_steps(step, step_input_file)

                        # get the right configuration files to run the step
                        local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                        # start the timer to compute the step running time
                        start_time = time.time()
                        msg = "Running pipeline step " + step
                        print(msg)
                        nptt_log.info(msg)
                        if local_pipe_cfg_path == "pipe_source_tree_code":
                            result = stp.call(step_input_file, bkg_list)
                        else:
                            result = stp.call(step_input_file, bkg_list, config_file=local_pipe_cfg_path+'/background.cfg')

                        if result is not None:
                            result.save(step_output_file)
                            # end the timer to compute the step running time
                            end_time = repr(time.time() - start_time)   # this is in seconds
                            msg = "Step "+step+" took "+end_time+" seconds to finish"
                            print(msg)
                            nptt_log.info(msg)
                            outhdr = fits.getheader(step_output_file)
                            step_completed = True

                            # add the running time for this step
                            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                            return outhdr, step_output_file, step_input_file, run_pytests, nptt_log

                else:
                    msg = "The input file does not exist. Skipping step."
                    print(msg)
                    nptt_log.info(msg)
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Skipping "+step+" because the input file does not exist.")

            else:
                msg = "Skipping running pipeline step "+step
                print(msg)
                nptt_log.info(msg)
                end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, output_directory)
                if os.path.isfile(step_output_file):
                    outhdr = fits.getheader(step_output_file)
                    step_completed = True
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return outhdr, step_output_file, step_input_file, run_pytests, nptt_log
                else:
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")

    else:
        msg = "Skipping "+step+" because data is BOTS."
        print(msg)
        nptt_log.info(msg)
        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
        pytest.skip(msg)


# FUNCTION FOR VALIDATION

# fixture to validate the subtraction works fine: re-run the step with the same file as msa_imprint file
@pytest.fixture(scope="module")
def check_output_is_zero(output_vars):
    """
    This test is a simple subtraction of the input file minus a copy of the input file (instead of the actual
    backfround file(s). The output is expected to be zero.
    Args:
        output_vars: list
    Returns:
        result: boolean
    """
    # output_vars = hdr, step_output_file, step_input_file, run_pytests, run_pipe_step, nptt_log
    step_input_file = output_vars[2]
    step_output_file = output_vars[1]
    # set specifics for the test
    bgfile_copy = copy.deepcopy(step_input_file)
    result_to_check = step_output_file.replace(".fits", "_zerotest.fits")
    # run the step with the specifics
    stp = BackgroundStep()
    res = stp.call(step_input_file, bgfile_copy)
    res.save(result_to_check)
    # check that the end product of image - image is zero
    c = fits.getdata(result_to_check)
    subtraction = sum(c.flatten())
    result = False
    if subtraction == 0.0:
        result = True
    # erase test output file
    subprocess.run(["rm", result_to_check])
    return result


# Unit tests

def test_s_bkdsub_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_hdu[3] = bkg_subtract_completion_tests, bkg_subtract_numerical_tests, bkg_subtract_validation_tests
    run_pytests = output_vars[3][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running completion pytest..."
        print(msg)
        nptt_log.info(msg)
        assert bkg_subtract_utils.s_bkdsub_exists(output_vars[0]), "The keyword S_BKDSUB was not added to the header " \
                                                                   "--> background step was not completed."


def test_check_output_is_zero(output_vars, request):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_hdu[3] = bkg_subtract_completion_tests, bkg_subtract_numerical_tests, bkg_subtract_validation_tests
    run_pytests = output_vars[3][1]
    if not run_pytests:
        msg = "Skipping pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running numerical accuracy pytest..."
        print(msg)
        nptt_log.info(msg)
        assert request.getfixturevalue('check_output_is_zero'), "Subtraction result is not equal to zero."

