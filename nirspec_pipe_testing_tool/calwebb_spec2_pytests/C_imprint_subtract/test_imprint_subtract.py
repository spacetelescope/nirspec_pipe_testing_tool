
"""
py.test module for unit testing the imprint_subtract step.
"""

import pytest
import os
import time
import copy
import subprocess
from glob import glob
from astropy.io import fits
from jwst.imprint.imprint_step import ImprintStep

from nirspec_pipe_testing_tool import core_utils
from . import imprint_subtract_utils



# HEADER
__author__ = "M. Pena-Guerrero"
__version__ = "1.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: separated completion from numerical tests
# Apr 2019 - Version 1.2: implemented nptt_log capability
# Apr 2023 - Version 1.3: Cleaned-up code


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files


# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "imprint_subtract"
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
    imprint_subtract_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    imprint_subtract_numerical_tests = config.getboolean("run_pytest", "_".join((step, "numerical", "tests")))
    #imprint_subtract_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [imprint_subtract_completion_tests, imprint_subtract_numerical_tests]#, imprint_subtract_validation_tests]

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # Only run step if data is IFU or MSA
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

    # Only run for IFU or MOS data
    if core_utils.check_IFU_true(inhdr) or core_utils.check_MOS_true(inhdr):
        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False
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
                # Create the pipeline step log
                stp_pipelog = "calspec2_" + step + "_" + detector + ".log"
                core_utils.mk_stpipe_log_cfg(output_dir, stp_pipelog)
                print("Pipeline step screen output will be logged in file: ", stp_pipelog)

                if os.path.isfile(step_input_file):
                    msg = " The input file "+step_input_file+" exists... will run step "+step
                    print(msg)
                    nptt_log.info(msg)
                    msa_imprint_structure = config.get("additional_arguments", "msa_imprint_structure")
                    msg = "msa_imprint_structure file: "+msa_imprint_structure
                    print(msg)
                    nptt_log.info(msg)

                    if not os.path.isfile(msa_imprint_structure):
                        msg = "Need msa_imprint_structure file to continue. Step will be skipped."
                        print(msg)
                        nptt_log.info(msg)
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                        msg = "Skipping "+step+" because msa_imprint_structure file in the configuration file does not exist."
                        print(msg)
                        nptt_log.info(msg)
                        pytest.skip(msg)
                    else:
                        msg = "*** Step "+step+" set to True"
                        print(msg)
                        nptt_log.info(msg)
                        stp = ImprintStep()

                        # check that previous pipeline steps were run up to this point
                        core_utils.check_completed_steps(step, step_input_file)

                        # get the right configuration files to run the step
                        local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                        # start the timer to compute the step running time
                        start_time = time.time()
                        if local_pipe_cfg_path == "pipe_source_tree_code":
                            result = stp.call(step_input_file, msa_imprint_structure)
                        else:
                            result = stp.call(step_input_file, msa_imprint_structure,
                                              config_file=local_pipe_cfg_path+'/imprint.cfg')
                        if result is not None:
                            result.save(step_output_file)
                            # end the timer to compute the step running time
                            end_time = repr(time.time() - start_time)   # this is in seconds
                            msg = "Step "+step+" took "+end_time+" seconds to finish"
                            print(msg)
                            nptt_log.info(msg)
                            outhdr = fits.getheader(step_output_file)
                            step_completed = True
                        else:
                            outhdr = fits.getheader(step_input_file)

                        # add the running time for this step
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                        return outhdr, step_output_file, step_input_file, run_pytests, nptt_log

                else:
                    msg = " The input file does not exist. Skipping step."
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
                    step_completed = False
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")

    else:
        pytest.skip("Skipping "+step+" because data is neither IFU or MOS.")


# VALIDATION FUNCTIONS

# fixture to validate the subtraction works fine: re-run the step with the same file as msa_imprint file
@pytest.fixture(scope="module")
def check_output_is_zero(output_vars):
    """
    This test is a simple subtraction of the input file minus a copy of the input file (instead of the actual
    MSA imprint file. The output is expected to be zero.
    Args:
        output_vars: list
    Returns:
        result: boolean
    """
    # output_vars = hdul, step_output_file, step_input_file, run_pytests
    step_input_file = output_vars[2]
    step_output_file = output_vars[1]
    # Only run test if data is IFU or MSA
    inhdr = fits.getheader(step_input_file)
    if core_utils.check_IFU_true(inhdr) or core_utils.check_MOS_true(inhdr):
        # set specifics for the test
        msa_imprint_structure = copy.deepcopy(step_input_file)
        result_to_check = step_output_file.replace(".fits", "_zerotest.fits")
        # run the step with the specifics
        stp = ImprintStep()
        res = stp.call(step_input_file, msa_imprint_structure)
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

def test_s_imprint_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[3] = imprint_subtract_completion_tests, imprint_subtract_numerical_tests,
    #                  imprint_subtract_validation_tests
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
        assert imprint_subtract_utils.s_imprint_exists(output_vars[0]), "The keyword S_IMPRINT was not added to the " \
                                                                        "header --> imprint_subtract step was " \
                                                                        "not completed."


def test_check_output_is_zero(output_vars, request):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[3] = imprint_subtract_completion_tests, imprint_subtract_numerical_tests,
    #                  imprint_subtract_validation_tests
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
