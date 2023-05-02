
"""
py.test module for unit testing the msa_flagging step.
"""

import pytest
import os
import time
from glob import glob
from astropy.io import fits
from jwst.msaflagopen.msaflagopen_step import MSAFlagOpenStep

from nirspec_pipe_testing_tool import core_utils
from . import msa_flagging_utils
from .. auxiliary_code import msa_flagging_testing


# HEADER
__author__ = "M. Pena-Guerrero"
__version__ = "1.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: separated completion tests from future tests
# Apr 2019 - Version 1.2: implemented nptt_log capability
# Apr 2023 - Version 1.3: Cleaned-up code


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

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # Only run step if data is MOS or IFU
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

    print("Is data MOS or IFU?", core_utils.check_MOS_true(inhdr), core_utils.check_IFU_true(inhdr))
    if core_utils.check_MOS_true(inhdr) or core_utils.check_IFU_true(inhdr):

        # determine which tests are to be run
        msa_flagging_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
        # msa_flagging_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
        msa_flagging_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
        run_pytests = [msa_flagging_completion_tests, msa_flagging_validation_tests]

        # get other info from the configuration file relevant only for this step
        msa_flagging_operability_ref = config.get("benchmark_intermediary_products", "msa_flagging_operability_ref")
        msa_flagging_threshold = config.get("additional_arguments", "msa_flagging_threshold")
        stellarity = config.get("additional_arguments", "stellarity")
        save_msa_flagging_figs = config.getboolean("additional_arguments", "save_msa_flagging_plots")
        if stellarity == "source_type":
            stellarity = None
        testing_arguments = [msa_flagging_operability_ref, msa_flagging_threshold, stellarity, save_msa_flagging_figs]

        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False
        if run_calwebb_spec2:
            if os.path.isfile(step_output_file):
                outhdr = fits.getheader(step_output_file)
            else:
                msg = "Skipping "+step+" because the output file does not exist."
                nptt_log.info(msg)
                pytest.skip(msg)
            return outhdr, step_output_file, step_input_file, run_pytests, testing_arguments, nptt_log
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
                    stp = MSAFlagOpenStep()

                    # check that previous pipeline steps were run up to this point
                    core_utils.check_completed_steps(step, step_input_file)

                    # get the right configuration files to run the step
                    # local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                    # start the timer to compute the step running time
                    start_time = time.time()
                    # if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                    # else:
                    #    result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/NOCONFIGFI.cfg')
                    result.save(step_output_file)
                    # end the timer to compute the step running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    msg = "Step "+step+" took "+end_time+" seconds to finish"
                    print(msg)
                    nptt_log.info(msg)

                    # add the running time for this step
                    step_completed = True
                    outhdr = fits.getheader(step_output_file)
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return outhdr, step_output_file, step_input_file, run_pytests, testing_arguments, nptt_log

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
                    return outhdr, step_output_file, step_input_file, run_pytests, testing_arguments, nptt_log
                else:
                    step_completed = False
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")

    else:
        pytest.skip("Skipping "+step+" because data is neither MOS or IFU.")


# fixture to validate the msa_flagging step
@pytest.fixture(scope="module")
def validate_msa_flagging(output_vars):
    step_output_file = output_vars[1]
    msa_flagging_operability_ref, msa_flagging_threshold, stellarity, save_msa_flagging_figs = output_vars[4]

    # get the logger instance and log result
    nptt_log = output_vars[-1]
    msg = " Performing MSA flagging validation test... "
    print(msg)
    nptt_log.info(msg)
    result, result_msg, log_msgs = msa_flagging_testing.run_msa_flagging_testing(
        step_output_file, msa_flagging_threshold=msa_flagging_threshold, stellarity=stellarity,
        operability_ref=msa_flagging_operability_ref, save_figs=save_msa_flagging_figs,
        show_figs=False, debug=False)
    for msg in log_msgs:
        nptt_log.info(msg)

    if result == "skip":
        pytest.skip("MSA flagging validation test will be skipped.")

    return result


# Unit tests

def test_msa_failed_open_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars = hdul, step_output_file, step_input_file, run_pytests, testing_arguments
    # output_vars[3] = msa_flagging_completion_tests, msa_flagging_validation_tests
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
        assert msa_flagging_utils.msa_failed_open_exists(output_vars[0]), "The keyword S_MSAFLG was not added to the " \
                                                                          "header --> msa_flagging step was not " \
                                                                          "completed."


def test_validate_msa_flagging(output_vars, request):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars = hdul, step_output_file, step_input_file, run_pytests, testing_arguments
    # output_vars[3] = msa_flagging_completion_tests, msa_flagging_validation_tests
    run_pytests = output_vars[3][1]
    if not run_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running validation pytest..."
        print(msg)
        nptt_log.info(msg)
        assert request.getfixturevalue("validate_msa_flagging")
