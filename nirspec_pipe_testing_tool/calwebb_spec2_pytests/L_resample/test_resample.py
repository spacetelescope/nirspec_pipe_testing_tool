
"""
py.test module for unit testing the resample_spec step.
"""

import os
import time
import pytest
import logging
from glob import glob
from astropy.io import fits
from jwst.resample import ResampleSpecStep

from nirspec_pipe_testing_tool.utils import change_filter_opaque2science
from . import resample_utils
from nirspec_pipe_testing_tool import core_utils
from .. import TESTSDIR



# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.2"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: separated completion from other tests
# Apr 2019 - Version 1.2: implemented logging capability


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "resample_spec"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    # determine if the pipeline is to be run in full, per steps, or skipped
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    if run_calwebb_spec2 == "skip":
        print('\n * PTT finished processing run_calwebb_spec2 is set to skip. \n')
        pytest.exit("Skipping pipeline run and tests for spec2, run_calwebb_spec2 is set to skip in PTT_config file.")
    elif "T" in run_calwebb_spec2:
        run_calwebb_spec2 = True
    else:
        run_calwebb_spec2 = False

    # get the general info
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, outstep_file_suffix = set_inandout_filenames_info
    run_pipe_step = config.getboolean("run_pipe_steps", step)
    # determine which tests are to be run
    resample_spec_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    #resample_spec_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    #resample_spec_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [resample_spec_completion_tests]#, resample_spec_reffile_tests, resample_spec_validation_tests]

    end_time = '0.0'
    # Only run step if data is not IFU or BOTS
    mode_used = config.get("calwebb_spec2_input_file", "mode_used").lower()
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(output_directory, initial_input_file)
    detector = fits.getval(initial_input_file, "DETECTOR", 0)
    if not os.path.isfile(initial_input_file):
        pytest.skip("Skipping "+step+" because the initial input file given in PTT_config.cfg does not exist.")

    if mode_used != "bots" and mode_used != "ifu":
        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False

        # check if the filter is to be changed
        change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
        if change_filter_opaque:
            is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file,
                                                                                                      step=step)
            if is_filter_opaque:
                filter_opaque_msg = "With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. " \
                                    "Resample pytest now set to Skip."
                print(filter_opaque_msg)
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because the input file does not exist.")

        if run_calwebb_spec2:
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            return hdul, step_output_file, run_pytests
        else:

            if run_pipe_step:
                # Create the logfile for PTT, but erase the previous one if it exists
                PTTcalspec2_log = os.path.join(output_directory, 'PTT_calspec2_'+detector+'_'+step+'.log')
                if os.path.isfile(PTTcalspec2_log):
                    os.remove(PTTcalspec2_log)
                print("Information outputed to screen from PTT will be logged in file: ", PTTcalspec2_log)
                for handler in logging.root.handlers[:]:
                    logging.root.removeHandler(handler)
                logging.basicConfig(filename=PTTcalspec2_log, level=logging.INFO)
                # print pipeline version
                import jwst
                pipeline_version = "\n *** Using jwst pipeline version: "+jwst.__version__+" *** \n"
                print(pipeline_version)
                logging.info(pipeline_version)
                if change_filter_opaque:
                    logging.info(filter_opaque_msg)

                if os.path.isfile(step_input_file):

                    msg = " *** Step "+step+" set to True"
                    print(msg)
                    logging.info(msg)
                    stp = ResampleSpecStep()

                    # check that previous pipeline steps were run up to this point
                    core_utils.check_completed_steps(step, step_input_file)

                    # get the right configuration files to run the step
                    local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                    # start the timer to compute the step running time
                    start_time = time.time()
                    if local_pipe_cfg_path == "pipe_source_tree_code":
                        result = stp.call(step_input_file)
                    else:
                        result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/resample_spec.cfg')
                    result.save(step_output_file)
                    # end the timer to compute the step running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    msg = "Step "+step+" took "+end_time+" seconds to finish"
                    print(msg)
                    logging.info(msg)
                    step_completed = True
                    hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)

                    # rename and move the pipeline log file
                    pipelog = "pipeline_" + detector + ".log"
                    try:
                        calspec2_pilelog = "calspec2_pipeline_" + step + "_" + detector + ".log"
                        pytest_workdir = TESTSDIR
                        logfile = glob(pytest_workdir + "/" + pipelog)[0]
                        os.rename(logfile, os.path.join(output_directory, calspec2_pilelog))
                    except IndexError:
                        print("\n* WARNING: Something went wrong. Could not find a ", pipelog, " file \n")

                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return hdul, step_output_file, run_pytests

                else:
                    msg = " The input file does not exist. Skipping step."
                    print(msg)
                    logging.info(msg)
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Skipping "+step+" because the input file does not exist.")

            else:
                msg = "Skipping running pipeline step "+step
                print(msg)
                logging.info(msg)
                end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, output_directory)
                if os.path.isfile(step_output_file):
                    hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                    step_completed = True
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return hdul, step_output_file, run_pytests
                else:
                    step_completed = False
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")

    else:
        pytest.skip("Skipping "+step+" because data is either IFU or BOTS.")


# Unit tests

def test_s_resample_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[2] = resample_spec_completion_tests, resample_spec_reffile_tests, resample_spec_validation_tests
    run_pytests = output_hdul[2][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running completion pytest...\n"
        print(msg)
        logging.info(msg)
        assert resample_utils.s_resamp_exists(output_hdul[0]), "The keyword S_RESAMP was not added to the header " \
                                                               "--> Resample step was not completed."

