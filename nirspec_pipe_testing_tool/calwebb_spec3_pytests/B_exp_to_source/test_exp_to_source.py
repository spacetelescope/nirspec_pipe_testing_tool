"""
py.test module for unit testing the assign_wcs step.
"""

import os
import time
import pytest
import logging
from astropy.io import fits
from glob import glob

import jwst
from jwst.datamodels import (MultiExposureModel, MultiSlitModel, ModelContainer)
from jwst.exp_to_source import exp_to_source, multislit_to_container

from . import exp_to_source_utils
from nirspec_pipe_testing_tool import core_utils
from nirspec_pipe_testing_tool.calwebb_spec3_pytests import TESTSDIR

# print pipeline version
pipeline_version = "\n *** Using jwst pipeline version: " + jwst.__version__ + " *** \n"
print(pipeline_version)

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"


# HISTORY
# Jul 2020 - Version 1.0: initial version


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "master_background"
    data_directory = config.get("calwebb_spec2_input_file", "data_directory")
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file_basename = config.get("calwebb_spec3", "s3_input_file")
    initial_input_file = os.path.join(data_directory, initial_input_file_basename)
    if os.path.isfile(initial_input_file):
        print("\n Taking initial input file from data_directory:")
    else:
        initial_input_file = os.path.join(output_directory, initial_input_file_basename)
        print("\n Taking initial file from output_directory: ")
    print(" Initial input file = ", initial_input_file, "\n")
    # Get the detector used
    detector = fits.getval(initial_input_file, "DETECTOR", 0)
    True_steps_suffix_map = "spec3_steps_suffix_map_" + detector + ".txt"
    pytests_directory = TESTSDIR
    True_steps_suffix_map = os.path.join(pytests_directory, True_steps_suffix_map)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, output_directory)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    # determine if the pipeline is to be run in full, per steps, or skipped
    run_calwebb_spec3 = config.get("calwebb_spec3", "run_calwebb_spec3")
    print("run_calwebb_spec3 = ", run_calwebb_spec3)
    if run_calwebb_spec3 == "skip":
        print('\n * PTT finished processing run_calwebb_spec3 is set to skip. \n')
        pytest.exit("Finished processing file, run_calwebb_spec3 is set to skip in configuration file.")
    else:
        run_calwebb_spec3 = bool(run_calwebb_spec3)

    # get the general info
    step, step_input_filename, output_file, in_file_suffix, outstep_file_suffix, True_steps_suffix_map = set_inandout_filenames
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    txt_name = os.path.join(output_directory, True_steps_suffix_map)
    step_input_file = os.path.join(output_directory, step_input_filename)
    step_output_file = os.path.join(output_directory, output_file)
    mode_used = config.get("calwebb_spec2_input_file", "mode_used").lower()

    # start the timer to compute the step running time of NPTT
    nptt_start_time = time.time()

    # determine if steps is to be run
    run_pipe_step = config.getboolean("run_pipe_steps", step)

    # determine which tests are to be run
    #completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    #reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    #validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    #run_pytests = [completion_tests, reffile_tests, validation_tests]
    # TODO - booleans for different tests have not been implemented in the configuration file yet
    run_pytests = [False]

    # Get the detector used
    detector = fits.getval(step_input_file, "DETECTOR", 0)

    # get main header from input file
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)

    # if run_calwebb_spec3 is True, calwebb_spec3 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # only do this step if data is NOT IFU
    if mode_used == 'ifu':
        pytest.skip("Skipping test for step "+step+" because data is IFU.")

    # run the pipeline
    if run_calwebb_spec3:

        hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
        return hdul, step_output_file, run_pytests, mode_used

    else:

        if run_pipe_step:

            # Create the logfile for NPTT, but erase the previous one if it exists
            npttcalspec3_log = os.path.join(output_directory, 'NPTT_calspec3_' + detector + '_' + step + '.log')
            if os.path.isfile(npttcalspec3_log):
                os.remove(npttcalspec3_log)
            print("Output information on screen from NPTT will be logged in file: ", npttcalspec3_log)
            for handler in logging.root.handlers[:]:
                logging.root.removeHandler(handler)
            logging.basicConfig(filename=npttcalspec3_log, level=logging.INFO)
            logging.info(pipeline_version)

            # check that previous pipeline steps were run up to this point
            core_utils.check_completed_steps(step, step_input_file)

            # TODO - need a function to gather all the same source files into a list to feed into this pipeline step

            if os.path.isfile(step_input_file):
                msg = " *** Step " + step + " set to True"
                print(msg)
                logging.info(msg)

                print("Not running this individual step at the moment.")

                """
                *** NEED THE COMMANDS TO RUN THIS STEP - requires a list of _cal.fits files
                
                
                stp = MasterBackgroundStep()

                # get the right configuration files to run the step
                local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")

                # start the timer to compute the step running time
                print("running pipeline...")
                start_time = time.time()
                if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                else:
                    result = stp.call(step_input_file, config_file=local_pipe_cfg_path + '/master_background.cfg')
                result.save(step_output_file)

                # end the timer to compute the step running time
                end_time = repr(time.time() - start_time)  # this is in seconds
                msg = "Step " + step + " took " + end_time + " seconds to finish"
                print(msg)
                logging.info(msg)

                # rename and move the pipeline log file
                pipelog = "pipeline_" + detector + ".log"
                try:
                    calspec3_pipelog = "calspec3_pipeline_" + step + "_" + detector + ".log"
                    pytest_workdir = TESTSDIR
                    logfile = glob(pytest_workdir + "/" + pipelog)[0]
                    os.rename(logfile, os.path.join(output_directory, calspec3_pipelog))
                except IndexError:
                    print("\n* WARNING: Something went wrong. Could not find a ", pipelog, " file \n")
                """

            else:
                msg = "Skipping step. Input file " + step_input_file + " does not exit."
                print(msg)
                logging.info(msg)
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping " + step + " because the input file does not exist.")

        else:
            print("Skipping running pipeline step ", step)
            # add the running time for this step
            end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, output_directory)

        if os.path.isfile(step_output_file):
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            step_completed = True
            # add the running time for this step
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            return hdul, step_output_file, run_pytests, mode_used
        else:
            step_completed = False
            # add the running time for this step
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")


# reference file tests

"""
def test_wavran_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest...\n"
        print(msg)
        logging.info(msg)
        result = master_background_utils.wavran_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]
"""

# keyword tests


def test_masterbg_exists(output_hdul):
    # want to run this pytest?
    # output_hdul = hdul, step_output_file, run_pytests, mode_used
    # run_pytests = completion_tests, reffile_tests, validation_tests
    run_pytests = output_hdul[2][0]
    if not run_pytests:
        msg = "Skipping pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest...\n"
        print(msg)
        logging.info(msg)
        assert exp_to_source_utils.masterbg_exists(output_hdul[1]), "The keyword MASTERBG was not added to " \
                                                                    "the header."


