"""
py.test module for unit testing the assign_wcs step.
"""

import os
import subprocess
import time
import pytest
from astropy.io import fits
from glob import glob

from nirspec_pipe_testing_tool import core_utils
from . import assign_wcs_utils
from .. import TESTSDIR
from ..auxiliary_code import compare_wcs_ifu
from ..auxiliary_code import compare_wcs_fs
from ..auxiliary_code import compare_wcs_mos

from jwst.assign_wcs.assign_wcs_step import AssignWcsStep


# HEADER
__author__ = "M. Pena-Guerrero & G. Kanarek"
__version__ = "2.8"


# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check
# Feb 2019 - Version 2.1: Maria made changes to be able to process 491 and 492 files in the same directory
# Apr 2019 - Version 2.2: implemented nptt_log capability
# Dec 2019 - Version 2.3: implemented image processing and text file name handling
# Jun 2020 - Version 2.4: Changed comparison file to be our own instead of ESA files
# Jul 2020 - Version 2.5: Added check to see if running spec2 is appropriate according to EXP_TYPE rules taken from CRDS
# Sep 2021 - Version 2.6: Adding print statement for NPTT version
# Oct 2022 - Version 2.7: Modified exit message to be generic and added or for BOTS test
# Apr 2023 - Version 2.8: Removed running the pipeline in full and cleaned-up code


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(config):
    step = "assign_wcs"
    data_directory = config.get("calwebb_spec2_input_file", "data_directory")
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file_basename = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(data_directory, initial_input_file_basename)
    if os.path.isfile(initial_input_file):
        print("\n Taking initial input file from data_directory:")
    else:
        initial_input_file = os.path.join(output_directory, initial_input_file_basename)
        print("\n Taking initial file from output_directory: ")
    print(" Initial input file = ", initial_input_file, "\n")
    # Get the detector used
    detector = fits.getval(initial_input_file, "DETECTOR", 0)
    True_steps_suffix_map = "spec2_step_run_map_" + detector + ".txt"
    True_steps_suffix_map = os.path.join(output_directory, True_steps_suffix_map)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, output_directory)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header and other variables
@pytest.fixture(scope="module")
def output_vars(set_inandout_filenames, config):
    # determine if the tests are to be run or skipped
    run_calwebb_spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    if run_calwebb_spec2 == "skip":
        print('\n * run_calwebb_spec2 is set to skip \n')
        pytest.exit("Skipping pipeline tests for spec2, run_calwebb_spec2 is set to skip in NPTT_config file.")

    # get the general info
    set_inandout_filenames_info = core_utils.read_info4output_vars(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, outstep_file_suffix = set_inandout_filenames_info
    mode_used =  config.get("calwebb_spec2_input_file", "mode_used")

    # start the timer to compute the step test running time
    stp_test_start_time = time.time()

    # determine which steps are to be run, if not run in full
    run_pipe_step = config.getboolean("run_spec2_steps", step)

    # determine which tests are to be run
    assign_wcs_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    assign_wcs_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    assign_wcs_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests]

    # get other relevant info from NPTT config file
    compare_assign_wcs_and_extract_2d_with_esa = config.getboolean("benchmark_intermediary_products",
                                                                   "compare_assign_wcs_and_extract_2d_with_esa")
    esa_files_path = config.get("benchmark_intermediary_products", "esa_files_path")
    data_directory = config.get("calwebb_spec2_input_file", "data_directory")
    truth_file = os.path.join(data_directory, config.get("benchmark_intermediary_products", "truth_file_assign_wcs"))
    if compare_assign_wcs_and_extract_2d_with_esa:
        truth_file = esa_files_path
    print("Will use this 'truth' file to compare result of assign_wcs: ")
    print(truth_file)
    msa_shutter_conf = config.get("benchmark_intermediary_products", "msa_conf_name")
    wcs_threshold_diff = config.get("additional_arguments", "wcs_threshold_diff")
    save_wcs_plots = config.getboolean("additional_arguments", "save_wcs_plots")
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")

    # Get the main header and the detector used
    inhdr = fits.getheader(step_input_file)
    detector = inhdr["DETECTOR"]

    # Get the logfile instance for NPTT created in the run.py script
    nptt_log = os.path.join(output_directory, 'NPTT_calspec2_' + detector + '.log')
    nptt_log = core_utils.mk_nptt_log(nptt_log, reset=False)

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    if not run_calwebb_spec2:
        # read the assign wcs fits file
        outhdr = fits.getheader(step_output_file)

        return outhdr, step_output_file, msa_shutter_conf, truth_file, wcs_threshold_diff, save_wcs_plots, \
               run_pytests, mode_used, compare_assign_wcs_and_extract_2d_with_esa, nptt_log
    else:
        if run_pipe_step:
            # run the pipeline step if the file exists
            if os.path.isfile(step_input_file):
                # Create the pipeline step log
                stp_pipelog = "calspec2_" + step + "_" + detector + ".log"
                core_utils.mk_stpipe_log_cfg(output_dir, stp_pipelog)
                print("Pipeline step screen output will be logged in file: ", stp_pipelog)

                # check that previous pipeline steps were run up to this point
                core_utils.check_completed_steps(step, step_input_file)

                msg = " *** Step " + step + " set to True"
                print(msg)
                nptt_log.info(msg)
                stp = AssignWcsStep()

                # get the right configuration files to run the step
                local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")

                # start the timer to compute the step running time
                msg = "Running pipeline step " + step
                print(msg)
                nptt_log.info(msg)
                start_time = time.time()

                if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                else:
                    result = stp.call(step_input_file, config_file=local_pipe_cfg_path + '/assign_wcs.cfg')
                result.save(step_output_file)

                # end the timer to compute the step running time
                end_time = repr(time.time() - start_time)  # this is in seconds
                msg = "Step " + step + " took " + end_time + " seconds to finish"
                print(msg)
                nptt_log.info(msg)

                # remove the copy of the MSA shutter configuration file, if needed
                if core_utils.check_MOS_true(inhdu):
                    if TESTSDIR == os.path.dirname(msametfl):
                        print("Removing MSA config file from: ", TESTSDIR)
                        subprocess.run(["rm", msametfl])

            else:
                msg = "Skipping step. Input file " + step_input_file + " does not exit."
                print(msg)
                nptt_log.info(msg)
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping " + step + " because the input file does not exist.")

        else:
            msg = "Skipping running pipeline step " + step
            print(msg)
            nptt_log.info(msg)
            # add the running time for this step
            end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, output_directory)

        if os.path.isfile(step_output_file):
            outhdr = fits.getheader(step_output_file)
            step_completed = True
            # add the running time for this step
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            return outhdr, step_output_file, msa_shutter_conf, truth_file, wcs_threshold_diff, save_wcs_plots, \
                   run_pytests, mode_used, compare_assign_wcs_and_extract_2d_with_esa, nptt_log
        else:
            step_completed = False
            # add the running time for this step
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Test skipped because input file " + step_output_file + " does not exist.")


# THESE FUNCTIONS ARE TO VALIDATE THE WCS STEP

# fixture to validate the WCS
@pytest.fixture(scope="module")
def validate_wcs(output_vars):
    # get the input information for the wcs routine
    hdr = output_vars[0]
    infile_name = output_vars[1]
    msa_conf_name = output_vars[2]
    truth_file = output_vars[3]
    mode_used = output_vars[7]
    compare_assign_wcs_and_extract_2d_with_esa = output_vars[8]
    esa_files_path = None
    if compare_assign_wcs_and_extract_2d_with_esa:
        esa_files_path = output_vars[3]
        truth_file = None

    # get the logger instance
    nptt_log = output_vars[9]

    # define the threshold difference between the pipeline output and the truth files for the pytest to pass or fail
    threshold_diff = float(output_vars[4])

    # save the output plots
    save_wcs_plots = output_vars[5]

    # show the figures
    show_figs = False

    msg = " Performing WCS validation test... "
    print(msg)
    nptt_log.info(msg)
    log_msgs = None
    if core_utils.check_FS_true(hdr) or core_utils.check_BOTS_true(hdr):
        result, log_msgs = compare_wcs_fs.compare_wcs(infile_name, truth_file=truth_file, esa_files_path=esa_files_path,
                                                      show_figs=show_figs, save_figs=save_wcs_plots,
                                                      threshold_diff=threshold_diff, raw_data_root_file=None,
                                                      output_directory=None,
                                                      debug=False)

    elif core_utils.check_MOS_true(hdr) and mode_used != "MOS_sim":
        result, log_msgs = compare_wcs_mos.compare_wcs(infile_name, msa_conf_name=msa_conf_name, truth_file=truth_file,
                                                       esa_files_path=esa_files_path,
                                                       show_figs=show_figs, save_figs=save_wcs_plots,
                                                       threshold_diff=threshold_diff,
                                                       raw_data_root_file=None,
                                                       output_directory=None, debug=False)

    elif core_utils.check_IFU_true(hdr):
        result, log_msgs = compare_wcs_ifu.compare_wcs(infile_name, truth_file=truth_file,
                                                       esa_files_path=esa_files_path,
                                                       show_figs=show_figs, save_figs=save_wcs_plots,
                                                       threshold_diff=threshold_diff, raw_data_root_file=None,
                                                       output_directory=None,
                                                       debug=False)

    else:
        pytest.skip("Skipping pytest: The fits file is not FS, MOS, or IFU.")

    if log_msgs is not None:
        for msg in log_msgs:
            nptt_log.info(msg)

    if "skip" in result:
        pytest.skip("Skipping assign_wcs pytest.")
    elif "PASS" in result:
        result = True
    else:
        result = False

    return result


# Unit tests

# reference files from running calwebb_spec1
def test_rmask_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        print(" * Running reference file pytest...")
        result = assign_wcs_utils.rmask_rfile_is_correct(output_vars)
        if result is None:
            pytest.skip()
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_saturation_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.saturation_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_superbias_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.superbias_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_linearity_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.linearity_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_dark_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.dark_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_readnoise_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        print(" * Running reference file pytest...")
        result = assign_wcs_utils.readnoise_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_gain_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.gain_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


# reference files specific to the WCS step
def test_camera_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.camera_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_colimator_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.colimator_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_disperser_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.disperser_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_fore_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.fore_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_fpa_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.fpa_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_ifufore_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.ifufore_rfile_is_correct(output_vars)
        if result:
            for log_msg in result[1]:
                print(log_msg)
                nptt_log.info(log_msg)
            assert not result[0], result[0]
        else:
            msg = "Skipping ref_file pytest because data is NOT IFU."
            print(msg)
            nptt_log.info(msg)
            pytest.skip(msg)


def test_ifupost_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.ifupost_rfile_is_correct(output_vars)
        if result:
            for log_msg in result[1]:
                print(log_msg)
                nptt_log.info(log_msg)
            assert not result[0], result[0]
        else:
            msg = "Skipping ref_file pytest because data is NOT IFU."
            print(msg)
            nptt_log.info(msg)
            pytest.skip(msg)


def test_ifuslicer_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.ifuslicer_rfile_is_correct(output_vars)
        if result:
            for log_msg in result[1]:
                print(log_msg)
                nptt_log.info(log_msg)
            assert not result[0], result[0]
        else:
            msg = "Skipping ref_file pytest because data is NOT IFU."
            print(msg)
            nptt_log.info(msg)
            pytest.skip(msg)


def test_msa_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.msa_rfile_is_correct(output_vars)
        if result:
            for log_msg in result[1]:
                print(log_msg)
                nptt_log.info(log_msg)
            assert not result[0], result[0]
        else:
            msg = "Skipping ref_file pytest because data is NOT MOS."
            print(msg)
            nptt_log.info(msg)
            pytest.skip(msg)


def test_ote_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.ote_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_wavran_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = assign_wcs_utils.wavran_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


# other tests specific to the WCS step
"""
def test_wavstart_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][0]
    if not run_pytests:
        msg = "Skipping pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        assert assign_wcs_utils.wavstart_exists(output_vars[1]), "The keyword WAVSTART was not added to the header."

def test_sporder_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][0]
    if not run_pytests:
        msg = "Skipping pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        assert assign_wcs_utils.sporder_exists(output_vars[1]), "The keyword SPORDER was not added to the header."
"""


def test_s_wcs_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running completion pytest..."
        print(msg)
        nptt_log.info(msg)
        assert assign_wcs_utils.s_wcs_exists(
                output_vars[0]), "The keyword S_WCS was not added to the header --> extract_2d step was not completed."


def test_validate_wcs(output_vars, request):
    # get the logger instance
    nptt_log = output_vars[9]
    # want to run this pytest?
    # output_vars[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_vars[6][2]
    if not run_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running validation pytest..."
        print(msg)
        nptt_log.info(msg)
        # print("\n", validate_wcs, "\n")
        assert request.getfixturevalue('validate_wcs'), "Output value from compare_wcs.py is greater than threshold."
