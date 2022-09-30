
"""
py.test module for unit testing the extract_2d step.
"""

import pytest
import os
import time
import logging
from glob import glob
from astropy.io import fits

from jwst.extract_2d.extract_2d_step import Extract2dStep

from nirspec_pipe_testing_tool import core_utils
from .. import TESTSDIR
from . import extract_2d_utils
from .. auxiliary_code import compare_wcs_mos
from .. auxiliary_code import check_corners_extract2d


# HEADER
__author__ = "M. A. Pena-Guerrero & G. Kanarek"
__version__ = "2.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Jan 2019 - Version 2.0: test separated from assign_wcs
# Mar 2019 - Version 2.1: separated completion from validation tests
# Apr 2019 - Version 2.2: implemented logging capability
# Jun 2020 - Version 2.3: Changed comparison file to be our own instead of ESA files


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "extract_2d"
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
    extract_2d_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    extract_2d_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    assign_wcs_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [extract_2d_completion_tests, extract_2d_validation_tests, assign_wcs_validation_tests]
    # get other relevant info from PTT config file
    compare_assign_wcs_and_extract_2d_with_esa = config.getboolean("benchmark_intermediary_products",
                                                                   "compare_assign_wcs_and_extract_2d_with_esa")
    esa_files_path = config.get("benchmark_intermediary_products", "esa_files_path")
    data_directory = config.get("calwebb_spec2_input_file", "data_directory")
    truth_file = os.path.join(data_directory, config.get("benchmark_intermediary_products", "truth_file_extract_2d"))
    if compare_assign_wcs_and_extract_2d_with_esa:
        truth_file = esa_files_path
    print("Will use this 'truth' file to compare result of extract_2d: ")
    print(truth_file)
    msa_conf_name = config.get("benchmark_intermediary_products", "msa_conf_name")
    extract_2d_threshold_diff = int(config.get("additional_arguments", "extract_2d_threshold_diff"))
    
    # Check if the mode used is MOS_sim and get the threshold for the assign_wcs test
    mode_used = config.get("calwebb_spec2_input_file", "mode_used").lower()
    wcs_threshold_diff = config.get("additional_arguments", "wcs_threshold_diff")
    save_wcs_plots = config.getboolean("additional_arguments", "save_wcs_plots")
    
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # only do this step if data is NOT IFU
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(output_directory, initial_input_file)
    if os.path.isfile(initial_input_file):
        inhdu = core_utils.read_hdrfits(initial_input_file, info=False, show_hdr=False)
        detector = fits.getval(initial_input_file, "DETECTOR", 0)
    else:
        pytest.skip("Skipping "+step+" because the initial input file given in PTT_config.cfg does not exist.")

    if not core_utils.check_IFU_true(inhdu):
        if run_calwebb_spec2:
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            return hdul, step_output_file, msa_conf_name, truth_file, run_pytests, mode_used, wcs_threshold_diff, \
                   save_wcs_plots, extract_2d_threshold_diff, compare_assign_wcs_and_extract_2d_with_esa
            
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

                if os.path.isfile(step_input_file):
                    msg = " The input file "+step_input_file+" exists... will run step "+step
                    print(msg)
                    logging.info(msg)
                    stp = Extract2dStep()
                        
                    # check that previous pipeline steps were run up to this point
                    core_utils.check_completed_steps(step, step_input_file)
                    
                    # get the right configuration files to run the step
                    local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")

                    # start the timer to compute the step running time
                    start_time = time.time()
                    if local_pipe_cfg_path == "pipe_source_tree_code":
                        result = stp.call(step_input_file)
                    else:
                        result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/extract_2d.cfg')
                    result.save(step_output_file)
                    step_completed = True
                    hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)

                    # end the timer to compute the step running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    msg = "Step "+step+" took "+end_time+" seconds to finish"
                    print(msg)
                    logging.info(msg)

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
                    return hdul, step_output_file, msa_conf_name, truth_file, run_pytests, mode_used, \
                           wcs_threshold_diff, save_wcs_plots, extract_2d_threshold_diff, compare_assign_wcs_and_extract_2d_with_esa

                else:
                    msg = " The input file does not exist. Skipping step."
                    print(msg)
                    logging.info(msg)
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Skiping "+step+" because the input file does not exist.")

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
                    return hdul, step_output_file, msa_conf_name, truth_file, run_pytests, mode_used, \
                           wcs_threshold_diff, save_wcs_plots, extract_2d_threshold_diff, compare_assign_wcs_and_extract_2d_with_esa
                else:
                    step_completed = False
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")

    else:
        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
        pytest.skip("Skipping "+step+" because data is IFU.")


# THESE FUNCTIONS ARE TO VALIDATE BOTH THE WCS AND THE 2D_EXTRACT STEPS

# fixture to validate the WCS and extract 2d steps only for MOS simulations
@pytest.fixture(scope="module")
def validate_MOSsim_wcs_extract2d(output_hdul):
    # get the input information for the wcs routine
    infile_name = output_hdul[1]
    msa_conf_name = output_hdul[2]
    truth_file = output_hdul[3]
    mode_used = output_hdul[5]
    compare_assign_wcs_and_extract_2d_with_esa = output_hdul[9]
    esa_files_path = None
    if compare_assign_wcs_and_extract_2d_with_esa:
        esa_files_path = output_hdul[3]
        truth_file = None

    # define the threshold difference between the pipeline output and the benchmark files for the pytest to pass or fail
    threshold_diff = float(output_hdul[6])
    
    # save the output plots
    save_wcs_plots = output_hdul[7]
    
    # show the figures
    show_figs = False
    
    msg = "\n Performing WCS validation test... "
    print(msg)
    logging.info(msg)
    if mode_used == "mos_sim":
        result, log_msgs = compare_wcs_mos.compare_wcs(infile_name, msa_conf_name=msa_conf_name, truth_file=truth_file,
                                                       esa_files_path=esa_files_path,
                                                       show_figs=show_figs, save_figs=save_wcs_plots,
                                                       threshold_diff=threshold_diff,
                                                       raw_data_root_file=None,
                                                       output_directory=None, debug=False)
        for msg in log_msgs:
            logging.info(msg)

    else:
        pytest.skip("Skipping pytest for WCS validation: The fits file is not MOS simulated data, the validation "
                    "test was done after the assign_wcs step.")
    
    if result == "skip":
        pytest.skip("Pytest for assign_wcs validation after extract_2d will be skipped.")
    
    return result


# fixture to validate extract 2d step
@pytest.fixture(scope="module")
def validate_extract2d(output_hdul):
    # get the input information for the wcs routine
    hdu = output_hdul[0]
    infile_name = output_hdul[1]
    msa_conf_name = output_hdul[2]
    truth_file = output_hdul[3]
    extract_2d_threshold_diff = output_hdul[8]
    compare_assign_wcs_and_extract_2d_with_esa = output_hdul[9]
    esa_files_path = None
    if compare_assign_wcs_and_extract_2d_with_esa:
        esa_files_path = output_hdul[3]
        truth_file = None
    print('Will be using this number of pixels as threshold for extract_2d test: ', extract_2d_threshold_diff)

    msg = "\n Performing extract_2d validation test... "
    print(msg)
    logging.info(msg)
    if core_utils.check_FS_true(hdu) or core_utils.check_BOTS_true(hdu):
        result, log_msgs = check_corners_extract2d.find_FSwindowcorners(infile_name, truth_file=truth_file,
                                                                        esa_files_path=esa_files_path,
                                                                        extract_2d_threshold_diff=
                                                                        extract_2d_threshold_diff)

    elif core_utils.check_MOS_true(hdu):
        result, log_msgs = check_corners_extract2d.find_MOSwindowcorners(infile_name, msa_conf_name,
                                                                         truth_file=truth_file,
                                                                         esa_files_path=esa_files_path,
                                                                         extract_2d_threshold_diff=
                                                                         extract_2d_threshold_diff)
        
    else:
        pytest.skip("Skipping pytest: The fits file is not FS or MOS.")
    
    for msg in log_msgs:
        logging.info(msg)

    if result == "skip":
        pytest.skip("Extract_2d validation will be skipped.")
        
    return result


# Unit tests

def test_s_ext2d_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[4] = extract_2d_completion_tests, extract_2d_validation_tests
    run_pytests = output_hdul[4][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running completion pytest...\n"
        print(msg)
        logging.info(msg)
        assert extract_2d_utils.s_ext2d_exists(output_hdul[0]), "The keyword S_EXTR2D was not added to the header " \
                                                                "--> extract_2d step was not completed."


def test_validate_MOSsim_wcs_extract2d(output_hdul, request):
    # want to run this pytest? For this particular case, check both for the extract_2d step and for assign_wcs
    # output_hdul[4] = extract_2d_completion_tests, extract_2d_validation_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[4][1]
    assign_wcs_pytests = output_hdul[4][2]
    if not run_pytests and not assign_wcs_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running validation pytest...\n"
        print(msg)
        logging.info(msg)
        assert request.getfixturevalue("validate_MOSsim_wcs_extract2d"), "Output value from compare_wcs.py is " \
                                                                         "greater than threshold."


def test_validate_extract2d(output_hdul, request):
    # want to run this pytest?
    # output_hdul[4] = extract_2d_completion_tests, extract_2d_validation_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[4][1]
    if not run_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running validation pytest...\n"
        print(msg)
        logging.info(msg)
        assert request.getfixturevalue("validate_extract2d")
