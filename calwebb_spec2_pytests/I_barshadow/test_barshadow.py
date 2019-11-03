
"""
py.test module for unit testing the barshadow step.
"""

import os
import time
import pytest
import logging
from glob import glob
from astropy.io import fits
from jwst.barshadow.barshadow_step import BarShadowStep

from . import barshadow_utils
from .. import core_utils
from .. auxiliary_code import change_filter_opaque2science
from .. auxiliary_code import barshadow_testing



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
    barshadow_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [barshadow_completion_tests, barshadow_validation_tests]
    barshadow_threshold_diff = config.get("additional_arguments", "barshadow_threshold_diff")
    save_barshadow_final_plot = config.getboolean("additional_arguments", "save_barshadow_final_plot")
    save_barshadow_intermediary_plots = config.getboolean("additional_arguments", "save_barshadow_intermediary_plots")
    write_barshadow_files = config.getboolean("additional_arguments", "write_barshadow_files")
    barshadow_switches = [barshadow_threshold_diff, save_barshadow_final_plot, save_barshadow_intermediary_plots, write_barshadow_files]

    end_time = '0.0'

    # Only run step if data is MOS
    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(working_directory, initial_input_file)
    if os.path.isfile(initial_input_file):
        inhdu = core_utils.read_hdrfits(initial_input_file, info=False, show_hdr=False)
        detector = fits.getval(initial_input_file, "DETECTOR", 0)
    else:
        pytest.skip("Skipping "+step+" because the initial input file given in PTT_config.cfg does not exist.")

    print("core_utils.check_MOS_true(inhdu)=", core_utils.check_MOS_true(inhdu))
    if core_utils.check_MOS_true(inhdu):
        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False

        # check if the filter is to be changed
        change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
        if change_filter_opaque:
            is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file, step=step)
            if is_filter_opaque:
                filter_opaque_msg = "With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Barshadow pytest now set to Skip."
                print(filter_opaque_msg)
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

        if run_calwebb_spec2:
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            return hdul, step_output_file, barshadow_switches, run_pytests

        else:
            if run_pipe_step:

                # Create the logfile for PTT, but erase the previous one if it exists
                PTTcalspec2_log = os.path.join(working_directory, 'PTT_calspec2_'+detector+'_'+step+'.log')
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
                    msg = " * calwebb_spec2 took "+end_time+" seconds to finish."
                    print(msg)
                    logging.info(msg)
                    step_completed = True
                    hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)

                    # rename and move the pipeline log file
                    try:
                        calspec2_pilelog = "calspec2_pipeline_"+step+"_"+detector+".log"
                        pytest_workdir = os.getcwd()
                        logfile = glob(pytest_workdir+"/pipeline.log")[0]
                        os.rename(logfile, os.path.join(working_directory, calspec2_pilelog))
                    except:
                        IndexError

                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return hdul, step_output_file, barshadow_switches, run_pytests

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
                end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, working_directory)
                if os.path.isfile(step_output_file):
                    hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                    step_completed = True
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return hdul, step_output_file, barshadow_switches, run_pytests
                else:
                    step_completed = False
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip()


    else:
        pytest.skip("Skipping "+step+" because data is not MOS.")



# fixture to validate the barshadow correction
@pytest.fixture(scope="module")
def validate_barshadow(output_hdul):
    hdu = output_hdul[0]
    log_msgs = None

    # show the figures
    show_figs = False

    # Determine if data is MOS
    if core_utils.check_MOS_true(hdu):
        # Determine if the source is point or extended. If extended, unknown, or not present, correction will be applied.
        if 'SRCTYPE' in hdu:
            if hdu['SRCTYPE'] == 'POINT':
                pytest.skip("Skipping pytest: The SRCTYPE keyword is set to POINT so barshadow correction is not applied.")

        plfile = output_hdul[1].replace('_barshadow', '_pathloss')
        bsfile = output_hdul[1]
        barshadow_threshold_diff, save_barshadow_final_plot, save_barshadow_intermediary_plots, write_barshadow_files = output_hdul[2]
        barshadow_testresult, result_msg, log_msgs = barshadow_testing.run_barshadow_tests(plfile, bsfile,
                                                        barshadow_threshold_diff=float(barshadow_threshold_diff),
                                                        save_final_figs=save_barshadow_final_plot,
                                                        show_final_figs=show_figs,
                                                        save_intermediary_figs=save_barshadow_intermediary_plots,
                                                        show_intermediary_figs=show_figs,
                                                        write_barshadow_files = write_barshadow_files,
                                                        debug=False)

    else:
        pytest.skip("Skipping pytest: The input fits file is not MOS.")

    if log_msgs is not None:
        for msg in log_msgs:
            logging.info(msg)

    if barshadow_testresult == "skip":
        logging.info(result_msg)
        pytest.skip(result_msg)
    else:
        print(result_msg)
        logging.info(result_msg)

    return barshadow_testresult



# Unit tests

def test_s_barsha_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[3] = barshadow_completion_tests, barshadow_validation_tests
    run_pytests = output_hdul[3][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running completion pytest...\n"
        print(msg)
        logging.info(msg)
        assert barshadow_utils.s_barsha_exists(output_hdul[0]), "The keyword S_BARSHA was not added to the header --> Barshadow step was not completed."


def test_validate_barshadow(output_hdul, request):
    # want to run this pytest?
    # output_hdul[3] = barshadow_completion_tests, barshadow_validation_tests
    run_pytests = output_hdul[3][1]
    if not run_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running validation pytest...\n"
        print(msg)
        logging.info(msg)
        assert request.getfixturevalue("validate_barshadow"), "Output value from barshadow_testing.py is greater than threshold."


