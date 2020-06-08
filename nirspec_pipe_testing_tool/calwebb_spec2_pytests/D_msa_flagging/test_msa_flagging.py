
"""
py.test module for unit testing the msa_flagging step.
"""

import pytest
import os
import time
import logging
from glob import glob
from astropy.io import fits
from jwst.msaflagopen.msaflagopen_step import MSAFlagOpenStep

from .. import core_utils
from .. import TESTSDIR
from . import msa_flagging_utils


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.2"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: separated completion tests from future tests
# Apr 2019 - Version 1.2: implemented logging capability


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
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    # determine which tests are to be run
    msa_flagging_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    #msa_flagging_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    #msa_flagging_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [msa_flagging_completion_tests]#, msa_flagging_reffile_tests, msa_flagging_validation_tests]
    # determine which steps are to be run, if not run in full
    run_pipe_step = config.getboolean("run_pipe_steps", step)

    end_time = '0.0'

    # Only run step if data is MOS or IFU
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(output_directory, initial_input_file)
    if os.path.isfile(initial_input_file):
        inhdu = core_utils.read_hdrfits(initial_input_file, info=False, show_hdr=False)
        detector = fits.getval(initial_input_file, "DETECTOR", 0)
    else:
        pytest.skip("Skipping "+step+" because the initial input file given in PTT_config.cfg does not exist.")

    print ("Is data MOS or IFU?", core_utils.check_MOS_true(inhdu), core_utils.check_IFU_true(inhdu))
    if core_utils.check_MOS_true(inhdu) or core_utils.check_IFU_true(inhdu):

        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False
        if run_calwebb_spec2:
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            return hdul, step_output_file, run_pytests
        else:
            if run_pipe_step:

                if os.path.isfile(step_input_file):

                    msg = " The input file "+step_input_file+" exists... will run step "+step
                    print(msg)
                    logging.info(msg)
                    stp = MSAFlagOpenStep()

                    # check that previous pipeline steps were run up to this point
                    core_utils.check_completed_steps(step, step_input_file)

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

                    # get the right configuration files to run the step
                    #local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                    # start the timer to compute the step running time
                    start_time = time.time()
                    #if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                    #else:
                    #    result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/NOCONFIGFI.cfg')
                    result.save(step_output_file)
                    # end the timer to compute the step running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    msg = "Step "+step+" took "+end_time+" seconds to finish"
                    print(msg)
                    logging.info(msg)

                    # rename and move the pipeline log file
                    try:
                        calspec2_pilelog = "calspec2_pipeline_" + step + "_" + detector + ".log"
                        pytest_workdir = TESTSDIR
                        logfile = glob(pytest_workdir + "/pipeline.log")[0]
                        os.rename(logfile, os.path.join(output_directory, calspec2_pilelog))
                    except:
                        IndexError

                    # add the running time for this step
                    step_completed = True
                    hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
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
                    return hdul, step_output_file, step_input_file, run_pytests
                else:
                    step_completed = False
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")

    else:
        pytest.skip("Skipping "+step+" because data is neither MOS or IFU.")


# Unit tests

def test_msa_failed_open_exists(output_hdul):
    # want to run this pytest?
    # output_hdu[2] = msa_flagging_completion_tests, msa_flagging_reffile_tests, msa_flagging_validation_tests
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
        assert msa_flagging_utils.msa_failed_open_exists(output_hdul[0]), "The keyword S_MSAFLG was not added to the header --> msa_flagging step was not completed."
