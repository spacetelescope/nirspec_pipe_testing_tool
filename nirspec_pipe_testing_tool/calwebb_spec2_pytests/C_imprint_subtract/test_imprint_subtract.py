
"""
py.test module for unit testing the imprint_subtract step.
"""

import pytest
import os
import time
import copy
import subprocess
import logging
from glob import glob
from astropy.io import fits
from jwst.imprint.imprint_step import ImprintStep

from nirspec_pipe_testing_tool import core_utils
from .. import TESTSDIR
from . import imprint_subtract_utils



# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.2"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: separated completion from numerical tests
# Apr 2019 - Version 1.2: implemented logging capability


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

    # determine which steps are to be run, if not run in full
    run_pipe_step = config.getboolean("run_pipe_steps", step)
    # determine which tests are to be run
    imprint_subtract_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    imprint_subtract_numerical_tests = config.getboolean("run_pytest", "_".join((step, "numerical", "tests")))
    #imprint_subtract_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [imprint_subtract_completion_tests, imprint_subtract_numerical_tests]#, imprint_subtract_validation_tests]

    end_time = '0.0'

    # Only run step if data is IFU or MSA
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(output_directory, initial_input_file)
    detector = fits.getval(initial_input_file, "DETECTOR", 0)
    calspec2_pilelog = "calspec2_pipeline_" + step + "_" + detector + ".log"
    pytest_workdir = TESTSDIR

    if os.path.isfile(initial_input_file):
        inhdu = core_utils.read_hdrfits(initial_input_file, info=False, show_hdr=False)
    else:
        pytest.skip("Skipping "+step+" because the initial input file given in PTT_config.cfg does not exist.")

    if core_utils.check_IFU_true(inhdu) or core_utils.check_MOS_true(inhdu):

        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False
        if run_calwebb_spec2:
            if os.path.isfile(step_output_file):
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            else:
                pytest.skip("Skipping "+step+" because the output file does not exist.")
            return hdul, step_output_file, run_pipe_step, step_input_file, run_pytests

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
                    msa_imprint_structure = config.get("additional_arguments", "msa_imprint_structure")
                    msg = "msa_imprint_structure file: "+msa_imprint_structure
                    print(msg)
                    logging.info(msg)

                    if not os.path.isfile(msa_imprint_structure):
                        print (" Need msa_imprint_structure file to continue. Step will be skipped.")
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                        pytest.skip("Skipping "+step+" because msa_imprint_structure file in the configuration file does not exist.")

                    else:

                        msg = "*** Step "+step+" set to True"
                        print(msg)
                        logging.info(msg)
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
                            logging.info(msg)
                            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                            step_completed = True
                        else:
                            hdul = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)

                        # rename and move the pipeline log file
                        try:
                            logfile = glob(pytest_workdir + "/pipeline.log")[0]
                            os.rename(logfile, os.path.join(output_directory, calspec2_pilelog))
                        except:
                            IndexError

                        # add the running time for this step
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                        return hdul, step_output_file, run_pipe_step, step_input_file, run_pytests

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
        pytest.skip("Skipping "+step+" because data is neither IFU or MOS.")


# VALIDATION FUNCTIONS

# fixture to validate the subtraction works fine: re-run the step with the same file as msa_imprint file
@pytest.fixture(scope="module")
def check_output_is_zero(output_hdul):
    run_step = output_hdul[2]
    step_input_file = output_hdul[3]
    step_output_file = output_hdul[1]
    # Only run test if data is IFU or MSA
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if core_utils.check_IFU_true(inhdu) or core_utils.check_MOS_true(inhdu):
        if run_step:
            # set specifics for the test
            msa_imprint_structure = copy.deepcopy(step_input_file)
            result_to_check = step_output_file.replace(".fits", "_zerotest.fits")
            # run the step with the specifics
            stp = ImprintStep()
            res = stp.call(step_input_file, msa_imprint_structure)
            res.save(result_to_check)
            # check that the end product of image - image is zero
            c = fits.getdata(result_to_check)
            substraction = sum(c.flatten())
            result = False
            if substraction == 0.0:
                result = True
            # erase test output file
            subprocess.run(["rm", result_to_check])
            return result
        pytest.skip("Skipping validation check for imprint_subtract. Step set to False in configuration file.")



### Unit tests

def test_s_imprint_exists(output_hdul):
    # want to run this pytest?
    # output_hdu[4] = imprint_subtract_completion_tests, imprint_subtract_numerical_tests, imprint_subtract_validation_tests
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
        assert imprint_subtract_utils.s_imprint_exists(output_hdul[0]), "The keyword S_IMPRINT was not added to the header --> imprint_subtract step was not completed."

def test_check_output_is_zero(output_hdul, request):
    # want to run this pytest?
    # output_hdu[4] = imprint_subtract_completion_tests, imprint_subtract_numerical_tests, imprint_subtract_validation_tests
    run_pytests = output_hdul[4][1]
    if not run_pytests:
        msg = "Skipping pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running numerical accuracy pytest...\n"
        print(msg)
        logging.info(msg)
        assert request.getfixturevalue('check_output_is_zero'), "Substraction result is not equal to zero."
