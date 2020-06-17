
"""
py.test module for unit testing the extract_1d step.
"""

import os
import time
import pytest
import logging
from astropy.io import fits
from glob import glob
from jwst.extract_1d.extract_1d_step import Extract1dStep

from nirspec_pipe_testing_tool.utils import change_filter_opaque2science
from . import extract_1d_utils
from nirspec_pipe_testing_tool import core_utils
from .. import TESTSDIR


# HEADER
__author__ = "M. A. Pena-Guerrero & Gray Kanarek"
__version__ = "2.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check
# Mar 2019 - Version 2.1: Maria added infrastructure to separate completion from other tests.
# Apr 2019 - Version 2.2: implemented logging capability
# Dec 2019 - Version 2.3: implemented imaging text file name handling capability

# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(config):
    step = "extract_1d"
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
    extract_1d_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    extract_1d_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    #extract_1d_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [extract_1d_completion_tests, extract_1d_reffile_tests]#, extract_1d_validation_tests]

    # Make sure at this point the MSA shutter configuration file is removed from the calwebb_spec2_pytests directory
    msa_shutter_conf = config.get("esa_intermediary_products", "msa_conf_name")
    msametfl = os.path.basename(msa_shutter_conf)
    if os.getcwd() != os.path.dirname(msa_shutter_conf):
        print("Removing MSA config file from: ", os.getcwd())
        os.remove(msametfl)

    # check if processing an image, then set proper variables
    imaging_mode = False
    mode_used = config.get("calwebb_spec2_input_file", "mode_used").lower()
    if mode_used in ('image', 'confirm', 'taconfirm', 'wata', 'msata', 'bota', 'focus', 'mimf'):
        run_calwebb_spec2 = True
        imaging_mode = True

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # Get the detector used
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(output_directory, initial_input_file)
    if os.path.isfile(initial_input_file):
        detector = fits.getval(initial_input_file, "DETECTOR", 0)
    else:
        pytest.skip("Skipping "+step+" because the initial input file given in PTT_config.cfg does not exist.")

    # check if the filter is to be changed
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    if change_filter_opaque:
        is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file,
                                                                                                  step=step)
        if is_filter_opaque:
            filter_opaque_msg = "With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. " \
                                "Extract_1d pytest now set to Skip."
            print(filter_opaque_msg)
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            #core_utils.convert_html2pdf()   # convert the html report into a pdf file
            # move the final reporting files to the working directory
            core_utils.move_txt_files_2workdir(config, detector)
            pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

    if run_calwebb_spec2:
        # read the output header
        hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)

        # end the timer to compute the step running time of PTT
        PTT_end_time = time.time()
        core_utils.start_end_PTT_time(txt_name, start_time=None, end_time=PTT_end_time)

        # move the final reporting files to the working directory
        core_utils.move_txt_files_2workdir(config, detector)

        return hdul, step_output_file, run_pytests

    else:

        if run_pipe_step:

            # Create the logfile for PTT, but erase the previous one if it exists
            PTTcalspec2_log = os.path.join(output_directory, 'PTT_calspec2_'+detector+'_'+step+'.log')
            if imaging_mode:
                PTTcalspec2_log = PTTcalspec2_log.replace('calspec2', 'calimage2')
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
                stp = Extract1dStep()

                # check that previous pipeline steps were run up to this point
                core_utils.check_completed_steps(step, step_input_file)

                # get the right configuration files to run the step
                local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                # start the timer to compute the step running time
                start_time = time.time()
                if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                else:
                    result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/extract_1d.cfg')
                result.save(step_output_file)
                # end the timer to compute the step running time
                end_time = repr(time.time() - start_time)   # this is in seconds
                msg = "Step "+step+" took "+end_time+" seconds to finish"
                print(msg)
                logging.info(msg)

                # rename and move the pipeline log file
                calspec2_pilelog = "calspec2_pipeline_"+step+"_"+detector+".log"
                if imaging_mode:
                    calspec2_pilelog = calspec2_pilelog.replace('calspec2', 'calimage2')
                pytest_workdir = TESTSDIR
                logfile = glob(pytest_workdir+"/pipeline.log")[0]
                os.rename(logfile, os.path.join(output_directory, calspec2_pilelog))

            else:
                msg = " The input file does not exist. Skipping step."
                print(msg)
                logging.info(msg)
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                #core_utils.convert_html2pdf()   # convert the html report into a pdf file
                # end the timer to compute the step running time of PTT
                PTT_end_time = time.time()
                core_utils.start_end_PTT_time(txt_name, start_time=None, end_time=PTT_end_time)
                # move the final reporting files to the working directory
                core_utils.move_txt_files_2workdir(config, detector)
                # skip the test if input file does not exist
                pytest.skip("Skipping "+step+" because the input file does not exist.")

        else:
            msg = "Skipping running pipeline step "+step
            print(msg)
            logging.info(msg)
            # get the running time for this step
            end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, output_directory)

        # add the running time for this step
        if os.path.isfile(step_output_file):
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            step_completed = True
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)

        else:
            step_completed = False
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            # move the final reporting files to the working directory
            core_utils.move_txt_files_2workdir(config, detector)
            pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")

        # get the total running time and print it in the file
        total_time = repr(core_utils.get_time_to_run_pipeline(txt_name))
        total_time_min = repr(round(float(total_time)/60.0, 2))
        msg = "\n\n **** The total time for the pipeline to run was "+total_time+" seconds."
        print(msg)
        logging.info(msg)
        line2write = "{:<20} {:<20} {:<20} {:<20}".format('', '', 'total_time  ', total_time+'  ='+total_time_min+'min')
        print(line2write)
        logging.info(line2write)
        with open(txt_name, "a") as tf:
            tf.write(line2write+"\n")

        # convert the html report into a pdf file
        #core_utils.convert_html2pdf()

        # end the timer to compute the step running time of PTT
        PTT_end_time = time.time()
        core_utils.start_end_PTT_time(txt_name, start_time=None, end_time=PTT_end_time)

        # move the final reporting text files to the working directory
        core_utils.move_txt_files_2workdir(config, detector)

        # rename and move the PTT log file
        try:
            PTT_calspec2_log = "PTT_calspec2_" + detector + ".log"
            if imaging_mode:
                PTT_calspec2_log = PTT_calspec2_log.replace('calspec2', 'calimage2')
            pytest_workdir = TESTSDIR
            logfile = glob(pytest_workdir + "/pipeline.log")[0]
            os.rename(logfile, os.path.join(output_directory, PTT_calspec2_log))
        except:
            IndexError

        return hdul, step_output_file, run_pytests


"""
# Move the files at the end of the pytests
@pytest.fixture(scope="function", autouse=True)
def move_output_files(request):
    def fin():
        #core_utils.convert_html2pdf()
        # move the final reporting files to the working directory
        core_utils.move_txt_files_2workdir(config, detector)
        msg = "Output report files have been moved to the output_directory path indicated in the PTT_config file."
        print(msg)
        logging.info(msg)
    request.addfinalizer(fin)

### This function is commented out because it only executes after the test, but not if it is skipped. The way to call
### the fixture is for instance:
###      def test_s_extr1d_exists(output_hdul, move_output_files):
### it will seem that move_output_files is not used, but it is called under the hood within the pytest framework.
"""


# Unit tests

def test_extract1d_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[2] = extract_1d_completion_tests, extract_1d_reffile_tests, extract_1d_validation_tests
    run_pytests = output_hdul[2][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest...\n"
        print(msg)
        logging.info(msg)
        result = extract_1d_utils.extract1d_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]


def test_s_extr1d_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[2] = extract_1d_completion_tests, extract_1d_reffile_tests, extract_1d_validation_tests
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
        assert extract_1d_utils.s_extr1d_exists(output_hdul[0]), "The keyword S_EXTR1D was not added to the " \
                                                                 "header --> Extract 1D step was not completed."


