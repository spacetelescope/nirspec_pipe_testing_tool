
"""
py.test module for unit testing the extract_1d step.
"""

import os
import time
import pytest
from astropy.io import fits
from jwst.extract_1d.extract_1d_step import Extract1dStep

from .. auxiliary_code import change_filter_opaque2science
from . import extract_1d_utils
from .. import core_utils


# HEADER
__author__ = "M. A. Pena-Guerrero & Gray Kanarek"
__version__ = "2.1"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check
# Mar 2019 - Version 2.1: Maria added infrastructure to separate completion from other tests.

# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "extract_1d"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
    run_pipe_step = config.getboolean("run_pipe_steps", step)
    # determine which tests are to be run
    extract_1d_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    extract_1d_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    #extract_1d_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [extract_1d_completion_tests, extract_1d_reffile_tests]#, extract_1d_validation_tests]

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # Get the detector used
    detector = fits.getval(step_input_file, "DETECTOR", 0)

    # check if the filter is to be changed
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    if change_filter_opaque:
        is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file, step=step)
        if is_filter_opaque:
            print ("With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Extract_1d pytest now set to Skip.")
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            #core_utils.convert_html2pdf()   # convert the html report into a pdf file
            # move the final reporting files to the working directory
            core_utils.move_latest_report_and_txt_2workdir(detector)
            pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

    if run_calwebb_spec2:
        # read the output header
        hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)

        # move the final reporting files to the working directory
        core_utils.move_latest_report_and_txt_2workdir(detector)

        # end the timer to compute the step running time of PTT
        PTT_end_time = time.time()
        core_utils.start_end_PTT_time(txt_name, start_time=None, end_time=PTT_end_time)

        return hdul, step_output_file, run_pytests

    else:
        if os.path.isfile(step_input_file):
            if run_pipe_step:
                print ("*** Step "+step+" set to True")
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
                print("Step "+step+" took "+end_time+" seconds to finish")

            else:
                print("Skipping running pipeline step ", step)
                # get the running time for this step
                end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, working_directory)

            # add the running time for this step
            step_completed = True
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)

            # get the total running time and print it in the file
            total_time = repr(core_utils.get_time_to_run_pipeline(txt_name))
            total_time_min = repr(round(float(total_time)/60.0, 2))
            print ("The total time for the pipeline to run was "+total_time+" seconds.")
            #print ("   ( = "+total_time_min+" minutes )")
            line2write = "{:<20} {:<20} {:<20} {:<20}".format('', '', 'total_time  ', total_time+'  ='+total_time_min+'min')
            print (line2write)
            with open(txt_name, "a") as tf:
                tf.write(line2write+"\n")

            # convert the html report into a pdf file
            #core_utils.convert_html2pdf()

            # end the timer to compute the step running time of PTT
            PTT_end_time = time.time()
            core_utils.start_end_PTT_time(txt_name, start_time=None, end_time=PTT_end_time)

            # move the final reporting files to the working directory
            core_utils.move_latest_report_and_txt_2workdir(detector)

            return hdul, step_output_file, run_pytests

        else:
            print (" The input file does not exist. Skipping step.")
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            #core_utils.convert_html2pdf()   # convert the html report into a pdf file
            # end the timer to compute the step running time of PTT
            PTT_end_time = time.time()
            core_utils.start_end_PTT_time(txt_name, start_time=None, end_time=PTT_end_time)
            # move the final reporting files to the working directory
            core_utils.move_latest_report_and_txt_2workdir(detector)
            # skip the test if input file does not exist
            pytest.skip("Skipping "+step+" because the input file does not exist.")


"""
# Move the files at the end of the pytests
@pytest.fixture(scope="function", autouse=True)
def move_output_files(request):
    def fin():
        #core_utils.convert_html2pdf()
        # move the final reporting files to the working directory
        core_utils.move_latest_report_and_txt_2workdir(detector)
        print("Output report files have been moved to the working_directory path indicated in the PTT_config file.")
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
        pytest.skip(msg)
    else:
        print("\n * Running reference file pytest...\n")
        result = extract_1d_utils.extract1d_rfile_is_correct(output_hdul)
        assert not result, result

def test_s_extr1d_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[2] = extract_1d_completion_tests, extract_1d_reffile_tests, extract_1d_validation_tests
    run_pytests = output_hdul[2][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running completion pytest...\n")
        assert extract_1d_utils.s_extr1d_exists(output_hdul[0]), "The keyword S_EXTR1D was not added to the header --> Extract 1D step was not completed."


