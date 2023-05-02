
"""
py.test module for unit testing the srctype step.
"""

import os
import time
import pytest
from glob import glob
from astropy.io import fits

from jwst.srctype.srctype_step import SourceTypeStep

from . import srctype_utils
from nirspec_pipe_testing_tool import core_utils
from nirspec_pipe_testing_tool.utils import change_filter_opaque2science


# HEADER
__author__ = "M. Pena-Guerrero"
__version__ = "1.4"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: introduced structure to separate completion from other tests if needed
# Apr 2019 - Version 1.2: implemented nptt_log capability
# Jun 2020 - Version 1.3: implemented APT user input vs pipeline logic
# Apr 2023 - Version 1.4: Cleaned-up code


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "srctype"
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
    run_pipe_step = config.getboolean("run_spec2_steps", step)
    # determine which tests are to be run
    run_pytests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # check if the filter is to be changed
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    if change_filter_opaque:
        is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file, step=step)
        if is_filter_opaque:
            filter_opaque_msg = "With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Source " \
                                "Type pytest now set to Skip."
            print(filter_opaque_msg)
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

    # only run this step if data is not BOTS
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

    if run_calwebb_spec2:
        outhdr = fits.getheader(step_output_file)
        scihdur = fits.getheader(step_output_file, 'SCI')
        return outhdr, step_output_file, run_pytests, scihdur, nptt_log

    else:
        if run_pipe_step:
            if os.path.isfile(step_input_file):
                # Create the pipeline step log
                stp_pipelog = "calspec2_" + step + "_" + detector + ".log"
                core_utils.mk_stpipe_log_cfg(output_dir, stp_pipelog)
                print("Pipeline step screen output will be logged in file: ", stp_pipelog)

                msg = " *** Step "+step+" set to True"
                print(msg)
                nptt_log.info(msg)
                stp = SourceTypeStep()

                # check that previous pipeline steps were run up to this point
                core_utils.check_completed_steps(step, step_input_file)

                # get the right configuration files to run the step
                local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                # start the timer to compute the step running time
                start_time = time.time()
                if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                else:
                    result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/srctype.cfg')
                result.save(step_output_file)
                # end the timer to compute the step running time
                end_time = repr(time.time() - start_time)   # this is in seconds
                msg = "Step "+step+" took "+end_time+" seconds to finish"
                print(msg)
                nptt_log.info(msg)

                step_completed = True
                outhdr = fits.getheader(step_output_file)
                scihdur = fits.getheader(step_output_file, 'SCI')

                # rename and move the pipeline log file
                pipelog = "pipeline_" + detector + ".log"
                try:
                    calspec2_pilelog = "calspec2_pipeline_" + step + "_" + detector + ".log"
                    pytest_workdir = TESTSDIR
                    logfile = glob(pytest_workdir + "/" + pipelog)[0]
                    os.rename(logfile, os.path.join(output_directory, calspec2_pilelog))
                except IndexError:
                    print("* WARNING: Something went wrong. Could not find a ", pipelog, " file ")

                # add the running time for this step
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                return outhdr, step_output_file, run_pytests, scihdur, nptt_log

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
                scihdur = fits.getheader(step_output_file, 'SCI')
                step_completed = True
                # add the running time for this step
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                return outhdr, step_output_file, run_pytests, scihdur, nptt_log
            else:
                step_completed = False
                # add the running time for this step
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")


# Unit tests

def test_s_srctyp_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    run_pytests = output_vars[2]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running completion pytest..."
        print(msg)
        nptt_log.info(msg)
        assert srctype_utils.s_srctyp_exists(output_vars[0]), "The keyword S_SRCTYP was not added to the header " \
                                                              "--> Srctype step was not completed."


def test_srctyp_logic(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    run_pytests = output_vars[2]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running source type logic pytest..."
        print(msg)
        nptt_log.info(msg)
        test_result, test_msg = srctype_utils.srctype_logic_correct(output_vars)
        assert test_result, test_msg

