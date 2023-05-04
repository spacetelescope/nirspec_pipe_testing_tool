
"""
py.test module for unit testing the pathloss step.
"""

import os
import time
import pytest
from glob import glob
from astropy.io import fits
import urllib.request

from jwst.pathloss.pathloss_step import PathLossStep

from . import pathloss_utils
from nirspec_pipe_testing_tool import core_utils
from nirspec_pipe_testing_tool.utils import change_filter_opaque2science
from .. auxiliary_code import pathloss_fs_ps
from .. auxiliary_code import pathloss_fs_uni
from .. auxiliary_code import pathloss_ifu_ps
from .. auxiliary_code import pathloss_ifu_uni
from .. auxiliary_code import pathloss_mos


# HEADER
__author__ = "M. Pena-Guerrero & G. Kanarek"
__version__ = "2.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check
# Mar 2019 - Version 2.1: Maria separated completion from validation tests
# Apr 2021 - Version 2.2: combined point-uniform pathloss test for MOS data
# Apr 2023 - Version 2.3: Cleaned-up code


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "pathloss"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_vars(set_inandout_filenames, config):
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
    set_inandout_filenames_info = core_utils.read_info4output_vars(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, outstep_file_suffix = set_inandout_filenames_info
    run_pipe_step = config.getboolean("run_spec2_steps", step)
    # determine which tests are to be run
    pathloss_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    pathloss_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    pathloss_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [pathloss_completion_tests, pathloss_reffile_tests, pathloss_validation_tests]
    pathloss_threshold_diff = config.get("additional_arguments", "pathloss_threshold_diff")
    save_flattest_plot = config.getboolean("additional_arguments", "save_pathloss_plot")
    write_flattest_files = config.getboolean("additional_arguments", "write_pathloss_files")
    pathloss_switches = [pathloss_threshold_diff, save_flattest_plot, write_flattest_files]

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # check if the filter is to be changed
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    if change_filter_opaque:
        is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file,
                                                                                                  step=step)
        if is_filter_opaque:
            filter_opaque_msg = "With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. " \
                                "Pathloss pytest now set to Skip."
            print(filter_opaque_msg)
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

    # only run this step if data is not BOTS
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(output_directory, initial_input_file)
    if os.path.isfile(initial_input_file):
        inhdr = fits.getheader(initial_input_file)
        detector = inhdr["DETECTOR"]
    else:
        msg = "Skipping "+step+" because the initial input file given in NPTT_config.cfg does not exist."
        pytest.skip(msg)

    # Get the logfile instance for NPTT created in the run.py script
    nptt_log = os.path.join(output_directory, 'NPTT_calspec2_' + detector + '.log')
    nptt_log = core_utils.mk_nptt_log(nptt_log, reset=False)

    if not core_utils.check_BOTS_true(inhdr):
        if run_calwebb_spec2:
            outhdr = fits.getheader(step_output_file)
            return outhdr, step_output_file, step_input_file, run_pytests, pathloss_switches, nptt_log
        else:
            if run_pipe_step:
                if os.path.isfile(step_input_file):
                    if change_filter_opaque:
                        nptt_log.info(filter_opaque_msg)

                    # Create the pipeline step log
                    stp_pipelog = "calspec2_" + step + "_" + detector + ".log"
                    core_utils.mk_stpipe_log_cfg(output_dir, stp_pipelog)
                    print("Pipeline step screen output will be logged in file: ", stp_pipelog)

                    msg = " *** Step "+step+" set to True"
                    print(msg)
                    nptt_log.info(msg)
                    stp = PathLossStep()

                    # check that previous pipeline steps were run up to this point
                    core_utils.check_completed_steps(step, step_input_file)

                    # get the right configuration files to run the step
                    local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")

                    # start the timer to compute the step running time
                    start_time = time.time()
                    if local_pipe_cfg_path == "pipe_source_tree_code":
                        result = stp.call(step_input_file)
                    else:
                        result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/pathloss.cfg')
                    result.save(step_output_file)

                    # For the moment, the pipeline is using the wrong reference file for slit 400A1, so the
                    # file needs to be re-processed with the right reference file BUT we are
                    # skipping this slit
                    """
                    if core_utils.check_FS_true(inhdr):
                        print("\n * WARNING: For the moment, the wrong reference file is being used for "
                              "processing slit 400A1. The file will be re-processed ")
                        # The file can manually be re-run with
                        # $ strun jwst.pathloss.PathLossStep final_output_caldet1_NRS1_srctype.fits
                        #                                                --override_pathloss=jwst-nirspec-a400.plrf.fits
                        
                        pathloss_400a1 = step_input_file.replace("srctype.fits", "pathloss_400A1.fits")
                        reffile_400a1 = "jwst-nirspec-a400.plrf.fits"
                        print("Re-processing slit with new reference file: ", reffile_400a1)
                        pl = PathLossStep()
                        pl.override_pathloss = reffile_400a1
                        pl_result = pl.run(step_input_file)
                        pl_result.save(pathloss_400a1)
                        #import subprocess
                        #subprocess.run(["mv", step_input_file.replace("srctype", "pathlossstep"), pathloss_400a1])
                        print("Saved pipeline re-processed file as: ", pathloss_400a1)
                    """

                    # end the timer to compute calwebb_spec2 running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    msg = " * calwebb_spec2 took "+end_time+" seconds to finish."
                    print(msg)
                    nptt_log.info(msg)
                    step_completed = True
                    outhdr = fits.getheader(step_output_file)

                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return outhdr, step_output_file, step_input_file, run_pytests, pathloss_switches, nptt_log

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
                    step_completed = True
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    return outhdr, step_output_file, step_input_file, run_pytests, pathloss_switches, nptt_log
                else:
                    step_completed = False
                    # add the running time for this step
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")

    else:
        pytest.skip("Skipping "+step+" because data is BOTS.")


# fixture to validate the pathloss step
@pytest.fixture(scope="module")
def validate_pathloss(output_vars):
    hdr = output_vars[0]
    step_input_filename = output_vars[2]
    comparison_filename = output_vars[1]
    threshold_diff, save_figs, writefile = output_vars[4]

    # other default variables
    show_figs = False
    log_msgs = None
    debug = False

    # get the logger instance
    nptt_log = output_vars[-1]

    # determine the type of source
    is_mos = True
    if not core_utils.check_MOS_true(hdr):
        is_mos = False
        source_type = fits.getval(comparison_filename, "SRCTYPE", "SCI", 1)
        msg = "Source type is: "+source_type
        print(msg)
        nptt_log.info(msg)

    # get the corresponding reference file
    reffile = hdr["R_PTHLOS"].replace("crds://", "")
    # download the file if necessary
    if not os.path.isfile(reffile):
        reffile_url = "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/" + reffile
        urllib.request.urlretrieve(reffile_url, reffile)
        print('  Pathloss reference file retreved from CRDS!')

    if is_mos:
        median_diff, result_msg, log_msgs = pathloss_mos.pathtest(step_input_filename, reffile,
                                                                  comparison_filename,
                                                                  writefile=writefile, show_figs=show_figs,
                                                                  save_figs=save_figs,
                                                                  threshold_diff=threshold_diff,
                                                                  debug=debug)

    else:
        if core_utils.check_FS_true(hdr) or core_utils.check_BOTS_true(hdr):
            if "point" in source_type.lower():
                median_diff, result_msg, log_msgs = pathloss_fs_ps.pathtest(step_input_filename, reffile,
                                                                            comparison_filename,
                                                                            writefile=writefile, show_figs=show_figs,
                                                                            save_figs=save_figs,
                                                                            threshold_diff=threshold_diff,
                                                                            debug=debug)
            elif "extend" in source_type.lower():
                median_diff, result_msg, log_msgs = pathloss_fs_uni.pathtest(step_input_filename, reffile,
                                                                             comparison_filename,
                                                                             writefile=writefile, show_figs=show_figs,
                                                                             save_figs=save_figs,
                                                                             threshold_diff=threshold_diff,
                                                                             debug=debug)
        elif core_utils.check_IFU_true(hdr):
            if "point" in source_type.lower():
                median_diff, result_msg, log_msgs = pathloss_ifu_ps.pathtest(step_input_filename, reffile,
                                                                             comparison_filename,
                                                                             writefile=writefile, show_figs=show_figs,
                                                                             save_figs=save_figs,
                                                                             threshold_diff=threshold_diff,
                                                                             debug=debug)
            elif "extend" in source_type.lower():
                median_diff, result_msg, log_msgs = pathloss_ifu_uni.pathtest(step_input_filename, reffile,
                                                                              comparison_filename,
                                                                              writefile=writefile, show_figs=show_figs,
                                                                              save_figs=save_figs,
                                                                              threshold_diff=threshold_diff,
                                                                              debug=debug)

        else:
            pytest.skip("Skipping pytest: The input fits file is not FS, MOS, or IFU. This tool does not include the "
                        "routine to verify this kind of file.")

    if log_msgs is not None:
        for msg in log_msgs:
            nptt_log.info(msg)

    if median_diff == "skip":
        nptt_log.info(result_msg)
        pytest.skip(result_msg)
    else:
        print(result_msg)
        nptt_log.info(result_msg)

    return median_diff


# Unit tests

def test_s_pthlos_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[2] = pathloss_completion_tests, pathloss_reffile_tests, pathloss_validation_tests
    run_pytests = output_vars[3][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running completion pytest..."
        print(msg)
        nptt_log.info(msg)
        assert pathloss_utils.s_pthlos_exists(output_vars[0]), "The keyword S_PTHLOS was not added to the header --> " \
                                                               "Pathloss step was not completed."


"""
def test_r_pthlos_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[2] = pathloss_completion_tests, pathloss_reffile_tests, pathloss_validation_tests
    run_pytests = output_vars[3][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        assert pathloss_utils.r_pthlos_exists(output_vars[0]), "The keyword R_PTHLOS was not added to the header --> " \
                                                                Not sure what reference file was used."
"""


def test_validate_pathloss(output_vars, request):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[2] = pathloss_completion_tests, pathloss_reffile_tests, pathloss_validation_tests
    run_pytests = output_vars[3][2]
    if not run_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running validation pytest..."
        print(msg)
        nptt_log.info(msg)
        assert request.getfixturevalue("validate_pathloss"), "Output value from pathloss is greater than threshold."


def test_pthlos_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[2] = pathloss_completion_tests, pathloss_reffile_tests, pathloss_validation_tests
    run_pytests = output_vars[3][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = " * Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = pathloss_utils.pthlos_rfile_is_correct(output_vars)
        if result is not None:
            for log_msg in result[1]:
                print(log_msg)
                nptt_log.info(log_msg)
            assert not result[0], result[0]
