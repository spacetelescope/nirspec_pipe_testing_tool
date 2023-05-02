
"""
py.test module for unit testing the flat field step.
"""

import os
import subprocess
import time
import pytest
from glob import glob
from astropy.io import fits

from jwst.flatfield.flat_field_step import FlatFieldStep

from . import flat_field_utils
from nirspec_pipe_testing_tool import core_utils
from .. auxiliary_code import flattest_fs
from .. auxiliary_code import flattest_ifu
from .. auxiliary_code import flattest_mos
from nirspec_pipe_testing_tool.utils import change_filter_opaque2science



# HEADER
__author__ = "M. Pena-Guerrero"
__version__ = "1.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Mar 2019 - Version 1.1: modified reference file tests and separated completion from validation tests
# Apr 2019 - Version 1.2: implemented nptt_log capability
# Apr 2023 - Version 1.3: Cleaned-up code


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "flat_field"
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
    msa_shutter_conf = config.get("benchmark_intermediary_products", "msa_conf_name")
    msametfl = os.path.basename(msa_shutter_conf)
    dflat_path = config.get("benchmark_intermediary_products", "dflat_path")
    sflat_path = config.get("benchmark_intermediary_products", "sflat_path")
    fflat_path = config.get("benchmark_intermediary_products", "fflat_path")
    flattest_threshold_diff = config.get("additional_arguments", "flattest_threshold_diff")
    save_flattest_plot = config.getboolean("additional_arguments", "save_flattest_plot")
    write_flattest_files = config.getboolean("additional_arguments", "write_flattest_files")
    flattest_paths = [step_output_file, msa_shutter_conf, dflat_path, sflat_path, fflat_path]
    flattest_switches = [flattest_threshold_diff, save_flattest_plot, write_flattest_files]
    run_pipe_step = config.getboolean("run_spec2_steps", step)
    # determine which tests are to be run
    flat_field_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    flat_field_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    flat_field_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [flat_field_completion_tests, flat_field_reffile_tests, flat_field_validation_tests]

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # check if the filter is to be changed
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    if change_filter_opaque:
        is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file,
                                                                                                  step=step)
        if is_filter_opaque:
            filter_opaque_msg = "With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Flat " \
                                "Field pytest now set to Skip."
            print(filter_opaque_msg)
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

    # get the MSA shutter configuration file full path only for MOS data
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
        flattest_paths = [step_output_file, msa_shutter_conf, dflat_path, sflat_path, fflat_path]
        return outhdr, step_output_file, flattest_paths, flattest_switches, run_pytests, nptt_log
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
                stp = FlatFieldStep()

                # check that previous pipeline steps were run up to this point
                core_utils.check_completed_steps(step, step_input_file)

                # start the timer to compute the step running time
                ontheflyflat = step_output_file.replace("flat_field.fits", "interpolatedflat.fits")
                msg1 = "Step product will be saved as: "+step_output_file
                msg2 = "on-the-fly flat will be saved as: "+ontheflyflat
                # get the right configuration files to run the step
                print(msg1)
                print(msg2)
                nptt_log.info(msg1)
                nptt_log.info(msg2)
                local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                # start the timer to compute the step running time
                start_time = time.time()
                if local_pipe_cfg_path == "pipe_source_tree_code":
                    stp.call(step_input_file, output_file=step_output_file, save_interpolated_flat=True)
                else:
                    stp.call(step_input_file, output_file=step_output_file, save_interpolated_flat=True,
                             config_file=local_pipe_cfg_path+'/flat_field.cfg')
                # end the timer to compute the step running time
                end_time = repr(time.time() - start_time)   # this is in seconds
                msg = "Step "+step+" took "+end_time+" seconds to finish"
                print(msg)
                nptt_log.info(msg)
                # move the on-the-fly flat to the working directory
                subprocess.run(["mv", os.path.basename(ontheflyflat), ontheflyflat])
                # raname and move the flat_field output
                subprocess.run(["mv", os.path.basename(step_output_file).replace("_flat_field.fits",
                                                                                 "_flatfieldstep.fits"),
                                step_output_file])

                step_completed = True
                outhdr = fits.getheader(step_output_file)

                # add the running time for this step
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                return outhdr, step_output_file, flattest_paths, flattest_switches, run_pytests, nptt_log

            else:
                msg = " The input file does not exist. Skipping step."
                print(msg)
                nptt_log.info(msg)
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping " + step + " because the input file does not exist.")

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
                return outhdr, step_output_file, flattest_paths, flattest_switches, run_pytests, nptt_log
            else:
                step_completed = False
                # add the running time for this step
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")


# fixture to validate the flat field step
@pytest.fixture(scope="module")
def validate_flat_field(output_vars):
    hdr = output_vars[0]
    step_output_file, msa_shutter_conf, dflat_path, sflat_path, fflat_path = output_vars[2]
    flattest_threshold_diff, save_flattest_plot, write_flattest_files = output_vars[3]

    # show the figures
    show_figs = False

    # get the logger instance
    nptt_log = output_vars[-1]
    log_msgs = None

    # determine the testing mode
    do_fs_flattest = False
    if core_utils.check_FS_true(hdr) or core_utils.check_BOTS_true(hdr):
        do_fs_flattest = True
    elif 'bots' in step_output_file.lower() or 'fs' in step_output_file.lower():
        do_fs_flattest = True
    do_mos_flattest = False
    if core_utils.check_MOS_true(hdr) or 'mos' in step_output_file.lower():
        do_mos_flattest = True
    do_ifu_flattest = False
    if core_utils.check_IFU_true(hdr) or 'ifu' in step_output_file.lower():
        do_ifu_flattest = True

    # run the test
    if do_fs_flattest:
        median_diff, result_msg, log_msgs = flattest_fs.flattest(step_output_file, dflat_path=dflat_path,
                                                                 sflat_path=sflat_path, fflat_path=fflat_path,
                                                                 writefile=write_flattest_files,
                                                                 show_figs=show_figs, save_figs=save_flattest_plot,
                                                                 interpolated_flat=None,
                                                                 threshold_diff=flattest_threshold_diff,
                                                                 output_directory=None, debug=False)

    elif do_mos_flattest:
        median_diff, result_msg, log_msgs = flattest_mos.flattest(step_output_file, dflat_path=dflat_path,
                                                                  sflat_path=sflat_path, fflat_path=fflat_path,
                                                                  msa_shutter_conf=msa_shutter_conf,
                                                                  writefile=write_flattest_files,
                                                                  show_figs=show_figs, save_figs=save_flattest_plot,
                                                                  interpolated_flat=None,
                                                                  threshold_diff=flattest_threshold_diff,
                                                                  debug=False)

    elif do_ifu_flattest:
        median_diff, result_msg, log_msgs = flattest_ifu.flattest(step_output_file, dflat_path=dflat_path,
                                                                  sflat_path=sflat_path, fflat_path=fflat_path,
                                                                  writefile=write_flattest_files,
                                                                  mk_all_slices_plt=False, show_figs=show_figs,
                                                                  save_figs=save_flattest_plot,
                                                                  interpolated_flat=None,
                                                                  threshold_diff=flattest_threshold_diff,
                                                                  debug=False)

    else:
        pytest.skip("Skipping pytest: The input fits file is not FS, MOS, or IFU. This tool does not yet include the "
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


def test_s_flat_exists(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[4] = flat_field_completion_tests, flat_field_reffile_tests, flat_field_validation_tests
    run_pytests = output_vars[4][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = "* Running completion pytest..."
        print(msg)
        nptt_log.info(msg)
        assert flat_field_utils.s_flat_exists(output_vars[0]), "The keyword S_FLAT was not added to the header --> " \
                                                               "flat_field step was not completed."


def test_validate_flat_field(output_vars, request):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[4] = flat_field_completion_tests, flat_field_reffile_tests, flat_field_validation_tests
    run_pytests = output_vars[4][2]
    if not run_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = "* Running validation pytest..."
        print(msg)
        nptt_log.info(msg)
        assert request.getfixturevalue("validate_flat_field"), "Output value from flat test is greater than threshold."


def test_fflat_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[4] = flat_field_completion_tests, flat_field_reffile_tests, flat_field_validation_tests
    run_pytests = output_vars[4][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = "* Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = flat_field_utils.fflat_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_sflat_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[4] = flat_field_completion_tests, flat_field_reffile_tests, flat_field_validation_tests
    run_pytests = output_vars[4][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = "* Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = flat_field_utils.sflat_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]


def test_dflat_rfile(output_vars):
    # get the logger instance
    nptt_log = output_vars[-1]
    # want to run this pytest?
    # output_vars[4] = flat_field_completion_tests, flat_field_reffile_tests, flat_field_validation_tests
    run_pytests = output_vars[4][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in NPTT_config.cfg file."
        print(msg)
        nptt_log.info(msg)
        pytest.skip(msg)
    else:
        msg = "* Running reference file pytest..."
        print(msg)
        nptt_log.info(msg)
        result = flat_field_utils.dflat_rfile_is_correct(output_vars)
        for log_msg in result[1]:
            print(log_msg)
            nptt_log.info(log_msg)
        assert not result[0], result[0]

