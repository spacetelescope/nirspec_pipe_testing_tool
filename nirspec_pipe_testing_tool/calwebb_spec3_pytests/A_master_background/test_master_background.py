"""
py.test module for unit testing the assign_wcs step.
"""

import os
import subprocess
import time
import pytest
import logging
from astropy.io import fits
from glob import glob
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline
from jwst.master_background.master_background_step import MasterBackgroundStep

from . import master_background_utils
from nirspec_pipe_testing_tool import core_utils
from nirspec_pipe_testing_tool.calwebb_spec3_pytests import TESTSDIR

# print pipeline version
import jwst

pipeline_version = "\n *** Using jwst pipeline version: " + jwst.__version__ + " *** \n"
print(pipeline_version)

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"


# HISTORY
# Jun 2020 - Version 1.0: initial version completed


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "master_background"
    data_directory = config.get("calwebb_spec2_input_file", "data_directory")
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file_basename = config.get("calwebb_spec3", "s3_input_file")
    initial_input_file = os.path.join(data_directory, initial_input_file_basename)
    if os.path.isfile(initial_input_file):
        print("\n Taking initial input file from data_directory:")
    else:
        initial_input_file = os.path.join(output_directory, initial_input_file_basename)
        print("\n Taking initial file from output_directory: ")
    print(" Initial input file = ", initial_input_file, "\n")
    # Get the detector used
    detector = fits.getval(initial_input_file, "DETECTOR", 0)
    True_steps_suffix_map = "spec3_steps_suffix_map_" + detector + ".txt"
    pytests_directory = TESTSDIR
    True_steps_suffix_map = os.path.join(pytests_directory, True_steps_suffix_map)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, output_directory)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    # determine if the pipeline is to be run in full, per steps, or skipped
    run_calwebb_spec3 = config.get("calwebb_spec3", "run_calwebb_spec3")
    print("run_calwebb_spec3 = ", run_calwebb_spec3)
    if run_calwebb_spec3 == "skip":
        print('\n * PTT finished processing run_calwebb_spec3 is set to skip. \n')
        pytest.exit("Finished processing file, run_calwebb_spec3 is set to skip in configuration file.")
    else:
        run_calwebb_spec3 = bool(run_calwebb_spec3)

    # get the general info
    step, step_input_filename, output_file, in_file_suffix, outstep_file_suffix, True_steps_suffix_map = set_inandout_filenames
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    txt_name = os.path.join(output_directory, True_steps_suffix_map)
    step_input_file = os.path.join(output_directory, step_input_filename)
    step_output_file = os.path.join(output_directory, output_file)
    mode_used = config.get("calwebb_spec2_input_file", "mode_used").lower()

    # start the timer to compute the step running time of NPTT
    nptt_start_time = time.time()

    # determine which steps are to be run, if not run in full
    run_pipe_step = config.getboolean("run_pipe_steps", step)

    # determine which tests are to be run
    master_background_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    master_background_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    master_background_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [master_background_completion_tests, master_background_reffile_tests,
                   master_background_validation_tests]

    # Get the detector used
    detector = fits.getval(step_input_file, "DETECTOR", 0)

    # get main header from input file
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)

    # if run_calwebb_spec3 is True, calwebb_spec3 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # get the shutter configuration file for MOS data only
    msa_shutter_conf = "No shutter configuration file will be used."
    if core_utils.check_MOS_true(inhdu):
        msa_shutter_conf = config.get("esa_intermediary_products", "msa_conf_name")

        # check if the configuration shutter file name is in the header of the fits file and if not add it
        msametfl = fits.getval(step_input_file, "MSAMETFL", 0)
        if os.path.basename(msa_shutter_conf) != msametfl:
            msametfl = os.path.basename(msa_shutter_conf)
            fits.setval(step_input_file, "MSAMETFL", 0, value=msametfl)

        # copy the MSA shutter configuration file to this directory if not the working directory
        if os.getcwd() != os.path.dirname(msa_shutter_conf):
            # copy the MSA shutter configuration file into the pytest directory
            print("Removing MSA config file from: ", os.getcwd())
            subprocess.run(["cp", msa_shutter_conf, "."])

    # check if processing an image, then set proper variables
    if mode_used in ('image', 'confirm', 'taconfirm', 'wata', 'msata', 'bota', 'focus', 'mimf'):
        run_calwebb_spec3 = True
        imaging_mode = True
        print('\n * Image processing stops after calwebb_spec2. ')
        # end script for imaging case
        if imaging_mode:
            print('\n * NPTT finished processing imaging mode. \n')
            pytest.exit("Imaging does not get processed through calwebb_spec3.")

    # get the name of the configuration file and run the pipeline
    calwebb_spec3_cfg = config.get("calwebb_spec3", "calwebb_spec3_cfg")

    # copy the configuration file to create the pipeline log
    if not os.path.isfile(os.path.join(output_directory, "stpipe-log.cfg")):
        stpipelogcfg = calwebb_spec3_cfg.replace("calwebb_spec3.cfg", "stpipe-log.cfg")
        subprocess.run(["cp", stpipelogcfg, os.getcwd()])

    # run the pipeline
    if run_calwebb_spec3:

        # Create the logfile for NPTT, but remove the previous log file
        npttcalspec3_log = os.path.join(output_directory, 'NPTT_calspec3_' + detector + '.log')
        if os.path.isfile(npttcalspec3_log):
            os.remove(npttcalspec3_log)
        print("Spec3 screen information output from NPTT will be logged in file: ", npttcalspec3_log)
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(filename=npttcalspec3_log, level=logging.INFO)
        logging.info(pipeline_version)

        run_calwebb_spec3_msg = " *** Will run pipeline in full ... "
        print(run_calwebb_spec3_msg)
        logging.info(run_calwebb_spec3_msg)

        # create the map
        txt_name = "spec3_full_run_map_" + detector + ".txt"
        if os.path.isfile(txt_name):
            os.remove(txt_name)
        master_background_utils.create_completed_steps_txtfile(txt_name, step_input_file)

        # start the timer to compute the step running time of NPTT
        core_utils.start_end_PTT_time(txt_name, start_time=nptt_start_time, end_time=None)

        if mode_used == "bots":
            calwebb_spec3_cfg = calwebb_spec3_cfg.replace("calwebb_spec3.cfg", "calwebb_tso-spec3.cfg")
            print('\nUsing the following configuration file to run TSO pipeline:')
        else:
            print('\nUsing the following configuration file to run spectroscopic pipeline:')
        print(calwebb_spec3_cfg, '\n')

        # start the timer to compute the step running time
        start_time = time.time()

        # run the pipeline
        print('Running pipeline... \n')
        Spec3Pipeline.call(step_input_file, config_file=calwebb_spec3_cfg)

        # end the timer to compute calwebb_spec3 running time
        end_time = repr(time.time() - start_time)  # this is in seconds
        calspec3_time = " * Pipeline took " + end_time + " seconds to finish.\n"
        print(calspec3_time)
        logging.info(calspec3_time)

        # add the detector string to the name of the files and move them to the working directory
        core_utils.add_detector2filename(output_directory, step_input_file)

        # state name of the final spec3 _cal file
        final_output_name = step_input_file.replace("spec2", "spec3")
        final_output_name_msg = "\nThe final pipeline product was saved in: " + final_output_name
        print(final_output_name_msg)
        logging.info(final_output_name_msg)

        # read the assign wcs fits file
        hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
        # scihdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False, ext=1)

        if core_utils.check_MOS_true(inhdu):
            if os.getcwd() != os.path.dirname(msa_shutter_conf):
                # remove the copy of the MSA shutter configuration file
                print("Removing MSA config file from: ", os.getcwd())
                subprocess.run(["rm", msametfl])

        # rename and move the pipeline log file
        calspec3_pipelog = "calspec3_pipeline_" + detector + ".log"
        try:
            path_where_pipeline_was_run = os.getcwd()
            logfile = glob(path_where_pipeline_was_run + "/pipeline.log")[0]
            print(logfile)
            os.rename(logfile, os.path.join(output_directory, calspec3_pipelog))
        except IndexError:
            print("\nWARNING: Something went wrong. Could not find a pipeline.log file \n")

        # make sure we are able to find calspec3_pipelog either in the calwebb_spec3 directory or in the working dir
        if not os.path.isfile(calspec3_pipelog):
            calspec3_pipelog = os.path.join(output_directory, calspec3_pipelog)

        # add the running time for all steps
        step_running_times = core_utils.calculate_step_run_time(calspec3_pipelog)
        end_time_list = []
        for stp in core_utils.step_string_dict:
            if stp in step_running_times:
                step_completed = True
                step_time = step_running_times[stp]["run_time"]
                out_suffix = core_utils.step_string_dict[stp]["suffix"]
                core_utils.add_completed_steps(txt_name, stp, out_suffix, step_completed, step_time)
                end_time_list.append(step_time)

        # print total running time in the text file and move it to the indicated directory
        string2print = "pipeline_total_time"
        if float(end_time) <= sum(end_time_list):
            tot_time = repr(sum(end_time_list))
        else:
            tot_time = end_time
        master_background_utils.print_time2file(txt_name, tot_time, string2print)
        nptt_runtimes_msg = "Pipeline and NPTT run times written in file: " + os.path.basename(
            txt_name) + " in working directory. \n"
        print(nptt_runtimes_msg)
        logging.info(nptt_runtimes_msg)

        # move the final reporting text files to the working directory
        if os.getcwd() != output_directory:
            core_utils.move_txt_files_2workdir(config, detector)

        return hdul, step_output_file, msa_shutter_conf, run_pytests, mode_used

    else:

        # create the map but remove a previous one if it exists
        if os.path.isfile(txt_name):
            os.remove(txt_name)
        master_background_utils.create_completed_steps_txtfile(txt_name, step_input_file)

        # start the timer to compute the step running time of NPTT
        core_utils.start_end_PTT_time(txt_name, start_time=nptt_start_time, end_time=None)
        msg = "\n Pipeline and NPTT run times will be written in file: " + os.path.basename(
            txt_name) + " in working directory. \n"
        print(msg)
        logging.info(msg)

        if run_pipe_step:

            # Create the logfile for NPTT, but erase the previous one if it exists
            npttcalspec3_log = os.path.join(output_directory, 'NPTT_calspec3_' + detector + '_' + step + '.log')
            if os.path.isfile(npttcalspec3_log):
                os.remove(npttcalspec3_log)
            print("Output information on screen from NPTT will be logged in file: ", npttcalspec3_log)
            for handler in logging.root.handlers[:]:
                logging.root.removeHandler(handler)
            logging.basicConfig(filename=npttcalspec3_log, level=logging.INFO)
            logging.info(pipeline_version)

            # check that previous pipeline steps were run up to this point
            core_utils.check_completed_steps(step, step_input_file)

            if os.path.isfile(step_input_file):
                msg = " *** Step " + step + " set to True"
                print(msg)
                logging.info(msg)
                stp = MasterBackgroundStep()

                # get the right configuration files to run the step
                local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")

                # start the timer to compute the step running time
                print("running pipeline...")
                start_time = time.time()
                if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                else:
                    result = stp.call(step_input_file, config_file=local_pipe_cfg_path + '/master_background.cfg')
                result.save(step_output_file)

                # end the timer to compute the step running time
                end_time = repr(time.time() - start_time)  # this is in seconds
                msg = "Step " + step + " took " + end_time + " seconds to finish"
                print(msg)
                logging.info(msg)

                if core_utils.check_MOS_true(inhdu):
                    # remove the copy of the MSA shutter configuration file
                    if os.getcwd() != os.path.dirname(msa_shutter_conf):
                        print("Removing MSA config file from: ", os.getcwd())
                        subprocess.run(["rm", msametfl])

                # rename and move the pipeline log file
                pipelog = "pipeline_" + detector + ".log"
                try:
                    calspec3_pipelog = "calspec3_pipeline_" + step + "_" + detector + ".log"
                    pytest_workdir = TESTSDIR
                    logfile = glob(pytest_workdir + "/" + pipelog)[0]
                    os.rename(logfile, os.path.join(output_directory, calspec3_pipelog))
                except IndexError:
                    print("\n* WARNING: Something went wrong. Could not find a ", pipelog, " file \n")

            else:
                msg = "Skipping step. Input file " + step_input_file + " does not exit."
                print(msg)
                logging.info(msg)
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping " + step + " because the input file does not exist.")

        else:
            print("Skipping running pipeline step ", step)
            # add the running time for this step
            end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, output_directory)

        if os.path.isfile(step_output_file):
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            step_completed = True
            # add the running time for this step
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            return hdul, step_output_file, msa_shutter_conf, run_pytests, mode_used
        else:
            step_completed = False
            # add the running time for this step
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Test skipped because input file "+step_output_file+" does not exist.")


"""
# fixture to validate 
@pytest.fixture(scope="module")
def validate_wcs(output_hdul):
    # get the input information for the wcs routine
    hdu = output_hdul[0]
    infile_name = output_hdul[1]
    msa_conf_name = output_hdul[2]
    esa_files_path = output_hdul[3]
    mode_used = output_hdul[7]

    # define the threshold difference between the pipeline output and the ESA files for the pytest to pass or fail
    threshold_diff = float(output_hdul[4])

    # save the output plots
    save_wcs_plots = output_hdul[5]

    # show the figures
    show_figs = False

    msg = "\n Performing WCS validation test... "
    print(msg)
    logging.info(msg)
    log_msgs = None
    if core_utils.check_FS_true(hdu):
        result, log_msgs = compare_wcs_fs.compare_wcs(infile_name, esa_files_path=esa_files_path, show_figs=show_figs,
                                                      save_figs=save_wcs_plots, threshold_diff=threshold_diff,
                                                      raw_data_root_file=None,
                                                      output_directory=None,
                                                      debug=False)

    elif core_utils.check_MOS_true(hdu) and mode_used != "MOS_sim":
        result, log_msgs = compare_wcs_mos.compare_wcs(infile_name, esa_files_path=esa_files_path,
                                                       msa_conf_name=msa_conf_name,
                                                       show_figs=show_figs, save_figs=save_wcs_plots,
                                                       threshold_diff=threshold_diff,
                                                       raw_data_root_file=None,
                                                       output_directory=None, debug=False)

    elif core_utils.check_IFU_true(hdu):
        result, log_msgs = compare_wcs_ifu.compare_wcs(infile_name, esa_files_path=esa_files_path, show_figs=show_figs,
                                                       save_figs=save_wcs_plots, threshold_diff=threshold_diff,
                                                       raw_data_root_file=None,
                                                       output_directory=None,
                                                       debug=False)

    else:
        # We do not have ESA data to compare with for BOTS
        pytest.skip("Skipping pytest: The fits file is not FS, MOS, or IFU.")

    if log_msgs is not None:
        for msg in log_msgs:
            logging.info(msg)

    if "skip" in result:
        pytest.skip("Skipping assign_wcs pytest.")
    elif "PASS" in result:
        result = True
    else:
        result = False

    return result
"""

# reference file tests

"""
def test_wavran_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest...\n"
        print(msg)
        logging.info(msg)
        result = master_background_utils.wavran_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]
"""

# keyword tests


def test_masterbg_exists(output_hdul):
    # want to run this pytest?
    # output_hdul = hdul, step_output_file, msa_shutter_conf, run_pytests, mode_used
    # run_pytests = master_background_completion_tests, master_background_reffile_tests,
    #               master_background_validation_tests
    run_pytests = output_hdul[3][0]
    if not run_pytests:
        msg = "Skipping pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest...\n"
        print(msg)
        logging.info(msg)
        assert master_background_utils.masterbg_exists(output_hdul[1]), "The keyword MASTERBG was not added to " \
                                                                        "the header."


# validation test

"""
def test_validate_wcs(output_hdul, request):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][2]
    if not run_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running validation pytest...\n"
        print(msg)
        logging.info(msg)
        # print("\n", validate_wcs, "\n")
        assert request.getfixturevalue('validate_wcs'), "Output value from compare_wcs.py is greater than threshold."
"""

