
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
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_image2 import Image2Pipeline

from . import assign_wcs_utils
from .. import core_utils
from .. auxiliary_code import change_filter_opaque2science
from .. auxiliary_code import compare_wcs_ifu
from .. auxiliary_code import compare_wcs_fs
from .. auxiliary_code import compare_wcs_mos


# print pipeline version
import jwst
pipeline_version = "\n *** Using jwst pipeline version: "+jwst.__version__+" *** \n"
print(pipeline_version)


# HEADER
__author__ = "M. A. Pena-Guerrero & Gray Kanarek"
__version__ = "2.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Gray added routine to generalize reference file check
# Feb 2019 - Version 2.1: Maria made changes to be able to process 491 and 492 files in the same directory
# Apr 2019 - Version 2.2: implemented logging capability
# Dec 2019 - Version 2.3: implemented image processing and text file name handling



# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(config):
    step = "assign_wcs"
    data_directory = config.get("calwebb_spec2_input_file", "data_directory")
    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
    initial_input_file_basename = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(data_directory, initial_input_file_basename)
    if os.path.isfile(initial_input_file):
        print("\n Taking initial input file from data_directory:")
    else:
        initial_input_file = os.path.join(working_directory, initial_input_file_basename)
        print("\n Taking initial file from working_directory: ")
    print(" Initial input file = ", initial_input_file , "\n")
    # Get the detector used
    detector = fits.getval(initial_input_file, "DETECTOR", 0)
    True_steps_suffix_map = "True_steps_suffix_map_"+detector+".txt"
    pytests_directory = os.getcwd()
    True_steps_suffix_map = os.path.join(pytests_directory, True_steps_suffix_map)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, working_directory)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info

    # start the timer to compute the step running time of PTT
    PTT_start_time = time.time()

    # check if the filter is to be changed
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    if change_filter_opaque:
        _, step_input_file = change_filter_opaque2science.change_filter_opaque(step_input_file)
        change_filter_opaque_msg = " * With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Further steps will be skipped. \n"
        print(change_filter_opaque_msg)

    # determine if the pipeline is to be run in full
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    # determine which steps are to be run, if not run in full
    run_pipe_step = config.getboolean("run_pipe_steps", step)
    # determine which tests are to be run
    assign_wcs_completion_tests = config.getboolean("run_pytest", "_".join((step, "completion", "tests")))
    assign_wcs_reffile_tests = config.getboolean("run_pytest", "_".join((step, "reffile", "tests")))
    assign_wcs_validation_tests = config.getboolean("run_pytest", "_".join((step, "validation", "tests")))
    run_pytests = [assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests]

    # get other relevant info from PTT config file
    esa_files_path = config.get("esa_intermediary_products", "esa_files_path")
    wcs_threshold_diff = config.get("additional_arguments", "wcs_threshold_diff")
    save_wcs_plots = config.getboolean("additional_arguments", "save_wcs_plots")
    working_directory = config.get("calwebb_spec2_input_file", "working_directory")

    # Get the detector used
    detector = fits.getval(step_input_file, "DETECTOR", 0)

    # get main header from input file
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # Check if data is IFU that the Image Model keyword is correct
    mode_used = config.get("calwebb_spec2_input_file", "mode_used")
    if mode_used == "IFU":
        DATAMODL = fits.getval(step_input_file, "DATAMODL", 0)
        if DATAMODL != "IFUImageModel":
            fits.setval(step_input_file, "DATAMODL", 0, value="IFUImageModel")
            print("DATAMODL keyword changed to IFUImageModel.")

    # get the shutter configuration file for MOS data only
    msa_shutter_conf = "No shutter configuration file will be used."
    if core_utils.check_MOS_true(inhdu):
        msa_shutter_conf = config.get("esa_intermediary_products", "msa_conf_name")

        # check if the configuration shutter file name is in the header of the fits file and if not add it
        msametfl = fits.getval(step_input_file, "MSAMETFL", 0)
        if os.path.basename(msa_shutter_conf) != msametfl:
            msametfl = os.path.basename(msa_shutter_conf)
            fits.setval(step_input_file, "MSAMETFL", 0, value=msametfl)

    # check if processing an image, then set proper variables
    imaging_mode = False
    if mode_used in ('image', 'confirm', 'taconfirm', 'wata', 'msata', 'bota', 'focus', 'mimf'):
        run_calwebb_spec2 = True
        imaging_mode = True
        print('\n * Image processing will only be run in full with PTT. All intermediary products will be saved.')
        print('     ->  For now, all pytests will be skipped since there are now image-specific routines yet. \n')
        #TODO: add imaging tests

    # run the pipeline
    if run_calwebb_spec2:

        # Create the logfile for PTT, but remove the previous log file
        PTTcalspec2_log = os.path.join(working_directory, 'PTT_calspec2_'+detector+'.log')
        if imaging_mode:
            PTTcalspec2_log = PTTcalspec2_log.replace('calspec2', 'calimage2')
        if os.path.isfile(PTTcalspec2_log):
            os.remove(PTTcalspec2_log)
        print("Information outputed to screen from PTT will be logged in file: ", PTTcalspec2_log)
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        logging.basicConfig(filename=PTTcalspec2_log, level=logging.INFO)
        logging.info(pipeline_version)
        if change_filter_opaque:
            logging.info(change_filter_opaque_msg)

        run_calwebb_spec2_msg = " *** Will run calwebb_spec2 in full ... "
        print(run_calwebb_spec2_msg)
        logging.info(run_calwebb_spec2_msg)

        # create the map
        txt_name = "full_run_map_"+detector+".txt"
        if os.path.isfile(txt_name):
            os.remove(txt_name)
        assign_wcs_utils.create_completed_steps_txtfile(txt_name, step_input_file)

        # start the timer to compute the step running time of PTT
        core_utils.start_end_PTT_time(txt_name, start_time=PTT_start_time, end_time=None)

        # get the name of the configuration file and run the pipeline
        calwebb_spec2_cfg = config.get("run_calwebb_spec2_in_full", "calwebb_spec2_cfg")
        calwebb_image2_cfg = False
        if mode_used == "BOTS":
            calwebb_spec2_cfg = calwebb_spec2_cfg.replace("calwebb_spec2.cfg", "calwebb_tso_spec2.cfg")
            print('\nUsing the following configuration file to run TSO pipeline:')
            print(calwebb_spec2_cfg, '\n')
        elif imaging_mode:
            calwebb_image2_cfg = calwebb_spec2_cfg.replace("calwebb_spec2.cfg", "calwebb_image2.cfg")
            print('\nUsing the following configuration file to run imaging pipeline:')
            print(calwebb_image2_cfg, '\n')
        else:
            print('\nUsing the following configuration file to run spectroscopic pipeline:')
            print(calwebb_spec2_cfg, '\n')

        input_file = config.get("calwebb_spec2_input_file", "input_file")
        if "_uncal_rate" in input_file:
            input_file = input_file.replace("_uncal_rate", "")
        if "_uncal" in input_file:
            input_file = input_file.replace("_uncal", "")
        final_output_name = input_file.replace(".fits", "_cal.fits")
        if core_utils.check_MOS_true(inhdu):
            # copy the MSA shutter configuration file into the pytest directory
            subprocess.run(["cp", msa_shutter_conf, "."])

        # start the timer to compute the step running time
        start_time = time.time()

        # run the pipeline
        print('Running pipeline... \n')
        if not calwebb_image2_cfg:
            Spec2Pipeline.call(step_input_file, config_file=calwebb_spec2_cfg)#, logcfg="stpipe-log.cfg")
        else:
            Image2Pipeline.call(step_input_file, config_file=calwebb_image2_cfg)

        # end the timer to compute calwebb_spec2 running time
        end_time = repr(time.time() - start_time)   # this is in seconds
        calspec2_time = " * Pipeline took "+end_time+" seconds to finish.\n"
        print(calspec2_time)
        logging.info(calspec2_time)

        # remove the copy of the MSA shutter configuration file
        if core_utils.check_MOS_true(inhdu):
            subprocess.run(["rm", msametfl])

        # add the detector string to the name of the files and move them to the working directory
        core_utils.add_detector2filename(working_directory, step_input_file)
        final_output_name_msg = "\nThe final pipeline product was saved in: "+final_output_name
        print(final_output_name_msg)
        logging.info(final_output_name_msg)

        # read the assign wcs fits file
        hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
        #scihdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False, ext=1)

        # rename and move the pipeline log file
        try:
            calspec2_pilelog = "calspec2_pipeline_"+detector+".log"
            if imaging_mode:
                calspec2_pilelog = calspec2_pilelog.replace('calspec2', 'calimage2')
            pytest_workdir = os.getcwd()
            logfile = glob(pytest_workdir+"/pipeline.log")[0]
            os.rename(logfile, os.path.join(working_directory, calspec2_pilelog))
        except:
            IndexError

        # make sure we are able to find calspec2_pilelog either in the calwebb_spec2 directory or in the working dir
        if not os.path.isfile(calspec2_pilelog):
            calspec2_pilelog = os.path.join(working_directory, calspec2_pilelog)

        # add the running time for all steps
        step_running_times = core_utils.calculate_step_run_time(calspec2_pilelog)
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
        assign_wcs_utils.print_time2file(txt_name, tot_time, string2print)
        PTT_runtimes_msg = "Pipeline and PTT run times written in file: "+os.path.basename(txt_name)+" in working directory. \n"
        print(PTT_runtimes_msg)
        logging.info(PTT_runtimes_msg)

        # move the final reporting text files to the working directory
        core_utils.move_txt_files_2workdir(detector)

        # end script for imaging case
        if imaging_mode:
            print('\nPTT finished processing imaging mode. \n')
            pytest.exit("Skipping pytests for now because they need to be written for imaging mode.")

        return hdul, step_output_file, msa_shutter_conf, esa_files_path, wcs_threshold_diff, save_wcs_plots, run_pytests, mode_used

    else:

        # create the map but remove a previous one if it exists
        if os.path.isfile(txt_name):
            os.remove(txt_name)
        assign_wcs_utils.create_completed_steps_txtfile(txt_name, step_input_file)

        # start the timer to compute the step running time of PTT
        core_utils.start_end_PTT_time(txt_name, start_time=PTT_start_time, end_time=None)
        msg = "\n Pipeline and PTT run times will be written in file: " + os.path.basename(
            txt_name) + " in working directory. \n"
        print(msg)
        logging.info(msg)

        if run_pipe_step:

            # Create the logfile for PTT, but erase the previous one if it exists
            PTTcalspec2_log = os.path.join(working_directory, 'PTT_calspec2_'+detector+'_'+step+'.log')
            if imaging_mode:
                PTTcalspec2_log = PTTcalspec2_log.replace('calspec2', 'calimage2')
            if os.path.isfile(PTTcalspec2_log):
                os.remove(PTTcalspec2_log)
            print("Information outputed to screen from PTT will be logged in file: ", PTTcalspec2_log)
            for handler in logging.root.handlers[:]:
                logging.root.removeHandler(handler)
            logging.basicConfig(filename=PTTcalspec2_log, level=logging.INFO)
            logging.info(pipeline_version)
            if change_filter_opaque:
                logging.info(change_filter_opaque_msg)

            # check that previous pipeline steps were run up to this point
            core_utils.check_completed_steps(step, step_input_file)

            if os.path.isfile(step_input_file):
                msg = " *** Step "+step+" set to True"
                print(msg)
                logging.info(msg)
                stp = AssignWcsStep()

                if core_utils.check_MOS_true(inhdu):
                    # copy the MSA shutter configuration file into the pytest directory
                    subprocess.run(["cp", msa_shutter_conf, "."])

                # get the right configuration files to run the step
                local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")

                # start the timer to compute the step running time
                print("running pipeline...")
                start_time = time.time()
                if local_pipe_cfg_path == "pipe_source_tree_code":
                    result = stp.call(step_input_file)
                else:
                    result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/assign_wcs.cfg')
                result.save(step_output_file)

                # end the timer to compute the step running time
                end_time = repr(time.time() - start_time)   # this is in seconds
                msg = "Step "+step+" took "+end_time+" seconds to finish"
                print(msg)
                logging.info(msg)

                if core_utils.check_MOS_true(inhdu):
                    # remove the copy of the MSA shutter configuration file
                    subprocess.run(["rm", msametfl])

                # rename and move the pipeline log file
                try:
                    calspec2_pilelog = "calspec2_pipeline_"+step+"_"+detector+".log"
                    if imaging_mode:
                        calspec2_pilelog = calspec2_pilelog.replace('calspec2', 'calimage2')
                    pytest_workdir = os.getcwd()
                    logfile = glob(pytest_workdir+"/pipeline.log")[0]
                    os.rename(logfile, os.path.join(working_directory, calspec2_pilelog))
                except:
                    IndexError

            else:
                msg = "Skipping step. Input file "+step_input_file+" does not exit."
                print(msg)
                logging.info(msg)
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because the input file does not exist.")


        else:
            print("Skipping running pipeline step ", step)
            # add the running time for this step
            end_time = core_utils.get_stp_run_time_from_screenfile(step, detector, working_directory)

        if os.path.isfile(step_output_file):
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            step_completed = True
            # add the running time for this step
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            return hdul, step_output_file, msa_shutter_conf, esa_files_path, wcs_threshold_diff, save_wcs_plots, run_pytests, mode_used
        else:
            step_completed = False
            # add the running time for this step
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip()




### THESE FUNCTIONS ARE TO VALIDATE THE WCS STEP

# fixture to validate the WCS
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
                                            save_figs=save_wcs_plots, threshold_diff=threshold_diff, debug=False)

    elif core_utils.check_MOS_true(hdu)  and  mode_used != "MOS_sim":
        result, log_msgs = compare_wcs_mos.compare_wcs(infile_name, esa_files_path=esa_files_path, msa_conf_name=msa_conf_name,
                                             show_figs=show_figs, save_figs=save_wcs_plots,
                                             threshold_diff=threshold_diff, debug=False)

    elif core_utils.check_IFU_true(hdu):
        result, log_msgs = compare_wcs_ifu.compare_wcs(infile_name, esa_files_path=esa_files_path, show_figs=show_figs,
                                            save_figs=save_wcs_plots, threshold_diff=threshold_diff, debug=False)

    else:#if core_utils.check_BOTS_true(hdu):
        #pytest.skip("Skipping pytest: BOTS files at the moment are not being compared against an ESA intermediary product.")
        pytest.skip("Skipping pytest: The fits file is not FS, MOS, or IFU.")

    if log_msgs is not None:
        for msg in log_msgs:
            logging.info(msg)

    if "skip" in result:
        pytest.skip("Skipping assign_wcs pytest.")

    return result


# Unit tests

# reference files from running calwebb_spec1
def test_rmask_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        print("\n * Running reference file pytest...")
        result = assign_wcs_utils.rmask_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_saturation_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest..."
        print(msg)
        logging.info(msg)
        result = assign_wcs_utils.saturation_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_superbias_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest..."
        print(msg)
        logging.info(msg)
        result = assign_wcs_utils.superbias_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_linearity_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest..."
        print(msg)
        logging.info(msg)
        result = assign_wcs_utils.linearity_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_dark_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest..."
        print(msg)
        logging.info(msg)
        result = assign_wcs_utils.dark_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_readnoise_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        print("\n * Running reference file pytest...")
        result = assign_wcs_utils.readnoise_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_gain_rfile(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][1]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest..."
        print(msg)
        logging.info(msg)
        result = assign_wcs_utils.gain_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]


# reference files specific to the WCS step
def test_camera_rfile(output_hdul):
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
        result = assign_wcs_utils.camera_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_colimator_rfile(output_hdul):
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
        result = assign_wcs_utils.colimator_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_disperser_rfile(output_hdul):
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
        result = assign_wcs_utils.disperser_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_fore_rfile(output_hdul):
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
        result = assign_wcs_utils.fore_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_fpa_rfile(output_hdul):
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
        result = assign_wcs_utils.fpa_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

def test_ifufore_rfile(output_hdul):
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
        result = assign_wcs_utils.ifufore_rfile_is_correct(output_hdul)
        if result:
            for log_msg in result[1]:
                print(log_msg)
                logging.info(log_msg)
            assert not result[0], result[0]
        else:
            msg = "Skipping ref_file pytest because data is NOT IFU.\n"
            print(msg)
            logging.info(msg)
            pytest.skip(msg)

def test_ifupost_rfile(output_hdul):
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
        result = assign_wcs_utils.ifupost_rfile_is_correct(output_hdul)
        if result:
            for log_msg in result[1]:
                print(log_msg)
                logging.info(log_msg)
            assert not result[0], result[0]
        else:
            msg = "Skipping ref_file pytest because data is NOT IFU.\n"
            print(msg)
            logging.info(msg)
            pytest.skip(msg)

def test_ifuslicer_rfile(output_hdul):
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
        result = assign_wcs_utils.ifuslicer_rfile_is_correct(output_hdul)
        if result:
            for log_msg in result[1]:
                print(log_msg)
                logging.info(log_msg)
            assert not result[0], result[0]
        else:
            msg = "Skipping ref_file pytest because data is NOT IFU.\n"
            print(msg)
            logging.info(msg)
            pytest.skip(msg)

def test_msa_rfile(output_hdul):
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
        result = assign_wcs_utils.msa_rfile_is_correct(output_hdul)
        if result:
            for log_msg in result[1]:
                print(log_msg)
                logging.info(log_msg)
            assert not result[0], result[0]
        else:
            msg = "Skipping ref_file pytest because data is NOT MOS.\n"
            print(msg)
            logging.info(msg)
            pytest.skip(msg)

def test_ote_rfile(output_hdul):
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
        result = assign_wcs_utils.ote_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]

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
        result = assign_wcs_utils.wavran_rfile_is_correct(output_hdul)
        for log_msg in result[1]:
            print(log_msg)
            logging.info(log_msg)
        assert not result[0], result[0]


# other tests specific to the WCS step
"""
def test_wavstart_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][0]
    if not run_pytests:
        msg = "Skipping pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest...\n"
        print(msg)
        logging.info(msg)
        assert assign_wcs_utils.wavstart_exists(output_hdul[1]), "The keyword WAVSTART was not added to the header."

def test_sporder_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][0]
    if not run_pytests:
        msg = "Skipping pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running reference file pytest...\n"
        print(msg)
        logging.info(msg)
        assert assign_wcs_utils.sporder_exists(output_hdul[1]), "The keyword SPORDER was not added to the header."
"""

def test_s_wcs_exists(output_hdul):
    # want to run this pytest?
    # output_hdul[6] = assign_wcs_completion_tests, assign_wcs_reffile_tests, assign_wcs_validation_tests
    run_pytests = output_hdul[6][0]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        logging.info(msg)
        pytest.skip(msg)
    else:
        msg = "\n * Running completion pytest...\n"
        print(msg)
        logging.info(msg)
        assert assign_wcs_utils.s_wcs_exists(output_hdul[0]), "The keyword S_WCS was not added to the header --> extract_2d step was not completed."

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
        #print("\n", validate_wcs, "\n")
        assert request.getfixturevalue('validate_wcs'), "Output value from compare_wcs.py is greater than threshold."


