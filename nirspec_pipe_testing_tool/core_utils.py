import collections
import os
import re
import glob
import time
import subprocess
import configparser
import logging
from logging.handlers import RotatingFileHandler
from astropy.io import fits
from nirspec_pipe_testing_tool.calwebb_spec2_pytests import TESTSDIR

import jwst
import nirspec_pipe_testing_tool as nptt

'''
This script contains functions frequently used in the test suite.
'''

# HEADER
__author__ = "M. Pena-Guerrero"
__version__ = "1.3"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Feb 2019 - Version 1.2: made changes to be able to process 491 and 492 files in the same directory
# Apr 2023 - Version 1.3: Cleaned-up code


# dictionary of the steps and corresponding strings to be added to the file name after the step has ran
step_string_dict = collections.OrderedDict()
# spec2
step_string_dict["assign_wcs"] = {"outfile": True, "suffix": "_assign_wcs"}
step_string_dict["bkg_subtract"] = {"outfile": True, "suffix": "_subtract_images"}
step_string_dict["imprint_subtract"] = {"outfile": True, "suffix": "_imprint"}
step_string_dict["msa_flagging"] = {"outfile": True, "suffix": "_msa_flagging"}
step_string_dict["extract_2d"] = {"outfile": True, "suffix": "_extract_2d"}
step_string_dict["srctype"] = {"outfile": True, "suffix": "_srctype"}
step_string_dict["wavecorr"] = {"outfile": True, "suffix": "_wavecorr"}
step_string_dict["flat_field"] = {"outfile": True, "suffix": "_flat_field"}
# step_string_dict["straylight"]       = {"outfile" : True, "suffix" : "_stray"}   # MIRI only
# step_string_dict["fringe"]           = {"outfile" : True, "suffix" : "_fringe"}   # MIRI only
step_string_dict["pathloss"] = {"outfile": True, "suffix": "_pathloss"}
step_string_dict["barshadow"] = {"outfile": True, "suffix": "_barshadow"}
step_string_dict["photom"] = {"outfile": True, "suffix": "_photom"}
step_string_dict["resample_spec"] = {"outfile": True, "suffix": "_s2d"}
step_string_dict["cube_build"] = {"outfile": True, "suffix": "_s3d"}
step_string_dict["extract_1d"] = {"outfile": True, "suffix": "_x1d"}

# spec3
spec3_step_dict = collections.OrderedDict()
spec3_step_dict["assign_mtwcs"] = {"outfile": True, "suffix": "assign_mtwcs"}
spec3_step_dict["master_background"] = {"outfile": True, "suffix": "_master_background"}
spec3_step_dict["exp_to_source"] = {"outfile": True, "suffix": "_cal"}
spec3_step_dict["outlier_detection"] = {"outfile": True, "suffix": "_crf"}
spec3_step_dict["resample_spec"] = {"outfile": True, "suffix": "_s2d"}
spec3_step_dict["cube_build"] = {"outfile": True, "suffix": "_s3d"}
spec3_step_dict["extract_1d"] = {"outfile": True, "suffix": "_x1d"}


def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default,
       split on a comma and strip whitespaces."""
    return [chunk.strip(chars) for chunk in option.split(sep)]


def mk_stpipe_log_cfg(outpur_dir, log_name):
    """
    Create a configuration file with the name log_name, where
    the pipeline will write all output.
    Args:
        outpur_dir: sstr: path of the output directory
        log_name: str, name of the log to record screen output
    Returns:
        nothing
    """
    config = configparser.ConfigParser()
    config.add_section("*")
    config.set("*", "handler", "file:" + log_name)
    config.set("*", "level", "INFO")
    pipe_log_config = os.path.join(outpur_dir, "pipeline-log.cfg")
    config.write(open(pipe_log_config, "w"))


def mk_nptt_log(logname, reset=True):
    """
    Wrapper for setting up a standard logger.

    Usage:
    lg = set_logger("myname", "/path/to/log/folder")
    lg.info("This is some info I want to give")
    lg.warning("This is a warning")
    lg.error("This is an error")

    Args:
        logname : str, name of the logger.
        reset: boolean, start a new file or not

    Returns:
        logger instance
    """
    # Start the logger and set the default log level
    logger = logging.getLogger(logname)
    logging_level = logging.INFO
    logger.setLevel(logging_level)
    # Restart the logger each time and set the report format
    if reset:
        file_handler = RotatingFileHandler(logname, mode='w')
    else:
        file_handler = logging.FileHandler(logname)
    formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(name)s : %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    print("\n * NPTT screen output will be logged in file: ", logname)
    return logger


def get_sci_extensions(fits_file_name):
    """
    This function obtains all the science extensions in the given file
    Args:
        fits_file_name: name of the fits file of interest

    Returns:
        sci_dicts: dictionary of the numbers of the science extensions
    """
    hdulist = fits.open(fits_file_name)
    sci_dicts = collections.OrderedDict()
    for ext, hdu in enumerate(hdulist):
        if hdu.name == "SCI":
            sci_dicts[hdu.header['SLTNAME']] = ext
    hdulist.close()
    return sci_dicts


def get_step_inandout_filename(step, initial_input_file, output_directory, debug=False):
    """
    This function determines the corresponding input file name for the step (i.e. the pipeline expects a specific
    format and name for each step). This is according to which steps where set to True  in the configuration file.
    Args:
        step: string, pipeline step to be ran
        initial_input_file: the base name of the input file for calwebb_spec2
        output_directory: string, path where the output files will be saved at

    Returns:
        step_input_filename: string, the base name of the input file for the specified step
        step_output_filename: string, the base name of the output file for the specified step
        in_file_suffix : string, the suffix added to the initial file name for the step input file
        out_file_suffix : string, the suffix added to the initial file name for the step output file
    """

    # make sure the input file name has the detector included in the name of the output files
    hdr = fits.getheader(initial_input_file)
    detector = hdr["DETECTOR"]
    exp_type = hdr["EXP_TYPE"]
    initial_input_file_basename = os.path.basename(initial_input_file)
    if "_uncal_rate" in initial_input_file_basename:
        initial_input_file_basename = initial_input_file_basename.replace("_uncal_rate", "")
    if detector.lower() not in initial_input_file_basename.lower():
        initial_input_file_basename = initial_input_file_basename.replace(".fits", "_" + detector + ".fits")

    # get the right input and output name according to the steps dictionary
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = "", "", "", ""
    for i, stp in enumerate(step_string_dict):
        if stp == step:
            if debug:
                print("Found step in dictionary: ", step)
            # find the input suffix according to the step
            if stp == "assign_wcs":
                step_input_filename = initial_input_file
            else:
                # recurrently go backwards in the step dictionary until an input file exists
                exit_while_loop = False
                counter = 0
                if debug:
                    print("Will enter while loop to find the appropriate input file name for this step.")
                while not exit_while_loop:
                    # look for the input file starting from the end of the dictionary
                    j = len(step_string_dict.items()) - 1
                    for s in reversed(step_string_dict.items()):
                        # This is to make sure we do not enter an infinite while loop
                        counter += 1
                        if counter > len(step_string_dict.items()):  # 13 is the number of steps of calwebb_spec2
                            exit_while_loop = True
                            print("Limiting number of iterations reached. No input file found. ")
                            print("NPTT will use initial input file for step ", step)
                            in_file_suffix = ""
                            step_input_filename = initial_input_file
                            break

                        if debug:
                            print("In the for loop of reversed steps at step: ", s[0], "   j=", j, " i=", i)
                        # only go through the loop if at the step before the goal step
                        if j < i:
                            outfile = step_string_dict[s[0]]["outfile"]
                            if outfile:
                                in_file_suffix = step_string_dict[s[0]]["suffix"]
                                step_input_basename = initial_input_file_basename.replace(".fits",
                                                                                          in_file_suffix + ".fits")
                                if "gain_scale" in step_input_basename:
                                    step_input_basename = step_input_basename.replace("_gain_scale", "")
                                step_input_filename = os.path.join(output_directory, step_input_basename)

                                # make sure the cube_build step only runs for IFU data
                                if "cube_build" in s and "IFU" not in exp_type:
                                    continue

                                # make sure the input file exists
                                if debug:
                                    print("Step creates output file, checking if it exists: ", step_input_filename)
                                if os.path.isfile(step_input_filename):
                                    # break and exit the while loop
                                    exit_while_loop = True
                                    break
                                else:
                                    # maybe the initial file had _rate and the output does not
                                    if '_rate' in step_input_filename:
                                        print("File does not exist, checking without '_rate' suffix... ")
                                        step_input_filename = step_input_filename.replace('_rate', '')
                                        # break and exit the while loop if the file exists
                                        if os.path.isfile(step_input_filename):
                                            exit_while_loop = True
                                            break
                                    if debug:
                                        print("Step did not run to create product file.")
                                    j -= 1
                                    continue
                            else:
                                if debug:
                                    print("Step does not produce an outfile, continuing in the while loop...")
                                continue
                        elif j == 0:
                            # make sure the while loop ends at this point
                            print("NPTT will use initial input file for step ", step)
                            in_file_suffix = ""
                            step_input_filename = initial_input_file
                            exit_while_loop = True
                            break
                        else:
                            j -= 1
                            continue

            # determine the output suffix according to the step
            if step_string_dict[stp]["outfile"]:
                out_file_suffix = step_string_dict[stp]["suffix"]
            else:
                out_file_suffix = in_file_suffix
            step_output_basename = initial_input_file_basename.replace(".fits", out_file_suffix + ".fits")
            if '_rate' in step_output_basename:
                step_output_basename = step_output_basename.replace('_rate', '')
            if "_gain_scale".lower() in step_output_basename.lower():
                step_output_basename = step_output_basename.replace("_gain_scale", "")
            step_output_filename = os.path.join(output_directory, step_output_basename)
            # now exit the for loop because step was reached
            break
    if debug:
        print("Step input file:  ", step_input_filename)
        print("Step output file: ", step_output_filename)
        input()

    return in_file_suffix, out_file_suffix, step_input_filename, step_output_filename


def add_detector2filename(output_directory, step_input_file):
    """
    This function is only used in the case when the pipeline is run in full. At the end of the run, before the files
    have been moved to the working directory, the names will be changed to include the detector name.
    Args:
        output_directory: string, path of working directory
        step_input_file: string, full name of initial input file

    Returns:
        nothing
    """
    # list of the steps suffix to be added to the file name after the step has ran
    step_strings = []
    for stp, stp_dict in step_string_dict.items():
        suffix = stp_dict['suffix']
        step_strings.append(suffix)

    # get the detector name and add it
    det = fits.getval(step_input_file, "DETECTOR", 0)
    # step_input_file_basename = os.path.basename(step_input_file).replace(".fits", "")
    nptt_directory = TESTSDIR  # directory where NPTT lives
    step_files = glob.glob(os.path.join(nptt_directory, "*.fits"))
    for sf in step_files:
        for stp_str in step_strings:
            if stp_str in sf.lower():  # look for the step string appearing in the name of the files
                # make sure to get the right name for the MSA flagging step output file
                if "msaflagopenstep" in sf:
                    sf = sf.replace("msaflagopenstep", "msa_flagging")
                if det.lower() not in sf.lower():  # only add the detector to the name if it is not part of the basename
                    new_name = os.path.basename(sf.replace(stp_str + ".fits", "_" + det + stp_str + ".fits"))
                else:
                    print("Detector ", det, " is already part of the file name ", os.path.basename(sf),
                          ", NPTT will not add it again.")
                    new_name = os.path.basename(sf)
                # move to the working directory
                new_name = os.path.join(output_directory, new_name)
                subprocess.run(["mv", sf, new_name])


def read_True_steps_suffix_map(txtfile_name_with_path):
    """
    This function reads the text file that contains all the steps in the configuration file, the
    corresponding suffix, and whether they were completed or not.

    Args:
        txtfile_name_with_path: string, full name and path of the text file

    Returns:
        steps_list: list, steps set to True in configuration file
        suffix_list: list, suffix for the output file corresponding to each step had it completed
        completion_list: list, strings of True or False depending on whether the step completed or not

    """
    steps_list, suffix_list, completion_list = [], [], []
    with open(txtfile_name_with_path, "r") as tf:
        for line in tf.readlines():
            if "#" not in line:
                info = line.split()
                steps_list.append(info[0])
                suffix_list.append(info[1])
                completion_list.append(info[2])
    return steps_list, suffix_list, completion_list


def read_completion_to_full_run_map(full_run_map, step):
    """
    This function adds the completion boolean to the full run map file.
    Args:
        full_run_map: string, path and name of the text file to be updated
        step: string, name of the step to be updated

    Returns:
        step_product = string, name of the intermediary fits product of the step
    """
    step_product = ""
    with open(full_run_map, "r") as tf:
        for line in tf.readlines():
            if "#" not in line:
                if step in line:
                    step_product = line.split()[1]
    return step_product


def calculate_step_run_time(screen_output_txt):
    """
    This function calculates the step run times from the screen_output_txt file.
    Args:
        screen_output_txt: string, path and name of the text file

    Returns:
        step_running_times: dictionary, with the following structure for all steps that where ran
                           step_running_times["assign_wcs"] = {start_time: float, end_time: float, run_time: float}
    """

    def get_timestamp(time_list):
        """
        This sub-function obtains a time stamp from the time strings read from the screen output text file.
        Args:
            time_list: list, contains 2 strings (date and time)

        Returns:
            timestamp: float, time stamp
        """
        time_tuple, secs_frac, timestamp = [], None, None
        stp_date = time_list[0].split("-")
        stp_time = time_list[1].split(":")
        # append the date
        for item in stp_date:
            check_for_letters = re.search('[a-zA-Z]', item)
            if check_for_letters is None:
                if "," not in item:  # make sure there are no fractional numbers
                    time_tuple.append(int(item))
            else:
                return timestamp
        # append the time
        for item in stp_time:
            if "," not in item:  # make sure there are no fractional numbers
                time_tuple.append(int(item))
            else:
                # separate the fractions of seconds from seconds and append them
                secs_frac = stp_time[-1].split(",")
                for sf in secs_frac:
                    time_tuple.append(int(sf))
        # convert the time tuple into a time stamp, which HAS to have nine integer numbers
        if len(time_tuple) < 9:
            while len(time_tuple) is not 9:
                time_tuple.append(0)
        time_tuple = tuple(time_tuple)
        timestamp = time.mktime(time_tuple)
        return timestamp

    # read the screen_output_txt file
    pipe_steps = ["assign_wcs", "bkg_subtract", "imprint_subtract", "msa_flagging", "extract_2d",
                  "srctype", "wavecorr", "flat_field", "pathloss", "barshadow", "photom",
                  "resample_spec", "cube_build", "extract_1d"]
    step_running_times = {}
    with open(screen_output_txt, "r") as sot:
        for line in sot.readlines():
            line = line.replace("\n", "")
            if "Ending" in line:  # make sure not to overwrite the dictionary
                continue
            for pstp in pipe_steps:
                if pstp in line and "running with args" in line:
                    start_time_list = line.split()[0:2]
                    start_time = get_timestamp(start_time_list)
                    if start_time is None:
                        continue
                if pstp in line and "done" in line:
                    end_time_list = line.split()[0:2]
                    end_time = get_timestamp(end_time_list)
                    if end_time is None:
                        continue
                    run_time = end_time - start_time
                    step_running_times[pstp] = {"start_time": start_time, "end_time": end_time, "run_time": run_time}
    return step_running_times


def get_stp_run_time_from_screenfile(step, det, output_directory):
    """
    This function calculates the running time for the given step from the screen output file. It is used when NPTT
    is told to skip running the pipeline step.
    Args:
        step: string, name of step
        det: string, detector used - DETECTOR keyword value
        output_directory: string, path to the working directory

    Returns:
        end_time: string, running time for this step

    """

    # add the running time for this step
    calspec2_pipelog = "calspec2_pipeline_" + det + ".log"
    calspec2_pipelog = os.path.join(output_directory, calspec2_pipelog)

    # if the pipelog is not found look for the step pipelog
    if not os.path.isfile(calspec2_pipelog):
        calspec2_pipelog = calspec2_pipelog.replace("pipeline", step)

    # if NPTT can find either log file
    end_time = None
    if os.path.isfile(calspec2_pipelog):
        step_running_times = calculate_step_run_time(calspec2_pipelog)

        for stp in step_string_dict:
            if stp in step_running_times:
                if stp == step:
                    step_time = step_running_times[stp]["run_time"]
                    end_time = step_time
                    break
            else:
                continue

    if end_time is None:
        print("\n * NPTT unable to calculate time from " + calspec2_pipelog + " for step ", step)
        print("   - This means that the step has not been run yet, or was set not to run in NPTT_config. \n")
        end_time = "0.0"

    return end_time


def add_completed_steps(True_steps_suffix_map, step, outstep_file_suffix, step_completed, end_time):
    """
    This function adds the completed steps along with the corresponding suffix of the output file name into a text file.
    Args:
        True_steps_suffix_map: string, full path of where the text file will be written into
        step: string, pipeline step just ran
        outstep_file_suffix: string, suffix added right before .fits to the input file
        step_completed: boolean, True if the step was completed and False if it was skipped
        end_time: string, time it took for the step to run (in seconds)

    Returns:
        nothing
    """
    # print ("Map saved at: ", True_steps_suffix_map)
    if (float(end_time)) > 60.0:
        end_time_min = float(end_time) / 60.  # this is in minutes
        if end_time_min > 60.0:
            end_time_hr = end_time_min / 60.  # this is in hours
            end_time = repr(end_time) + "  =" + repr(round(end_time_hr, 1)) + "hr"
        else:
            end_time = repr(end_time) + "  =" + repr(round(end_time_min, 1)) + "min"
    line2write = "{:<20} {:<20} {:<20} {:<20}".format(step, outstep_file_suffix, str(step_completed), end_time)
    with open(True_steps_suffix_map, "a") as tf:
        tf.write(line2write + "\n")


def start_end_nptt_time(txt_name, start_time=None, end_time=None):
    """
    This function calculates and prints the starting/ending NPTT running time in the True_steps_suffix_map.txt or
    full_run_map.txt file.
    Args:
        txt_name: string, path and name of the text file
        start_time: float, starting time
        end_time: float, ending time

    Returns:
        Nothing.
    """
    if start_time is not None:
        line1 = "# {:<20} {:<20}".format('Starting NPTT running time: ', repr(start_time))
        line2 = "# {:<17} {:<20} {:<20} {:<20}".format("Step", "Added suffix", "Step complition", "Time to run [s]")
        # start the timer to compute the step running time of NPTT
        print(line2)
        with open(txt_name, "w") as tf:
            tf.write(line1 + "\n")
            tf.write(line2 + "\n")

    if end_time is not None:
        # get the start time from the file
        with open(txt_name, "r") as tf:
            for line in tf.readlines():
                if "Starting NPTT running time" in line:
                    nptt_start_time = float(line.split(":")[-1])
                    break
        # compute end the timer to compute NPTT running time
        nptt_total_time = end_time - nptt_start_time  # this is in seconds
        if nptt_total_time > 60.0:
            nptt_total_time_min = round(nptt_total_time / 60.0, 1)  # in minutes
            nptt_total_run_time = repr(nptt_total_time_min) + "min"
            if nptt_total_time_min > 60.0:
                nptt_total_time_hr = round(nptt_total_time_min / 60.0, 1)  # in hrs
                nptt_total_run_time = repr(nptt_total_time_hr) + "hr"
        else:
            nptt_total_run_time = repr(round(nptt_total_time, 1)) + "sec"

        print("# The total time for NPTT to run (including pipeline) was " + repr(nptt_total_time) + " seconds.")
        line2write = "{:<20} {:<20} {:<20} {:<20}".format('', '', 'NPTT+pipe_total_time ',
                                                          repr(nptt_total_time) + '  =' + nptt_total_run_time)
        print(line2write)
        with open(txt_name, "a") as tf:
            tf.write(line2write + "\n")
        print("writen to: ", txt_name)


def get_correct_input_step_filename(initial_input_file, steps_list, suffix_list, completion_list):
    """
    This function gets the name of the input step file depending on the info stored in the steps mapping
    text file.

    Args:
        initial_input_file: string, name of the input file to initialize calwebb_spec2
        steps_list: list, steps set to True in configuration file
        suffix_list: list, suffix corresponding to the output file for the step
        completion_list: list, strings of booleans depending if the step ran to completion or not

    Returns:
        step_input_filename: string, base name of the input file for the step
    """
    step_input_filename = initial_input_file
    for i, stp in enumerate(steps_list):
        if str_to_bool(completion_list[i]):
            step_input_filename = step_input_filename.replace(".fits", suffix_list[i] + ".fits")
    return step_input_filename


def str_to_bool(s):
    """
    This function converts a string of True or False into a boolean value.
    Args:
        s: string

    Returns:
        boolean value of s
    """
    if s == 'True':
        return True
    elif s == 'False':
        return False


def set_inandout_filenames(step, config):
    """
    Set the step input and output file names, and add the step suffix to the step map.
    Args:
        step: string, name of the step as in the pipeline
        config: object, this is the configuration file object

    Returns:
        step_input_filename = string
        step_output_filename = string
        in_file_suffix = string, suffix of the input file
        out_file_suffix = string, suffix of the output file
        True_steps_suffix_map = string, path to the suffix map
    """
    data_directory = config.get("calwebb_spec2_input_file", "data_directory")
    output_directory = config.get("calwebb_spec2_input_file", "output_directory")
    initial_input_file_basename = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(data_directory, initial_input_file_basename)
    # Get the detector used
    detector = fits.getval(initial_input_file, "DETECTOR", 0)
    True_steps_suffix_map = "spec2_full_run_map_" + detector + ".txt"
    if os.path.isfile(initial_input_file):
        print("\n Taking initial input file from data_directory:")
    else:
        initial_input_file = os.path.join(output_directory, initial_input_file_basename)
        print("\n Taking initial file from output_directory: ")
    print(" Initial input file = ", initial_input_file, "\n")
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")

    if not run_calwebb_spec2:
        True_steps_suffix_map = os.path.join(output_directory, "spec2_step_run_map_" + detector + ".txt")
        print("Pipeline was set to run step by step. Steps and suffix map file named: ", True_steps_suffix_map,
              ", located in output directory.")
    else:
        print("Pipeline was set to run in full. Suffix map named: spec2_full_run_map_DETECTOR.txt, "
              "located in output directory.")
    suffix_and_filenames = get_step_inandout_filename(step, initial_input_file, output_directory)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    print("step_input_filename = ", step_input_filename)
    print("step_output_filename = ", step_output_filename)
    return step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


def read_info4output_vars(config, step_info):
    """
    Unfold the variables from the return of the function set_inandout_filenames, and get variables from the
    configuration file.
    Args:
        config: object, this is the configuration file object
        step_info: list, this is the output from the function set_inandout_filenames

    Returns:
        set_inandout_filenames_info = list, variables from set_inandout_filenames function and configuration file
    """
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    output_directory = config.get(initiate_calwebb_spc2, "output_directory")
    step, step_input_filename, output_file, in_file_suffix, outstep_file_suffix, True_steps_suffix_map = step_info
    txt_name = os.path.join(output_directory, True_steps_suffix_map)
    step_input_file = os.path.join(output_directory, step_input_filename)
    step_output_file = os.path.join(output_directory, output_file)
    set_inandout_filenames_info = [step, txt_name, step_input_file, step_output_file, outstep_file_suffix]

    return set_inandout_filenames_info


def check_FS_true(hdr):
    """
    This function checks if the fits file is a Fixed Slit.
    Args:
        hdr: dict, header keywords

    Returns:
        result: boolean, if true, the file is assumed to be Fixed Slit
    """
    result = False
    if "EXP_TYPE" in hdr:
        if hdr["EXP_TYPE"] == "NRS_FIXEDSLIT":
            result = True
    return result


def check_BOTS_true(hdr):
    """
    This function checks if the fits file is a Bright Object Time Series, BOTS.
    Args:
        hdr: dict, header keywords

    Returns:
        result: boolean, if true, the file is assumed to be BOTS
    """
    result = False
    if "EXP_TYPE" in hdr:
        if hdr["EXP_TYPE"] == "NRS_BRIGHTOBJ":
            result = True
    return result


def check_MOS_true(hdr):
    """
    This function checks if the fits file is Multi-Object Spectroscopy (MOS).
    Args:
        hdr: dict, header keywords

    Returns:
        result: boolean, if true, the file is assumed to be MOS data
    """
    result = False
    if "EXP_TYPE" in hdr:
        if "MSA" in hdr["EXP_TYPE"]:
            result = True
    return result


def check_IFU_true(hdr):
    """
    This function checks if the fits file is IFU data.
    Args:
        hdr: dict, header keywords

    Returns:
        result: boolean, if true, the file is assumed to be IFU data
    """
    result = False
    if "EXP_TYPE" in hdr:
        if "IFU" in hdr["EXP_TYPE"]:
            result = True
    return result


def check_completed_steps(step, step_input_file, caldet1=False):
    """
    This function checks the header of the file to make sure that steps up to this point were run.
    Args:
        step: string, calwebb_spec2 step
        step_input_file: string, full path and name of step input file
        caldet1: boolean, if True this function will also check for calwebb_detector1 completed steps

    Returns:
        nothing but prints warnings
    """
    # get the header of the file
    hdr = fits.getheader(step_input_file)

    # get all the completed steps
    steps_caldet1 = ["GRPSCL", "DQINIT", "SATURA", "SUPERB", "REFPIX", "LINEAR", "DARK", "JUMP", "RAMP", "GANSCL"]
    steps_calwebbspec2 = collections.OrderedDict()
    steps_calwebbspec2["assign_wcs"] = "WCS"
    steps_calwebbspec2["msa_flagging"] = "MSAFLG"
    steps_calwebbspec2["extract_2d"] = "EXTR2D"
    steps_calwebbspec2["srctype"] = "SRCTYP"
    steps_calwebbspec2["wavecorr"] = "WAVCOR"
    steps_calwebbspec2["flat_field"] = "FLAT"
    steps_calwebbspec2["pathloss"] = "PTHLOS"
    steps_calwebbspec2["barshadow"] = "BARSHA"
    steps_calwebbspec2["photom"] = "PHOTOM"
    steps_calwebbspec2["resample_spec"] = "RESAMP"
    steps_calwebbspec2["cube_build"] = "IFUCUB"
    steps_calwebbspec2["extract_1d"] = "EXTR1D"
    for key, val in hdr.items():
        # check that calwebb_detector1 was run
        if caldet1:
            for pstp in steps_caldet1:
                if "S_" + pstp in hdr:
                    if val == "COMPLETE" or val == "SKIPPED":
                        continue
                else:
                    print("*** WARNING: calwebb_detector1 step ", pstp,
                          " was not run. Please make sure this is intentional. \n")
        # check that calwebb_spec2 steps up to relevant step were ran
        step_list = []
        for pstp in steps_calwebbspec2:
            if step == pstp:
                break
            step_list.append(pstp)

    for s in step_list:
        keyword = "S_" + steps_calwebbspec2[s]
        if keyword not in hdr:
            print("*** WARNING: Keyword ", keyword,
                  " for step completion does not appear in the header of the input file. This means the step "
                  "may not have been ran, ")
            print(
                "             please make sure this is intentional or output products beyond this step "
                "may be incorrect. \n")


def find_which_slit(hdr):
    """
    This function determines which Fixed Slit was used
    Args:
        hdr: dict, header keywords

    Returns:
        s: string, slit used or is None if not in the list
    """
    # the order of this list corresponds to the
    slits = ["S200A1", "S200A2", "S200B1", "S400A1", "S1600A1"]
    if "FXD_SLIT" in hdr:
        for i, s in enumerate(slits):
            print("slit: ", s, i)
            if s in hdr["FXD_SLIT"]:
                return i + 1, s


def find_DETECTOR(hdr):
    """
    This function determines which detector was used
    Args:
        hdr: dict, header keywords

    Returns:
        det: string, either NRS1 or NRS2
    """
    if "DETECTOR" in hdr:
        det = hdr["DETECTOR"]
        return det


def get_time_to_run_pipeline(True_steps_suffix_map):
    """
    This function calculates the total running time by adding the individual time per step.
    Args:
        True_steps_suffix_map: string, full path of where the text file will be written into

    Returns:
        total_time: string, total calculate the total time by reading the time of each step from the map

    """
    times_per_step = []
    with open(True_steps_suffix_map, "r") as tf:
        for line in tf.readlines():
            if '#' in line:
                continue
            # this is in case there are = signs in the line
            if "=" in line:
                line = line.split("=")[0]
            t = line.split()[3]
            if "'" in t:
                t = t.replace("'", "")
            t = float(t)
            times_per_step.append(t)
            if "extract_1d" in line:
                break
    total_time = sum(times_per_step)
    return total_time


def get_latest_file(filetype, detector=None, disregard_known_files=False):
    """
    This function gets the latest modified file of filetype. This function will look into the calwebb_spec2_pytests
    directory for the given file.
    Args:
        filetype: string, name/type of file type, e.g. *.txt, *.html, full_run_map_DETECTOR.txt
        detector: string, name in header keyword DETECTOR
        disregard_known_files: boolean, if True function will not look for True_steps_suffix_map_DETECTOR.txt
        or full_run_map_DETECTOR.txt

    Returns:
        latest_filetypefile: string, name of the latest file modified of type filetype
    """
    # get a list of all the filetype files in the calwebb_spec2_pytests dir
    list_of_filetypefiles = glob.glob(filetype)
    # find the latest of the filetype files but exclude known file names
    if disregard_known_files:
        if detector is None:
            print("get_latest_file: DETECTOR not defined.")
            exit()
        if "spec2_suffix_map_" + detector + ".txt" in list_of_filetypefiles:
            idx = list_of_filetypefiles.index("spec2_suffix_map_" + detector + ".txt")
            list_of_filetypefiles.pop(idx)
        if "spec2_full_run_map_" + detector + ".txt" in list_of_filetypefiles:
            idx = list_of_filetypefiles.index("spec2_full_run_map_" + detector + ".txt")
            list_of_filetypefiles.pop(idx)
    latest_filetypefile = "File not found."
    if len(list_of_filetypefiles) > 0:
        latest_filetypefile = max(list_of_filetypefiles, key=os.path.getctime)
    return latest_filetypefile


def convert_html2pdf(detector):
    """
    This function converts the latest html file into a pdf.
    In order to work it needs this plugin (https://pypi.org/project/pdfkit/#description):
    pip install pdfkit
    and the OSX 64-bit download from: https://wkhtmltopdf.org/downloads.html
    """
    # get a list of all the html files in the calwebb_spec2_pytests dir
    latest_htmlfile = get_latest_file("*" + detector + "*.html")
    # create the pdf output name
    pdf_file = latest_htmlfile.replace(".html", ".pdf")
    # convert the html report into a pdf file
    # import pdfkit
    # options = {'dpi': 96}
    # pdfkit.from_file(latest_htmlfile, pdf_file, options=options)
    print("\n Converted ", latest_htmlfile, " to ", pdf_file, ". Both files are available in current directory. \n")


def move_txt_files_2outdir(output_dir, detector):
    """
    This function moves the NPTT output reporting files into the output directory.
    Args:
        output_dir: string, path to the output directory
        detector: string, value of the DETECTOR hearder keyword

    Returns:
        Nothing
    """
    # get a list of all the txt files in the calwebb_spec2_pytests dir
    latest_screenoutputtxtfile = get_latest_file("*screen*" + detector + "*.txt", detector,
                                                 disregard_known_files=True)  # this picks up the output_screen file
    latest_suffixmaptxtfile = get_latest_file("spec2_suffix_map_" + detector + ".txt")
    latest_fullrunmaptxtfile = get_latest_file("spec2_full_run_map_" + detector + ".txt")
    spec3_suffixmaptxtfile = get_latest_file("spec3_suffix_map_" + detector + ".txt")
    spec3_fullrunmaptxtfile = get_latest_file("spec3_full_run_map_" + detector + ".txt")
    # move these files into the working directory
    files2move = [latest_screenoutputtxtfile,
                  latest_suffixmaptxtfile, latest_fullrunmaptxtfile,
                  spec3_suffixmaptxtfile, spec3_fullrunmaptxtfile]
    for f in files2move:
        if f != "File not found." and os.path.isfile(f):
            subprocess.run(["mv", f, output_dir])
            print(f, " moved to ", output_dir)
