import os
import re
import time
import argparse
import subprocess
import configparser
import logging
from glob import glob
from astropy.io import fits

from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
from jwst.group_scale.group_scale_step import GroupScaleStep
from jwst.dq_init.dq_init_step import DQInitStep
from jwst.saturation.saturation_step import SaturationStep
from jwst.superbias.superbias_step import SuperBiasStep
from jwst.refpix.refpix_step import RefPixStep
from jwst.rscd.rscd_step import RSCD_Step
from jwst.lastframe.lastframe_step import LastFrameStep
from jwst.linearity.linearity_step import LinearityStep
from jwst.dark_current.dark_current_step import DarkCurrentStep
from jwst.jump.jump_step import JumpStep
from jwst.ramp_fitting.ramp_fit_step import RampFitStep
from jwst.gain_scale.gain_scale_step import GainScaleStep


"""
This script will perform calwebb_detector1 in one single run, outputing intermediary files named:
 intermediary_step_ran.fits
 These files will be saved in the working directory.

 Usage:
    In a terminal, in the directory where the testing tool lives, and within the pipeline environment, type
    > python run_cal_detector1.py uncal_file.fits
"""


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.4"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Jul 2018 - Version 1.1: Added functions to calculate running times from screen log
# Feb 2019 - Version 1.2: made changes to be able to process 491 and 492 files in the same directory
# Mar 2019 - Version 1.3: added logging capability
# May 2019 - Version 1.4: added capability to process darks


def get_caldet1cfg_and_workingdir():
    """
    This function gets the path to where the pipeline testing tool lives
    Returns:
        calwebb_detector1_cfg = string, path to where the config file lives
    """
    config = configparser.ConfigParser()
    config.read(['../calwebb_spec2_pytests/PTT_config.cfg'])
    pipe_testing_tool_path = config.get("calwebb_spec2_input_file", "pipe_testing_tool_path")
    calwebb_detector1_cfg = os.path.join(pipe_testing_tool_path, "utils/data/calwebb_detector1.cfg")
    calwebb_tso1_cfg = os.path.join(pipe_testing_tool_path, "utils/data/calwebb_tso1.cfg")
    calwebb_dark_cfg = os.path.join(pipe_testing_tool_path, "utils/data/calwebb_dark.cfg")
    working_dir = config.get("calwebb_spec2_input_file", "working_directory")
    mode_used = config.get("calwebb_spec2_input_file", "mode_used")
    raw_data_root_file = config.get("calwebb_spec2_input_file", "raw_data_root_file")
    return calwebb_detector1_cfg, calwebb_tso1_cfg, calwebb_dark_cfg, working_dir, mode_used, raw_data_root_file


def calculate_step_run_time(pipelog_file, pipe_steps):
    """
    This function calculates the step run times from the pipelog_file.
    Args:
        pipelog_file: string, path and name of the log file
        pipe_steps: list, steps ran in calwebb_detector_1

    Returns:
        step_running_times: dictionary, with the following structure for all steps that where ran
                           step_running_times["assign_wcs"] = {start_time: float, end_time: float, run_time: float}
    """
    def get_timestamp(time_list):
        """
        This sub-function obtains a time stamp from the time strings read from the pipelog_file.
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
            while len(time_tuple) != 9:
                time_tuple.append(0)
        time_tuple = tuple(time_tuple)
        timestamp = time.mktime(time_tuple)
        return timestamp

    # read the pipelog_file
    step_running_times = {}
    with open(pipelog_file, "r") as sot:
        for line in sot.readlines():
            line = line.replace("\n", "")
            if "Ending" in line:   # make sure not to overwrite the dictionary
                break
            for pstp in pipe_steps:
                pstp = pstp.replace(".fits", "")
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
                    step_running_times[pstp] = {"start_time" : start_time, "end_time" : end_time, "run_time" : run_time}
    return step_running_times



# Get arguments to run script
parser = argparse.ArgumentParser(description='')
parser.add_argument("fits_input_uncal_file",
                    action='store',
                    default=None,
                    help='Name of fits file, i.e. blah.fits')
parser.add_argument("-sbs",
                    dest="step_by_step",
                    action='store_true',
                    default=False,
                    help='By default, calwebb detecetor 1 will be executed in one single run, if -sxs is used, execution will go step by step.')
args = parser.parse_args()

# Set the variables inputed from the command line
fits_input_uncal_file = args.fits_input_uncal_file
step_by_step = args.step_by_step

# Get the detector used
detector = fits.getval(fits_input_uncal_file, "DETECTOR", 0)

# Get the cfg file
calwebb_detector1_cfg, calwebb_tso1_cfg, calwebb_dark_cfg, working_dir, mode_used, rawdatrt = get_caldet1cfg_and_workingdir()
if mode_used != "BOTS" and mode_used.lower() != "dark":
    cfg_file = calwebb_detector1_cfg
elif mode_used == "BOTS":
    cfg_file = calwebb_tso1_cfg
elif mode_used.lower() == "dark":
    cfg_file = calwebb_dark_cfg
configfile_used = "Using this configuration file: "+cfg_file

# Initiate the PTT log file
PTTcaldetector1_log = os.path.join(working_dir, 'PTT_caldetector1_'+detector+'.log')
print("Information outputed to screen from this script will be logged in file: ", PTTcaldetector1_log)
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.basicConfig(filename=PTTcaldetector1_log, level=logging.DEBUG)

print(configfile_used)
logging.info(configfile_used)

# create the text file to record the names of the output files and the time the pipeline took to run
txt_outputs_summary = "cal_detector1_outputs_and_times_"+detector+".txt"
end_time_total = []
line0 = "# {:<20}".format("Input file: "+fits_input_uncal_file)
line1 = "# {:<16} {:<19} {:<20}".format("Step", "Output file", "Time to run [s]")
with open(txt_outputs_summary, "w+") as tf:
    tf.write(line0+"\n")
    tf.write(line1+"\n")

# Name of the file containing all the pipeline output
caldetector1_pipeline_log = "pipeline.log"

#final_output_caldet1 = "gain_scale.fits"
final_output_caldet1 = "final_output_caldet1_"+detector+".fits"
output_names = ["group_scale.fits", "dq_init.fits", "saturation.fits", "superbias.fits", "refpix.fits",
                "lastframe.fits", "linearity.fits", "dark_current.fits", "jump.fits", "ramp_fit.fits",
                final_output_caldet1]

if not step_by_step:
    print("Got arguments and will run the calwebb_detector1 pipeline in full. This may take a while...")

    # start the timer to compute the step running time
    start_time = time.time()
    result = Detector1Pipeline.call(fits_input_uncal_file, config_file=cfg_file)
    result.save(final_output_caldet1)

    # end the timer to compute pipeline running time
    end_time = time.time() - start_time  # this is in seconds
    time2finish_string = " * calwebb_detector1 took "+repr(end_time)+" seconds to finish *"
    print(time2finish_string)
    logging.info(time2finish_string)
    if end_time > 60.0:
        end_time_min = round(end_time / 60.0, 1)   # in minutes
        tot_time = repr(end_time_min)+"min"
        if end_time_min > 60.0:
            end_time_hr = round(end_time_min / 60.0, 1)   # in hrs
            tot_time = repr(end_time_hr)+"hr"
    else:
        tot_time = str(round(end_time, 1))+"sec"
    total_time = "{:<18} {:<20} {:<20}".format("", "total_time = ", repr(end_time)+"  ="+tot_time)

    # get the running time for the individual steps
    step_running_times = calculate_step_run_time(caldetector1_pipeline_log, output_names)

    # write step running times in the text file
    end_time_list = []
    for outnm in output_names:
        stp = outnm.replace(".fits", "")
        if stp in step_running_times:
            step_time = step_running_times[stp]["run_time"]
        else:
            step_time = "N/A"
        end_time_list.append(step_time)
        line2write = "{:<18} {:<20} {:<20}".format(stp, outnm, step_time)
        with open(txt_outputs_summary, "a") as tf:
            tf.write(line2write+"\n")

else:
    print("Got arguments and will run the calwebb_detector1 pipeline step by step.")

    # steps to be ran, in order
    steps_to_run = [GroupScaleStep(), DQInitStep(), SaturationStep(), SuperBiasStep(), RefPixStep(),
                    LastFrameStep(), LinearityStep(), DarkCurrentStep(), JumpStep(), RampFitStep(), GainScaleStep()]
    comp_keys = ["S_GRPSCL", "S_DQINIT", "S_SATURA", "S_SUPERB", "S_REFPIX", "S_LASTFR", "S_LINEAR",
                 "S_DARK", "S_JUMP", "S_RAMP", "S_GANSCL"]

    # run the pipeline step by step
    for i, stp_instance in enumerate(steps_to_run):
        stp = stp_instance
        if i == 0:
            step_input_file = fits_input_uncal_file
        else:
            # check if step was completed and find the appropriate input file
            j = 1
            continue_while = True
            while continue_while:
                step_input_file = output_names[i-j]
                if (i-j == 0):
                    step_input_file = fits_input_uncal_file
                    break
                if i == len(output_names)-1:
                    step_input_file = glob("*ramp*fit.fits")[0]
                    break
                if os.path.isfile(step_input_file):
                    completion_key_val = fits.getval(step_input_file, comp_keys[i-j])
                    msg = "Checking for next step... "
                    completion_keywd_msg = " * Completion keyword: "+comp_keys[i-j]+"   and value: "+completion_key_val
                    print(msg)
                    print(completion_keywd_msg)
                    logging.info(msg)
                    logging.info(completion_keywd_msg)
                    if "SKIPPED" in completion_key_val:
                        j += 1
                    elif "COMPLETE" in completion_key_val:
                        continue_while = False
        running_stp_msg = "\n-> Running step: "+str(stp)+"   with input file: "+step_input_file
        output_msg = "   output will be saved as: "+output_names[i]
        logging.info(running_stp_msg)
        logging.info(output_msg)

        # start the timer to compute the step running time
        start_time = time.time()

        if "ramp" not in output_names[i]:
            result = stp.call(step_input_file)
            result.save(output_names[i])
        else:
            # the pipeline works differently for the ramp_fit step because it has more than one output
            # this step is also hanging from the command line
            #subprocess.call(["strun", "jwst.ramp_fitting.RampFitStep", "jump.fits"])
            (out_slope, int_slope) = stp.call(step_input_file)
            out_slope.save(output_names[i])
            try:
                int_slope.save("ramp_fit_int.fits")
            except:
                AttributeError
                msg = "File has only 1 integration."
                print(msg)
                logging.info(msg)

        # end the timer to compute cal_detector1 running time
        et = time.time() - start_time   # this is in seconds
        end_time = repr(et)
        end_time_total.append(et)
        step = output_names[i].replace(".fits", "")
        msg = " * calwebb_detector1 step "+step+" took "+end_time+" seconds to finish * \n"
        print(msg)
        logging.info(msg)
        if et > 60.0:
            end_time_min = round(et / 60.0, 1)   # in minutes
            end_time = repr(end_time_min)+"min"
            if end_time_min > 60.0:
                end_time_hr = round(end_time_min / 60.0, 1)   # in hrs
                end_time = repr(end_time_hr)+"hr"
        else:
            end_time = repr(round(et, 1))+"sec"

        # record results in text file
        line2write = "{:<18} {:<20} {:<20}".format(step, output_names[i], end_time)
        with open(txt_outputs_summary, "a") as tf:
            tf.write(line2write+"\n")


    # record total time in text file
    tot_time_sec = sum(end_time_total)
    if tot_time_sec > 60.0:
        tot_time_min = round((tot_time_sec/60.0), 1)
        tot_time = repr(tot_time_min)+"min"
        if tot_time_min > 60.0:
            tot_time_hr = round((tot_time_min/60.0), 1)
            tot_time = repr(tot_time_hr)+"hr"
    else:
        tot_time = round((tot_time_sec/60.0), 1)
    total_time = "{:<18} {:>20} {:>20}".format("", "total_time  ", repr(tot_time_sec)+"  ="+tot_time)


# record total time in text file
with open(txt_outputs_summary, "a") as tf:
    tf.write(total_time+"\n")
msg = "\n ** Calwebb_detector 1 took "+repr(tot_time)+" to complete **"
print(msg)
logging.info(msg)


# Move products to working dir

# rename and move the pipeline log file
pytest_workdir = os.getcwd()
logfile = glob(pytest_workdir+"/*pipeline*.log")[0]
new_name = "caldetector1_pipeline_"+detector+".log"
os.rename(logfile, os.path.join(working_dir, new_name))

# move the PTT log file and the fits intermediary product files
fits_list = glob("*.fits")
if len(fits_list) >= 1:
    msg = "Output fits files are located at: "+working_dir
    print(msg)
    logging.info(msg)
    for ff in fits_list:
        if "step_" in ff:
            ff_newname = os.path.join(working_dir, ff.replace("step_", ""))
        else:
            ff_newname = os.path.join(working_dir, ff)
        if detector.lower() not in ff.lower():
            ff_newname = ff_newname.replace(".fits", "_"+detector+".fits")
        subprocess.run(["mv", ff, ff_newname])
    # move text files too
    subprocess.run(["mv", txt_outputs_summary, working_dir])

else:
    msg = "No fits files detected after calwbb_detector1 finished. Exiting script."
    print(msg)
    logging.info(msg)

msg = "Script  run_cal_detector1.py  finished."
print(msg)
logging.info(msg)
