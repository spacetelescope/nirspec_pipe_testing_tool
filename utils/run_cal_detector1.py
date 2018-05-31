import os
import time
import argparse
import subprocess
import configparser
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
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed


def get_caldet1cfg_and_workingdir():
    """
    This function gets the path to where the pipeline testing tool lives
    Returns:
        calwebb_detector1_cfg = string, path to where the config file lives
    """
    config = configparser.ConfigParser()
    config.read(['../calwebb_spec2_pytests/PTT_config.cfg'])
    pipe_testing_tool_path = config.get("calwebb_spec2_input_file", "pipe_testing_tool_path")
    calwebb_detector1_cfg = os.path.join(pipe_testing_tool_path, "utils/calwebb_detector1.cfg")
    working_dir = config.get("calwebb_spec2_input_file", "working_directory")
    mode_used = config.get("calwebb_spec2_input_file", "mode_used")
    raw_data_root_file = config.get("calwebb_spec2_input_file", "raw_data_root_file")
    return calwebb_detector1_cfg, working_dir, mode_used, raw_data_root_file



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

# Get the calwebb_detector1.cfg file
#try:
calwebb_detector1_cfg, working_dir, mode_used, rawdatrt = get_caldet1cfg_and_workingdir()
print ("Using this configuration file: ", calwebb_detector1_cfg)
#except:
    #KeyError
    #print ("This script expects to find the calwebb_detector1.cfg file in current directory.")
    #calwebb_detector1_cfg = "calwebb_detector1.cfg"
    #working_dir = "."
    #exp_type = fits.getval(fits_input_uncal_file, 'EXP_TYPE', 0)
    #if "IFU" in exp_type:
    #    mode_used = "IFU"
    #elif "MSA" in exp_type:
    #    mode_used = "MOS"
    #elif "FIXED" in exp_type:
    #    mode_used = "FS"
    #raw_data_root_file = ""

# create the text file to record the names of the output files and the time the pipeline took to run
txt_outputs_summary = "cal_detector1_outputs_and_times.txt"
end_time_total = []
line0 = "# {:<20}".format("Input file: "+fits_input_uncal_file)
line1 = "# {:<17} {:<20} {:<20}".format("Step", "Output file", "Time to run [s]")
with open(txt_outputs_summary, "w+") as tf:
    tf.write(line0+"\n")
    tf.write(line1+"\n")

final_out = "gain_scale.fits"
output_names = ["group_scale.fits", "dq_init.fits", "saturation.fits", "superbias.fits", "refpix.fits",
                "lastframe.fits", "linearity.fits", "dark_current.fits", "jump.fits", "ramp_fit.fits", final_out]

# Get and save the value of the raw data root name to add at the end of calwebb_detector1
#rawdatrt = fits.getval(fits_input_uncal_file, 'rawdatrt', 0)

if not step_by_step:
    # start the timer to compute the step running time
    start_time = time.time()
    Detector1Pipeline.call(fits_input_uncal_file, config_file=calwebb_detector1_cfg)
    # end the timer to compute calwebb_spec2 running time
    end_time = time.time() - start_time  # this is in seconds
    print(" * calwebb_detector1 took "+repr(end_time)+" seconds to finish *")
    total_time = "{:<18} {:<20} {:<20}".format("", "total time = ", repr(end_time))
    tot_time_min = round((end_time/60.0), 2)
    for outnm in output_names:
        line2write = "{:<18} {:<20} {:<20}".format(outnm.replace(".fits", ""), outnm, "")
        with open(txt_outputs_summary, "a") as tf:
            tf.write(line2write+"\n")
else:
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
                    print("Checking for next step... ")
                    print (" * Completion keyword: ", comp_keys[i-j], completion_key_val)
                    if "SKIPPED" in completion_key_val:
                        j += 1
                    elif "COMPLETE" in completion_key_val:
                        continue_while = False
        print ("\n-> Running step: ", stp, "   with input file: ", step_input_file)
        print ("   output will be saved as: ", output_names[i])

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
                print ("File has only 1 integration.")

        # end the timer to compute cal_detector1 running time
        end_time = repr(time.time() - start_time)   # this is in seconds
        end_time_total.append(time.time() - start_time)
        step = output_names[i].replace(".fits", "")
        print(" * calwebb_detector1 step ", step, " took "+end_time+" seconds to finish * \n")
        if (time.time() - start_time) > 60.0:
            end_time_min = (time.time() - start_time)/60.  # this is in minutes
            print ("    ( = "+repr(end_time_min)+" minutes )")
            end_time = end_time+"  ="+repr(round(end_time_min, 2))+"min"

        # record results in text file
        line2write = "{:<18} {:<20} {:<20}".format(step, output_names[i], end_time)
        with open(txt_outputs_summary, "a") as tf:
            tf.write(line2write+"\n")

    # record total time in text file
    tot_time_sec = sum(end_time_total)
    tot_time_min = round((tot_time_sec/60.0), 2)
    total_time = "{:<18} {:>20} {:<20}".format("", "total_time  ", repr(tot_time_sec)+"  ="+repr(tot_time_min)+"min")

# record total time in text file
with open(txt_outputs_summary, "a") as tf:
    tf.write(total_time+"\n")
print ("\n ** Calwebb_detector 1 took "+repr(tot_time_min)+" minutes to complete **")

# Move fits products to working dir
fits_list = glob("*.fits")
if len(fits_list) >= 1:
    print("Output fits files are located at: ", working_dir)
    for ff in fits_list:
        if "step_" in ff:
            ff_newname = os.path.join(working_dir, ff.replace("step_", ""))
        else:
            ff_newname = working_dir
        subprocess.run(["mv", ff, ff_newname])
    # move the text file too
    subprocess.run(["mv", txt_outputs_summary, working_dir])
else:
    print("No fits files detected after calwbb_detector1 finished. Exiting script.")

print ("Script  run_cal_detector1.py  finished.")
