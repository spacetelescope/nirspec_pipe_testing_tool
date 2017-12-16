import os
import time
import argparse
import subprocess
import configparser
from glob import glob

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


def get_caldet1cfg_and_workingdir():
    """
    This function gets the path to where the pipeline testing tool lives
    Returns:
        calwebb_detector1_cfg = string, path to where the config file lives
    """
    config = configparser.ConfigParser()
    config.read(['../calwebb_spec2_pytests/cwspec2_config.cfg'])
    pipe_testing_tool_path = config.get("calwebb_spec2_input_file", "pipe_testing_tool_path")
    calwebb_detector1_cfg = os.path.join(pipe_testing_tool_path, "utils/calwebb_detector1.cfg")
    working_dir = config.get("calwebb_spec2_input_file", "working_directory")
    return calwebb_detector1_cfg, working_dir



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
calwebb_detector1_cfg, working_dir = get_caldet1cfg_and_workingdir()
print ("Using this configuration file: ", calwebb_detector1_cfg)

if not step_by_step:
    # start the timer to compute the step running time
    start_time = time.time()
    Detector1Pipeline.call(fits_input_uncal_file, config_file=calwebb_detector1_cfg)
    # end the timer to compute calwebb_spec2 running time
    end_time = repr(time.time() - start_time)   # this is in seconds
    print(" * calwebb_detector1 took "+end_time+" seconds to finish *")
else:
    # steps to be ran, in order
    steps_to_run = [GroupScaleStep(), DQInitStep(), SaturationStep(), SuperBiasStep(), RefPixStep(), RSCD_Step(),
                    LastFrameStep(), LinearityStep(), DarkCurrentStep(), JumpStep(), RampFitStep(), GainScaleStep()]
    final_out = fits_input_uncal_file.replace(".fits", "_caldet1.fits")
    output_names = ["group_scale.fits", "dq_init.fits", "saturation.fits", "superbias.fits", "refpix.fits", "rscd.fits",
                    "lastframe.fits", "linearity.fits", "dark_current.fits", "jump.fits", "ramp_fit.fits", final_out]

    # start the timer to compute the step running time
    start_time = time.time()

    # run the pipeline step by step
    for i, stp_instance in enumerate(steps_to_run):
        stp = stp_instance
        if i == 0:
            step_input_file = fits_input_uncal_file
        else:
            step_input_file = output_names[i-1]
        result = stp.call(step_input_file)
        result.save(output_names[i])

    # end the timer to compute calwebb_spec2 running time
    end_time = repr(time.time() - start_time)   # this is in seconds
    print(" * calwebb_detector1 step by step took "+end_time+" seconds to finish *")


# Move fits products to working dir
fits_list = glob("*.fits")
if len(fits_list) >= 1:
    print("Output fits files are located at: ", working_dir)
    #subprocess.run(["mv", "*.fits", working_dir])
else:
    print("No fits files detected after calwbb_detector1 finished. Exiting script.")

print ("Script  run_cal_detector1.py  finished.")