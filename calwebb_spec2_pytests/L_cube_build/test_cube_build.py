
"""
py.test module for unit testing the cube_build step.
"""

import os
import subprocess
import time
from glob import glob

import pytest
from astropy.io import fits
from jwst.cube_build.cube_build_step import CubeBuildStep

from .. auxiliary_code import change_filter_opaque2science
from . import cube_build_utils
from .. import core_utils



# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "cube_build"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    run_pipe_step = config.getboolean("run_pipe_steps", step)
    run_pytests = config.getboolean("run_pytest", "_".join((step, "tests")))
    # Only run step if data is IFU
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    end_time = '0.0'
    if core_utils.check_IFU_true(inhdu):
        stp = CubeBuildStep()
        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False

        # check if the filter is to be changed
        change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
        if change_filter_opaque:
            is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file, step=step)
            if is_filter_opaque:
                print ("With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Flat Field pytest now set to Skip.")
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

        if run_calwebb_spec2:
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            return hdul, step_output_file, run_pytests
        else:
            if config.getboolean("steps", step):
                print ("*** Step "+step+" set to True")
                if os.path.isfile(step_input_file):
                    if run_pipe_step:
                        # check that previous pipeline steps were run up to this point
                        core_utils.check_completed_steps(step, step_input_file)

                        # get the right configuration files to run the step
                        local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                        # start the timer to compute the step running time
                        start_time = time.time()
                        if local_pipe_cfg_path == "pipe_source_tree_code":
                            result = stp.call(step_input_file)
                        else:
                            result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/cube_build.cfg')
                        result.save(step_output_file)
                        # end the timer to compute the step running time
                        end_time = repr(time.time() - start_time)   # this is in seconds
                        print("Step "+step+" took "+end_time+" seconds to finish")
                    # determine the specific output of the cube step
                    filt = fits.getval(step_input_file, 'filter')
                    grat = fits.getval(step_input_file, 'grating')
                    gratfilt = grat+"-"+filt+"_s3d"
                    specific_output_file = glob(step_output_file.replace('cube.fits', (gratfilt+'*.fits').lower()))[0]
                    cube_suffix = specific_output_file.split('photom_')[-1].replace('.fits', '')
                    # record info
                    step_completed = True
                    core_utils.add_completed_steps(txt_name, step, "_"+cube_suffix, step_completed, end_time)
                    hdul = core_utils.read_hdrfits(specific_output_file, info=False, show_hdr=False)
                    return hdul, step_output_file, run_pytests
                else:
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Skipping "+step+" because the input file does not exist.")
            else:
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+". Step set to False in configuration file.")
    else:
        pytest.skip("Skipping "+step+" because data is not IFU.")



# Unit tests

def test_s_ifucub_exists(output_hdul):
    # want to run this pytest?
    run_pytests = output_hdul[2]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        assert cube_build_utils.s_ifucub_exists(output_hdul[0]), "The keyword S_IFUCUB was not added to the header --> IFU cube_build step was not completed."

