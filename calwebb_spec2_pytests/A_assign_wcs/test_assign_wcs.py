
"""
py.test module for unit testing the assign_wcs step.
"""

import pytest
import os
import time
import subprocess
from astropy.io import fits
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep

from .. import core_utils
from . import assign_wcs_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(config):
    step = "assign_wcs"
    step_dict = dict(config.items("steps"))
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    True_steps_suffix_map = config.get("calwebb_spec2_input_file", "True_steps_suffix_map")
    pytests_directory = os.getcwd()
    True_steps_suffix_map = os.path.join(pytests_directory, True_steps_suffix_map)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, step_dict)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    stp = AssignWcsStep()
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    skip_runing_pipe_step = config.getboolean("tests_only", "_".join((step, "tests")))
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'
    # get the MSA shutter configuration file full path only for MOS data
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if core_utils.check_MOS_true(inhdu):
        msa_conf_root = config.get("esa_intermediary_products", "msa_conf_root")
        msametfl = fits.getval(step_input_file, "MSAMETFL", 0)
        msa_shutter_conf = os.path.join(msa_conf_root, msametfl)
    # run the pipeline
    if run_calwebb_spec2:
        print ("*** Will run calwebb_spec2... ")
        # create the map
        full_run_map = "full_run_map.txt"
        assign_wcs_utils.create_map_from_full_run(full_run_map, step_input_file)
        # get the name of the configuration file and run the pipeline
        calwebb_spec2_cfg = config.get("run_calwebb_spec2_in_full", "calwebb_spec2_cfg")
        input_file_basename_list = step_input_file.split("_")[:4]
        input_file_basename = "_".join((input_file_basename_list))
        final_output_name = "_".join((input_file_basename, "cal.fits"))
        final_output_name_basename = os.path.basename(final_output_name)
        if core_utils.check_MOS_true(inhdu):
            # copy the MSA shutter configuration file into the pytest directory
            subprocess.run(["cp", msa_shutter_conf, "."])
        # start the timer to compute the step running time
        start_time = time.time()
        Spec2Pipeline.call(step_input_file, config_file=calwebb_spec2_cfg)
        # end the timer to compute calwebb_spec2 running time
        end_time = repr(time.time() - start_time)   # this is in seconds
        print(" * calwebb_spec2 took "+end_time+" seconds to finish.")
        if core_utils.check_MOS_true(inhdu):
            # remove the copy of the MSA shutter configuration file
            subprocess.run(["rm", msametfl])
        # move the output file into the working directory
        print ("The final calwebb_spec2 product was saved in: ", final_output_name)
        subprocess.run(["mv", final_output_name_basename, final_output_name])
        # read the assign wcs fits file
        step_output_file = core_utils.read_completion_to_full_run_map(full_run_map, step)
        hdul = core_utils.read_hdrfits(step_output_file, info=True, show_hdr=True)
        scihdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False, ext=1)
        return hdul, scihdul
    else:
        # create the Map of file names
        assign_wcs_utils.create_completed_steps_txtfile(txt_name, step_input_file)
        if config.getboolean("steps", step):
            print ("*** Step "+step+" set to True")
            if os.path.isfile(step_input_file):
                if not skip_runing_pipe_step:
                    if core_utils.check_MOS_true(inhdu):
                        # copy the MSA shutter configuration file into the pytest directory
                        subprocess.run(["cp", msa_shutter_conf, "."])
                    # start the timer to compute the step running time
                    start_time = time.time()
                    result = stp.call(step_input_file)
                    result.save(step_output_file)
                    # end the timer to compute the step running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    print("Step "+step+" took "+end_time+" seconds to finish")
                    if core_utils.check_MOS_true(inhdu):
                        # remove the copy of the MSA shutter configuration file
                        subprocess.run(["rm", msametfl])
                step_completed = True
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False, ext=0)
                scihdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False, ext=1)
                return hdul, scihdul
            else:
                print("Skipping step. Intput file "+step_input_file+" does not exit.")
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because the input file does not exist.")
        else:
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+". Step set to False in configuration file.")



# Unit tests

def test_wavstart_exists(output_hdul):
    assert assign_wcs_utils.wavstart_exists(output_hdul[1]), "The keyword WAVSTART was not added to the header."

def test_sporder_exists(output_hdul):
    assert assign_wcs_utils.sporder_exists(output_hdul[1]), "The keyword SPORDER was not added to the header."

def test_s_wcs_exists(output_hdul):
    assert assign_wcs_utils.s_wcs_exists(output_hdul[0]), "The keyword S_WCS was not added to the header --> extract_2d step was not completed."

