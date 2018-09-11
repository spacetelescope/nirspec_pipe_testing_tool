
"""
py.test module for unit testing the extract_2d step.
"""

import pytest
import os
import time
import subprocess
from astropy.io import fits

from jwst.extract_2d.extract_2d_step import Extract2dStep
from .. import core_utils
from . import extract_2d_utils
from .. auxiliary_code import compare_wcs_mos


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed



# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "extract_2d"
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
    run_pytests_assign_wcs = config.getboolean("run_pytest", "_".join(("assign_wcs", "tests")))
    esa_files_path = config.get("esa_intermediary_products", "esa_files_path")
    msa_conf_name = config.get("esa_intermediary_products", "msa_conf_name")

    # Check if the mode used is MOS_sim and get the threshold for the assign_wcs test
    mode_used = config.get("calwebb_spec2_input_file", "mode_used")
    wcs_threshold_diff = config.get("additional_arguments", "wcs_threshold_diff")
    save_wcs_plots = config.getboolean("additional_arguments", "save_wcs_plots")

    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'
    # only do this step if data is NOT IFU
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if not core_utils.check_IFU_true(inhdu):
        if run_calwebb_spec2:
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            return hdul, step_output_file, msa_conf_name, esa_files_path, run_pytests, mode_used, wcs_threshold_diff, save_wcs_plots, run_pytests_assign_wcs
        else:
            if os.path.isfile(step_input_file):
                if run_pipe_step:
                    print ("*** Step "+step+" set to True")
                    stp = Extract2dStep()

                    # check that previous pipeline steps were run up to this point
                    core_utils.check_completed_steps(step, step_input_file)

                    # get the right configuration files to run the step
                    local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                    # start the timer to compute the step running time
                    start_time = time.time()
                    if local_pipe_cfg_path == "pipe_source_tree_code":
                        result = stp.call(step_input_file)
                    else:
                        result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/extract_2d.cfg')
                    result.save(step_output_file)
                    # end the timer to compute the step running time
                    end_time = repr(time.time() - start_time)   # this is in seconds
                    print("Step "+step+" took "+end_time+" seconds to finish")
                else:
                    print("Skipping running pipeline step ", step)
                    # add the running time for this step
                    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
                    end_time = core_utils.get_stp_run_time_from_screenfile(step, working_directory)
                step_completed = True
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                return hdul, step_output_file, msa_conf_name, esa_files_path, run_pytests, mode_used, wcs_threshold_diff, save_wcs_plots, run_pytests_assign_wcs

            else:
                print (" The input file does not exist. Skipping step.")
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skiping "+step+" because the input file does not exist.")

    else:
        pytest.skip("Skipping "+step+" because data is IFU.")



### THESE FUNCTIONS ARE TO VALIDATE BOTH THE WCS AND THE 2D_EXTRACT STEPS

# fixture to validate the WCS and extract 2d steps, only for MOS simulations
@pytest.fixture(scope="module")
def validate_wcs_extract2d(output_hdul):
    # get the input information for the wcs routine
    infile_name = output_hdul[1]
    msa_conf_name = output_hdul[2]
    esa_files_path = output_hdul[3]
    mode_used = output_hdul[5]

    # define the threshold difference between the pipeline output and the ESA files for the pytest to pass or fail
    threshold_diff = float(output_hdul[6])

    # save the output plots
    save_wcs_plots = output_hdul[7]

    # show the figures
    show_figs = False

    if mode_used == "MOS_sim":
        result = compare_wcs_mos.compare_wcs(infile_name, esa_files_path=esa_files_path, msa_conf_name=msa_conf_name,
                                             show_figs=show_figs, save_figs=save_wcs_plots,
                                             threshold_diff=threshold_diff, mode_used=mode_used, debug=False)

    else:
        pytest.skip("Skipping pytest for WCS validation: The fits file is not MOS simulated data, the validation test was done after the assign_wcs step.")

    if result == "skip":
        pytest.skip("Pytest for assign_wcs validation after extract_2d will be skipped.")

    return result



### Unit tests

def test_s_ext2d_exists(output_hdul):
    # want to run this pytest?
    run_pytests = output_hdul[4]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running completion pytest...\n")
        assert extract_2d_utils.s_ext2d_exists(output_hdul[0]), "The keyword S_EXTR2D was not added to the header --> extract_2d step was not completed."


def test_validate_wcs_extract2d(output_hdul):
    # want to run this pytest? For this particular case, check both for the extract_2d step and for assign_wcs
    run_pytests = output_hdul[4]
    assign_wcs_pytests = output_hdul[8]
    if not run_pytests and not assign_wcs_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
       print("\n * Running validation pytest...\n")
       assert validate_wcs_extract2d(output_hdul), "Output value from compare_wcs.py is greater than threshold."
