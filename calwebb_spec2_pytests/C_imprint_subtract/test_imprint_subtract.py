
"""
py.test module for unit testing the imprint_subtract step.
"""

import pytest
import os
import time
import copy
import subprocess
from astropy.io import fits
from jwst.imprint.imprint_step import ImprintStep

from .. import core_utils
from . import imprint_subtract_utils


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files


# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "imprint_subtract"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    skip_runing_pipe_step = config.getboolean("tests_only", "_".join((step, "tests")))
    # Only run step if data is IFU or MSA
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    end_time = '0.0'
    run_step = config.getboolean("steps", step)
    if core_utils.check_IFU_true(inhdu) or core_utils.check_MOS_true(inhdu):
        stp = ImprintStep()
        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False
        if run_calwebb_spec2:
            # read the assign wcs fits file
            input_file = config.get("calwebb_spec2_input_file", "input_file")
            local_step_output_file = input_file.replace(".fits", "_imprint_substract.fits")
            if os.path.isfile(local_step_output_file):
                hdul = core_utils.read_hdrfits(local_step_output_file, info=False, show_hdr=False)
            else:
                pytest.skip("Skipping "+step+" because the output file does not exist.")
            # move the output file into the working directory
            working_directory = config.get("calwebb_spec2_input_file", "working_directory")
            step_output_file = os.path.join(working_directory, local_step_output_file)
            print ("Step product was saved as: ", step_output_file)
            subprocess.run(["mv", local_step_output_file, step_output_file])
            return hdul, run_step, step_input_file, step_output_file
        else:
            if run_step:
                print ("*** Step "+step+" set to True")
                if os.path.isfile(step_input_file):
                    print(" The input file ", step_input_file,"exists... will run step "+step)
                    msa_imprint_structure = config.get("additional_arguments", "msa_imprint_structure")
                    print("msa_imprint_structure file: ", msa_imprint_structure)
                    if not os.path.isfile(msa_imprint_structure):
                        print (" Need msa_imprint_structure file to continue. Step will be skipped.")
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                        pytest.skip("Skipping "+step+" because msa_imprint_structure file in the configuration file does not exist.")
                    else:
                        if not skip_runing_pipe_step:
                            # get the right configuration files to run the step
                            local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                            # start the timer to compute the step running time
                            start_time = time.time()
                            if local_pipe_cfg_path == "pipe_source_tree_code":
                                result = stp.call(step_input_file, msa_imprint_structure)
                            else:
                                result = stp.call(step_input_file, msa_imprint_structure,
                                                  config_file=local_pipe_cfg_path+'/imprint.cfg')
                            if result is not None:
                                result.save(step_output_file)
                                # end the timer to compute the step running time
                                end_time = repr(time.time() - start_time)   # this is in seconds
                                print("Step "+step+" took "+end_time+" seconds to finish")
                                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                                step_completed = True
                            else:
                                hdul = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
                        else:
                            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                            step_completed = True
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                        return hdul, run_step, step_input_file, step_output_file
                else:
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                    pytest.skip("Skipping "+step+" because the input file does not exist.")
            else:
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+". Step set to False in configuration file.")
    else:
        pytest.skip("Skipping "+step+" because data is neither IFU or MOS.")



### VALIDATION FUNCTIONS

# fixture to validate the subtraction works fine: re-run the step with the same file as msa_imprint file
@pytest.fixture(scope="module")
def check_output_is_zero(output_hdul):
    run_step = output_hdul[1]
    step_input_file = output_hdul[2]
    step_output_file = output_hdul[3]
    # Only run test if data is IFU or MSA
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if core_utils.check_IFU_true(inhdu) or core_utils.check_MOS_true(inhdu):
        if run_step:
            # set specifics for the test
            msa_imprint_structure = copy.deepcopy(step_input_file)
            result_to_check = step_output_file.replace(".fits", "_zerotest.fits")
            # run the step with the specifics
            stp = ImprintStep()
            res = stp.call(step_input_file, msa_imprint_structure)
            res.save(result_to_check)
            # check that the end product of image - image is zero
            c = fits.getdata(result_to_check)
            substraction = sum(c.flatten())
            result = False
            if substraction == 0.0:
                result = True
            # erase test output file
            subprocess.run(["rm", result_to_check])
            return result
        pytest.skip("Skipping validation check for imprint_subtract. Step set to False in configuration file.")



### Unit tests

def test_s_imprint_exists(output_hdul):
    assert imprint_subtract_utils.s_imprint_exists(output_hdul[0]), "The keyword S_IMPRINT was not added to the header --> imprint_subtract step was not completed."

def test_check_output_is_zero(output_hdul):
    assert check_output_is_zero(output_hdul), "Substraction result is not equal to zero."