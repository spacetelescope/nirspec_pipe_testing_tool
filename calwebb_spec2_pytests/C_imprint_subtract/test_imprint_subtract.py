
"""
py.test module for unit testing the imprint_subtract step.
"""

import pytest
import os
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
    # Only run step if data is IFU
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if core_utils.check_IFU_true(inhdu) or core_utils.check_MOS_true(inhdu):
        stp = ImprintStep()
        # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
        step_completed = False
        if not run_calwebb_spec2:
            if config.getboolean("steps", step):
                print ("*** Step "+step+" set to True")
                if os.path.isfile(step_input_file):
                    print(" The input file ", step_input_file,"exists... will run step "+step)
                    msa_imprint_structure = config.get("additional_arguments", "msa_imprint_structure")
                    if not os.path.isfile(msa_imprint_structure):
                        print (" Need msa_imprint_structure file to continue. Step will be skipped.")
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                        pytest.skip("Skipping "+step+" because msa_imprint_structure file in the configuration file does not exist.")
                    else:
                        if not skip_runing_pipe_step:
                            result = stp.call(step_input_file, msa_imprint_structure)
                            if result is not None:
                                result.save(step_output_file)
                                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                                step_completed = True
                            else:
                                hdul = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
                        else:
                            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                            step_completed = True
                        core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                        return hdul
                else:
                    core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                    pytest.skip("Skipping "+step+" because the input file does not exist.")
            else:
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                pytest.skip("Skipping "+step+". Step set to False in configuration file.")
    else:
        pytest.skip("Skipping "+step+" because data is neither IFU or MOS.")


# Unit tests

def test_s_imprint_exists(output_hdul):
    assert imprint_subtract_utils.s_imprint_exists(output_hdul), "The keyword S_IMPRINT was not added to the header --> imprint_subtract step was not completed."
