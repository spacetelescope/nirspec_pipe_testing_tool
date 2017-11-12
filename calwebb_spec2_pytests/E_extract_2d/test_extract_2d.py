from __future__ import print_function, division

"""
py.test module for unit testing the extract_2d step.
"""

import pytest
import os
from jwst.extract_2d.extract_2d_step import Extract2dStep

from .. import core_utils
from . import extract_2d_utils
from ..auxiliary_code import compare_wcs_fs
from ..auxiliary_code import compare_wcs_mos
from ..auxiliary_code import compare_wcs_ifu


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "extract_2d"
    step_dict = dict(config.items("steps"))
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
    True_steps_suffix_map = config.get("calwebb_spec2_input_file", "True_steps_suffix_map")
    txt_name = os.path.join(working_directory, True_steps_suffix_map)
    steps_list, suffix_list, completion_list = core_utils.read_True_steps_suffix_map(txt_name)
    step_input_filename = core_utils.get_correct_input_step_filename(initial_input_file, steps_list,
                                                                     suffix_list, completion_list)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, step_dict)
    in_file_suffix, out_file_suffix, _, _ = suffix_and_filenames
    step_output_filename = step_input_filename.replace(".fits", out_file_suffix+".fits")
    print ("step_input_filename = ", step_input_filename)
    print ("step_output_filename = ", step_output_filename)
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    step = set_inandout_filenames[0]
    step_input_filename = set_inandout_filenames[1]
    output_file = set_inandout_filenames[2]
    outstep_file_suffix = set_inandout_filenames[4]
    True_steps_suffix_map = set_inandout_filenames[5]
    txt_name = os.path.join(working_directory, True_steps_suffix_map)
    step_input_file = os.path.join(working_directory, step_input_filename)
    step_output_file = os.path.join(working_directory, output_file)
    stp = Extract2dStep()
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    esa_files_path = config.get("esa_intermediary_products", "esa_files_path")
    wcs_threshold_diff = config.get("additional_arguments", "wcs_threshold_diff")
    save_wcs_plots = config.getboolean("additional_arguments", "save_wcs_plots")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    if not run_calwebb_spec2:
        if config.getboolean("steps", step):
            print ("*** Step "+step+" set to True")
            if os.path.isfile(step_input_file):
                #result = stp.call(step_input_file)
                #result.save(step_output_file)
                step_completed = True
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                return hdul, step_output_file, esa_files_path, wcs_threshold_diff, save_wcs_plots
            else:
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
                pytest.skip("Skiping "+step+" because the input file does not exist.")
        else:
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed)
            pytest.skip("Skiping "+step+". Step set to False in configuration file.")



### THESE FUNCTIONS ARE TO VALIDATE BOTH THE WCS AND THE 2D_EXTRACT STEPS

# fixture to validate the WCS and extract 2d steps
@pytest.fixture(scope="module")
def validate_wcs_extract2d(output_hdul):
    # get the input information for the wcs routine
    hdu = output_hdul[0]
    infile_name = output_hdul[1]
    esa_files_path = output_hdul[2]

    # define the threshold difference between the pipeline output and the ESA files for the pytest to pass or fail
    threshold_diff = float(output_hdul[3])

    # save the output plots
    save_wcs_plots = output_hdul[4]

    # show the figures
    show_figs = False

    if core_utils.check_FS_true(hdu):
        # Find what slit the data corresponds to
        ext, slit = core_utils.find_which_slit(hdu)
        if (slit is not None) or (slit != "NULL"):
            median_diff = compare_wcs_fs.compare_wcs(infile_name, esa_files_path=esa_files_path,
                                                     auxiliary_code_path=None, plot_names=None,
                                                     show_figs=show_figs, save_figs=save_wcs_plots,
                                                  threshold_diff=threshold_diff)

    elif core_utils.check_MOS_true(hdu):
        msa_conf_root = output_hdul.msa_conf_root
        median_diff = compare_wcs_mos.compare_wcs(infile_name, msa_conf_root=msa_conf_root,
                                                  esa_files_path=esa_files_path, auxiliary_code_path=None,
                                                  plot_names=None, show_figs=show_figs, save_figs=save_wcs_plots,
                                                  threshold_diff=threshold_diff)

    elif core_utils.check_IFU_true(hdu):
        median_diff = compare_wcs_ifu.compare_wcs(infile_name, esa_files_path=esa_files_path, auxiliary_code_path=None,
                                                  plot_names=None, show_figs=show_figs, save_figs=save_wcs_plots,
                                                  threshold_diff=threshold_diff)
    else:
        pytest.skip("Skipping pytest: The fits file is not FS, MOS, or IFU. \n "
                    "This pytest tool does not yet include the routine to verify this kind of file.")
    return median_diff



### Unit tests

def test_s_ext2d_exists(output_hdul):
    assert extract_2d_utils.s_ext2d_exists(output_hdul[0])


def test_validate_wcs_extract2d(output_hdul):
    assert validate_wcs_extract2d(output_hdul)
