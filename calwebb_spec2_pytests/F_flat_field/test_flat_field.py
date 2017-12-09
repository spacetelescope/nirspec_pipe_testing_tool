
"""
py.test module for unit testing the flat field step.
"""

import pytest
import os
import time

from jwst.flatfield.flat_field_step import FlatFieldStep
from .. import core_utils
from . import flat_field_utils
from .. auxiliary_code import flattest_mos
from .. auxiliary_code import flattest_fs
from .. auxiliary_code import flattest_ifu


# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(request, config):
    step = "flat_field"
    step_info = core_utils.set_inandout_filenames(step, config)
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map = step_info
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info
    stp = FlatFieldStep()
    msa_conf_root = config.get("esa_intermediary_products", "msa_conf_root")
    dflat_path = config.get("esa_intermediary_products", "dflat_path")
    sflat_path = config.get("esa_intermediary_products", "sflat_path")
    fflat_path = config.get("esa_intermediary_products", "fflat_path")
    flattest_threshold_diff = config.get("additional_arguments", "flattest_threshold_diff")
    save_flattest_plot = config.getboolean("additional_arguments", "save_flattest_plot")
    write_flattest_files = config.getboolean("additional_arguments", "write_flattest_files")
    flattest_paths = [step_output_file, msa_conf_root, dflat_path, sflat_path, fflat_path]
    flattest_switches = [flattest_threshold_diff, save_flattest_plot, write_flattest_files]
    skip_runing_pipe_step = config.getboolean("tests_only", "_".join((step, "tests")))
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'
    if not run_calwebb_spec2:
        if config.getboolean("steps", step):
            print ("*** Step "+step+" set to True")
            if os.path.isfile(step_input_file):
                if not skip_runing_pipe_step:
                    # start the timer to compute the step running time
                    start_time = time.time()
                    result = stp.call(step_input_file)
                    result.save(step_output_file)
                    # end the timer to compute the step running time
                    end_time = time.time() - start_time   # this is in seconds
                    print("Step "+step+" took "+end_time+" seconds to finish")
                step_completed = True
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
                return hdul, flattest_paths, flattest_switches
            else:
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because the input file does not exist.")
        else:
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+". Step set to False in configuration file.")



### THESE FUNCTIONS ARE TO VALIDATE BOTH THE WCS AND THE 2D_EXTRACT STEPS

# fixture to validate the WCS and extract 2d steps
@pytest.fixture(scope="module")
def validate_flat_field(output_hdul):
    # get the input information for the wcs routine
    hdu = output_hdul[0]
    step_output_file, msa_conf_root, dflatref_path, sfile_path, fflat_path = output_hdul[1]
    flattest_threshold_diff, save_flattest_plot, write_flattest_files = output_hdul[2]

    # show the figures
    show_figs = False

    if core_utils.check_FS_true(hdu):
        # Find what slit the data corresponds to
        ext, slit = core_utils.find_which_slit(hdu)
        if (slit is not None) or (slit != "NULL"):
            median_diff = flattest_fs.flattest(step_output_file, dflatref_path=dflatref_path, sfile_path=sfile_path,
                                                fflat_path=fflat_path, writefile=write_flattest_files,
                                                show_figs=show_figs, save_figs=save_flattest_plot, plot_name=None,
                                                threshold_diff=flattest_threshold_diff, debug=False)

    elif core_utils.check_MOS_true(hdu):
        median_diff = flattest_mos.flattest(step_output_file, dflatref_path=dflatref_path, sfile_path=sfile_path,
                                           fflat_path=fflat_path, msa_conf_root=msa_conf_root,
                                           writefile=write_flattest_files,
                                           show_figs=show_figs, save_figs=save_flattest_plot, plot_name=None,
                                           threshold_diff=flattest_threshold_diff, debug=False)

    elif core_utils.check_IFU_true(hdu):
        median_diff = flattest_ifu.flattest(step_output_file, dflatref_path=dflatref_path, sfile_path=sfile_path,
                                            fflat_path=fflat_path, writefile=write_flattest_files,
                                            mk_all_slices_plt=False, show_figs=show_figs,
                                            save_figs=save_flattest_plot, plot_name=None,
                                            threshold_diff=flattest_threshold_diff, debug=False)

    else:
        pytest.skip("Skipping pytest: The input fits file is not FS, MOS, or IFU. This tool does not yet include the "
                    "routine to verify this kind of file.")
    return median_diff



# Unit tests

def test_s_flat_exists(output_hdul):
    assert flat_field_utils.s_flat_exists(output_hdul[0]), "The keyword S_FLAT was not added to the header --> flat_field step was not completed."

def test_validate_flat_field(output_hdul):
    assert validate_flat_field(output_hdul), "Output value from flattest.py is greater than threshold."
