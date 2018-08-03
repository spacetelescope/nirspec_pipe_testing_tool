
"""
py.test module for unit testing the flat field step.
"""

import os
import subprocess
import time

import pytest
from jwst.flatfield.flat_field_step import FlatFieldStep

from . import flat_field_utils
from .. import core_utils
from .. auxiliary_code import flattest_fs
from .. auxiliary_code import flattest_ifu
from .. auxiliary_code import flattest_mos
from .. auxiliary_code import change_filter_opaque2science



# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed


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
    msa_shutter_conf = config.get("esa_intermediary_products", "msa_conf_name")
    msametfl = os.path.basename(msa_shutter_conf)
    dflat_path = config.get("esa_intermediary_products", "dflat_path")
    sflat_path = config.get("esa_intermediary_products", "sflat_path")
    fflat_path = config.get("esa_intermediary_products", "fflat_path")
    flattest_threshold_diff = config.get("additional_arguments", "flattest_threshold_diff")
    save_flattest_plot = config.getboolean("additional_arguments", "save_flattest_plot")
    write_flattest_files = config.getboolean("additional_arguments", "write_flattest_files")
    flattest_paths = [step_output_file, msa_shutter_conf, dflat_path, sflat_path, fflat_path]
    flattest_switches = [flattest_threshold_diff, save_flattest_plot, write_flattest_files]
    run_pipe_step = config.getboolean("run_pipe_steps", step)
    run_pytests = config.getboolean("run_pytest", "_".join((step, "tests")))
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'

    # check if the filter is to be changed
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    if change_filter_opaque:
        is_filter_opaque, step_input_filename = change_filter_opaque2science.change_filter_opaque(step_input_file, step=step)
        if is_filter_opaque:
            print ("With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Flat Field pytest now set to Skip.")
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+" because FILTER=OPAQUE.")

    # get the MSA shutter configuration file full path only for MOS data
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if run_calwebb_spec2:
        hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
        flattest_paths = [step_output_file, msa_shutter_conf, dflat_path, sflat_path, fflat_path]
        return hdul, step_output_file, flattest_paths, flattest_switches, run_pytests
    else:
        if os.path.isfile(step_input_file):
            if run_pipe_step:
                print ("*** Step "+step+" set to True")
                stp = FlatFieldStep()

                # check that previous pipeline steps were run up to this point
                core_utils.check_completed_steps(step, step_input_file)

                flat_suffix = "intflat"
                if core_utils.check_MOS_true(inhdu):
                    # copy the MSA shutter configuration file into the pytest directory
                    subprocess.run(["cp", msa_shutter_conf, "."])
                # start the timer to compute the step running time
                ontheflyflat = step_output_file.replace("flat_field.fits", "intflat.fits")
                print ("Step product will be saved as: ", step_output_file)
                print ("on-the-fly flat will be saved as: ", ontheflyflat)
                # get the right configuration files to run the step
                local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                if core_utils.check_IFU_true(inhdu):
                    flat_suffix = "intflat"
                # start the timer to compute the step running time
                start_time = time.time()
                if local_pipe_cfg_path == "pipe_source_tree_code":
                    stp.call(step_input_file, output_file=step_output_file, flat_suffix=flat_suffix)
                              #override_dflat="/grp/crds/jwst/references/jwst/jwst_nirspec_dflat_0001.fits",
                              #override_fflat="/grp/crds/jwst/references/jwst/jwst_nirspec_fflat_0015.fits",
                              #override_sflat="/grp/crds/jwst/references/jwst/jwst_nirspec_sflat_0034.fits")
                else:
                    stp.call(step_input_file, output_file=step_output_file, flat_suffix=flat_suffix,
                             config_file=local_pipe_cfg_path+'/flat_field.cfg')
                # end the timer to compute the step running time
                end_time = repr(time.time() - start_time)   # this is in seconds
                print("Step "+step+" took "+end_time+" seconds to finish")
                # move the on-the-fly flat to the working directory
                subprocess.run(["mv", os.path.basename(ontheflyflat), ontheflyflat])
                # raname and move the flat_field output
                subprocess.run(["mv", os.path.basename(step_output_file).replace("_flat_field.fits", "_flatfieldstep.fits"),
                                step_output_file])
                if core_utils.check_MOS_true(inhdu):
                    # remove the copy of the MSA shutter configuration file
                    subprocess.run(["rm", msametfl])
            else:
                print("Skipping running pipeline step ", step)
                # add the running time for this step
                working_directory = config.get("calwebb_spec2_input_file", "working_directory")
                end_time = core_utils.get_stp_run_time_from_screenfile(step, working_directory)
            step_completed = True
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            hdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False)
            return hdul, step_output_file, flattest_paths, flattest_switches, run_pytests

        else:
            print (" The input file does not exist. Skipping step.")
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+" because the input file does not exist.")



### THESE FUNCTIONS ARE TO VALIDATE BOTH THE WCS AND THE 2D_EXTRACT STEPS

# fixture to validate the WCS and extract 2d steps
@pytest.fixture(scope="module")
def validate_flat_field(output_hdul):
    # get the input information for the wcs routine
    hdu = output_hdul[0]
    step_output_file, msa_shutter_conf, dflatref_path, sfile_path, fflat_path = output_hdul[2]
    flattest_threshold_diff, save_flattest_plot, write_flattest_files = output_hdul[3]

    # show the figures
    show_figs = False

    if core_utils.check_FS_true(hdu) or core_utils.check_BOTS_true(hdu):
        median_diff, msg = flattest_fs.flattest(step_output_file, dflatref_path=dflatref_path, sfile_path=sfile_path,
                                                fflat_path=fflat_path, writefile=write_flattest_files,
                                                show_figs=show_figs, save_figs=save_flattest_plot, plot_name=None,
                                                threshold_diff=flattest_threshold_diff, debug=False)

    elif core_utils.check_MOS_true(hdu):
        median_diff, msg = flattest_mos.flattest(step_output_file, dflatref_path=dflatref_path, sfile_path=sfile_path,
                                               fflat_path=fflat_path, msa_shutter_conf=msa_shutter_conf,
                                               writefile=write_flattest_files,
                                               show_figs=show_figs, save_figs=save_flattest_plot, plot_name=None,
                                               threshold_diff=flattest_threshold_diff, debug=False)

    elif core_utils.check_IFU_true(hdu):
        median_diff, msg = flattest_ifu.flattest(step_output_file, dflatref_path=dflatref_path, sfile_path=sfile_path,
                                                fflat_path=fflat_path, writefile=write_flattest_files,
                                                mk_all_slices_plt=False, show_figs=show_figs,
                                                save_figs=save_flattest_plot, plot_name=None,
                                                threshold_diff=flattest_threshold_diff, debug=False)

    else:
        pytest.skip("Skipping pytest: The input fits file is not FS, MOS, or IFU. This tool does not yet include the "
                    "routine to verify this kind of file.")

    if median_diff == "skip":
        pytest.skip(msg)
    else:
        print(msg)

    return median_diff



# Unit tests

def test_s_flat_exists(output_hdul):
    # want to run this pytest?
    run_pytests = output_hdul[4]
    if not run_pytests:
        msg = "Skipping completion pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running completion pytest...\n")
        assert flat_field_utils.s_flat_exists(output_hdul[0]), "The keyword S_FLAT was not added to the header --> flat_field step was not completed."

def test_validate_flat_field(output_hdul):
    # want to run this pytest?
    run_pytests = output_hdul[4]
    if not run_pytests:
        msg = "Skipping validation pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running validation pytest...\n")
        assert validate_flat_field(output_hdul), "Output value from flattest.py is greater than threshold."

def test_fflat_rfile(output_hdul):
    # want to run this pytest?
    run_pytests = output_hdul[4]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running reference file pytest...\n")
        result = flat_field_utils.fflat_rfile_is_correct(output_hdul)
        assert not result, result

def test_sflat_sfile(output_hdul):
    # want to run this pytest?
    run_pytests = output_hdul[4]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running reference file pytest...\n")
        result = flat_field_utils.sflat_rfile_is_correct(output_hdul)
        assert not result, result

def test_dflat_dfile(output_hdul):
    # want to run this pytest?
    run_pytests = output_hdul[4]
    if not run_pytests:
        msg = "Skipping ref_file pytest: option to run Pytest is set to False in PTT_config.cfg file.\n"
        print(msg)
        pytest.skip(msg)
    else:
        print("\n * Running reference file pytest...\n")
        result = flat_field_utils.dflat_rfile_is_correct(output_hdul)
        assert not result, result

