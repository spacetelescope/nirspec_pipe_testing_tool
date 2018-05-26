
"""
py.test module for unit testing the assign_wcs step.
"""

import os
import subprocess
import time
import pytest
from astropy.io import fits
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline

from . import assign_wcs_utils
from .. import core_utils
from .. auxiliary_code import change_filter_opaque2science
from .. auxiliary_code import compare_wcs_ifu
from .. auxiliary_code import compare_wcs_fs
from .. auxiliary_code import compare_wcs_mos


# print pipeline version
import jwst
print("\n*** Using jwst pipeline version: ", jwst.__version__, " *** \n")

# Set up the fixtures needed for all of the tests, i.e. open up all of the FITS files

# Default names of pipeline input and output files
@pytest.fixture(scope="module")
def set_inandout_filenames(config):
    step = "assign_wcs"
    step_dict = dict(config.items("steps"))
    True_steps_suffix_map = config.get("calwebb_spec2_input_file", "True_steps_suffix_map")
    data_directory = config.get("calwebb_spec2_input_file", "data_directory")
    working_directory = config.get("calwebb_spec2_input_file", "working_directory")
    initial_input_file_basename = config.get("calwebb_spec2_input_file", "input_file")
    initial_input_file = os.path.join(data_directory, initial_input_file_basename)
    print("\n Using initial input file: ", initial_input_file , "\n")
    pytests_directory = os.getcwd()
    True_steps_suffix_map = os.path.join(pytests_directory, True_steps_suffix_map)
    suffix_and_filenames = core_utils.get_step_inandout_filename(step, initial_input_file, step_dict, working_directory)
    in_file_suffix, out_file_suffix, step_input_filename, step_output_filename = suffix_and_filenames
    return step, step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


# fixture to read the output file header
@pytest.fixture(scope="module")
def output_hdul(set_inandout_filenames, config):
    set_inandout_filenames_info = core_utils.read_info4outputhdul(config, set_inandout_filenames)
    step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix = set_inandout_filenames_info

    # check if the filter is to be changed
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    if change_filter_opaque:
        _, step_input_file = change_filter_opaque2science.change_filter_opaque(step_input_file)
        print (" * With FILTER=OPAQUE, the calwebb_spec2 will run up to the extract_2d step. Further steps will be skipped. \n")

    stp = AssignWcsStep()
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    skip_runing_pipe_step = config.getboolean("tests_only", "_".join((step, "tests")))
    esa_files_path = config.get("esa_intermediary_products", "esa_files_path")
    msa_conf_name = config.get("esa_intermediary_products", "msa_conf_name")
    wcs_threshold_diff = config.get("additional_arguments", "wcs_threshold_diff")
    save_wcs_plots = config.getboolean("additional_arguments", "save_wcs_plots")
    # if run_calwebb_spec2 is True calwebb_spec2 will be called, else individual steps will be ran
    step_completed = False
    end_time = '0.0'
    # Check if data is IFU that the Image Model keyword is correct
    mode_used = config.get("calwebb_spec2_input_file", "mode_used")
    if mode_used == "IFU":
        DATAMODL = fits.getval(step_input_file, "DATAMODL", 0)
        if DATAMODL != "IFUImageModel":
            fits.setval(step_input_file, "DATAMODL", 0, value="IFUImageModel")
            print("DATAMODL keyword changed to IFUImageModel.")
    # get the MSA shutter configuration file full path only for MOS data
    inhdu = core_utils.read_hdrfits(step_input_file, info=False, show_hdr=False)
    if core_utils.check_MOS_true(inhdu):
        msa_shutter_conf = config.get("esa_intermediary_products", "msa_conf_name")
        # check if the configuration shutter file name is in the header of the fits file and if not add it
        msametfl = fits.getval(step_input_file, "MSAMETFL", 0)
        if os.path.basename(msa_shutter_conf) != msametfl:
            msametfl = os.path.basename(msa_shutter_conf)
            fits.setval(step_input_file, "MSAMETFL", 0, value=msametfl)
    # run the pipeline
    if run_calwebb_spec2:
        print ("*** Will run calwebb_spec2... ")
        # create the map
        full_run_map = "full_run_map.txt"
        assign_wcs_utils.create_map_from_full_run(full_run_map, step_input_file)
        # get the name of the configuration file and run the pipeline
        calwebb_spec2_cfg = config.get("run_calwebb_spec2_in_full", "calwebb_spec2_cfg")
        input_file = config.get("calwebb_spec2_input_file", "input_file")
        final_output_name = input_file.replace(".fits", "_cal.fits")
        working_directory = config.get("calwebb_spec2_input_file", "working_directory")
        workingdir_final_output = os.path.join(working_directory, final_output_name)
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
        subprocess.run(["mv", final_output_name, workingdir_final_output])
        local_step_output_file = input_file.replace(".fits", "_assign_wcs.fits")
        step_output_file = os.path.join(working_directory, local_step_output_file)
        subprocess.run(["mv", local_step_output_file, step_output_file])
        # read the assign wcs fits file
        hdul = core_utils.read_hdrfits(step_output_file, info=True, show_hdr=True)
        scihdul = core_utils.read_hdrfits(step_output_file, info=False, show_hdr=False, ext=1)
        return hdul, step_output_file, msa_conf_name, esa_files_path, wcs_threshold_diff, save_wcs_plots
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
                    # get the right configuration files to run the step
                    local_pipe_cfg_path = config.get("calwebb_spec2_input_file", "local_pipe_cfg_path")
                    # start the timer to compute the step running time
                    start_time = time.time()
                    if local_pipe_cfg_path == "pipe_source_tree_code":
                        result = stp.call(step_input_file)
                    else:
                        result = stp.call(step_input_file, config_file=local_pipe_cfg_path+'/assign_wcs.cfg')
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
                return hdul, step_output_file, msa_conf_name, esa_files_path, wcs_threshold_diff, save_wcs_plots
            else:
                print("Skipping step. Intput file "+step_input_file+" does not exit.")
                core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
                pytest.skip("Skipping "+step+" because the input file does not exist.")
        else:
            core_utils.add_completed_steps(txt_name, step, outstep_file_suffix, step_completed, end_time)
            pytest.skip("Skipping "+step+". Step set to False in configuration file.")



### THESE FUNCTIONS ARE TO VALIDATE THE WCS STEP

# fixture to validate the WCS
@pytest.fixture(scope="module")
def validate_wcs(output_hdul):
    # get the input information for the wcs routine
    hdu = output_hdul[0]
    infile_name = output_hdul[1]
    msa_conf_name = output_hdul[2]
    esa_files_path = output_hdul[3]

    # define the threshold difference between the pipeline output and the ESA files for the pytest to pass or fail
    threshold_diff = float(output_hdul[4])

    # save the output plots
    save_wcs_plots = output_hdul[5]

    # show the figures
    show_figs = False

    if core_utils.check_FS_true(hdu) or core_utils.check_BOTS_true(hdu):
        result = compare_wcs_fs.compare_wcs(infile_name, esa_files_path=esa_files_path, show_figs=show_figs,
                                            save_figs=save_wcs_plots, threshold_diff=threshold_diff, debug=False)

    elif core_utils.check_MOS_true(hdu):
        result = compare_wcs_mos.compare_wcs(infile_name, esa_files_path=esa_files_path, msa_conf_name=msa_conf_name,
                                             show_figs=show_figs, save_figs=save_wcs_plots,
                                             threshold_diff=threshold_diff, debug=False)

    elif core_utils.check_IFU_true(hdu):
        result = compare_wcs_ifu.compare_wcs(infile_name, esa_files_path=esa_files_path, show_figs=show_figs,
                                            save_figs=save_wcs_plots, threshold_diff=threshold_diff, debug=False)
    else:
        pytest.skip("Skipping pytest: The fits file is not FS, MOS, or IFU. PTT does not yet include the routine to verify this kind of file.")

    return result


# Unit tests

# reference files from running calwebb_spec1
def test_rmask_rfile(output_hdul):
    result = assign_wcs_utils.rmask_rfile_is_correct(output_hdul)
    assert not result, result

def test_saturation_rfile(output_hdul):
    result = assign_wcs_utils.saturation_rfile_is_correct(output_hdul)
    assert not result, result

def test_superbias_rfile(output_hdul):
    result = assign_wcs_utils.superbias_rfile_is_correct(output_hdul)
    assert not result, result

def test_linearity_rfile(output_hdul):
    result = assign_wcs_utils.linearity_rfile_is_correct(output_hdul)
    assert not result, result

def test_dark_rfile(output_hdul):
    result = assign_wcs_utils.dark_rfile_is_correct(output_hdul)
    assert not result, result

def test_readnoise_rfile(output_hdul):
    result = assign_wcs_utils.readnoise_rfile_is_correct(output_hdul)
    assert not result, result

def test_gain_rfile(output_hdul):
    result = assign_wcs_utils.gain_rfile_is_correct(output_hdul)
    assert not result, result



# reference files specific to the WCS step
def test_camera_rfile(output_hdul):
    result = assign_wcs_utils.camera_rfile_is_correct(output_hdul)
    assert not result, result

def test_colimator_rfile(output_hdul):
    result = assign_wcs_utils.colimator_rfile_is_correct(output_hdul)
    assert not result, result

def test_disperser_rfile(output_hdul):
    result = assign_wcs_utils.disperser_rfile_is_correct(output_hdul)
    assert not result, result

def test_fore_rfile(output_hdul):
    result = assign_wcs_utils.fore_rfile_is_correct(output_hdul)
    assert not result, result

def test_fpa_rfile(output_hdul):
    result = assign_wcs_utils.fpa_rfile_is_correct(output_hdul)
    assert not result, result

def test_ifufore_rfile(output_hdul):
    result = assign_wcs_utils.ifufore_rfile_is_correct(output_hdul)
    assert not result, result

def test_ifupost_rfile(output_hdul):
    result = assign_wcs_utils.ifupost_rfile_is_correct(output_hdul)
    assert not result, result

def test_ifuslicer_rfile(output_hdul):
    result = assign_wcs_utils.ifuslicer_rfile_is_correct(output_hdul)
    assert not result, result

def test_msa_rfile(output_hdul):
    result = assign_wcs_utils.msa_rfile_is_correct(output_hdul)
    assert not result, result

def test_ote_rfile(output_hdul):
    result = assign_wcs_utils.ote_rfile_is_correct(output_hdul)
    assert not result, result

def test_wavran_rfile(output_hdul):
    result = assign_wcs_utils.wavran_rfile_is_correct(output_hdul)
    assert not result, result


# other tests specific to the WCS step
def test_wavstart_exists(output_hdul):
    assert assign_wcs_utils.wavstart_exists(output_hdul[1]), "The keyword WAVSTART was not added to the header."

def test_sporder_exists(output_hdul):
    assert assign_wcs_utils.sporder_exists(output_hdul[1]), "The keyword SPORDER was not added to the header."

def test_s_wcs_exists(output_hdul):
    assert assign_wcs_utils.s_wcs_exists(output_hdul[0]), "The keyword S_WCS was not added to the header --> extract_2d step was not completed."

def test_validate_wcs(output_hdul):
    assert validate_wcs(output_hdul), "Output value from compare_wcs.py is greater than threshold."


