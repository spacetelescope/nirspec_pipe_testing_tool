"""
This script creates the NPTT input configuration file.

Example usage:
    The code works from the terminal or called as a module.

    Terminal
        To create a the NPTT configuration file go to the output directory, then type:
        $ nptt_mk_npttconfig_file output_directory input_file mode_used raw_data_root_file

        This will create a NPTT config file in the output directory with many of the default values. Please open the
        config file and make sure that all variables are properly set. The variables are:
        output_directory = path where NPTT will place all output from the pipeline and tests
        input_file = basename of the count rate file
        mode_used = FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf
        raw_data_root_file = basename of the file that generated the count rate file

        There are several optional parameters, please see the help.

    As a module
        # import the tool
        import nirspec_pipe_testing_tool as nptt

        # set the required variables
        output_directory = string,
        input_file = string, basename of the count rate file
        mode_used = string
        raw_data_root_file = string, basename of the data where the count rate file was generated from

        # create the NPTT config file with the minimum information and all default values
        nptt.utils.mk_pnttconfig_file.mk_nptt_cfg(output_directory, input_file, mode_used, raw_data_root_file)

        # OR create the NPTT config file with specified variables

        # set the optional variables
        data_directory = string
        local_pipe_cfg_path = string, local path to the pipeline configuration files
        comparison_file_path = string, path and name of the file to compare with for assign_wcs
        msa_conf_name = string, basename of the shutter configuration file
        dflat_path = string, path and name of the D-flat
        sflat_path = string, path and name of the S-flat for this filter/grating configuration
        fflat_path = string, path and name of the F-flat for this filter configuration
        run_calwebb_spec2 = string, name of the step to run in spec2
        wcs_threshold_diff = float, acceptable difference between the comparison and the pipeline file
        save_plots = boolean
        change_filter_opaque = boolean
        extract_2d_threshold_diff = float, acceptable difference between the comparison and the pipeline file
        flattest_threshold_diff = float, acceptable difference between the comparison and the pipeline file

        # create the NPTT config file
        mk_npttconfig_file.mk_nptt_cfg(output_directory, input_file, mode_used, raw_data_root_file,
                                       data_directory=data_directory,
                                       comparison_file_path=comparison_file_path, msa_conf_name=msa_conf_name,
                                       dflat_path=dflat_path, sflat_path=sflat_path, fflat_path=fflat_path,
                                       run_calwebb_spec2=run_calwebb_spec2, wcs_threshold_diff=wcs_threshold_diff,
                                       save_plots=save_plots, change_filter_opaque=change_filter_opaque,
                                       extract_2d_threshold_diff=extract_2d_threshold_diff,
                                       flattest_threshold_diff=flattest_threshold_diff)


"""

import os
import sys
import configparser
import argparse
from astropy.io import fits


def write_nptt_cfg(calwebb_spec2_input_file, benchmark_intermediary_products, run_calwebb_spec2_in_full, run_pipe_steps,
                  run_pytest, spec3_args, additional_arguments):

    config = configparser.ConfigParser(allow_no_value=True)
    config.add_section("calwebb_spec2_input_file")
    config.set("calwebb_spec2_input_file", "output_directory", calwebb_spec2_input_file[0])
    config.set("calwebb_spec2_input_file", "data_directory", calwebb_spec2_input_file[1])
    config.set("calwebb_spec2_input_file", "input_file", calwebb_spec2_input_file[2])
    config.set("calwebb_spec2_input_file", "mode_used", calwebb_spec2_input_file[3])
    config.set("calwebb_spec2_input_file", "change_filter_opaque", calwebb_spec2_input_file[4])
    config.set("calwebb_spec2_input_file", "raw_data_root_file", calwebb_spec2_input_file[5])
    config.set("calwebb_spec2_input_file", "local_pipe_cfg_path", calwebb_spec2_input_file[6])

    config.add_section("benchmark_intermediary_products")
    esa_files_path, msa_conf_name, dflat_path, sflat_path, fflat_path, truth_awcs, truth_e2d, msa_flag_opref = benchmark_intermediary_products
    config.set("benchmark_intermediary_products", "compare_assign_wcs_and_extract_2d_with_esa", "True")
    config.set("benchmark_intermediary_products", "esa_files_path", esa_files_path)
    config.set("benchmark_intermediary_products", "# the 'truths' files (or benchmark file to compare to) is expected "
                                                  "to be in the data_directory", None)
    config.set("benchmark_intermediary_products", "truth_file_assign_wcs", truth_awcs)
    config.set("benchmark_intermediary_products", "truth_file_extract_2d", truth_e2d)
    config.set("benchmark_intermediary_products", "# other necessary paths or files", None)
    config.set("benchmark_intermediary_products", "msa_conf_name", msa_conf_name)
    config.set("benchmark_intermediary_products", "msa_flagging_operability_ref", msa_flag_opref)
    config.set("benchmark_intermediary_products", "dflat_path", dflat_path)
    config.set("benchmark_intermediary_products", "sflat_path", sflat_path)
    config.set("benchmark_intermediary_products", "fflat_path", fflat_path)

    config.add_section("run_calwebb_spec2_in_full")
    config.set("run_calwebb_spec2_in_full", "# options for run_calwebb_spec2: True (will run in full), "
                                            "False (run individual steps), skip (go to spec3)", None)
    config.set("run_calwebb_spec2_in_full", "run_calwebb_spec2", run_calwebb_spec2_in_full)

    config.add_section("calwebb_spec3")
    config.set("calwebb_spec3", "# options for run_calwebb_spec3: True (will run in full), "
                                "False (run individual steps), skip (only do spec2)", None)
    config.set("calwebb_spec3", "run_calwebb_spec3", spec3_args[0])
    config.set("calwebb_spec3", "s3_input_file", spec3_args[1])

    config.add_section("run_spec2_steps")
    config.set("run_spec2_steps", "# Spec2 steps", None)
    config.set("run_spec2_steps", "assign_wcs", run_pipe_steps[0])
    config.set("run_spec2_steps", "bkg_subtract", run_pipe_steps[1])
    config.set("run_spec2_steps", "imprint_subtract", run_pipe_steps[2])
    config.set("run_spec2_steps", "msa_flagging", run_pipe_steps[3])
    config.set("run_spec2_steps", "extract_2d", run_pipe_steps[4])
    config.set("run_spec2_steps", "srctype", run_pipe_steps[5])
    config.set("run_spec2_steps", "wavecorr", run_pipe_steps[6])
    config.set("run_spec2_steps", "flat_field", run_pipe_steps[7])
    config.set("run_spec2_steps", "pathloss", run_pipe_steps[8])
    config.set("run_spec2_steps", "barshadow", run_pipe_steps[9])
    config.set("run_spec2_steps", "photom", run_pipe_steps[10])
    config.set("run_spec2_steps", "resample_spec", run_pipe_steps[11])
    config.set("run_spec2_steps", "cube_build", run_pipe_steps[12])
    config.set("run_spec2_steps", "extract_1d", run_pipe_steps[13])

    config.add_section("run_spec3_steps")
    config.set("run_spec3_steps", "# Spec3 steps", None)
    config.set("run_spec3_steps", "assign_mtwcs", run_pipe_steps[14])
    config.set("run_spec3_steps", "master_background", run_pipe_steps[15])
    config.set("run_spec3_steps", "exp_to_source", run_pipe_steps[16])
    config.set("run_spec3_steps", "outlier_detection", run_pipe_steps[17])
    config.set("run_spec3_steps", "resample_spec", run_pipe_steps[18])
    config.set("run_spec3_steps", "cube_build", run_pipe_steps[19])
    config.set("run_spec3_steps", "extract_1d", run_pipe_steps[20])

    config.add_section("run_pytest")
    config.set("run_pytest", "# Spec2 tests", None)
    config.set("run_pytest", "assign_wcs_completion_tests", run_pytest[0])
    config.set("run_pytest", "assign_wcs_reffile_tests", run_pytest[1])
    config.set("run_pytest", "assign_wcs_validation_tests", run_pytest[2])
    config.set("run_pytest", "bkg_subtract_completion_tests", run_pytest[3])
    config.set("run_pytest", "bkg_subtract_numerical_tests", run_pytest[4])
    config.set("run_pytest", "imprint_subtract_completion_tests", run_pytest[5])
    config.set("run_pytest", "imprint_subtract_numerical_tests", run_pytest[6])
    config.set("run_pytest", "msa_flagging_completion_tests", run_pytest[7])
    config.set("run_pytest", "msa_flagging_validation_tests", run_pytest[8])
    config.set("run_pytest", "extract_2d_completion_tests", run_pytest[9])
    config.set("run_pytest", "extract_2d_validation_tests", run_pytest[10])
    config.set("run_pytest", "srctype_completion_tests", run_pytest[11])
    config.set("run_pytest", "wavecorr_completion_tests", run_pytest[12])
    config.set("run_pytest", "wavecorr_reffile_tests", run_pytest[13])
    config.set("run_pytest", "flat_field_completion_tests", run_pytest[14])
    config.set("run_pytest", "flat_field_reffile_tests", run_pytest[15])
    config.set("run_pytest", "flat_field_validation_tests", run_pytest[16])
    config.set("run_pytest", "pathloss_completion_tests", run_pytest[17])
    config.set("run_pytest", "pathloss_reffile_tests", run_pytest[18])
    config.set("run_pytest", "pathloss_validation_tests", run_pytest[19])
    config.set("run_pytest", "barshadow_completion_tests", run_pytest[20])
    config.set("run_pytest", "barshadow_validation_tests", run_pytest[21])
    config.set("run_pytest", "photom_completion_tests", run_pytest[22])
    config.set("run_pytest", "resample_spec_completion_tests", run_pytest[23])
    config.set("run_pytest", "cube_build_completion_tests", run_pytest[24])
    config.set("run_pytest", "extract_1d_completion_tests", run_pytest[25])
    config.set("run_pytest", "extract_1d_reffile_tests", run_pytest[26])
    config.set("run_pytest", "# Spec3 tests", None)
    config.set("run_pytest", "master_background_completion_tests", run_pytest[27])
    config.set("run_pytest", "master_background_reffile_tests", run_pytest[28])
    config.set("run_pytest", "master_background_validation_tests", run_pytest[29])

    config.add_section("additional_arguments")
    config.set("additional_arguments", "wcs_threshold_diff", additional_arguments[0])
    config.set("additional_arguments", "save_wcs_plots", additional_arguments[1])
    config.set("additional_arguments", "bkg_list", additional_arguments[2])
    config.set("additional_arguments", "msa_imprint_structure", additional_arguments[3])
    config.set("additional_arguments", "msa_flagging_threshold", additional_arguments[4])
    config.set("additional_arguments", "stellarity", additional_arguments[5])
    config.set("additional_arguments", "save_msa_flagging_plots", additional_arguments[6])
    config.set("additional_arguments", "extract_2d_threshold_diff", additional_arguments[7])
    config.set("additional_arguments", "flattest_threshold_diff", additional_arguments[8])
    config.set("additional_arguments", "save_flattest_plot", additional_arguments[9])
    config.set("additional_arguments", "write_flattest_files", additional_arguments[10])
    config.set("additional_arguments", "pathloss_threshold_diff", additional_arguments[11])
    config.set("additional_arguments", "save_pathloss_plot", additional_arguments[12])
    config.set("additional_arguments", "write_pathloss_files", additional_arguments[13])
    config.set("additional_arguments", "barshadow_threshold_diff", additional_arguments[14])
    config.set("additional_arguments", "save_barshadow_final_plot", additional_arguments[15])
    config.set("additional_arguments", "save_barshadow_intermediary_plots", additional_arguments[16])
    config.set("additional_arguments", "write_barshadow_files", additional_arguments[17])

    detector = fits.getval(calwebb_spec2_input_file[5], "DETECTOR")
    nptt_config = os.path.join(calwebb_spec2_input_file[0], "NPTT_config_"+detector+".cfg")
    config.write(open(nptt_config, "w"))


def set_nptt_immutable_paths():
    wit4_path = os.environ.get('WIT4_PATH')
    if wit4_path is None:
        print("(msa_flagging_testing): The environment variable WIT4_PATH is not defined. To set it, follow the "
              "instructions at: \n"
              "                        https://github.com/spacetelescope/nirspec_pipe_testing_tool")
    test_data_suite = "nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/"
    esa_files_path = os.path.join(wit4_path, test_data_suite)
    return esa_files_path


def list2_allstrings_list(a_list):
    for i, item in enumerate(a_list):
        if not isinstance(item, str):
            a_list[i] = str(item)
    return a_list


def prepare_variables(output_directory, rate_input_file, mode_used, raw_data_root_file, data_directory=None,
                      comparison_file_path=None, msa_conf_name=None, dflat_path=None,
                      sflat_path=None, fflat_path=None, msa_flag_opref=None,
                      run_calwebb_spec2=None, wcs_threshold_diff=None,
                      save_plots=True, change_filter_opaque=False,
                      msa_flagging_threshold=None, stellarity=None, extract_2d_threshold_diff=None,
                      flattest_threshold_diff=None):
    """
    This function prepares all the input variables for the ConfigParser to write the NPTT configuration file.
    Args:
        output_directory: string
        rate_input_file: string
        mode_used: string
        raw_data_root_file: string, basename of the raw data
        data_directory: string
        comparison_file_path: string, path and name of the file to compare with for assign_wcs and extract_2d
        msa_conf_name: string, basename of the shutter configuration file
        dflat_path: string, path and name of the D-flat
        sflat_path: string, path and name of the S-flat for this filter/grating configuration
        fflat_path: string, path and name of the F-flat for this filter configuration
        msa_flag_opref: string, path and name of the msa_flagging operability reference file
        run_calwebb_spec2: string, name of the step to run in spec2
        wcs_threshold_diff: float, acceptable difference between the comparison and the pipeline file
        save_plots: boolean
        change_filter_opaque: boolean
        msa_flagging_threshold: string, acceptable percentage of overlap of values found in index_opens and
                                   index_trace for all slits with more than 100 pixels
        stellarity: string, value for new desired stellarity
        extract_2d_threshold_diff: float, acceptable difference between the comparison and the pipeline file
        flattest_threshold_diff: float, acceptable difference between the comparison and the pipeline file
    Returns:
        variables: list, contains lists of the variables included in each section of the NPTT config file
    """
    if data_directory is None:
        data_directory = output_directory

    local_pipe_cfg_path = 'pipe_source_tree_code'

    if msa_conf_name is None:
        msa_conf_name = '/path_to_corresponding_MSA_shutter_configuration_file/MSA_shutter_config.fits'

    wit4_path = os.environ.get('WIT4_PATH')
    if wit4_path is None:
        print("(mk_npttconfig_file): The environment variable CRDS_PATH is not defined. To set it, follow the "
              "instructions at: \n"
              "                     https://github.com/spacetelescope/nirspec_pipe_testing_tool")
        exit()
    if dflat_path is None:
        dflat_path = os.path.join(wit4_path, 'jwst_nirspec_dflat_0001.fits')
    mode_used = mode_used.upper()
    mu = mode_used
    if 'bots' in mode_used.lower():
        mu = 'FS'
    elif 'mos' in mode_used or 'msa' in mode_used:
        mu = 'MOS'
    elif 'ifu' in mode_used or 'ifs' in mode_used:
        mu = 'IFU'
    if sflat_path is None:
        sflat = os.path.join(wit4_path, 'jwst_nirspec_sflat_0007.fits')
    if fflat_path is None:
        fflat = os.path.join(wit4_path, 'jwst_nirspec_fflat_0004.fits')

    truth_assign_wcs = rate_input_file.replace(".fits", "_assign_wcs_truth.fits")
    truth_extract_2d = rate_input_file.replace(".fits", "_extract_2d_truth.fits")
    if 'rate' in truth_assign_wcs:
        truth_assign_wcs = truth_assign_wcs.replace("_rate", "")
        truth_extract_2d = truth_extract_2d.replace("_rate", "")

    crds_path = os.environ.get('CRDS_PATH')
    if wit4_path is None:
        print("(mk_npttconfig_file): The environment variable CRDS_PATH is not defined. To set it, follow the "
              "instructions at: \n"
              "                     https://github.com/spacetelescope/nirspec_pipe_testing_tool")
        exit()
    if msa_flag_opref is None:
        msa_flag_opref = os.path.join(crds_path, 'references/jwst/jwst_nirspec_msaoper_0001.json')

    #              spec2
    pipe_steps = ['assign_wcs', 'bkg_subtract', 'imprint_subtract', 'msa_flagging', 'extract_2d', 'srctype',
                  'wavecorr', 'flat_field', 'pathloss', 'barshadow', 'photom', 'resample_spec', 'cube_build',
                  'extract_1d',
                  # spec3
                  'assign_mtwcs', 'master_background', 'exp_to_source', 'outlier_detection', 'resample_spec',
                  'cube_build', 'extract_1d']

    nptt_pytests = ['assign_wcs_completion_tests', 'assign_wcs_reffile_tests', 'assign_wcs_validation_tests',
                   'bkg_subtract_completion_tests', 'bkg_subtract_numerical_tests',
                   'imprint_subtract_completion_tests', 'imprint_subtract_numerical_tests',
                   'msa_flagging_completion_tests', 'msa_flagging_validation_tests',
                   'extract_2d_completion_tests', 'extract_2d_validation_tests',
                   'srctype_completion_tests',
                   'wavecorr_completion_tests', 'wavecorr_reffile_tests',
                   'flat_field_completion_tests', 'flat_field_reffile_tests', 'flat_field_validation_tests',
                   'pathloss_completion_tests', 'pathloss_reffile_tests', 'pathloss_validation_tests',
                   'barshadow_completion_tests', 'barshadow_validation_tests',
                   'photom_completion_tests',
                   'resample_spec_completion_tests',
                   'cube_build_completion_tests',
                   'extract_1d_completion_tests', 'extract_1d_reffile_tests',   # spec2 steps up to here
                   # spec3 steps from here on
                   'master_background_completion_tests', 'master_background_reffile_tests',
                   'master_background_validation_tests']

    run_pipe_steps, run_pytests = [], []
    if run_calwebb_spec2 is None:
        run_calwebb_spec2 = True
        print("\n * The NPTT configuration file will be created with running the spec2 pipeline in full and ALL pytests \n"
              "   will be set to True. \n "
              "    -> If you need to skip any tests, please open the NPTT config file and change the pytest values to \n"
              "        False for the tests you are interested in skipping. \n"
              "    -> The spec3 pipeline steps will also be set to True. \n")
        # set individual steps to False and all NPTT tests to True
        for _ in pipe_steps:
            run_pipe_steps.append(False)
        for _ in nptt_pytests:
            run_pytests.append(True)
    else:
        step2run = run_calwebb_spec2
        run_calwebb_spec2 = False
        for step in pipe_steps:
            # set individual steps to True or False
            if step in step2run:
                run_pipe_steps.append(True)
            else:
                run_pipe_steps.append(False)
        # set individual NPTT tests to True or False
        for ptest in nptt_pytests:
            if step2run in ptest:
                run_pytests.append(True)
            else:
                run_pytests.append(False)

    # get the immutable paths
    esa_files_path = set_nptt_immutable_paths()

    # set the full ESA path to compare the data
    if comparison_file_path is None:
        esa_files_full_path = "".join([esa_files_path, mode_used, "_CV3/ESA_Int_products"])
        cutout = False
        if "fs" in mode_used.lower():
            esa_subarray = fits.getval(raw_data_root_file, "SUBARRAY")
            if isinstance(esa_subarray, bool):
                if esa_subarray:
                    cutout = True
            elif 'full' not in esa_subarray.lower():
                cutout = True
            if cutout:
                esa_files_full_path = "".join([esa_files_path, mode_used, "_CV3_cutouts/ESA_Int_products"])
    else:
        esa_files_full_path = comparison_file_path

    # spec3 variables
    run_calwebb_spec3 = 'skip'
    s3_input_file = rate_input_file.lower().split("_nrs")[0] + "_asn.json"

    # set the additional parameters section
    if wcs_threshold_diff is None:
        wcs_threshold_diff = 1.0e-7
    if save_plots:
        print(" * NPTT will save all test output plots in the output directory.")
        save_wcs_plots = True
        save_msa_flagging_plots = True
        save_flattest_plot = True
        save_pathloss_plot = True
        save_barshadow_final_plot = True
        save_barshadow_intermediary_plots = False
    else:
        print(" * NPTT will NOT save any test output plots.")
        save_wcs_plots = False
        save_msa_flagging_plots = False
        save_flattest_plot = False
        save_pathloss_plot = False
        save_barshadow_final_plot = False
        save_barshadow_intermediary_plots = False
    bkg_list = "/path_to_this_file/bkg_example1.fits, /path_to_this_file/bkg_example2.fits"
    msa_imprint_structure = "/path_to_this_file/msa_structure_imprint.fits"
    if msa_flagging_threshold is None:
        msa_flagging_threshold = 99.5
    if stellarity is None:
        stellarity = "source_type"
    if extract_2d_threshold_diff is None:
        extract_2d_threshold_diff = 4
    if flattest_threshold_diff is None:
        flattest_threshold_diff = 9.999e-5
    write_flattest_files = True
    pathloss_threshold_diff = 9.999e-5
    write_pathloss_files = True
    barshadow_threshold_diff = 0.0025
    write_barshadow_files = True

    # set the config file list sections
    calwebb_spec2_input_file = [output_directory, data_directory, rate_input_file, mode_used,
                                change_filter_opaque, raw_data_root_file, local_pipe_cfg_path]
    benchmark_intermediary_products = [esa_files_full_path, msa_conf_name, dflat_path, sflat_path, fflat_path,
                                       truth_assign_wcs, truth_extract_2d, msa_flag_opref]
    run_calwebb_spec2_in_full = str(run_calwebb_spec2)
    spec3_args = [run_calwebb_spec3, s3_input_file]
    additional_arguments = [wcs_threshold_diff, save_wcs_plots, bkg_list, msa_imprint_structure,
                            msa_flagging_threshold, stellarity, save_msa_flagging_plots,
                            extract_2d_threshold_diff, flattest_threshold_diff, save_flattest_plot,
                            write_flattest_files, pathloss_threshold_diff, save_pathloss_plot, write_pathloss_files,
                            barshadow_threshold_diff, save_barshadow_final_plot, save_barshadow_intermediary_plots,
                            write_barshadow_files]

    # make sure all variables are strings for creating the configuration file
    calwebb_spec2_input_file = list2_allstrings_list(calwebb_spec2_input_file)
    benchmark_intermediary_products = list2_allstrings_list(benchmark_intermediary_products)
    run_pipe_steps = list2_allstrings_list(run_pipe_steps)
    run_pytests = list2_allstrings_list(run_pytests)
    spec3_args = list2_allstrings_list(spec3_args)
    additional_arguments = list2_allstrings_list(additional_arguments)

    variables = [calwebb_spec2_input_file, benchmark_intermediary_products, run_calwebb_spec2_in_full, run_pipe_steps,
                 run_pytests, spec3_args, additional_arguments]
    return variables


def mk_nptt_cfg(output_directory, input_file, mode_used, raw_data_root_file, data_directory=None,
               comparison_file_path=None, msa_conf_name=None, dflat_path=None,
               sflat_path=None, fflat_path=None, msa_flag_opref=None, run_calwebb_spec2=None, wcs_threshold_diff=None,
               save_plots=True, change_filter_opaque=False, msa_flagging_threshold=None, stellarity=None,
               extract_2d_threshold_diff=None, flattest_threshold_diff=None, association=None):
    """
    This function makes the NPTT configuration file.
    Args:
        output_directory: string
        input_file: string
        mode_used: string
        raw_data_root_file: string, basename of the raw data
        data_directory: string
        comparison_file_path: string, path and name of the file to compare with for assign_wcs
        msa_conf_name: string, basename of the shutter configuration file
        dflat_path: string, path and name of the D-flat
        sflat_path: string, path and name of the S-flat for this filter/grating configuration
        fflat_path: string, path and name of the F-flat for this filter configuration
        msa_flag_opref: string, path and name of the msa_flagging operability reference file
        run_calwebb_spec2: string, name of the step to run in spec2
        wcs_threshold_diff: string, acceptable difference between the comparison and the pipeline file
        save_plots: boolean
        change_filter_opaque: boolean
        msa_flagging_threshold: string, acceptable percentage of overlap of values found in index_opens and
                                   index_trace for all slits with more than 100 pixels
        stellarity: string, value for new desired stellarity
        extract_2d_threshold_diff: string, acceptable difference between the comparison and the pipeline file
        flattest_threshold_diff: string, acceptable difference between the comparison and the pipeline file
    Returns:
        nothing
    """
    # prepare all variables to create the NPTT config file
    variables = prepare_variables(output_directory, input_file, mode_used, raw_data_root_file,
                                  data_directory=data_directory,
                                  comparison_file_path=comparison_file_path, msa_conf_name=msa_conf_name,
                                  dflat_path=dflat_path, sflat_path=sflat_path, fflat_path=fflat_path,
                                  msa_flag_opref=msa_flag_opref,
                                  run_calwebb_spec2=run_calwebb_spec2, wcs_threshold_diff=wcs_threshold_diff,
                                  save_plots=save_plots, change_filter_opaque=change_filter_opaque,
                                  msa_flagging_threshold=msa_flagging_threshold, stellarity=stellarity,
                                  extract_2d_threshold_diff=extract_2d_threshold_diff,
                                  flattest_threshold_diff=flattest_threshold_diff)

    (calwebb_spec2_input_file, benchmark_intermediary_products, run_calwebb_spec2_in_full,
    run_pipe_steps, run_pytests, spec3_args, additional_arguments) = variables

    # Create the NPTT config file
    write_nptt_cfg(calwebb_spec2_input_file, benchmark_intermediary_products, run_calwebb_spec2_in_full, run_pipe_steps,
                  run_pytests, spec3_args, additional_arguments)
    print('\n * Script  mk_npttconfig_file.py  finished * \n')


def main():

    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    # get the required variables
    parser.add_argument("output_directory",
                        action='store',
                        default=None,
                        help='Required variable. Path to directory where all the NPTT output will be saved.')
    parser.add_argument("input_file",
                        action='store',
                        default=None,
                        help='Required variable. Basename of input count rate fits file, i.e. blah_rate.fits')
    parser.add_argument("mode_used",
                        action='store',
                        default=None,
                        help='Required variable. Observation mode used: FS, MOS, IFU, BOTS, dark, image, confirm, '
                             'taconfirm, wata, msata, focus, mimf - The code is not case-sensitive.')
    parser.add_argument("raw_data_root_file",
                        action='store',
                        default=None,
                        help='Required variable. Basename of raw fits file, e.g. '
                             'NRSV84600010001P0000000002101_4_491_SE_2016-01-17T17h34m08.fits')

    # set the optional variables
    parser.add_argument("-d",
                        dest="data_directory",
                        action='store',
                        default=None,
                        help='Use the -d flag to change the path where NPTT will look for the input file.')
    parser.add_argument("-c",
                        dest="comparison_file_path",
                        action='store',
                        default=None,
                        help='Use the -c flag to change the path and name of the file NPTT will use as comparison '
                             'for assign_wcs test.')
    parser.add_argument("-m",
                        dest="msa_conf_name",
                        action='store',
                        default=None,
                        help='Use the -m flag to change the path and name of the MSA shutter configuration file, '
                             'e.g. blah_msa.fits')
    parser.add_argument("-D",
                        dest="dflat_path",
                        action='store',
                        default=None,
                        help='Use the -D flag to change the path and name of the D-flat file, e.g. /path_to_file/'
                             'jwst_CRDS_file_for_dflat.fits')
    parser.add_argument("-S",
                        dest="sflat_path",
                        action='store',
                        default=None,
                        help='Use the -S flag to change the path and name of the S-flat file, e.g. /path_to_file/'
                             'jwst_CRDS_file_for_sflat.fits')
    parser.add_argument("-F",
                        dest="fflat_path",
                        action='store',
                        default=None,
                        help='Use the -F flag to change the path and name of the F-flat file, '
                             'e.g. /path_to_file/jwst_CRDS_file_for_fflat.fits')
    parser.add_argument("-O",
                        dest="msa_flag_opref",
                        action='store',
                        default=None,
                        help='Use the -O flag to change the path and name of the msa_flagging operability reference'
                             'file , e.g. /path_to_file/jwst_nirspec_msaoper_0001.json')
    parser.add_argument("-r",
                        dest="run_calwebb_spec2",
                        action='store',
                        default=None,
                        help='Use the -r flag to only run the indicated step instead the full spec2 pipeline.')
    parser.add_argument("-w",
                        dest="wcs_threshold_diff",
                        action='store',
                        default=None,
                        help='Use the -w flag to change the float value that NPTT will use as acceptable difference '
                             'between the pipeline product and the comparison file for the assign_wcs test.')
    parser.add_argument("-p",
                        dest="save_plots",
                        action='store_false',
                        default=True,
                        help='Use the -p flag to skip saving the plots produced by the NPTT tests.')
    parser.add_argument("-o",
                        dest="change_filter_opaque",
                        action='store_true',
                        default=False,
                        help='Use the -o flag if the filter needs to be changed from OPAQUE to science or '
                             'vice versa.')
    parser.add_argument("-e",
                        dest="extract_2d_threshold_diff",
                        action='store',
                        default=None,
                        help='Use the -e flag to change the float value that NPTT will use as acceptable difference '
                             'between the pipeline product and the comparison file for the extract_2d test.')
    parser.add_argument("-f",
                        dest="flattest_threshold_diff",
                        action='store',
                        default=None,
                        help='Use the -f flag to change the float value that NPTT will use as acceptable difference '
                             'between the pipeline product and the comparison file for the flat_field test.')
    parser.add_argument("-t",
                        dest="msa_flagging_threshold",
                        action='store',
                        default=None,
                        help='Use the -t flag to change the float value that NPTT will use as acceptable percentage '
                             'of overlap of values found in index_opens and index_trace for all slits with more '
                             'than 100 pixels.')
    parser.add_argument("-s",
                        dest="stellarity",
                        action='store',
                        default=None,
                        help='Use the -s flag if to provide a specific value of stellarity, e.g. -s=0.8')
    args = parser.parse_args()

    # Set the required variables
    output_directory = args.output_directory
    input_file = args.input_file
    mode_used = args.mode_used
    raw_data_root_file = args.raw_data_root_file

    # set the other variables
    data_directory = args.data_directory
    change_filter_opaque = args.change_filter_opaque
    comparison_file_path = args.comparison_file_path
    msa_conf_name = args.msa_conf_name
    dflat_path = args.dflat_path
    sflat_path = args.sflat_path
    fflat_path = args.fflat_path
    msa_flag_opref = args.msa_flag_opref
    run_calwebb_spec2 = args.run_calwebb_spec2
    wcs_threshold_diff = args.wcs_threshold_diff
    save_plots = args.save_plots
    extract_2d_threshold_diff = args.extract_2d_threshold_diff
    flattest_threshold_diff = args.flattest_threshold_diff
    msa_flagging_threshold = args.msa_flagging_threshold
    stellarity = args.stellarity

    # create the NPTT config file
    mk_nptt_cfg(output_directory, input_file, mode_used, raw_data_root_file,
               data_directory=data_directory,
               comparison_file_path=comparison_file_path,  msa_conf_name=msa_conf_name,
               dflat_path=dflat_path, sflat_path=sflat_path, fflat_path=fflat_path, msa_flag_opref=msa_flag_opref,
               run_calwebb_spec2=run_calwebb_spec2, wcs_threshold_diff=wcs_threshold_diff,
               save_plots=save_plots, change_filter_opaque=change_filter_opaque,
               msa_flagging_threshold=msa_flagging_threshold, stellarity=stellarity,
               extract_2d_threshold_diff=extract_2d_threshold_diff,
               flattest_threshold_diff=flattest_threshold_diff)


if __name__ == '__main__':
    sys.exit(main())

