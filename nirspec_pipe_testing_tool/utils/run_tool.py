"""
This script is a wrapper to run NPTT and then move the report to the right directory.

Example usage:
    The code works from the terminal or called as a module.

    Terminal
        Simply type the command
        $ nptt_run report_name NPTT_config.cfg

        Optional arguments:
        -d1=blah_uncal.fits  this will run calwebb_detector1 on the given file
        -q  this is a boolean, if used several on-screen messages will be avoided
        
    As a module
        # import the script
        import nirspec_pipe_testing_tool as nptt

        # set the variables
        report_name = 'my_report'
        config_file = 'NPTT_config_NRS1.cfg'
        uncal_file = 'blah_uncal.fits'
        quiet = False

        # run the module
        nptt.utils.run.run_nptt(report_name=report_name, config_path=config,
                                uncal_file=uncal_file, verbose=quiet)
"""

import os
import subprocess
import sys
import configparser
import argparse
import time
from astropy.io import fits

import jwst
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline

import nirspec_pipe_testing_tool as nptt
from nirspec_pipe_testing_tool.utils import change_filter_opaque2science
from nirspec_pipe_testing_tool.utils.data import DATADIR
from nirspec_pipe_testing_tool import core_utils
from .. import calwebb_spec2_pytests
from .. import calwebb_spec3_pytests


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.2"

# HISTORY
# Sep 2019 - Version 1.0: initial version completed
# Jun 2022 - Version 1.1: Updated file name output and documentation of usage
# Apr 2023 - Version 1.2: Updated running pipelines code and fixed format


def run_det1(uncal_file, output_dir, intermediary_products=True):
    """
    Run the Detector1 pipeline on the given uncal_file.
    Args:
        uncal_file: string, full path of the uncal_file to run
        output_dir: string, path to place pipeline output
        intermediary_products: boolean, save all intermediary step files or not
    Returns:
        end_time: string, pretty print of how long the pipeline took
    """
    # start the timer
    start_time = time.time()
    # run detector1
    det1 = Detector1Pipeline()
    parameter_dict = {"group_scale": {"save_results": intermediary_products},
                      "dq_init": {"save_results": intermediary_products},
                      "saturation": {"save_results": intermediary_products},
                      "superbias": {"save_results": intermediary_products},
                      "refpix": {"save_results": intermediary_products},
                      "lastframe": {"save_results": intermediary_products},
                      "linearity": {"save_results": intermediary_products},
                      "dark_current": {"save_results": intermediary_products},
                      "jump": {"save_results": intermediary_products},
                      "ramp_fit": {"save_results": intermediary_products},
                      "gain_scale": {"save_results": intermediary_products}
                     }
    det1.call(uncal_file, output_dir=output_dir, steps=parameter_dict,
              logcfg="pipeline-log.cfg", save_results=True)
    # end the timer and report
    end_time, end_time_string = calc_run_time(start_time)
    print(' * The detector1 pipeline finished in ', end_time_string)
    return end_time, end_time_string


def run_spec2(rate_file, output_dir, pipe_stps2run, intermediary_products=True):
    """
    Run the Spec2 pipeline on the given rate_file.
    Args:
        rate_file: string, full path of the file to run
        output_dir: string, path to place pipeline output
        pipe_stps2run: dictionary, steps to be run
        intermediary_products: boolean, save all intermediary step files or not
    Returns:
        end_time: string, pretty print of how long the pipeline took
    """
    # start the timer
    start_time = time.time()
    # run spec2
    spec2 = Spec2Pipeline()
    parameter_dict = {"assign_wcs": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['assign_wcs']},
                      "barshadow": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['barshadow']},
                      "bkg_subtract": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['bkg_subtract']},
                      "extract_2d": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['extract_2d']},
                      "srctype": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['srctype']},
                      "wavecorr": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['wavecorr']},
                      "flat_field": {"save_interpolated_flat": intermediary_products,
                                    "save_results": intermediary_products,
                                    "skip": pipe_stps2run['flat_field']},
                      "msa_flagging": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['msa_flagging']},
                      "pathloss": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['pathloss']},
                      "barshadow": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['barshadow']},
                      "photom": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['photom']},
                      "resample_spec": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['resample_spec']},
                      "cube_build": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['cube_build']},
                      "extract_1d": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['extract_1d']}
                     }
    spec2.call(rate_file, output_dir=output_dir, steps=parameter_dict,
               logcfg="pipeline-log.cfg", save_results=True)
    # end the timer and report
    end_time, end_time_string = calc_run_time(start_time)
    print(' * The spec2 pipeline finished in ', end_time_string)
    return end_time, end_time_string


def run_spec3(asn_json, output_dir, pipe_stps2run, intermediary_products=True):
    """
    Run the Spec3 pipeline on the given rate_file.
    Args:
        asn_json: string, full path of the association file to run
        output_dir: string, path to place pipeline output
        pipe_stps2run: dictionary, steps to be run
        intermediary_products: boolean, save all intermediary step files or not
    Returns:
        end_time: string, pretty print of how long the pipeline took
    """
    # start the timer
    start_time = time.time()
    spec3 = Spec3Pipeline()
    parameter_dict = {"assign_mtwcs": {"save_results": intermediary_products,
                                      "skip": pipe_stps2run['assign_mtwcs']},
                      "master_background": {"save_results": intermediary_products,
                                           "skip": pipe_stps2run['master_background']},
                      "exp_to_source": {"save_results": intermediary_products,
                                       "skip": pipe_stps2run['exp_to_source']},
                      "outlier_detection": {"save_results": intermediary_products,
                                           "skip": pipe_stps2run['outlier_detection']},
                      "resample_spec": {"save_results": intermediary_products,
                                       "skip": pipe_stps2run['resample_spec']},
                      "cube_build": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['cube_build']},
                      "extract_1d": {"save_results": intermediary_products,
                                    "skip": pipe_stps2run['extract_1d']}
                     }
    spec3.call(asn_json, output_dir=output_dir, steps=parameter_dict,
               logcfg="pipeline-log.cfg", save_results=True)
    # end the timer and report
    end_time, end_time_string = calc_run_time(start_time)
    print(' * The spec3 pipeline finished in ', end_time_string)
    return end_time, end_time_string


def calc_run_time(start_time):
    """
    Calculate the total running time.
    Args:
        start_time: string, statting time - UTC in seconds
    Returns:
        end_time: string, end time in seconds and corresponding min or hr conversion
    """
    end_time_sec = time.time() - start_time
    if end_time_sec > 60.0:
        end_time_min = end_time_sec / 60.  # this is in minutes
        if end_time_min > 60.0:
            end_time_hr = end_time_min / 60.  # this is in hours
            end_time = repr(end_time_sec) + "sec = " + repr(round(end_time_hr, 1)) + "hr"
        else:
            end_time = repr(end_time_sec) + "sec = " + repr(round(end_time_min, 1)) + "min"
    else:
        end_time = repr(end_time_sec) + "sec"
    return end_time_sec, end_time


def print_time2file(txt_name, end_time, string2print):
    """
    This function is only used in the case of running the pipeline in full. It
    changes the total running time to the appropriate units and returns a string
    of the right spacing.
    Args:
        txt_name: string, name of the suffix map file
        end_time: string, time it took to run
        string2print: string, phrase to be printed in the file before the time

    Returns:
        nothing
    """
    if float(end_time) > 60.0:
        total_time_min = repr(round(float(end_time)/60.0, 1))   # in minutes
        total_time = total_time_min+'min'
        if float(total_time_min) > 60.0:
            total_time_hr = repr(round(float(end_time)/60.0, 1))   # in hours
            total_time = total_time_hr+'hr'
        end_time = end_time+'  ='+total_time
    line2write = "{:<20} {:<20} {:<20} {:<20}".format('', '', string2print, end_time)
    print(line2write)
    with open(txt_name, "a") as tf:
        tf.write(line2write+"\n")


def read_nptt_config_file(config_path):
    """
    This function reads the NPTT configuration file to get needed info.
    Args:
        config_path: string, path of the NPTT configuration file to use
    Returns:
        cfg_info: list of read info
    """

    if not os.path.exists(config_path):
        raise FileNotFoundError(config_path)

    config = configparser.ConfigParser()
    config.read([config_path])
    output_dir = config.get("calwebb_spec2_input_file", "output_directory")
    data_dir = config.get("calwebb_spec2_input_file", "data_directory")
    input_file = config.get("calwebb_spec2_input_file", "input_file")
    input_file = os.path.join(data_dir, input_file)
    mode_used =  config.get("calwebb_spec2_input_file", "mode_used")
    change_filter_opaque = config.getboolean("calwebb_spec2_input_file", "change_filter_opaque")
    msa_shutter_conf = config.get("benchmark_intermediary_products", "msa_conf_name")
    spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    spec3 = config.get("calwebb_spec3", "run_calwebb_spec3")
    spec3_asn = config.get("calwebb_spec3", "s3_input_file")
    awcs = config.getboolean("run_spec2_steps", "assign_wcs")
    bds = config.getboolean("run_spec2_steps", "bkg_subtract")
    ims = config.getboolean("run_spec2_steps", "imprint_subtract")
    msaf = config.getboolean("run_spec2_steps", "msa_flagging")
    e2d = config.getboolean("run_spec2_steps", "extract_2d")
    sct = config.getboolean("run_spec2_steps", "srctype")
    wco = config.getboolean("run_spec2_steps", "wavecorr")
    ffd = config.getboolean("run_spec2_steps", "flat_field")
    ptl = config.getboolean("run_spec2_steps", "pathloss")
    bsh = config.getboolean("run_spec2_steps", "barshadow")
    pho = config.getboolean("run_spec2_steps", "photom")
    rsp = config.getboolean("run_spec2_steps", "resample_spec")
    cbd = config.getboolean("run_spec2_steps", "cube_build")
    e1d = config.getboolean("run_spec2_steps", "extract_1d")
    spec2_stps2run = {'assign_wcs': awcs, 'bkg_subtract': bds,
                      'imprint_subtract': ims, 'msa_flagging': msaf,
                      'extract_2d': e2d, 'srctype': sct,
                      'wavecorr': wco, 'flat_field': ffd,
                      'pathloss': ptl, 'barshadow': bsh,
                      'photom': pho, 'resample_spec': rsp,
                      'cube_build': cbd, 'extract_1d': e1d}
    wcs = config.getboolean("run_spec3_steps", "assign_mtwcs")
    mbg = config.getboolean("run_spec3_steps", "master_background")
    e2s = config.getboolean("run_spec3_steps", "exp_to_source")
    old = config.getboolean("run_spec3_steps", "outlier_detection")
    rsp = config.getboolean("run_spec3_steps", "resample_spec")
    cbd = config.getboolean("run_spec3_steps", "cube_build")
    e1d = config.getboolean("run_spec3_steps", "extract_1d")
    spec3_stps2run = {'assign_mtwcs': wcs, 'master_background': mbg,
                      'exp_to_source': e2s, 'outlier_detection': old,
                      'resample_spec': rsp, 'cube_build': cbd, 'extract_1d': e1d}
    cfg_info = [output_dir, data_dir, input_file, mode_used, change_filter_opaque,
                msa_shutter_conf, spec2, spec3, spec3_asn, spec2_stps2run,
                spec3_stps2run]
    return cfg_info


def create_completed_steps_txtfile(txt_suffix_map, step_input_file):
    """
    This function creates the completed steps along with the corresponding suffix of
    the output file name into a text file.
    Args:
        txt_suffix_map: string, full path of where the text file will be written into
        step_input_file: string, name of the input file for the pipeline step

    Returns:
        Nothing. A text file will be created in the pytests directory where all steps will be added
    """
    # name of the text file to collect the step name and suffix
    line0 = "# {:<20}".format("Input file: "+step_input_file)
    line1 = "# {:<17} {:<20} {:<20} {:<20}".format("Step", "Added suffix",
            "Step complition", "Time to run [s]")
    print(line1)
    with open(txt_suffix_map, "w+") as tf:
        tf.write(line0+"\n")
        tf.write(line1+"\n")


def run_nptt(report_name=None, config_path=None, uncal_file=None, verbose=False):
    """
    This function runs NPTT and then moves the html report into the working directory specified
    in the NPTT configuration file.
    Args:
        report_name: string, name of the html report
        config_path: string, path of the NPTT configuration file to use
        uncal_file: string, name of the fits file to use as input for the calwebb_detector1 pipeline
        verbose: boolean
    Returns:
        nothing
    """
    # start the timer for NPTT
    nptt_start_time = time.time()
    print('\n *** Running NPTT *** \n')

    if config_path is None:
        # if there is no NPTT config file, request to make one
        print("No NPTT config file was provided. To create a NPTT config file, the following variables are required: \n"
              "     output_directory = full path where NPTT will place all output from the pipeline and tests \n "
              "     input_file = basename of the count rate file \n"
              "     mode_used = FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf \n"
              "     raw_data_root_file = basename of the file that generated the count rate file \n"
              "Keep in mind that this NPTT config file (in the output directory) has multiple default values. Include \n"
              "the following line in your script: \n "
              "nptt.utils.run_tool.run_nptt(report_name NPTT_config_NRS1.cfg)")
        exit()

    # get the html report and the info from the PTT config file
    cfg_info = read_nptt_config_file(config_path)
    (output_dir, data_dir, input_file, mode_used, change_filter_opaque,
     msa_shutter_conf, spec2, spec3, spec3_asn, spec2_stps2run, spec3_stps2run) = cfg_info

    # detemine if spec2 and spec3 are booleans
    if not "skip" in spec2.lower():
        if "t" in spec2.lower():
            spec2 = True
        else:
            spec2 = False
    if not "skip" in spec3.lower():
        if "t" in spec3.lower():
            spec3 = True
        else:
            spec3 = False

    # make sure you are in the data directory to run the pipeline
    cwd = os.getcwd()
    if cwd != output_dir:
        os.chdir(output_dir)

    # run the detector1 pipeline
    detector = None
    if uncal_file is not None:
        print("-> Running the detector1 pipeline - no intermediary files will be saved for this pipeline   \n")
        print("  * NPTT does not have tests for detector1, with the exception of reference file checks. \n")
        if os.path.isfile(uncal_file):
            detector = fits.getval(uncal_file, 'DETECTOR')
            pipelog = "detector1_pipeline_" + detector + ".log"
            core_utils.mk_stpipe_log_cfg(output_dir, pipelog)
            run_det1(uncal_file, output_dir, intermediary_products=False)
            print(' \n Detector1 pipeline finished! \n')
        else:
            print("\n * OH NO! Unable to run the calwebb_detector1 pipeline because the input file does not exist: ",
                          uncal_file)
            exit()

    # get the detector and make sure it is in the name of the output html report
    hdr = fits.getheader(input_file)
    if detector is None:
        detector = hdr['DETECTOR']
    if report_name is None:
        report_name = 'report'
    if 'html' not in report_name:
        report_name = report_name+'.html'
    if detector not in report_name:
        report_name_list = report_name.split(".html")
        report_name = report_name_list[0]+'_'+detector+".html"
        print('\n -> The test pytest results will be recorded in the html report: ', report_name)

    # create the ST pipeline configuration file in the current working directory
    pipelog = "calspec2_pipeline_" + detector + ".log"
    core_utils.mk_stpipe_log_cfg(output_dir, pipelog)
    print('\n calspec2_pipeline log file created! \n')

    # Create the logfile for NPTT, but remove the previous log file
    pipe_ver = "* Using jwst pipeline version: " + jwst.__version__ + " *"
    nptt_ver = "* Using NPTT version: " + nptt.__version__ + " *"
    nptt_log = os.path.join(output_dir, 'NPTT_calspec2_' + detector + '.log')
    nptt_log = core_utils.mk_nptt_log(nptt_log, reset=True)
    nptt_log.info(pipe_ver)
    nptt_log.info(nptt_ver)

    # run tests for spec2
    if not isinstance(spec2, str):
        # check if processing an image
        imaging_mode = False
        if mode_used in ('image', 'confirm', 'taconfirm', 'wata', 'msata', 'bota', 'focus', 'mimf'):
            imaging_mode = True

        # the pipeline is set to run in full
        if spec2:
            print("\n Determining if Spec2 should run or not for this data set")

            # determine if the spec2 pipeline should be run or not - obtained from CRDS repo:
            # https://github.com/spacetelescope/crds/blob/7fc216e73bd81d334b5d98d165be2caa3b6175e1/crds/jwst/jwst_system_crdscfg_b7.5.yaml
            exp_type = hdr['EXP_TYPE']
            allowed_exp_type_values = ['NRS_AUTOFLAT', 'NRS_AUTOWAVE', 'NRS_BRIGHTOBJ',
                                       'NRS_CONFIRM', 'NRS_FIXEDSLIT', 'NRS_FOCUS',
                                       'NRS_IFU', 'NRS_IMAGE', 'NRS_LAMP', 'NRS_MIMF',
                                       'NRS_MSASPEC', 'NRS_MSATA', 'NRS_TACONFIRM',
                                       'NRS_TACQ', 'NRS_TASLIT', 'NRS_WATA']
            if exp_type not in allowed_exp_type_values:
                print("\n OH NO! EXP_TYPE=", exp_type, " is not in exposure types expected for running the stage 2 pipeline: ",
                      allowed_exp_type_values)
                print("\n * Exiting code * \n")
                exit()
            else:
                print("   And the answer is YES.")

            # check if the filter is to be changed
            if change_filter_opaque:
                _, input_file = change_filter_opaque2science.change_filter_opaque(input_file)
                change_filter_opaque_msg = "WARNING: * With FILTER=OPAQUE, the calwebb_spec2 will run up-to" \
                                           "the extract_2d step. Further steps will be skipped. \n"
                print(change_filter_opaque_msg)
                nptt_log.info(change_filter_opaque_msg)

            # get the shutter configuration file for MOS data only
            if 'MSA' in exp_type:
                # check if the configuration shutter file name is in the header of the fits file and if not add it
                msametfl = hdr['MSAMETFL']
                if os.path.basename(msa_shutter_conf) != msametfl:
                    msametfl = os.path.basename(msa_shutter_conf)
                    fits.setval(input_file, "MSAMETFL", 0, value=msametfl)

            # Check if data is IFU that the Image Model keyword is correct
            if mode_used.lower() == "ifu":
                datamdl = hdr["DATAMODL"]
                if datamdl != "IFUImageModel":
                    fits.setval(input_file, "DATAMODL", 0, value="IFUImageModel")
                    print("DATAMODL keyword changed to IFUImageModel.")

            # create the map to record what steps did run (replace previous file)
            txt_name = os.path.join(output_dir, "spec2_full_run_map_" + detector + ".txt")

            # start the timer to compute the step running time of NPTT
            core_utils.start_end_nptt_time(txt_name, start_time=nptt_start_time, end_time=None)

            # now run the pipeline in full
            msg = "Running spec2 in full -> all intermediary files will be saved... "
            print(msg)
            nptt_log.info(msg)
            end_time, end_time_string = run_spec2(input_file, output_dir, spec2_stps2run,
                                                  intermediary_products=True)
            nptt_log.info(end_time_string)

            # add the running time for all steps to the map file from the pipeline log
            step_running_times = core_utils.calculate_step_run_time(pipelog)
            end_time_list = []
            for stp in core_utils.step_string_dict:
                if stp in step_running_times:
                    step_completed = True
                    step_time = step_running_times[stp]["run_time"]
                    out_suffix = core_utils.step_string_dict[stp]["suffix"]
                    core_utils.add_completed_steps(txt_name, stp, out_suffix, step_completed, step_time)
                    end_time_list.append(step_time)

            # print total running time in the text file and move it to the indicated directory
            string2print = "pipeline_total_time"
            if float(end_time) <= sum(end_time_list):
                tot_time = repr(sum(end_time_list))
            else:
                tot_time = end_time
            print_time2file(txt_name, tot_time, string2print)
            nptt_runtimes_msg = "Pipeline and NPTT run times written in file: " + os.path.basename(
                txt_name) + " in working directory. \n"
            print(nptt_runtimes_msg)
            nptt_log.info(nptt_runtimes_msg)

            # move the final reporting text files to the output directory
            core_utils.move_txt_files_2outdir(output_dir, detector)

            # check if processing an image
            if imaging_mode:
                print('\n * Image processing: All intermediary products will be saved, but NO tests for this mode.')
                print('     ->  For now, all pytest will be skipped since there are no image-specific tests yet. \n')
                exit()

        else:
            # create the map to record what steps did run (replace previous file)
            txt_name = os.path.join(output_dir, "spec2_step_run_map_" + detector + ".txt")
            if os.path.isfile(txt_name):
                os.remove(txt_name)
            create_completed_steps_txtfile(txt_name, input_file)
            msg = "Map of steps ran created at: " + txt_name
            print(msg)
            nptt_log.info(msg)

            # start the timer to compute the step running time of NPTT
            core_utils.start_end_nptt_time(txt_name, start_time=nptt_start_time, end_time=None)

            # now run the pipeline step-by-step, if the step is true
            msg = "Running Spec2 step-by-step... "
            print(msg)
            nptt_log.info(msg)
            # individual pipeline steps will be run in the corresponging script

        # run the tests
        print("\n Running Spec2 tests... ")
        report_name = report_name.replace(".html", "_spec2.html")
        args = ['pytest', '-s', '--config_file='+config_path, '--html='+report_name,
                '--self-contained-html', calwebb_spec2_pytests.TESTSDIR]
        if not verbose:
            args.pop(1)
        subprocess.run(args)
    else:
        print("\n NPTT will skip running the Spec2 pipeline and the corresponding tests ")

    # run tests for spec3
    if spec3 != "skip":
        if spec3:
            print("\n Running spec3 - all intermediary files will be saved... ")
            core_utils.mk_stpipe_log_cfg(output_dir, "calspec3_pipeline.log")
            run_spec3(spec3_asn, output_dir, spec3_stps2run, intermediary_products=True)
        print("\n Running spec3 tests... ")
        if "spec2" in report_name:
            report_name = report_name.replace("_spec2", "")
        report_name = report_name.replace(".html", "_spec3.html")
        args = ['pytest', '-s', '--config_file='+config_path, '--html='+report_name,
                '--self-contained-html', calwebb_spec3_pytests.TESTSDIR]
        if not verbose:
            args.pop(1)
        subprocess.run(args)
    else:
        print("\n NPTT will skip running the Spec3 pipeline and the corresponding tests ")

    print('\n * NPTT run finished * \n')


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("report_name",
                        action='store',
                        default=None,
                        help='Name of the html report, e.g. report_NRS2_v2')
    parser.add_argument('config',
                        action='store',
                        default=None,
                        help="Name of NPTT configuration file. To create a NPTT config file has been created, the"
                             "following variables will be needed: "
                             "  output_directory = full path where NPTT will place all output from the pipeline "
                             "and tests; "
                             "  input_file = basename of the count rate file; "
                             "  mode_used = FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, "
                             "focus, mimf; "
                             "  raw_data_root_file = basename of the file that generated the count rate file. "
                             "Keep in mind that this NPTT config file (in the output directory) has multiple "
                             "default values. Create the NPTT config from the terminal with the command: "
                             "$ nptt_mk_npttconfig_file output_directory input_file mode_used raw_data_root_file")
    parser.add_argument("-d1",
                        dest="uncal_file",
                        action='store',
                        default=None,
                        help='Use the -d1 flag to tell NPTT to run the calwebb_detector_1 pipeline as well, '
                             'e.g. -d1=jwdata0010010_11010_0001_NRS1_uncal.fits')
    parser.add_argument("-q",
                        dest="quiet",
                        action='store_false',
                        default=True,
                        help='Use the -q flag to set verbose to False, and not show progress messages on-screen.')

    args = parser.parse_args()
                        
    # Set the variables
    report_name = args.report_name
    config = args.config
    uncal_file = args.uncal_file
    quiet = args.quiet

    # Run NPTT
    run_nptt(report_name=report_name, config_path=config, uncal_file=uncal_file, verbose=quiet)


if __name__ == '__main__':
    sys.exit(main())
