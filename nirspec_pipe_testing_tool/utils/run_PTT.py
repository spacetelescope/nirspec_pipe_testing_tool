"""
This script is a wrapper to run PTT and then move the report to the right directory.

Example usage:
    The code works from the terminal or called as a module.

    Terminal
        Simply type the command
        $ nptt_run_PTT report_name PTT_config.cfg
        
    As a module
        # import the script
        import nirspec_pipe_testing_tool as nptt

        # set the variables
        report_name = 'my_report'
        config_file = 'PTT_config_NRS1.cfg'
        quiet = False

        # run the module
        nptt.utils.run_PTT.run_PTT(report_name, config_file, quiet)
"""

import os
import subprocess
import sys
import configparser
import argparse
from astropy.io import fits

from . import run_cal_detector1
from .. import calwebb_spec2_pytests
from .. import calwebb_spec3_pytests


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Sep 2019 - Version 1.0: initial version completed


def mk_stpipe_log_cfg(detector):
    """
    Create a configuration file with the name pipeline_detector.log
    :param detector: str, name of the detector to use
    :return: nothing
    """
    config = configparser.ConfigParser()
    config.add_section("*")
    config.set("*", "handler", "file:pipeline_"+detector+".log")
    config.set("*", "level", "INFO")
    pipe_log_config = os.path.join(os.getcwd(), "stpipe-log.cfg")
    config.write(open(pipe_log_config, "w"))


def read_PTTconfig_file(config_path):
    """
    This function reads the PTT configuration file to get needed info.
    :param config_path: string, path of the PTT configuration file to use
    :return: cfg_info: list of read info
    """

    if not os.path.exists(config_path):
        raise FileNotFoundError(config_path)

    config = configparser.ConfigParser()
    config.read([config_path])
    output_dir = config.get("calwebb_spec2_input_file", "output_directory")
    input_file = config.get("calwebb_spec2_input_file", "input_file")
    input_file = os.path.join(output_dir, input_file)
    spec2 = config.get("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    spec3 = config.get("calwebb_spec3", "run_calwebb_spec3")
    cfg_info = [output_dir, input_file, spec2, spec3]
    return cfg_info


def run_PTT(report_name=None, config_path=None, run_detector1=None, verbose=False):
    """
    This function runs PTT and then moves the html report into the working directory specified
    in the PTT configuration file.
    :param report_name: string, name of the html report
    :param config_path: string, path of the PTT configuration file to use
    :param run_detector1: string, name of the fits file to use as input for the calwebb_detector1 pipeline
    :param verbose: boolean
    """
    print('Running PTT. This may take a while...')

    if config_path is None:
        # if there is no PTT config file, request to make one
        print("No PTT config file was provided. To create a PTT config file, the following variables are required: \n"
              "     output_directory = full path where PTT will place all output from the pipeline and tests \n "
              "     input_file = basename of the count rate file \n"
              "     mode_used = FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf \n"
              "     raw_data_root_file = basename of the file that generated the count rate file \n"
              "Keep in mind that this PTT config file (in the output directory) has multiple default values. Include \n"
              "the following line in your script: \n "
              "nptt.utils.run_PTT.run_PTT(report_name PTT_config.cfg)")
        exit()

    # get the html report and the info from the PTT config file
    cfg_info = read_PTTconfig_file(config_path)
    output_dir, input_file, spec2, spec3 = cfg_info
    skip_spec2 = True
    if spec2 != "skip":
        skip_spec2 = False
    skip_spec3 = True
    if spec3 != "skip":
        skip_spec3 = False

    # make sure you are in the data directory to run the pipeline
    cwd = os.getcwd()
    if cwd != output_dir:
        os.chdir(output_dir)

    # run the detector_1 pipeline
    if run_detector1 is not None:
        print("Running the detector1 pipeline ")
        if os.path.isfile(run_detector1):
            run_cal_detector1.run_caldet1(run_detector1, step_by_step=False)
        else:
            print("Unable to run the calwebb_detector1 pipeline because the input file does not exist: ",
                  run_detector1)

    # get the detector and make sure it is in the name of the output html report
    detector = fits.getval(input_file, "DETECTOR", 0)
    if report_name is None:
        report_name = 'report'
    if 'html' not in report_name:
        report_name = report_name+'.html'
    if detector not in report_name:
        report_name_list = report_name.split(".html")
        report_name = report_name_list[0]+'_'+detector+".html"
        print('-> The detector used added to the html report name: ', report_name)

    # create the ST pipeline configuration file in the current working directory
    mk_stpipe_log_cfg(detector)

    # run tests for spec2
    if not skip_spec2:
        print("Running spec2 and tests")
        report_name = report_name.replace(".html", "_spec2.html")
        args = ['pytest', '-s', '--config_file='+config_path, '--html='+report_name,
                '--self-contained-html', calwebb_spec2_pytests.TESTSDIR]
        if not verbose:
            args.pop(1)
        subprocess.run(args)

    # run tests for spec3
    if not skip_spec3:
        print("Running spec3 and tests")
        if "spec2" in report_name:
            report_name = report_name.replace("_spec2", "")
        report_name = report_name.replace(".html", "_spec3.html")
        args = ['pytest', '-s', '--config_file='+config_path, '--html='+report_name,
                '--self-contained-html', calwebb_spec3_pytests.TESTSDIR]
        if not verbose:
            args.pop(1)
        subprocess.run(args)

    # move the html report
    if os.path.isfile(report_name):
        print('Moving PTT html report to working directory')
        os.rename(report_name, os.path.join(output_dir, report_name))
        print('Done.')
    else:
        print('WARNING: The html report was not created, something went wrong!')

    print('\n * Script  run_PTT.py  finished * \n')


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
                        help="Name of PTT configuration file. To create a PTT config file has been created, the"
                             "following variables will be needed: "
                             "  output_directory = full path where PTT will place all output from the pipeline "
                             "and tests; "
                             "  input_file = basename of the count rate file; "
                             "  mode_used = FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, "
                             "focus, mimf; "
                             "  raw_data_root_file = basename of the file that generated the count rate file. "
                             "Keep in mind that this PTT config file (in the output directory) has multiple "
                             "default values. Create the PTT config from the terminal with the command: "
                             "$ nptt_mk_pttconfig_file output_directory input_file mode_used raw_data_root_file")
    parser.add_argument("-d1",
                        dest="run_detector1",
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
    run_detector1 = args.run_detector1
    quiet = args.quiet

    # Run NPTT
    run_PTT(report_name=report_name, config_path=config, run_detector1=run_detector1, verbose=quiet)


if __name__ == '__main__':
    sys.exit(main())
