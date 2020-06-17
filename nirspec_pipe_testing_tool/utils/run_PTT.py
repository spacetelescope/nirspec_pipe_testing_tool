"""
This script is a wrapper to run PTT and then move the report to the right directory. The script
is to be run from the calwebb_spec2_pytests directory.

Example usage:
    The code works from the terminal or called as a module.

    Terminal
        Simply type the command
        $ nptt_run_PTT report_name PTT_config.cfg
        
    As a module
        # import the script
        import nirspec_pipe_testing_tool as nptt
        nptt.utils.run_PTT.run_PTT(report_name, PTT_config.cfg)
"""

import os
import subprocess
import sys
import configparser
import argparse
from astropy.io import fits
from .. import calwebb_spec2_pytests
from .. import calwebb_spec3_pytests


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Sep 2019 - Version 1.0: initial version completed


def read_PTTconfig_file(config_path):
    """
    This function reads the PTT configuration file to get needed info.
    
    Returns:
        cfg_info = list of read info
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


def run_PTT(report_name, config_path=None):
    """
    This function runs PTT and then moves the html report into the working directory specified
    in the PTT configuration file.
    Args:
    report_name: string, name of the html report
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

    # get the detector and make sure it is in the name of the output html report
    detector = fits.getval(input_file, "DETECTOR", 0)
    if 'html' not in report_name:
        report_name = report_name+'.html'
    if detector not in report_name:
        report_name_list = report_name.split(".html")
        report_name = report_name_list[0]+'_'+detector+".html"
        print('-> The detector used added to the html report name: ', report_name)

    # run tests for spec2
    if not skip_spec2:
        print("Running spec2 and tests")
        report_name = report_name.replace(".html", "_spec2.html")
        args = ['pytest', '-s', '--config_file='+config_path, '--html='+report_name,
                '--self-contained-html', calwebb_spec2_pytests.TESTSDIR]
        subprocess.run(args)

    # run tests for spec3
    if not skip_spec3:
        print("Running spec3 and tests")
        if "spec2" in report_name:
            report_name = report_name.replace("_spec2", "")
        report_name = report_name.replace(".html", "_spec3.html")
        args = ['pytest', '-s', '--config_file='+config_path, '--html='+report_name,
                '--self-contained-html', calwebb_spec3_pytests.TESTSDIR]
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

    args = parser.parse_args()
                        
    # Set the variables
    report_name = args.report_name
    config_path = args.config

    # Perform data move to the science extension and the keyword check on the file with the right number of extensions
    run_PTT(report_name, config_path)


if __name__ == '__main__':
    sys.exit(main())
