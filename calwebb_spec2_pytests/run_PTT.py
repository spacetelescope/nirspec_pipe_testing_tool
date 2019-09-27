"""
This script is a wrapper to run PTT and then move the report to the right directory. The script
is to be run from the calwebb_spec2_pytests directory.

Example usage:
    The code works from the terminal or called as a module.

    Terminal
        Simply type the command
        > python run_PTT report_name
        
    As a module
        import run_PTT
        run_PTT.run_PTT(report_name)
"""

import os
import re
import subprocess
import configparser
import argparse
from astropy.io import fits


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Sep 2019 - Version 1.0: initial version completed


def read_PTTconfig_file():
    """
    This function reads the PTT configuration file to get needed info.
    
    Returns:
        cfg_info = list of read info
    """
    config = configparser.ConfigParser()
    config.read(['../calwebb_spec2_pytests/PTT_config.cfg'])
    pipe_testing_tool_path = config.get("calwebb_spec2_input_file", "pipe_testing_tool_path")
    working_dir = config.get("calwebb_spec2_input_file", "working_directory")
    input_file = config.get("calwebb_spec2_input_file", "input_file")
    input_file = os.path.join(working_dir, input_file)
    cfg_info = [working_dir, input_file]
    return cfg_info


def run_PTT(report_name):
    """
    This function runs PTT and then moves the html report into the working directory specified
    in the PTT configuration file.
    Args:
    report_name: string, name of the html report
    """
    print('Running PTT. This may take a while...')
    
    # get the html report and the info from the PTT config file
    cfg_info = read_PTTconfig_file()
    
    # get the detector and make sure it is in the name of the output html report
    detector = fits.getval(cfg_info[1], "DETECTOR", 0)
    if detector not in report_name:
        print('-> The detector used is being added to the html report name.')
        report_name = report_name+'_'+detector
    if 'html' not in report_name:
        report_name = report_name+'.html'

    # run PTT
    cmd = ['pytest', '-s', '--config_file=PTT_config.cfg', '--html='+report_name,
           '--self-contained-html']
    subprocess.call(cmd)

    # move the html report
    if os.path.isfile(report_name):
        print('Moving PTT html report to working directory')
        os.rename(report_name, os.path.join(cfg_info[0], report_name))
        print('Done.')
    else:
        print('WARNING: The html report was not created, something went wrong!')



if __name__ == '__main__':
    
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("report_name",
                        action='store',
                        default=None,
                        help='Name of the html report, e.g. report_NRS2_v2')
    args = parser.parse_args()
                        
    # Set the variables
    report_name = args.report_name
        
    # Perform data move to the science extension and the keyword check on the file with the right number of extensions
    run_PTT(report_name)


    print ('\n * Script  run_PTT.py  finished * \n')


