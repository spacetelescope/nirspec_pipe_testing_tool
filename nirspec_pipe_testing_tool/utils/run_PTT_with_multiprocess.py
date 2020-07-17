"""
This script runs NPTT fro multiple data sets with multiprocessing.

Example usage:
    The code works from the terminal or called as a module.

    Terminal
        Simply type the command
        $ nptt_run_PTT_with_multiprocess multiprocessing_PTT_config.cfg

    As a module
        # import the script
        import nirspec_pipe_testing_tool as nptt
        nptt.utils.run_PTT_with_multiprocess.run_PTT_with_multiprocess(multiprocessing_PTT_config.cfg)
"""

import os
import multiprocessing
import sys
import configparser
import argparse
import time
from glob import glob

from . import run_PTT

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"


# HISTORY
# Jul 2020 - Version 1.0: initial version completed


def read_multiprocess_config_file(config_path):
    """
    This function reads the multiprocess PTT configuration file to get needed info.
    :param config_path: str, full path and name of the NPTT multiprocessing configuration file
    :return: cfg_info: list, read information
    """

    if not os.path.exists(config_path):
        raise FileNotFoundError(config_path)

    config = configparser.ConfigParser()
    config.read([config_path])
    common_path = config.get("data_sets_to_run", "common_path")
    data_sets = config.get("data_sets_to_run", "data_sets")
    cal_det1 = config.get("data_sets_to_run", "cal_det1")
    if cal_det1 == "all" or cal_det1 == "skip":
        cal_det1_input_file_list = cal_det1
    else:
        cal_det1_input_file_list = cal_det1.split(",")
    cores2use = config.get("data_sets_to_run", "cores2use")
    cfg_info = [common_path, data_sets, cal_det1_input_file_list, cores2use]
    return cfg_info


def get_cfg_files2run_NPTT(common_path, data_sets):
    """
    Get all the configuration files that will be used as input to run NPTT
    :param common_path: str, common part of the path
    :param data_sets: str, all the data sets separated by a comma
    :return: cfg_files2run_NPTT: list, configuration files to run NPTT
    """
    cfg_files2run_NPTT = []
    data_sets = data_sets.split(",")
    for ds in data_sets:
        data_dir = os.path.join(common_path, ds)
        cfg_files_in_data_dir = glob(data_dir+"/*PTT*.cfg")
        for cfg_file in cfg_files_in_data_dir:
            cfg_files2run_NPTT.append(cfg_file)
    return cfg_files2run_NPTT


def run_PTT_with_multiprocess(multiprocessing_PTT_config_file):
    """
    This function runs PTT and then moves the html report into the working directory specified
    in the PTT configuration file.
    :param multiprocessing_PTT_config_file: string, name and path of the multiprocessing configuration file to use
    """
    print('Running NPTT with multiprocessing. This may take a while...')

    # start timer
    start_time = time.time()

    # get the information from the multiprocessing configuration file
    cfg_info = read_multiprocess_config_file(multiprocessing_PTT_config_file)
    common_path, data_sets, cal_det1_input_file_list, cores2use = cfg_info

    # get all the NPTT configuration files to use - i.e. the times NPTT needs to be run
    cfg_files2run_NPTT = get_cfg_files2run_NPTT(common_path, data_sets)

    # get all the report names to use
    report_names = []
    for ptt_cfg in cfg_files2run_NPTT:
        d = os.path.dirname(ptt_cfg)
        r = os.path.join(d, "report")
        report_names.append(r)

    # gather the list of input files for the stage 1 pipeline if option set to all
    if cal_det1_input_file_list == "all" or cal_det1_input_file_list == "skip":
        det1_input_file_list = []
    else:
        det1_input_file_list = cal_det1_input_file_list
    for ptt_cfg in cfg_files2run_NPTT:
        if cal_det1_input_file_list == "all":
            if "NRS1" in ptt_cfg:
                det = "NRS1"
            if "NRS2" in ptt_cfg:
                det = "NRS2"
            det1_input_file_list.append("jwdata0010010_11010_0001_"+det+"_uncal.fits")
        if cal_det1_input_file_list == "skip":
            det1_input_file_list.append(None)

    # set the cores to use for NPTT
    if cores2use == "all":
        cores2use = os.cpu_count()
    else:
        cores2use = int(cores2use)

    # set the pool of cores to use for NPTT
    p = multiprocessing.Pool(cores2use)
    p.starmap(run_PTT.run_PTT, zip(report_names, cfg_files2run_NPTT, det1_input_file_list))
    p.close()
    p.join()

    # end timer
    end_time = time.time() - start_time

    # pretty print the total time
    if end_time > 60.0:
        total_time = round(end_time / 60.0, 1)  # in minutes
        total_time_str = repr(total_time) + " min"
        if total_time > 60.0:
            total_time = round(total_time / 60.0, 1)  # in hrs
            total_time_str = repr(total_time) + " hr"
    else:
        total_time_str = repr(round(end_time, 1)) + " sec"

    print(f'\n * It took {total_time_str} to process {len(cfg_files2run_NPTT)} NPTT runs with {cores2use} cores * \n')


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('config',
                        action='store',
                        default=None,
                        help="Path and name of multiprocessing configuration file to use.")

    args = parser.parse_args()

    # Set the variables
    config_path = args.config

    # Run NPTT with multiprocessing
    run_PTT_with_multiprocess(config_path)


if __name__ == '__main__':
    sys.exit(main())
