"""
This script creates a sample configuration file needed to run the run_PTT_with_multiprocessing.py script. A
configuration file will be created in the working directory.
"""

import os
import sys
import configparser


def write_ptt_multiprocessing_cfg():
    config = configparser.ConfigParser(allow_no_value=True)

    config.add_section("data_sets_to_run")
    config.set("data_sets_to_run", "", None)

    config.set("data_sets_to_run", "# path shared by all data sets", None)
    config.set("data_sets_to_run", "common_path", "path_to_testing_data/FS_FULL_FRAME/multiprocess_testing \n")

    config.set("data_sets_to_run", "# specific data sets", None)
    config.set("data_sets_to_run", "data_sets", "FS_G140H_opaque,MOS_G140M_opaque \n")
    config.set("data_sets_to_run", "", None)

    config.set("data_sets_to_run", "# cores to use for multiprocessing - to use all set variable as cores2use=all",
               None)
    config.set("data_sets_to_run", "cores2use", "8 \n")

    ptt_m_config = os.path.join(os.getcwd(), "multiprocessing_PTT_config.cfg")
    config.write(open(ptt_m_config, "w"))


def main():
    # create the sample multiprocessing configuration file in the current working directory
    write_ptt_multiprocessing_cfg()


if __name__ == '__main__':
    sys.exit(main())

