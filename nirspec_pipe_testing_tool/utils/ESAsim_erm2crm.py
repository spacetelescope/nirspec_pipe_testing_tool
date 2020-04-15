#!/usr/bin/env python
# IDENT         pthl_simulations_method.py
# LANGUAGE      Python 3.X
# AUTHOR        E.PUGA, M. PENA-GUERRERO
# PURPOSE       Wrapper Script to post-process IPS simulations
#               from erm --> crm
#               Ingredients: Input simulations, config file
# VERSION
# 1.0.0  12/04/2019 EP Creation

import json
import argparse
import sys
from pprint import pprint
from .ESAsim_post_processing import f_post_processing

def main():
    # =======================================================================
    # Preparing the list of input arguments
    # =======================================================================
    parser = argparse.ArgumentParser()
    # Mandatory input arguments.
    parser.add_argument('configPath', type=str, help='Path of the configuration file.')
    args = parser.parse_args()

    _config_file = args.configPath
    data = json.load(open(_config_file))
    print("# Using configuration file: {:s}".format(_config_file))


    input_inpath, input_outpath = (data['paths'].values())

    input_fwa, input_gwa, mode = (data['spectrograph_configuration'].values())

    detector_settings = data['detector_settings']
    detector_configuration = detector_settings['configuration']
    nint, input_ng, nf = (detector_settings['multiaccum'].values())
    input_nexp = detector_settings['nexp']

    crm_parameters = data['crm_parameters']
    sampling_grid, rep, components, apertures, input_obj_suffix, input_bckg_suffix, input_ffld_suffix, input_factor = (crm_parameters.values())

    archive_structure = data['archive_structure']
    input_base_obs_id, input_datetime_start, input_datetime_end, input_datetime_file, input_daily_folder = (archive_structure.values())

    if 'crm_input_list' in data:
        input_list = data['crm_input_list'] #dictionary with 2-dimensional list and other 1-dimensional lists
    else:
        input_list = None

    pprint(data)

    f_post_processing(input_inpath, input_fwa, input_gwa, input_nexp, input_ng, input_outpath, input_base_obs_id,
                      input_mode=mode,
                      detector_config=detector_configuration, nint=nint, nf=nf, sampling_grid=sampling_grid, rep=rep, components=components, apertures=apertures,
                      input_list=input_list,
                      input_seed=None, input_background=input_bckg_suffix, input_unity=input_ffld_suffix, input_object=input_obj_suffix, input_factor=input_factor,
                      input_datetime_start=input_datetime_start, input_datetime_end=input_datetime_end,
                      input_datetime_file=input_datetime_file, input_daily_folder=input_daily_folder)

if __name__ == '__main__':
    sys.exit(main())
