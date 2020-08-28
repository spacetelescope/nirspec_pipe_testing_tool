#!/usr/bin/env python
# IDENT         pthl_simulations_method.py
# LANGUAGE      Python 3.X
# AUTHOR        E.PUGA
# PURPOSE       Script with function f_post_processing (top-down structure).
#               This function can be scripted or directly
#               called from command line.
#               Note: This implementation makes use of
#               standard synthetic FPA106 model from
#               nrspydet.misc.fpa106_toolbox
# VERSION
# 1.0.0  12/04/2019 EP Creation

# ===============================================================
# Imports
# ===============================================================
import argparse
import datetime
import os.path
import sys

from astropy.io import fits
import numpy as np
from pprint import pprint
from collections import namedtuple

try:
    import nrspydet
except ImportError:
    print("Missing optional package: 'nrspydet':\n"
          "    This script requires the NIRSpec Instrument Performance Simulator"
          " (IPS) environment to run.", file=sys.stderr)
    exit(1)

from nrspydet.model import configurationfpa as c_configurationfpa
from nrspydet.model import multiaccum as c_multiaccum
from nrspydet.misc import fpa106_toolbox as fpa106_toolbox
from nrspysim.products import electronratemap as c_electronratemap
from nips.archives import toolbox as archive_toolbox
from nrspysim.misc import erm_to_sci as erm_to_sci
from nrspylib.geometry import tiltpoly
import subarray_dict as subdict


# ===============================================================
# Module-wide variables
# ===============================================================
_name = "ESAsim_post_processing.py"
_version = "1.0.0"

#These are parameters that are too big to be a default
# Filenames of real darks used to set the detectors quality flags
_base_path = '/Users/epuga/JWST_NIRSpec/NIRSpec_Calibration/IPS/IPSScenesPreparation/DataPackage/data/'
_dark491fname = _base_path+'darks/NRSDET-DARK-IRS2-7245080538_3_491_SE_2017-09-02T10h49m27.cts.fits'
_dark492fname = _base_path+'darks/NRSDET-DARK-IRS2-7245080538_3_492_SE_2017-09-02T10h49m27.cts.fits'

#TODO: make this function parameters and set these values as defaults
#These are only needed to set the GWA keywords - should no longer be needed in later version of IPS (IPS> 4.0)
#Hard coding model path!!! For the simulation we always use CV3!!!!
_model_path = '/Users/epuga/Library/Software/JWST_Python/data/IQLAC/'
_model_name = 'NIRS_FM2_05_CV3_FIT1'

def f_post_processing(input_inpath, input_fwa, input_gwa, input_nexp, input_ng, input_outpath, input_base_obs_id,
                      input_mode,
                      detector_config={"readout_mode": "IRS2","exposure_type": "full-frame", "subarray_id": None}, nint=1, nf=1, sampling_grid=None, rep=0, components=['object'],
                      apertures={"FS": ["slit"], "MOS": ["mos", "slit"], "IFS": ["ifs", "mos", "slit"]},
                      input_list=None,
                      input_seed=None, input_background=None, input_unity=None, input_object=None, input_factor=1.,
                      input_datetime_start='20180303T120000.000', input_datetime_end='20180303T120000.000',
                      input_datetime_file='2018-03-03T12h15m00', input_daily_folder='Day2018062'):

    '''Main function to be called from
        - a wrapper script that uses a configuration file
        - command line specifying all parameters (inherited from original PF scripts).
        The original parameters coming from this approach have the prefix input_,
        I kept the mandatory parameters as few as possible and adopted them from PF
        scripts (except "mode"), the ones I added new are optional + default value.

        Note: observing mode is handled only at this level, rest of parameters are handled in lower function levels'''

    input_path = input_inpath
    print("# Input path: {:s}".format(input_path))
    output_path = input_outpath
    print("# Output path: {:s}".format(output_path))
    # ===============================================================
    # Printout: Spectrograph config
    # ===============================================================
    mode = input_mode
    print("# Observing mode: {:s}".format(mode))
    print("# =====================================")
    FWA = input_fwa
    GWA = input_gwa
    print("# Spectral configuration: {:s}/{:s}".format(FWA, GWA))
    if (GWA == 'PRISM'):
        filetype = 'prism'
    elif (GWA == 'MIRROR'):
        filetype = 'image'
    else:
        filetype = 'grating'
    print("# Type of input file: {:s}".format(filetype))
    print("# =====================================")

    # ===============================================================
    # Printout: Exposure parameters
    # ===============================================================
    print("# Detector configuration:")
    nexp = input_nexp
    print("# Number of exposures per scene: {:d}".format(nexp))
    ng = input_ng
    nint = nint
    nf = nf
    print(" Number of groups {:d} - number of integrations {:d} - number of frames {:d}".format(ng, nint, nf))
    print("# =====================================")

    # ===============================================================
    # Printout: Output archive structure
    # ===============================================================
    print("# Output archive structure:")
    base_nid = 1000
    print("# Base NID number: {:d}".format(base_nid))
    env = 'IPS'
    pipeline_id = archive_toolbox._dic_env[env][0]
    jlab_id = archive_toolbox._dic_env[env][1]
    print("# Environment: {:s} ({:s} , {:s})".format(env, pipeline_id, jlab_id))
    base_obs_id = input_base_obs_id
    print("# Base observation ID: {:s}".format(base_obs_id))
    datetime_start = input_datetime_start
    datetime_end = input_datetime_end
    print("# Datetime start and end: {} : {}".format(datetime_start, datetime_end))
    datetime_file = input_datetime_file
    print("# Datetime file: {}".format(datetime_file))
    daily_folder = input_daily_folder
    print("# Daily folder: {:s}".format(daily_folder))
    print("# =====================================")

    # ===============================================================
    # Printout: CRM building blocks
    # ===============================================================
    print("# Optional arguments:")
    background_suffix = input_background
    if input_background is not None:
        print("# Suffix for the background electron-rate map: {:s}".format(background_suffix))
    unity_suffix = input_unity
    if unity_suffix is not None:
        print("# Suffix for the unity electron-rate map: {:s}".format(unity_suffix))
    object_suffix = input_object
    if object_suffix is not None:
        print("# Suffix for the object electron-rate map: {:s}".format(object_suffix))
    seed = input_seed
    if (seed is None):
        print("# No random generator seed provided.")
    else:
        print("# Random generator seed: {:d}".format(seed))
        np.random.seed(seed=seed)
    factor = input_factor
    print("# Normalisation factor: {:8.4e}".format(factor))
    print("# =====================================")

    # =======================================================================
    # Generating FPA configuration and multiaccum parameters
    # =======================================================================
    readout_mode = detector_config['readout_mode']
    exposure_type = detector_config['exposure_type']
    subarray_id = detector_config['subarray_id']
    multiaccum = f_detector(readout_mode, exposure_type, ng, nint=nint, nf=nf, GWA=GWA, subarray_id=subarray_id)

    # =======================================================================
    # Generating the instances of the reference files (nrspydet.references)
    # Note: when using real darks, one needs less calls
    # =======================================================================
    references = f_references(readout_mode, exposure_type, seed=seed)

    #TODO: create the logic for sampling grid

    #If we have not provided an explicit erm filename dictionary
    if input_list is None:
        #We need to generate the erm filename list with the logic
        input_list = f_erm_filename(mode, FWA, GWA, sampling_grid, object_suffix, background_suffix=background_suffix, unity_suffix=input_unity, erm_rep=rep, apertures=None)

    #If necessary, massage the erm filename dictionary to make it compatible
    input_list = f_modify_dict(input_list, mode, apertures, sampling_grid['dither_pts'], verbose=True)

    # =======================================================================
    # CRM:
    # - Loops in dither and nod positions
    # - Creates a background if it is provided as one of components, adds it
    # to the object, and writes it out if specifically requested
    # - Creates the object (plus background) crm exposures
    # - Dark flag and GWA keyword population
    # - Writes out the object crm
    # =======================================================================
    sel_apertures = apertures[mode]
    f_crm(input_path, input_list, FWA, GWA, sampling_grid, components, sel_apertures, nexp, multiaccum, references,
          output_path, base_nid, pipeline_id, jlab_id, base_obs_id, datetime_start, datetime_end, datetime_file, daily_folder,
          save_background=False)

    return None

def f_detector(readout_mode, exposure_type, ng, nint=1, nf=5, GWA=None, subarray_id=None):
    '''
        Basic function to parse the fpa configuration and multiaccum settings
        :param readout_mode:
        :param exposure_type:
        :param ng: number of groups
        :param nint: number of integrations. default for comissioning nint=1
        :param nf: number of frames. default is IRS2 nf=5
        :param GWA:
        :param subarray_id:
        :return:
        '''
    if subarray_id is not None:
        if GWA is None:
            print("ERROR - In subarray mode, the GWA is necessary.")
            print("ERROR - {:s} ".format(subarray_id))
            print("ERROR - {:s} ".format(GWA))
    configfpa = c_configurationfpa.ConfigurationFPA()
    configfpa.m_set_mode(readout_mode, exposure_type)

    if subarray_id is None or subarray_id == 'FULL-FRAME':
        configfpa.m_set_array_parameters()
    else:
        start_i = (subdict.subarray_dict[subarray_id]["substrt1"],)*2
        start_j = subdict.subarray_dict[subarray_id]["substrt2"][GWA]
        size_i = subdict.subarray_dict[subarray_id]["subsize1"]
        size_j = subdict.subarray_dict[subarray_id]["subsize2"]
        configfpa.m_set_array_parameters(start_i=start_i, start_j=start_j, size_i=size_i, size_j=size_j)
        configfpa.m_set_header('none', 'none', 'none')

    multiaccum = c_multiaccum.Multiaccum()
    multiaccum.m_set_configuration(configfpa)
    multiaccum.m_set_header('none', 'none', 'none')
    multiaccum.m_set(nint, ng, nf)
    multiaccum.m_info()

    return multiaccum

def f_references(readout_mode, exposure_type, seed=-1):
    '''
        Basic function to generate the necessary reference files which
        include the noise model

        Note: uses namedtuple package

        :param readout_mode:
        :param exposure_type:
        :param seed:
        :return:
        '''
    #
    if readout_mode == 'traditional': readout_mode='TRAD'

    dark = fpa106_toolbox.f_generate_dark(uniform=False, seed=seed)
    gain = fpa106_toolbox.f_generate_gain(exposure_type.upper())
    ctm = fpa106_toolbox.f_generate_ctm()
    readout = fpa106_toolbox.f_generate_readout_noise(readout_mode)
    noise_model = fpa106_toolbox.f_generate_noise_model(readout_mode)
    references = namedtuple("references", ['dark', 'gain', 'ctm', 'readout', 'noise_model'])

    return references(dark, gain, ctm, readout, noise_model)

def f_erm_filename(mode, FWA, GWA, sampling_grid, object_suffix, background_suffix=None, unity_suffix=None, erm_rep=0, apertures=None):
    '''
        Not tested.
        Utility function for components other than object, we are broadcasting in number of dither_pts
        We are not using components, assuming that object_suffix wont be empty
        :param mode:
        :param FWA:
        :param GWA:
        :param sampling_grid:
        :param object_suffix:
        :param background_suffix:
        :param unity_suffix:
        :param int erm_rep: IPS simulation repetition, defaulted to 0
        :param dict apertures:
        :return: erm_dict
        '''
    # =======================================================================
    # Generating the names of the electron-rate maps in a dictionary of keys
    # determined by the number of components.
    # The apertures depend on the obs mode:
    #   - object[dither_idx][nod_idx]
    #   - background[aper][dither_idx] (optional)
    #   - flatfield[aper][dither_idx] (unity) (optional)
    # TODO: they should be 2-dimensional lists
    # =======================================================================
    print("# Generating the names of the electron-rate maps.")
    erm_dict = {}
    list_object = []
    dither_pts = sampling_grid['dither_pts']
    for dither_index in range(dither_pts):
        if not object_suffix:
            name = '{:s}_{:s}_{:s}_{:02s}_{:03d}.erm'.format(FWA, GWA, mode, dither_index, erm_rep)
        else:
            name = '{:s}_{:s}_{:s}_{:s}_{:02}_{:03d}.erm'.format(FWA, GWA, mode, object_suffix, dither_index, erm_rep)
        list_object.append(name)
    print("# Object ERMs: {}".format(list_object))
    erm_dict['target_erm'] = list_object

    #TODO: the loop on dither_idx seems to be only necessary in IFS, handle this in IFS mode
    if background_suffix is not None:
        list_background = []
        for dither_index in range(dither_pts):
            list_names = []
            if apertures is None:
                name = '{:s}_{:s}_{:s}_{:s}.erm'.format(FWA, GWA, mode, background_suffix)
                list_names.append(name)
            else:
                for aperture in apertures[mode]:
                    name = '{:s}_{:s}_{:s}_{:s}.erm'.format(FWA, GWA, aperture, background_suffix)
                    list_names.append(name)
            list_background.append(list_names)
        print("# Background ERMs: {}".format(list_background))
        erm_dict['bckg_erm'] = list_background

    if unity_suffix is not None:
        list_unity = []
        for dither_index in range(dither_pts):
            list_names = []
            if apertures is None:
                name = '{:s}_{:s}_{:s}_{:s}.erm'.format(FWA, GWA, mode, unity_suffix)
                list_names.append(name)
            else:
                for aperture in apertures[mode]:
                    name = '{:s}_{:s}_{:s}_{:s}.erm'.format(FWA, GWA, aperture, unity_suffix)
                    list_names.append(name)
            list_unity.append(list_names)
        print("# Background ERMs: {}".format(list_unity))
        erm_dict['ff_erm'] = list_unity

    return erm_dict

def  f_output_archive_structure(output_path, dither_idx, nod_idx,
                                base_nid, pipeline_id, jlab_id, base_obs_id, datetime_start, datetime_end, datetime_file, daily_folder,
                                is_background=False):
    '''
        Basic (NIPS specific) function. This function generates for a given dither and nod position (i) nid, (iI) the folder name
        (creates it too if it not existing) and (iiI) the filenames (oldcts and one filename per sca).
        Note: The design may be adapted to accept FWA, GWA, like PF's example:
        #IFS example going through different configurations
        #obs_id = '{:s}-B-{:s}-{:s}'.format(base_obs_id, FWA, GWA)

        :param output_path:
        :param dither_idx:
        :param nod_idx:
        :param base_nid:
        :param pipeline_id:
        :param jlab_id:
        :param base_obs_id:
        :param datetime_start:
        :param datetime_end:
        :param datetime_file:
        :param daily_folder:
        :param is_background:
        :return:
        '''
    #TODO: bundle together output archive parameters: base_nid, pipeline_id, jlab_id, base_obs_id, datetime_start, datetime_end, datetime_file, daily_folder
    base_nid = 1000
    h_nid = 100 * dither_idx
    d_nid = nod_idx #* nexp + exp_idx  # this breaks if this number is larger than 100
    nid = base_nid + h_nid + d_nid + 1

    if is_background:
        obs_id = '{:s}-B-{:02d}'.format(base_obs_id, d_nid + 1)
    else:
        obs_id = '{:s}-{:02d}'.format(base_obs_id, d_nid + 1)


    # Generating folder name (following GG example)
    # NRSDEEP-PRM-02_1_1202_JW1_IPS_20180518T120000.000_20180518T124800.000
    folder_name = 'NRS{:s}_1_{:d}_{:s}_{:s}_{:s}_{:s}'.format(obs_id, nid, jlab_id, pipeline_id,
                                                              datetime_start, datetime_end)
                                                              # Generating folder
    work = os.path.join(output_path, daily_folder, folder_name)
    try:
        os.makedirs(work)
    except Exception as error:
        if (error.args[0] != 17):
            print('ERROR - Encountered an error when trying to create the exposure folder:')
            print('ERROR - {:s}'.format(work))
            print('ERROR - {}'.format(error.args))
            raise ValueError

    # Generating filename
    #TODO: make a namedtuple to output the filenames
    filenames = []

    # if I am not using old NIPS, there is no need to generate that format
    # filename = 'NRS_{:s}.oldcts.fits'.format(obs_id)
    # filenames.append(filename)

    for sca_index in range(2):
        filename = 'NRS{:s}_1_{:d}_SE_{:s}.cts.fits'.format(obs_id, 491 + sca_index, datetime_file)
        filenames.append(filename)

    return nid, obs_id, folder_name, filenames


def f_crm(input_path, input_list, FWA, GWA, sampling_grid, components, sel_apertures, nexp, multiaccum, references,
          output_path, base_nid, pipeline_id, jlab_id, base_obs_id, datetime_start, datetime_end, datetime_file, daily_folder,
          save_background=True):
    '''
        Basic (for now NIPS centric and commissioning specific (nexp~=1)) function to take all basic ingredients and
        perform
        - transform to count-rate-maps
        - deal with dither, nods and (potentially) exposures
        - create the archive structure, count-rate filenames and NID identifier
        - generate background exposure (save_background=True/False is hardcoded)
        - Include realistic dark flags
        - Populate GWA information (current workaround IPS>4.0 should take care of this)
        - Write the .fits files

        Note: Currently NIPS-centric functions are embedded in code #TODO: extract into functions

        :param input_path:
        :param input_list:
        :param FWA:
        :param GWA:
        :param sampling_grid:
        :param components:
        :param sel_apertures:
        :param nexp:
        :param multiaccum:
        :param references:
        :param output_path:
        :param base_nid:
        :param pipeline_id:
        :param jlab_id:
        :param base_obs_id:
        :param datetime_start:
        :param datetime_end:
        :param datetime_file:
        :param daily_folder:
        :param save_background:
        :return:
        '''

    #TODO: Idea, maybe promote this function upward to f_post_processing body, as it needs yet more parameters
    # =======================================================================
    # NIPS Archive structure: daily folder generation
    # =======================================================================
    work = os.path.join(output_path, daily_folder)
    try:
        os.makedirs(work)
    except Exception as error:
        if (error.args[0] != 17):
            print('ERROR - Encountered an error when trying to create the daily folder:')
            print('ERROR - {:s}'.format(work))
            print('ERROR - {}'.format(error.args))
            raise ValueError

    # =======================================================================
    # Generating the count-rate maps
    # =======================================================================
    print("# Generating the count-rate maps.")

    dither_pts = sampling_grid['dither_pts']
    nod_pts = sampling_grid['nod_pts']


    target_list = input_list['target_erm']
    if 'background' in components: bckg_list = input_list['bckg_erm']

    for dither_idx in range(dither_pts):
        if 'background' in components:
            # -------------------------------------------------------------------
            # Background ERM: we need to have a background to add, even if it is the same
            # -------------------------------------------------------------------
            b_erm = f_add_aper_erm(input_path, bckg_list, sel_apertures, dither_idx=dither_idx, verbose=True)

        for nod_idx in range(nod_pts):

            o_name = target_list[dither_idx][nod_idx]
            # -------------------------------------------------------------------
            # Target/object ERM
            # -------------------------------------------------------------------
            o_erm = c_electronratemap.ElectronRateMap()
            o_erm.m_read_from_fits(os.path.join(input_path, o_name))

            if 'background' in components:
                # -------------------------------------------------------------------
                # Obj ERM + Background ERM
                # -------------------------------------------------------------------
                o_erm.m_add(b_erm)
                # -------------------------------------------------------------------
                # Save first background ERM (mimmic off-source exposure for master background in IFS)
                # -------------------------------------------------------------------
                #TODO: boolean parameter save_background deeply buried in this function. Should be promoted to config file
                if save_background and (dither_idx == 0) and (nod_idx == 0):
                    b_ctms = erm_to_sci.f_standard(b_erm, multiaccum, noise_category='extended', noise_model=references.noise_model,
                                                   scale=1.0, nexp=nexp, planes=None, bias=None, use_superbias=False,
                                                   dark=references.dark, dark_sub=True, readout=references.readout,
                                                   gain=references.gain, ctm=references.ctm,
                                                   verbose=True, seed=-1, force_trimming=False, old=True)
                # -------------------------------------------------------------------
                # Generating archive folders
                # -------------------------------------------------------------------
                nid, obs_id, folder_name, filenames = f_output_archive_structure(output_path, dither_idx, nod_idx,
                                                                                 base_nid, pipeline_id, jlab_id,
                                                                                 base_obs_id, datetime_start, datetime_end,
                                                                                 datetime_file, daily_folder, is_background=True)

                # -------------------------------------------------------------------
                # Write out to archive folders
                # -------------------------------------------------------------------
                f_write_to_fits(b_ctms, output_path, daily_folder, folder_name, filenames, obs_id)

            # -------------------------------------------------------------------
            # Note: Just for NIPS2.0 and commissioning; otherwise introduce the
            # following for loop and change the f_standard parameter to nexp=1
            # GG: Generating the exposure folders and the associated exposures - 4 exposures for pointing
            # -------------------------------------------------------------------
            #for iexp in range(nexp):
            # GG: Note nexp in this function is only used to scale noise -> nexp=1

            o_ctms = erm_to_sci.f_standard(o_erm, multiaccum, noise_category='extended', noise_model=references.noise_model,
                                           scale=1.0, nexp=nexp, planes=None, bias=None, use_superbias=False,
                                           dark=references.dark, dark_sub=True, readout=references.readout, gain=references.gain,
                                           ctm=references.ctm, verbose=True, seed=-1,force_trimming=False, old=True)
            # -------------------------------------------------------------------
            # Generating archive folders
            # -------------------------------------------------------------------
            nid, obs_id, folder_name, filenames = f_output_archive_structure(output_path, dither_idx, nod_idx, base_nid,
                                                                             pipeline_id, jlab_id, base_obs_id, datetime_start,
                                                                             datetime_end, datetime_file, daily_folder,
                                                                             is_background=False)

            #TODO: Extract all below code into a NIPS specific function and include switch parameter in config file to branch off
            # -------------------------------------------------------------------
            # Populate header keywords: (specific to NIPS2.0)
            # - dark flags
            # - GWA information
            # -------------------------------------------------------------------
            flag_arrays = []

            # Reading in quality flags from one of our darks
            hdu = fits.open(_dark491fname)
            flag_arrays.append(hdu[3].data)
            print('# Reading in FLAG array 1')

            hdu = fits.open(_dark492fname)
            flag_arrays.append(hdu[3].data)
            print('# Reading in FLAG array 2')
            hdu.close()

            # Before writing out the cont-rate maps (new format for NIPS2.0) we need to update the GWA keywords because
            # the version of the IPS with which they were created did not populated the GWA keyword with the correct value
            # Therefore opening the model gtp file:

            #TODO: insert a try catch to skip this if IPS already populates it (not sure on which attribute)
            filename = _model_path + _model_name + '/Description/disperser_' + GWA + '_TiltY.gtp'
            poly = tiltpoly.TiltPoly()
            poly.m_readFromFile(filename)
            GWA_XTIL = poly.zeroreadings[0]  # dispersion direction

            filename = _model_path + _model_name + '/Description/disperser_' + GWA + '_TiltX.gtp'
            poly = tiltpoly.TiltPoly()
            poly.m_readFromFile(filename)
            GWA_YTIL = poly.zeroreadings[0]  # spatial direction


            mode_gg = 'unknown'
            slit_gg = 'unknown'
            aperture_gg = 'unknown'
            sca_ids = ['NRS1', 'NRS2']
            read_out = 'NRSRAPID'
            #TODO: Hardcoding warning: create a function of (readout_mode, nf) that spits out the cooked readout
            for o_idx in range(1, len(filenames)+1):
                o_ctms[o_idx].fits_header['GWA_XTIL'] = GWA_XTIL
                o_ctms[o_idx].fits_header['GWA_YTIL'] = GWA_YTIL
                o_ctms[o_idx].m_set_keywords(nid, 'IPS', sca_ids[o_idx-1], mode_gg, slit_gg, aperture_gg,
                                            FWA, GWA, read_out, False)
                #o_ctms[o_idx].quality = flag_arrays[o_idx-1]


            # -------------------------------------------------------------------
            # Write out to archive folders
            # We keep the old count-rate for NIPS 2.7 compatibility and plotting
            # -------------------------------------------------------------------
            f_write_to_fits(o_ctms, output_path, daily_folder, folder_name, filenames, obs_id)

    return None

def f_add_aper_erm(input_path, input_list, sel_apertures, dither_idx=0, verbose=False):
    '''
        Utility function to add the apertures (e.g. slit, mos) of a component (e.g. background)
        :param input_path:
        :param list input_list: 2-dimensional list of filenames with indices [aperture][dither]
        :param sel_apertures:
        :param dither_idx:
        :param verbose:
        :return:
        '''
    #this function is only applicable to background
    # and flatfield components, which depend on apertures
    #add all aperture components to the erm file
    #TODO: introduce logic to have always a full input component list and select by string, by checking the components in the global list
    for aperture_idx in range(len(sel_apertures)):
        filename = input_list[dither_idx][aperture_idx]
        if aperture_idx == 0:
            comp_erm = c_electronratemap.ElectronRateMap()
            comp_erm.m_read_from_fits(os.path.join(input_path, filename))
            if verbose: print('Standard component: {:s}'.format(sel_apertures[aperture_idx]))
        else:
            work_erm = c_electronratemap.ElectronRateMap()
            work_erm.m_read_from_fits(os.path.join(input_path, filename))
            comp_erm.m_add(work_erm)
            if verbose: print('Adding :{:s}'.format(input_list[dither_idx][aperture_idx]))

    return comp_erm

def f_list_dim(testlist, dim=0):
    """
    Utility function: tests if testlist is a list and how many dimensions it has
    returns -1 if it is no list at all, 0 if list is empty
    and otherwise the dimensions of it. All elements must be equal"""
    if isinstance(testlist, list):
        if testlist == []:
            return dim
        dim = dim + 1
        dim = f_list_dim(testlist[0], dim)
        return dim
    else:
        if dim == 0:
            return -1
        else:
            return dim


def f_print_dict(input_dict):
    for keys,values in input_dict.items():
        print(keys)
        pprint(values)

    return None


def f_modify_dict(input_list, mode, apertures, dither_pts, verbose=True):
    '''Utility function'''
    #TODO: check if dimensions of target_erm[0] is the same as other component, if not, broadcast them to match
    if verbose:
        print('Start with this dictionary....')
        f_print_dict(input_list)
    new_list = {}
    for field in input_list.keys():
        my_list = input_list[field]
        # check each keyword list is 2-dimensional, if it is 1-dimensional, we
        # assume that its is due to only one dither pos and make it 2-dimensional
        dim = f_list_dim(my_list)
        if dim == 1:
            if verbose: print('Reformating input list to make it 2-dimensional')
            my_list = [[file] for file in my_list]
        if field != 'target_erm':
            #1 broadcast to all dither positions
            if len(my_list[0]) == 1:
                my_list = [item * dither_pts for item in my_list]
            #2 Mask to check if all apertures are in one row (I changed the order of the for loops)
            aper_0 = apertures[mode][0]
            matching = [[1 if aper_0.lower() in aper[i].lower() else 0 for i in range(len(my_list[0]))] for aper in my_list]
            mymask = np.array(matching, dtype=bool)
            if (mymask.all(axis=1).any(axis=0)) and (len(my_list[0]) != 1):  # if it is true, we need to transpose
                # transpose
                my_list = [[aper[i] for aper in my_list] for i in range(len(my_list[0]))]
            else:
                # leave as it is
                pass
        new_list[field] = my_list
    if verbose:
        print('Finish with this dictionary....')
        f_print_dict(new_list)

    return new_list


def f_write_to_fits(ctms, output_path, daily_folder, folder_name, filenames, obs_id):
    '''
        Basic function to write count-rate maps compatible with NIPS2.0 (not with old NIPS)
        :param ctms:
        :param output_path:
        :param daily_folder:
        :param folder_name:
        :param filenames:
        :param obs_id:
        :return:
        '''
    #TODO (EP): Like for old NIPS, it would be good to create the display
    for idx in range(1,len(filenames)+1):
        # if idx == 0:
        #     ctms[idx].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filenames[idx]))
        #     ctms[idx].m_display_full(title=obs_id, filename=os.path.join(output_path, daily_folder, folder_name, '{:s}.pdf'.format(filenames[idx])),
        #                                  display=False, log_scale=True, category='data', comment=folder_name, figsize=(10.8, 8.5),
        #                                  dpi=100, axisbg='k', close=True, contour=False, cbar='horizontal', transparent=False,
        #                                  axis_font_size=None, time_stamp=True)
        # else:
        #     ctms[idx].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filenames[idx]))
        print('Writing to {:s}'.format(os.path.join(output_path, daily_folder, folder_name, filenames[idx-1])))
        ctms[idx].m_write_to_fits(os.path.join(output_path, daily_folder, folder_name, filenames[idx-1]))

    return None

def main():
    # =======================================================================
    # Parsing the input arguments
    # =======================================================================
    #TODO: This part has not been tested and some parameters are lagging behind, like subarray
    parser = argparse.ArgumentParser()
    parser.add_argument('inpath', type=str, help='Name of the folder where the input files are present.')
    listOfValidPositionsFWA = ['CLEAR', 'F070LP', 'F100LP', 'F170LP', 'F290LP', 'F140X', 'F110W']
    parser.add_argument('FWA', type=str, help='Filter wheel position string.', choices=listOfValidPositionsFWA)
    listOfValidPositionsGWA = ['PRISM', 'G140M', 'G235M', 'G395M', 'G140H', 'G235H', 'G395H', 'MIRROR']
    parser.add_argument('GWA', type=str, help='Filter wheel position string.', choices=listOfValidPositionsGWA)
    parser.add_argument('nexp', type=int, help='Number of exposures.')
    parser.add_argument('ng', type=int, help='Number of groups (NRSIRS2 readout mode).')
    parser.add_argument('outpath', type=str,
                        help='Name of the folder where the output files and folders will be generated.')
    parser.add_argument('obsid', type=str, help='Base for the generation of the observation ID.')
    listOfValidPositionsMode = ['FS', 'MOS', 'IFS']
    parser.add_argument('mode', type=str, help='Observing mode.', choices=listOfValidPositionsMode)
    #optional parameters
    parser.add_argument('-reado', '--readout', type=str,
                        help='Readout mode', choices= ['traditional', 'IRS2'], default='traditional')
    parser.add_argument('-expt', '--exposure_type', type=str,
                        help='Exposure type', choices= ['full-frame', 'subarray'], default='full-frame')
    parser.add_argument('-subid', '--subarray_id', type=str,
                        help='Optional subarray id', choices= ['FULL-FRAME', 'ALLSLITS', 'S200A1', 'S200A2', 'S200B1', 'S400A1', 'SUB1024A', 'SUB1024B', 'SUB2048', 'SUB32', 'SUB512', 'SUB512S']
                        , default=None)
    parser.add_argument('-s', '--seed', type=int,
                        help='Optional seed (integer) for the random number generator.', default=None)
    parser.add_argument('-b', '--background', type=str,
                        help='Suffix of the background electron-rate map.', default='background')
    parser.add_argument('-u', '--unity', type=str, help='Suffix of the unity electron-rate map.', default='unity')
    parser.add_argument('-o', '--object', type=str, help='Suffix of the object electron-rate map.', default='galaxy')
    parser.add_argument('-f', '--factor', type=float,
                        help='Optional scaling factor to be applied to the object cube.', default=1.0)

    parser.add_argument('-dts', '--datetime_start', type=str,
                        help='Start date and time.', default='20180303T120000.000')
    parser.add_argument('-dte', '--datetime_end', type=str,
                        help='End date and time.', default='20180303T120000.000')
    parser.add_argument('-dtf', '--datetime_file', type=str,
                        help='End date and time.', default='2018-03-03T12h15m00')
    parser.add_argument('-day', '--daily_folder', type=str,
                        help='Daily folder.', default='Day2018062')

    # =======================================================================
    # Loading the input arguments
    # =======================================================================
    args = parser.parse_args()
    argv = sys.argv
    narg = len(argv) - 1
    print("# =====================================")
    print("# Running :{:s}".format(argv[0]))
    print("# _name   : {:s}".format(_name))
    print("# _version: {:s}".format(_version))
    print("# Date    : {:s}".format((datetime.datetime.now()).isoformat()))
    print("# =====================================")
    # =======================================================================
    # Paths and filenames
    # =======================================================================
    input_path = args.inpath
    print("# Input path: {:s}".format(input_path))
    output_path = args.outpath
    print("# Output path: {:s}".format(output_path))
    FWA = args.FWA
    GWA = args.GWA
    if (GWA == 'PRISM'):
        filetype = 'prism'
    elif (GWA == 'MIRROR'):
        filetype = 'image'
    else:
        filetype = 'grating'
    mode = args.mode
    print("# Observing mode: {:s}".format(input_mode))
    print("# Type of input file: {:s}".format(filetype))
    print("# Spectral configuration: {:s}/{:s}".format(FWA, GWA))
    print("")
    print("# =====================================")
    base_obs_id = args.obsid
    print("# Base observation ID: {:s}".format(base_obs_id))
    datetime_start = args.datetime_start
    datetime_end = args.datetime_end
    print("# Datetime start and end: {} : {}".format(datetime_start, datetime_end))
    datetime_file = args.datetime_file
    print("# Datetime file: {}".format(datetime_file))
    daily_folder = args.daily_folder
    print("# Daily folder: {:s}".format(daily_folder))
    print("# =====================================")
    print("# Detector configuration:")
    nexp = args.nexp
    print("# Number of exposures: {:d}".format(nexp))
    print("")
    ng = args.ng
    print("# =====================================")
    print("# Optional arguments:")
    background_suffix = args.background
    print("# Suffix for the background electron-rate map: {:s}".format(background_suffix))
    unity_suffix = args.unity
    print("# Suffix for the unity electron-rate map: {:s}".format(unity_suffix))
    object_suffix = args.object
    print("# Suffix for the object electron-rate map: {:s}".format(unity_suffix))
    seed = args.seed
    if (seed is None):
        print("# No random generator seed provided.")
    else:
        print("# Random generator seed: {:d}".format(seed))
        np.random.seed(seed=seed)
    factor = args.factor
    print("# Normalisation factor: {:8.4e}".format(factor))
    print("# =====================================")

    f_post_processing(input_path, FWA, GWA, nexp, ng, output_path, base_obs_id, mode,
                      #several parameters missing, using the default in the function definition
                      seed=seed, input_background=background_suffix, input_unity=unity_suffix,
                      input_object=object_suffix, input_factor=factor, input_datetime_start=datetime_start,
                      input_datetime_end=datetime_end, input_datetime_file=datetime_file, input_daily_folder=daily_folder)


if __name__ == "__main__":
    sys.exit(main())

