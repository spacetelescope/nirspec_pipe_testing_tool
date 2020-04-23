import numpy as np
import os
import subprocess
import argparse
import sys
from collections import OrderedDict
from astropy.io import fits
from astropy import wcs

from jwst.assign_wcs import nirspec
from jwst import datamodels

from . import auxiliary_functions as auxfunc


"""
This script compares pipeline WCS info with ESA results for Multi-Object Spectroscopy (MOS) data.

"""


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "2.1"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Completely changed script to use the datamodel instead of the compute_world_coordinates
#                         script, and added new routines for plot making and statistics calculations.
# Aug 2018 - Version 2.1: Modified slit-y differences to be reported in absolute numbers rather than relative


def compare_wcs(infile_name, esa_files_path, msa_conf_name, show_figs=True, save_figs=False, threshold_diff=1.0e-7,
                mode_used=None, raw_data_root_file=None, output_directory=None, debug=False):
    """
    This function does the WCS comparison from the world coordinates calculated using the pipeline
    data model with the ESA intermediary files.

    Args:
        infile_name: str, name of the output fits file from the assign_wcs step (with full path)
        esa_files_path: str, full path of where to find all ESA intermediary products to make comparisons for the tests
        msa_conf_name: str, full path where to find the shutter configuration file
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots or not
        threshold_diff: float, threshold difference between pipeline output and ESA file
        mode_used: string, mode used in the PTT configuration file
        raw_data_root_file: None or string, name of the raw file that produced the _uncal.fits file for caldetector1
        output_directory: None or string, path to the output_directory where to save the plots
        debug: boolean, if true a series of print statements will show on-screen

    Returns:
        - plots, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to threshold
        - log_msgs: list, all print statements captured in this variable

    """

    log_msgs = []

    # get grating and filter info from the rate file header
    if isinstance(infile_name, str):
        # get the datamodel from the assign_wcs output file
        img = datamodels.ImageModel(infile_name)
        msg = 'infile_name='+infile_name
        print(msg)
        log_msgs.append(msg)
        # check that shutter configuration file in header is the same as given in PTT_config file
        msametfl = fits.getval(infile_name, "MSAMETFL", 0)
        if msametfl != os.path.basename(msa_conf_name):
            msg = "* WARNING! MSA config file name given in PTT_config file does not match the MSAMETFL " \
                  "keyword in main header.\n"
            print(msg)
            log_msgs.append(msg)
        basenameinfile_name = os.path.basename(infile_name)

    else:
        img = infile_name
        basenameinfile_name = ""

    det = img.meta.instrument.detector   # fits.getval(infile_name, "DETECTOR", 0)
    lamp = img.meta.instrument.lamp_state   # fits.getval(infile_name, "LAMP", 0)
    grat = img.meta.instrument.grating   # fits.getval(infile_name, "GRATING", 0)
    filt = img.meta.instrument.filter   # fits.getval(infile_name, "FILTER", 0)
    msg = "from assign_wcs file  -->     Detector: " + det + "   Grating: " + grat + "   Filter: " + \
          filt + "   Lamp: " + lamp
    print(msg)
    log_msgs.append(msg)

    # copy the MSA shutter configuration file into the pytest directory
    try:
        subprocess.run(["cp", msa_conf_name, "."])
    except FileNotFoundError:
        msg1 = " * PTT is not able to locate the MSA shutter configuration file. Please make sure that the " \
               "msa_conf_name variable in"
        msg2 = "   the PTT_config.cfg file is pointing exactly to where the fits file exists (i.e. full " \
               "path and name). "
        msg3 = "   -> The WCS test is now set to skip and no plots will be generated. "
        print(msg1)
        print(msg2)
        print(msg3)
        log_msgs.append(msg1)
        log_msgs.append(msg2)
        log_msgs.append(msg3)
        FINAL_TEST_RESULT = "skip"
        return FINAL_TEST_RESULT, log_msgs

    # get shutter info from metadata
    shutter_info = fits.getdata(msa_conf_name, extname="SHUTTER_INFO") # this is generally ext=2
    pslit = shutter_info.field("slitlet_id")
    quad = shutter_info.field("shutter_quadrant")
    row = shutter_info.field("shutter_row")
    col = shutter_info.field("shutter_column")
    msg = 'Using this MSA shutter configuration file: '+msa_conf_name
    print(msg)
    log_msgs.append(msg)

    # get the datamodel from the assign_wcs output file
    if mode_used is None or mode_used != "MOS_sim":
        img = datamodels.ImageModel(infile_name)
        # these commands only work for the assign_wcs output file loop over the slits
        # slits_list = nirspec.get_open_slits(img)
        # this function returns all open slitlets as defined in msa meta file, however, some of them may
        # not be projected on the detector, and those are later removed from the list of open slitlets. To
        # get the open and projected on the detector slitlets we use the following:
        slits_list = img.meta.wcs.get_transform('gwa', 'slit_frame').slits
        #print ('Open slits: ', slits_list, '\n')

        if debug:
            print("Instrument Configuration")
            print("Detector: {}".format(img.meta.instrument.detector))
            print("GWA: {}".format(img.meta.instrument.grating))
            print("Filter: {}".format(img.meta.instrument.filter))
            print("Lamp: {}".format(img.meta.instrument.lamp_state))
            print("GWA_XTILT: {}".format(img.meta.instrument.gwa_xtilt))
            print("GWA_YTILT: {}".format(img.meta.instrument.gwa_ytilt))
            print("GWA_TTILT: {}".format(img.meta.instrument.gwa_tilt))

    elif mode_used == "MOS_sim":
        if isinstance(infile_name, str):
            # this command works for the extract_2d and flat_field output files
            model = datamodels.MultiSlitModel(infile_name)
            slits_list = model.slits
        else:
            slits_list = infile_name.slits

    # list to determine if pytest is passed or not
    total_test_result = OrderedDict()

    # loop over the slices
    for slit in slits_list:
        name = slit.name
        msg = "\nWorking with slit: "+str(name)
        print(msg)
        log_msgs.append(msg)

        # get the right index in the list of open shutters
        pslit_list = pslit.tolist()
        slitlet_idx = pslit_list.index(int(name))

        # Get the ESA trace
        if raw_data_root_file is None:
            _, raw_data_root_file = auxfunc.get_modeused_and_rawdatrt_PTT_cfg_file(infile_name)
        msg = "Using this raw data file to find the corresponding ESA file: "+raw_data_root_file
        print(msg)
        log_msgs.append(msg)
        q, r, c = quad[slitlet_idx], row[slitlet_idx], col[slitlet_idx]
        msg = "Pipeline shutter info:   quadrant= "+str(q)+"   row= "+str(r)+"   col="+str(c)
        print(msg)
        log_msgs.append(msg)
        specifics = [q, r, c]
        esafile = auxfunc.get_esafile(esa_files_path, raw_data_root_file, "MOS", specifics)

        # skip the test if the esafile was not found
        if "ESA file not found" in esafile:
            msg1 = " * compare_wcs_mos.py is exiting because the corresponding ESA file was not found."
            msg2 = "   -> The WCS test is now set to skip and no plots will be generated. "
            print(msg1)
            print(msg2)
            log_msgs.append(msg1)
            log_msgs.append(msg2)
            FINAL_TEST_RESULT = "skip"
            return FINAL_TEST_RESULT, log_msgs

        # Open the trace in the esafile
        if len(esafile) == 2:
            print(len(esafile[-1]))
            if len(esafile[-1]) == 0:
                esafile = esafile[0]
        msg = "Using this ESA file: \n"+str(esafile)
        print(msg)
        log_msgs.append(msg)
        with fits.open(esafile) as esahdulist:
            print("* ESA file contents ")
            esahdulist.info()
            esa_shutter_i = esahdulist[0].header['SHUTTERI']
            esa_shutter_j = esahdulist[0].header['SHUTTERJ']
            esa_quadrant = esahdulist[0].header['QUADRANT']
            if debug:
                msg = "ESA shutter info:   quadrant=" + esa_quadrant + "   shutter_i="+esa_shutter_i + \
                      "   shutter_j="+esa_shutter_j
                print(msg)
                log_msgs.append(msg)
            # first check if ESA shutter info is the same as pipeline
            msg = "For slitlet"+str(name)
            print(msg)
            log_msgs.append(msg)
            if q == esa_quadrant:
                msg = "\n -> Same quadrant for pipeline and ESA data: "+str(q)
                print(msg)
                log_msgs.append(msg)
            else:
                msg = "\n -> Missmatch of quadrant for pipeline and ESA data: "+str(q)+esa_quadrant
                print(msg)
                log_msgs.append(msg)
            if r == esa_shutter_i:
                msg = "\n -> Same row for pipeline and ESA data: "+str(r)
                print(msg)
                log_msgs.append(msg)
            else:
                msg = "\n -> Missmatch of row for pipeline and ESA data: "+str(r)+esa_shutter_i
                print(msg)
                log_msgs.append(msg)
            if c == esa_shutter_j:
                msg = "\n -> Same column for pipeline and ESA data: "+str(c)+"\n"
                print(msg)
                log_msgs.append(msg)
            else:
                msg = "\n -> Missmatch of column for pipeline and ESA data: "+str(c)+esa_shutter_j+"\n"
                print(msg)
                log_msgs.append(msg)

            # Assign variables according to detector
            skipv2v3test = True
            if det == "NRS1":
                try:
                    esa_flux = fits.getdata(esafile, "DATA1")
                    esa_wave = fits.getdata(esafile, "LAMBDA1")
                    esa_slity = fits.getdata(esafile, "SLITY1")
                    esa_msax = fits.getdata(esafile, "MSAX1")
                    esa_msay = fits.getdata(esafile, "MSAY1")
                    pyw = wcs.WCS(esahdulist['LAMBDA1'].header)
                    try:
                        esa_v2v3x = fits.getdata(esafile, "V2V3X1")
                        esa_v2v3y = fits.getdata(esafile, "V2V3Y1")
                        skipv2v3test = False
                    except KeyError:
                        msg = "Skipping tests for V2 and V3 because ESA file does not contain corresponding extensions."
                        print(msg)
                        log_msgs.append(msg)
                except KeyError:
                    msg = "PTT did not find ESA extensions that match detector NRS1, skipping test for this slitlet..."
                    print(msg)
                    log_msgs.append(msg)
                    continue

            if det == "NRS2":
                try:
                    esa_flux = fits.getdata(esafile, "DATA2")
                    esa_wave = fits.getdata(esafile, "LAMBDA2")
                    esa_slity = fits.getdata(esafile, "SLITY2")
                    esa_msax = fits.getdata(esafile, "MSAX2")
                    esa_msay = fits.getdata(esafile, "MSAY2")
                    pyw = wcs.WCS(esahdulist['LAMBDA2'].header)
                    try:
                        esa_v2v3x = fits.getdata(esafile, "V2V3X2")
                        esa_v2v3y = fits.getdata(esafile, "V2V3Y2")
                        skipv2v3test = False
                    except KeyError:
                        msg = "Skipping tests for V2 and V3 because ESA file does not contain corresponding extensions."
                        print(msg)
                        log_msgs.append(msg)
                except KeyError:
                    msg = "PTT did not find ESA extensions that match detector NRS2, skipping test for this slitlet..."
                    print(msg)
                    log_msgs.append(msg)
                    continue

        # get the WCS object for this particular slit
        if mode_used is None  or  mode_used != "MOS_sim":
            try:
                wcs_slice = nirspec.nrs_wcs_set_input(img, name)
            except:
                ValueError
                msg = "* WARNING: Slitlet "+name+" was not found in the model. Skipping test for this slitlet."
                print(msg)
                log_msgs.append(msg)
                continue
        elif mode_used == "MOS_sim":
            wcs_slice = model.slits[0].wcs

        # if we want to print all available transforms, uncomment line below
        #print(wcs_slice)

        # The WCS object attribute bounding_box shows all valid inputs, i.e. the actual area of the data according
        # to the slice. Inputs outside of the bounding_box return NaN values.
        #bbox = wcs_slice.bounding_box
        #print('wcs_slice.bounding_box: ', wcs_slice.bounding_box)

        # In different observing modes the WCS may have different coordinate frames. To see available frames
        # uncomment line below.
        #print("Avalable frames: ", wcs_slice.available_frames)

        if mode_used is None  or  mode_used != "MOS_sim":
            if debug:
                # To get specific pixel values use following syntax:
                det2slit = wcs_slice.get_transform('detector', 'slit_frame')
                slitx, slity, lam = det2slit(700, 1080)
                print("slitx: " , slitx)
                print("slity: " , slity)
                print("lambda: " , lam)

            if debug:
                # The number of inputs and outputs in each frame can vary. This can be checked with:
                print('Number on inputs: ', det2slit.n_inputs)
                print('Number on outputs: ', det2slit.n_outputs)

        # Create x, y indices using the Trace WCS
        pipey, pipex = np.mgrid[:esa_wave.shape[0], : esa_wave.shape[1]]
        esax, esay = pyw.all_pix2world(pipex, pipey, 0)

        if det == "NRS2":
            msg = "NRS2 needs a flip"
            print(msg)
            log_msgs.append(msg)
            esax = 2049-esax
            esay = 2049-esay

        # Compute pipeline RA, DEC, and lambda
        slitlet_test_result_list = []
        pra, pdec, pwave = wcs_slice(esax-1, esay-1)
        # => RETURNS: RA, DEC, LAMBDA (lam *= 10**-6 to convert to microns)
        pwave *= 10**-6
        # calculate and print statistics for slit-y and x relative differences
        slitlet_name = repr(r)+"_"+repr(c)
        tested_quantity = "Wavelength Difference"
        rel_diff_pwave_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_wave, pwave,
                                                               tested_quantity)
        rel_diff_pwave_img, notnan_rel_diff_pwave, notnan_rel_diff_pwave_stats, print_stats_strings = rel_diff_pwave_data
        for msg in print_stats_strings:
            log_msgs.append(msg)
        result = auxfunc.does_median_pass_tes(notnan_rel_diff_pwave_stats[1], threshold_diff)
        msg = 'Result for test of ' + tested_quantity + ': ' + result
        print(msg)
        log_msgs.append(msg)
        slitlet_test_result_list.append({tested_quantity: result})

        # get the transforms for pipeline slit-y
        det2slit = wcs_slice.get_transform('detector', 'slit_frame')
        slitx, slity, _ = det2slit(esax-1, esay-1)
        tested_quantity = "Slit-Y Difference"
        # calculate and print statistics for slit-y and x relative differences
        rel_diff_pslity_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_slity, slity,
                                                                tested_quantity, abs=False)
        # calculate and print statistics for slit-y and x absolute differences
        #rel_diff_pslity_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_slity, slity, tested_quantity, abs=True)
        rel_diff_pslity_img, notnan_rel_diff_pslity, notnan_rel_diff_pslity_stats, print_stats_strings = rel_diff_pslity_data
        for msg in print_stats_strings:
            log_msgs.append(msg)
        result = auxfunc.does_median_pass_tes(notnan_rel_diff_pslity_stats[1], threshold_diff)
        msg = 'Result for test of ' + tested_quantity + ': ' + result
        print(msg)
        log_msgs.append(msg)
        slitlet_test_result_list.append({tested_quantity: result})

        # do the same for MSA x, y and V2, V3
        detector2msa = wcs_slice.get_transform("detector", "msa_frame")
        pmsax, pmsay, _ = detector2msa(esax-1, esay-1)
        # => RETURNS: msaX, msaY, LAMBDA (lam *= 10**-6 to convert to microns)
        # MSA-x
        tested_quantity = "MSA_X Difference"
        reldiffpmsax_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_msax, pmsax,
                                                             tested_quantity)
        reldiffpmsax_img, notnan_reldiffpmsax, notnan_reldiffpmsax_stats, print_stats_strings = reldiffpmsax_data
        for msg in print_stats_strings:
            log_msgs.append(msg)
        result = auxfunc.does_median_pass_tes(notnan_reldiffpmsax_stats[1], threshold_diff)
        msg = 'Result for test of ' + tested_quantity + ': ' + result
        print(msg)
        log_msgs.append(msg)
        slitlet_test_result_list.append({tested_quantity: result})
        # MSA-y
        tested_quantity = "MSA_Y Difference"
        reldiffpmsay_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_msay, pmsay,
                                                             tested_quantity)
        reldiffpmsay_img, notnan_reldiffpmsay, notnan_reldiffpmsay_stats, print_stats_strings = reldiffpmsay_data
        for msg in print_stats_strings:
            log_msgs.append(msg)
        result = auxfunc.does_median_pass_tes(notnan_reldiffpmsay_stats[1], threshold_diff)
        msg = 'Result for test of ' + tested_quantity + ': ' + result
        print(msg)
        log_msgs.append(msg)
        slitlet_test_result_list.append({tested_quantity: result})

        # V2 and V3
        if not skipv2v3test:
            detector2v2v3 = wcs_slice.get_transform("detector", "v2v3")
            pv2, pv3, _ = detector2v2v3(esax-1, esay-1)
            # => RETURNS: v2, v3, LAMBDA (lam *= 10**-6 to convert to microns)
            tested_quantity = "V2 difference"
            # converting to degrees to compare with ESA
            reldiffpv2_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_v2v3x, pv2,
                                                               tested_quantity)
            if reldiffpv2_data[-2][0] > 0.0:
                print("\nConverting pipeline results to degrees to compare with ESA")
                pv2 = pv2/3600.
                reldiffpv2_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_v2v3x, pv2,
                                                                   tested_quantity)
            reldiffpv2_img, notnan_reldiffpv2, notnan_reldiffpv2_stats, print_stats_strings = reldiffpv2_data
            for msg in print_stats_strings:
                log_msgs.append(msg)
            result = auxfunc.does_median_pass_tes(notnan_reldiffpv2_stats[1], threshold_diff)
            msg = 'Result for test of '+tested_quantity+': '+result
            print(msg)
            log_msgs.append(msg)
            slitlet_test_result_list.append({tested_quantity: result})
            tested_quantity = "V3 difference"
            # converting to degrees to compare with ESA
            reldiffpv3_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_v2v3y, pv3,
                                                               tested_quantity)
            if reldiffpv3_data[-2][0] > 0.0:
                print("\nConverting pipeline results to degrees to compare with ESA")
                pv3 = pv3/3600.
                reldiffpv3_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_v2v3y, pv3,
                                                                   tested_quantity)
            reldiffpv3_img, notnan_reldiffpv3, notnan_reldiffpv3_stats, print_stats_strings = reldiffpv3_data
            for msg in print_stats_strings:
                log_msgs.append(msg)
            result = auxfunc.does_median_pass_tes(notnan_reldiffpv3_stats[1], threshold_diff)
            msg = 'Result for test of '+tested_quantity+': '+result
            print(msg)
            log_msgs.append(msg)
            slitlet_test_result_list.append({tested_quantity: result})

        total_test_result[slitlet_name] = slitlet_test_result_list

        # PLOTS
        if show_figs or save_figs:
            # set the common variables
            main_title = filt+"   "+grat+"   SLITLET="+slitlet_name+"\n"
            bins = 15   # binning for the histograms, if None the function will automatically calculate them
            #             lolim_x, uplim_x, lolim_y, uplim_y
            plt_origin = None

            # Wavelength
            title = main_title+r"Relative wavelength difference = $\Delta \lambda$"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta \lambda$ = ($\lambda_{pipe} - \lambda_{ESA}) / \lambda_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_rel_diff_pwave_stats]
            if notnan_rel_diff_pwave_stats[1] is np.nan:
                msg = "Unable to create plot of relative wavelength difference."
                print(msg)
                log_msgs.append(msg)
            else:
                specific_plt_name = "_rel_wave_diffs.jpg"
                if isinstance(infile_name, str):
                    plt_name = infile_name.replace(basenameinfile_name, slitlet_name+"_"+det+specific_plt_name)
                else:
                    if output_directory is not None:
                        plt_name = os.path.join(output_directory, slitlet_name+"_"+det+specific_plt_name)
                    else:
                        plt_name = None
                        save_figs = False
                        print("No output_directory was provided. Figures will NOT be saved.")
                auxfunc.plt_two_2Dimgandhist(rel_diff_pwave_img, notnan_rel_diff_pwave, info_img, info_hist,
                                             plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs,
                                             save_figs=save_figs)

            # Slit-y
            title = main_title+r"Relative slit position = $\Delta$slit_y"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$slit_y = (slit_y$_{pipe}$ - slit_y$_{ESA}$)/slit_y$_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_rel_diff_pslity_stats]
            if notnan_rel_diff_pslity_stats[1] is np.nan:
                msg = "Unable to create plot of relative slit-y difference."
                print(msg)
                log_msgs.append(msg)
            else:
                specific_plt_name = "_rel_slitY_diffs.jpg"
                if isinstance(infile_name, str):
                    plt_name = infile_name.replace(basenameinfile_name, slitlet_name+"_"+det+specific_plt_name)
                else:
                    if output_directory is not None:
                        plt_name = os.path.join(output_directory, slitlet_name+"_"+det+specific_plt_name)
                    else:
                        plt_name = None
                        save_figs = False
                        print("No output_directory was provided. Figures will NOT be saved.")
                auxfunc.plt_two_2Dimgandhist(rel_diff_pslity_img, notnan_rel_diff_pslity, info_img, info_hist,
                                             plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs,
                                             save_figs=save_figs)

            # MSA-x
            title = main_title+r"Relative MSA-x Difference = $\Delta$MSA_x"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$MSA_x = (MSA_x$_{pipe}$ - MSA_x$_{ESA}$)/MSA_x$_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_reldiffpmsax_stats]
            if notnan_reldiffpmsax_stats[1] is np.nan:
                msg = "Unable to create plot of relative MSA-x difference."
                print(msg)
                log_msgs.append(msg)
            else:
                specific_plt_name = "_rel_MSAx_diffs.jpg"
                if isinstance(infile_name, str):
                    plt_name = infile_name.replace(basenameinfile_name, slitlet_name+"_"+det+specific_plt_name)
                else:
                    if output_directory is not None:
                        plt_name = os.path.join(output_directory, slitlet_name+"_"+det+specific_plt_name)
                    else:
                        plt_name = None
                        save_figs = False
                        print("No output_directory was provided. Figures will NOT be saved.")
                auxfunc.plt_two_2Dimgandhist(reldiffpmsax_img, notnan_reldiffpmsax, info_img, info_hist,
                                             plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs,
                                             save_figs=save_figs)

            # MSA-y
            title = main_title+r"Relative MSA-y Difference = $\Delta$MSA_y"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$MSA_y = (MSA_y$_{pipe}$ - MSA_y$_{ESA}$)/MSA_y$_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_reldiffpmsay_stats]
            if notnan_reldiffpmsay_stats[1] is np.nan:
                msg = "Unable to create plot of relative MSA-y difference."
                print(msg)
                log_msgs.append(msg)
            else:
                specific_plt_name = "_rel_MSAy_diffs.jpg"
                if isinstance(infile_name, str):
                    plt_name = infile_name.replace(basenameinfile_name, slitlet_name+"_"+det+specific_plt_name)
                else:
                    if output_directory is not None:
                        plt_name = os.path.join(output_directory, slitlet_name+"_"+det+specific_plt_name)
                    else:
                        plt_name = None
                        save_figs = False
                        print("No output_directory was provided. Figures will NOT be saved.")
                auxfunc.plt_two_2Dimgandhist(reldiffpmsay_img, notnan_reldiffpmsay, info_img, info_hist,
                                             plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs,
                                             save_figs=save_figs)

            if not skipv2v3test:
                # V2
                title = main_title+r"Relative V2 Difference = $\Delta$V2"+"\n"
                info_img = [title, "x (pixels)", "y (pixels)"]
                xlabel, ylabel = r"Relative $\Delta$V2 = (V2$_{pipe}$ - V2$_{ESA}$)/V2$_{ESA}$", "N"
                hist_data = notnan_reldiffpv2
                info_hist = [xlabel, ylabel, bins, notnan_reldiffpv2_stats]
                if notnan_reldiffpv2_stats[1] is np.nan:
                    msg = "Unable to create plot of relative V2 difference."
                    print(msg)
                    log_msgs.append(msg)
                else:
                    specific_plt_name = "_rel_V2_diffs.jpg"
                    if isinstance(infile_name, str):
                        plt_name = infile_name.replace(basenameinfile_name,
                                                       slitlet_name + "_" + det + specific_plt_name)
                    else:
                        if output_directory is not None:
                            plt_name = os.path.join(output_directory, slitlet_name + "_" + det + specific_plt_name)
                        else:
                            plt_name = None
                            save_figs = False
                            print("No output_directory was provided. Figures will NOT be saved.")
                    auxfunc.plt_two_2Dimgandhist(reldiffpv2_img, hist_data, info_img, info_hist,
                                                 plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs,
                                                 save_figs=save_figs)

                # V3
                title = main_title+r"Relative V3 Difference = $\Delta$V3"+"\n"
                info_img = [title, "x (pixels)", "y (pixels)"]
                xlabel, ylabel = r"Relative $\Delta$V3 = (V3$_{pipe}$ - V3$_{ESA}$)/V3$_{ESA}$", "N"
                hist_data = notnan_reldiffpv3
                info_hist = [xlabel, ylabel, bins, notnan_reldiffpv3_stats]
                if notnan_reldiffpv3_stats[1] is np.nan:
                    msg = "Unable to create plot of relative V3 difference."
                    print(msg)
                    log_msgs.append(msg)
                else:
                    specific_plt_name = "_rel_V3_diffs.jpg"
                    if isinstance(infile_name, str):
                        plt_name = infile_name.replace(basenameinfile_name,
                                                       slitlet_name + "_" + det + specific_plt_name)
                    else:
                        if output_directory is not None:
                            plt_name = os.path.join(output_directory, slitlet_name + "_" + det + specific_plt_name)
                        else:
                            plt_name = None
                            save_figs = False
                            print("No output_directory was provided. Figures will NOT be saved.")
                    auxfunc.plt_two_2Dimgandhist(reldiffpv3_img, hist_data, info_img, info_hist,
                                                 plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs,
                                                 save_figs=save_figs)

        else:
            msg = "NO plots were made because show_figs and save_figs were both set to False. \n"
            print(msg)
            log_msgs.append(msg)

    # remove the copy of the MSA shutter configuration file
    if isinstance(infile_name, str):
        subprocess.run(["rm", msametfl])

    # If all tests passed then pytest will be marked as PASSED, else it will be FAILED
    FINAL_TEST_RESULT = "FAILED"
    for sl, testlist in total_test_result.items():
        for tdict in testlist:
            for t, tr in tdict.items():
                if tr == "FAILED":
                    FINAL_TEST_RESULT = "FAILED"
                    msg = "\n * The test of "+t+" for slitlet "+sl+"  FAILED."
                    print(msg)
                    log_msgs.append(msg)
                else:
                    FINAL_TEST_RESULT = "PASSED"
                    msg = "\n * The test of "+t+" for slitlet "+sl+ "  PASSED."
                    print(msg)
                    log_msgs.append(msg)

    if FINAL_TEST_RESULT == "PASSED":
        msg = "\n *** Final result for assign_wcs test will be reported as PASSED *** \n"
        print(msg)
        log_msgs.append(msg)
    else:
        msg = "\n *** Final result for assign_wcs test will be reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)

    return FINAL_TEST_RESULT, log_msgs


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("infile_name",
                        action='store',
                        default=None,
                        help='Name of input fits file prior to assign_wcs step, i.e. blah_rate.fits')
    parser.add_argument("esa_files_path",
                        action='store',
                        default=None,
                        help='Path were to locate the benchmark files for comparison')
    parser.add_argument("msa_conf_name",
                        action='store',
                        default=None,
                        help='Name of the MSA shutter configuration file in pipeline format, e.g. blah_msa.fits')
    parser.add_argument("-f",
                        dest="show_figs",
                        action='store_true',
                        default=False,
                        help='Use flag -s to show on figures.')
    parser.add_argument("-s",
                        dest="save_figs",
                        action='store_false',
                        default=True,
                        help='Use flag -s to save on figures.')
    parser.add_argument("-t",
                        dest="threshold_diff",
                        action='store',
                        default=1.0e-7,
                        type=float,
                        help='Use flag -t change the default threshold (currently set to 1.0e-7).')
    parser.add_argument("-m",
                        dest="MOS_sim",
                        action='store_true',
                        default=False,
                        help='Use flag -m to turn on mode for MOS simulations.')
    parser.add_argument("-r",
                        dest="raw_data_root_file",
                        action='store',
                        default=None,
                        help='Use flag -r to give a raw_data_root_file.')
    parser.add_argument("-o",
                        dest="output_directory",
                        action='store',
                        default=None,
                        help='Use flag -o to provide the output_directory to save the figs.')
    parser.add_argument("-d",
                        dest="debug",
                        action='store_true',
                        default=False,
                        help='Use flag -d to turn on debug mode.')
    args = parser.parse_args()

    # Set the variables input from the command line
    infile_name = args.infile_name
    esa_files_path = args.esa_files_path
    msa_conf_name = args.msa_conf_name
    show_figs = args.show_figs
    save_figs = args.save_figs
    threshold_diff = args.threshold_diff
    MOS_sim = args.MOS_sim
    raw_data_root_file = args.raw_data_root_file
    output_directory = args.output_directory
    debug = args.debug

    # choose None or MOS_sim, only for MOS simulations
    if MOS_sim:
        mode_used = "MOS_sim"
    else:
        mode_used = "MOS"

    # Run the principal function of the script
    compare_wcs(infile_name, esa_files_path=esa_files_path, msa_conf_name=msa_conf_name,
                show_figs=show_figs, save_figs=save_figs, threshold_diff=threshold_diff,
                mode_used=mode_used, raw_data_root_file=raw_data_root_file,
                output_directory=output_directory, debug=debug)


if __name__ == '__main__':
    sys.exit(main())

