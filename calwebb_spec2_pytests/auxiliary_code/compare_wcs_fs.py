import numpy as np
import os
from astropy.io import fits
from astropy import wcs
from collections import OrderedDict

from jwst.assign_wcs import nirspec
from jwst import datamodels
from . import auxiliary_functions as auxfunc


"""
This script compares pipeline WCS info with ESA results for FIXED SLIT.

"""


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "2.4"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Completely changed script to use the datamodel instead of the compute_world_coordinates
#                         script, and added new routines for plot making and statistics calculations.
# Jun 2018 - Version 2.1: Fixed code to work for subarrays too.
# Aug 2018 - Version 2.2: Modified slit-y differences to be reported in absolute numbers rather than relative
# Apr 2019 - Version 2.3: Added capability to log on-screen messages
# May 2019 - Version 2.4: Fixed indexing missmatch with pipeline and ESA ALLSLITS subarray files


def compare_wcs(infile_name, esa_files_path=None, show_figs=True, save_figs=False, threshold_diff=1.0e-7, debug=False):
    """
    This function does the WCS comparison from the world coordinates calculated using the pipeline
    data model with the ESA intermediary files.

    Args:
        infile_name: str, name of the output fits file from the assign_wcs step (with full path)
        esa_files_path: str, full path of where to find all ESA intermediary products to make comparisons for the tests
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots or not
        threshold_diff: float, threshold difference between pipeline output and ESA file
        debug: boolean, if true a series of print statements will show on-screen

    Returns:
        - plots, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to threshold
        - log_msgs: list, all print statements captured in this variable

    """

    log_msgs = []

    # get grating and filter info from the rate file header
    det = fits.getval(infile_name, "DETECTOR", 0)
    msg = 'infile_name='+infile_name
    print(msg)
    log_msgs.append(msg)
    lamp = fits.getval(infile_name, "LAMP", 0)
    grat = fits.getval(infile_name, "GRATING", 0)
    filt = fits.getval(infile_name, "FILTER", 0)
    msg = "from assign_wcs file  -->     Detector: "+det+"   Grating: "+grat+"   Filter: "+filt+"   Lamp: "+lamp
    print(msg)
    log_msgs.append(msg)

    # mapping the ESA slit names to pipeline names
    map_slit_names = {'SLIT_A_1600' : 'S1600A1',
                      'SLIT_A_200_1': 'S200A1',
                      'SLIT_A_200_2': 'S200A2',
                      'SLIT_A_400':   'S400A1',
                      'SLIT_B_200':   'S200B1',
                      }

    # list to determine if pytest is passed or not
    total_test_result = OrderedDict()

    # get the datamodel from the assign_wcs output file
    img = datamodels.ImageModel(infile_name)

    # To get the open and projected on the detector
    open_slits = img.meta.wcs.get_transform('gwa', 'slit_frame').slits
    for opslit in open_slits:
        pipeslit = opslit.name
        msg = "\nWorking with slit: "+pipeslit
        print(msg)
        log_msgs.append(msg)

        # Get the ESA trace
        #raw_data_root_file = "NRSV84600010001P0000000002101_4_491_SE_2016-01-17T17h34m08.fits"  # for testing with G140M FULLFRAME
        #raw_data_root_file = "NRSSMOS-MOD-G1H-02-5344031756_1_491_SE_2015-12-10T03h25m56.fits"  # for testing with G140H FULLFRAME
        #raw_data_root_file = "NRSSDRK-ALLSLITS-5345150216_1_491_SE_2015-12-11T15h40m25.fits"  # for testing with G140H ALLSLITS
        #raw_data_root_file = "NRSV84600002001P0000000002101_1_491_SE_2016-01-17T15h09m16.fits"  # for testing with G140M ALLSLITS
        #raw_data_root_file = "NRSV84600004001P0000000002101_1_491_SE_2016-01-17T15h41m16.fits"  # for testing with G235H ALLSLITS
        _, raw_data_root_file = auxfunc.get_modeused_and_rawdatrt_PTT_cfg_file()
        specifics = [pipeslit]

        # check if ESA data is not in the regular directory tree, these files are exceptions
        NIDs = ["30055", "30055", "30205", "30133", "30133"]
        special_cutout_files = ["NRSSMOS-MOD-G1H-02-5344031756_1_491_SE_2015-12-10T03h25m56.fits",
                                "NRSSMOS-MOD-G1H-02-5344031756_1_492_SE_2015-12-10T03h25m56.fits",
                                "NRSSMOS-MOD-G2M-01-5344191938_1_491_SE_2015-12-10T19h29m26.fits",
                                "NRSSMOS-MOD-G3H-02-5344120942_1_491_SE_2015-12-10T12h18m25.fits",
                                "NRSSMOS-MOD-G3H-02-5344120942_1_492_SE_2015-12-10T12h18m25.fits"]
        if raw_data_root_file in special_cutout_files:
            nid = NIDs[special_cutout_files.index(raw_data_root_file)]
            msg = "Using NID = "+nid
            print(msg)
            log_msgs.append(msg)
        else:
            nid = None

        #esafile= "/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/FS_CV3/ESA_Int_products/Trace_SLIT_A_1600_V84600004001P0000000002101_39530_JLAB88_000001.fits"
        esafile, esafile_log_msgs = auxfunc.get_esafile(esa_files_path, raw_data_root_file, "FS", specifics, nid=nid)
        for msg in esafile_log_msgs:
            log_msgs.append(msg)

        # skip the test if the esafile was not found
        if esafile == "ESA file not found":
            msg1 = " * compare_wcs_fs.py is exiting because the corresponding ESA file was not found."
            msg2 = "   -> The WCS test is now set to skip and no plots will be generated. "
            print(msg1)
            print(msg2)
            log_msgs.append(msg1)
            log_msgs.append(msg2)
            FINAL_TEST_RESULT = "skip"
            return FINAL_TEST_RESULT, log_msgs

        """
        # comparison of filter, grating and x and y tilt
        gwa_xtil = fits.getval(infile_name, "gwa_xtil", 0)
        gwa_ytil = fits.getval(infile_name, "gwa_ytil", 0)
        esagrat = fits.getval(esafile, "GWA_POS", 0)
        esafilt = fits.getval(esafile, "FWA_POS", 0)
        esa_xtil = fits.getval(esafile, "gwa_xtil", 0)
        esa_ytil = fits.getval(esafile, "gwa_ytil", 0)
        print("pipeline: ")
        print("grating=", grat, " Filter=", filt, " gwa_xtil=", gwa_xtil, " gwa_ytil=", gwa_ytil)
        print("ESA:")
        print("grating=", esagrat, " Filter=", esafilt, " gwa_xtil=", esa_xtil, " gwa_ytil=", esa_ytil)
        """

        # Open the trace in the esafile
        msg = "Using this ESA file: \n"+esafile
        print(msg)
        log_msgs.append(msg)
        with fits.open(esafile) as esahdulist:
            print("* ESA file contents ")
            esahdulist.info()
            esa_slit_id = map_slit_names[esahdulist[0].header['SLITID']]
            # first check is esa_slit == to pipe_slit?
            if pipeslit == esa_slit_id:
                msg = "\n -> Same slit found for pipeline and ESA data: "+pipeslit+"\n"
                print(msg)
                log_msgs.append(msg)
            else:
                msg = "\n -> Missmatch of slits for pipeline and ESA data: "+pipeslit, esa_slit_id+"\n"
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
                    except:
                        KeyError
                        msg = "Skipping tests for V2 and V3 because ESA file does not contain corresponding extensions."
                        print(msg)
                        log_msgs.append(msg)
                except:
                    KeyError
                    msg = "This file does not contain data for detector NRS1. Skipping test for this slit."
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
                    except:
                        KeyError
                        msg = "Skipping tests for V2 and V3 because ESA file does not contain corresponding extensions."
                        print(msg)
                        log_msgs.append(msg)
                except:
                    KeyError
                    msg1 = "\n * compare_wcs_fs.py is exiting because there are no extensions that match detector NRS2 in the ESA file."
                    msg2 = "   -> The WCS test is now set to skip and no plots will be generated. \n"
                    print(msg1)
                    print(msg2)
                    log_msgs.append(msg1)
                    log_msgs.append(msg2)
                    FINAL_TEST_RESULT = "skip"
                    return FINAL_TEST_RESULT, log_msgs


        # get the WCS object for this particular slit
        wcs_slit = nirspec.nrs_wcs_set_input(img, pipeslit)

        # if we want to print all available transforms, uncomment line below
        #print(wcs_slit)

        # The WCS object attribute bounding_box shows all valid inputs, i.e. the actual area of the data according
        # to the slit. Inputs outside of the bounding_box return NaN values.
        #bbox = wcs_slit.bounding_box
        #print('bounding_box: ', bbox)

        # In different observing modes the WCS may have different coordinate frames. To see available frames
        # uncomment line below.
        #print("Avalable frames: ", wcs_slit.available_frames)

        if debug:
            # To get specific pixel values use following syntax:
            det2slit = wcs_slit.get_transform('detector', 'slit_frame')
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
            esax = 2049-esax
            esay = 2049-esay
            msg = "Flipped ESA data for detector NRS2 comparison with pipeline."
            print(msg)
            log_msgs.append(msg)

        # check if subarray is not FULL FRAME
        subarray = fits.getval(infile_name, "SUBARRAY", 0)

        if "FULL" not in subarray:
            # In subarray coordinates
            # subtract xstart and ystart values in order to get subarray coords instead of full frame
            # wcs_slit.x(y)start are 1-based, turn them to 0-based for extraction
            xstart, ystart = img.meta.subarray.xstart, img.meta.subarray.ystart
            esay = esay - (ystart - 1)
            esax = esax - (xstart - 1)

            #print("img.meta.subarray._instance = ", img.meta.subarray._instance)

            bounding_box = False

            # In full frame coordinates
            #pipey, pipex = np.mgrid[:esa_wave.shape[0], : esa_wave.shape[1]]
            #esax, esay = pyw.all_pix2world(pipex, pipey, 0)
            #sca2world = wcs_slit.get_transform('sca', 'world')
            #pra, pdec, pwave = sca2world(esax - 1, esay - 1)
        else:
            bounding_box = True

        # Compute pipeline RA, DEC, and lambda
        pra, pdec, pwave = wcs_slit(esax - 1, esay - 1, with_bounding_box=bounding_box)# => RETURNS: RA, DEC, LAMBDA
        pwave *= 10 ** -6    # (lam *= 10**-6 to convert to microns)

        """
        # checking that both ESA and pipeline have non NAN values
        no_nansp, no_nanse = [], []
        for vp, ve in zip(pwave, esa_wave):
            if np.nan not in vp:
                print(vp, ve)
                no_nansp.append(vp)
                no_nanse.append(ve)
        print(len(no_nansp), len(no_nanse))
        """

        # calculate and print statistics for slit-y and x relative differences
        tested_quantity = "Wavelength Difference"
        rel_diff_pwave_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_wave, pwave, tested_quantity)
        rel_diff_pwave_img, notnan_rel_diff_pwave, notnan_rel_diff_pwave_stats, print_stats = rel_diff_pwave_data
        for msg in print_stats:
            log_msgs.append(msg)
        test_result = auxfunc.does_median_pass_tes(notnan_rel_diff_pwave_stats[1], threshold_diff)
        msg = "\n * Result of the test for "+tested_quantity+":  "+test_result+"\n"
        print(msg)
        log_msgs.append(msg)
        total_test_result[pipeslit] = {tested_quantity : test_result}

        # get the transforms for pipeline slit-y
        det2slit = wcs_slit.get_transform('detector', 'slit_frame')
        slitx, slity, _ = det2slit(esax-1, esay-1, with_bounding_box=bounding_box)
        tested_quantity = "Slit-Y Difference"
        # calculate and print statistics for slit-y and x relative differences
        rel_diff_pslity_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_slity, slity, tested_quantity)
        # calculate and print statistics for slit-y and x absolute differences
        #rel_diff_pslity_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_slity, slity, tested_quantity, abs=True)
        rel_diff_pslity_img, notnan_rel_diff_pslity, notnan_rel_diff_pslity_stats, print_stats = rel_diff_pslity_data
        for msg in print_stats:
            log_msgs.append(msg)
        test_result = auxfunc.does_median_pass_tes(notnan_rel_diff_pslity_stats[1], threshold_diff)
        msg = "\n * Result of the test for "+tested_quantity+":  "+test_result+"\n"
        print(msg)
        log_msgs.append(msg)
        total_test_result[pipeslit] = {tested_quantity : test_result}

        # do the same for MSA x, y and V2, V3
        detector2msa = wcs_slit.get_transform("detector", "msa_frame")
        pmsax, pmsay, _ = detector2msa(esax-1, esay-1, with_bounding_box=bounding_box)   # => RETURNS: msaX, msaY, LAMBDA (lam *= 10**-6 to convert to microns)
        # MSA-x
        tested_quantity = "MSA_X Difference"
        reldiffpmsax_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_msax, pmsax, tested_quantity)
        reldiffpmsax_img, notnan_reldiffpmsax, notnan_reldiffpmsax_stats, print_stats = reldiffpmsax_data
        for msg in print_stats:
            log_msgs.append(msg)
        test_result = auxfunc.does_median_pass_tes(notnan_reldiffpmsax_stats[1], threshold_diff)
        msg = "\n * Result of the test for "+tested_quantity+":  "+test_result+"\n"
        print(msg)
        log_msgs.append(msg)
        total_test_result[pipeslit] = {tested_quantity : test_result}
        # MSA-y
        tested_quantity = "MSA_Y Difference"
        reldiffpmsay_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_msay, pmsay, tested_quantity)
        reldiffpmsay_img, notnan_reldiffpmsay, notnan_reldiffpmsay_stats, print_stats = reldiffpmsay_data
        for msg in print_stats:
            log_msgs.append(msg)
        test_result = auxfunc.does_median_pass_tes(notnan_reldiffpmsay_stats[1], threshold_diff)
        msg = "\n * Result of the test for "+tested_quantity+":  "+test_result+"\n"
        print(msg)
        log_msgs.append(msg)
        total_test_result[pipeslit] = {tested_quantity : test_result}

        # V2 and V3
        if not skipv2v3test:
            detector2v2v3 = wcs_slit.get_transform("detector", "v2v3")
            pv2, pv3, _ = detector2v2v3(esax-1, esay-1, with_bounding_box=bounding_box)   # => RETURNS: v2, v3, LAMBDA (lam *= 10**-6 to convert to microns)
            tested_quantity = "V2 difference"
            # converting to degrees to compare with ESA, pipeline is in arcsec
            reldiffpv2_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_v2v3x, pv2, tested_quantity)
            if reldiffpv2_data[-2][0] > 0.0:
                print("\nConverting pipeline results to degrees to compare with ESA")
                pv2 = pv2/3600.
                reldiffpv2_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_v2v3x, pv2, tested_quantity)
            reldiffpv2_img, notnan_reldiffpv2, notnan_reldiffpv2_stats, print_stats = reldiffpv2_data
            for msg in print_stats:
                log_msgs.append(msg)
            test_result = auxfunc.does_median_pass_tes(notnan_reldiffpv2_stats[1], threshold_diff)
            msg = "\n * Result of the test for " + tested_quantity + ":  " + test_result + "\n"
            print(msg)
            log_msgs.append(msg)
            total_test_result[pipeslit] = {tested_quantity : test_result}

            tested_quantity = "V3 difference"
            # converting to degrees to compare with ESA
            reldiffpv3_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_v2v3y, pv3, tested_quantity)
            if reldiffpv3_data[-2][0] > 0.0:
                print("\nConverting pipeline results to degrees to compare with ESA")
                pv3 = pv3/3600.
                reldiffpv3_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_v2v3y, pv3, tested_quantity)
            reldiffpv3_img, notnan_reldiffpv3, notnan_reldiffpv3_stats, print_stats = reldiffpv3_data
            for msg in print_stats:
                log_msgs.append(msg)
            test_result = auxfunc.does_median_pass_tes(notnan_reldiffpv3_stats[1], threshold_diff)
            msg = "\n * Result of the test for "+tested_quantity+":  "+test_result+"\n"
            print(msg)
            log_msgs.append(msg)
            total_test_result[pipeslit] = {tested_quantity : test_result}

        # PLOTS
        if show_figs or save_figs:
            # set the common variables
            basenameinfile_name = os.path.basename(infile_name)
            main_title = filt+"   "+grat+"   SLIT="+pipeslit+"\n"
            bins = 15   # binning for the histograms, if None the function will automatically calculate number
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
                plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_"+det+"_rel_wave_diffs.jpg")
                auxfunc.plt_two_2Dimgandhist(rel_diff_pwave_img, notnan_rel_diff_pwave, info_img, info_hist,
                                             plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

            # Slit-y
            title = main_title+r"Relative slit position = $\Delta$slit_y"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$slit_y = (slit_y$_{pipe}$ - slit_y$_{ESA}$)/slit_y$_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_rel_diff_pslity_stats]
            if notnan_rel_diff_pslity_stats[1] is np.nan:
                msg = "Unable to create plot of relative slit position."
                print(msg)
                log_msgs.append(msg)
            else:
                plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_"+det+"_rel_slitY_diffs.jpg")
                auxfunc.plt_two_2Dimgandhist(rel_diff_pslity_img, notnan_rel_diff_pslity, info_img, info_hist,
                                             plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

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
                plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_"+det+"_rel_MSAx_diffs.jpg")
                auxfunc.plt_two_2Dimgandhist(reldiffpmsax_img, notnan_reldiffpmsax, info_img, info_hist,
                                             plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

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
                plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_"+det+"_rel_MSAy_diffs.jpg")
                auxfunc.plt_two_2Dimgandhist(reldiffpmsay_img, notnan_reldiffpmsay, info_img, info_hist,
                                             plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

            if not skipv2v3test:
                # V2
                title = main_title+r"Relative V2 Difference = $\Delta$V2"+"\n"
                info_img = [title, "x (pixels)", "y (pixels)"]
                xlabel, ylabel = r"Relative $\Delta$V2 = (V2$_{pipe}$ - V2$_{ESA}$)/V2$_{ESA}$", "N"
                info_hist = [xlabel, ylabel, bins, notnan_reldiffpv2_stats]
                if notnan_reldiffpv2_stats[1] is np.nan:
                    msg = "Unable to create plot of relative V2 difference."
                    print(msg)
                    log_msgs.append(msg)
                else:
                    plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_"+det+"_rel_V2_diffs.jpg")
                    auxfunc.plt_two_2Dimgandhist(reldiffpv2_img, notnan_reldiffpv2_stats, info_img, info_hist,
                                                 plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

                # V3
                title = main_title+r"Relative V3 Difference = $\Delta$V3"+"\n"
                info_img = [title, "x (pixels)", "y (pixels)"]
                xlabel, ylabel = r"Relative $\Delta$V3 = (V3$_{pipe}$ - V3$_{ESA}$)/V3$_{ESA}$", "N"
                info_hist = [xlabel, ylabel, bins, notnan_reldiffpv3_stats]
                if notnan_reldiffpv3_stats[1] is np.nan:
                    msg = "Unable to create plot of relative V3 difference."
                    print(msg)
                    log_msgs.append(msg)
                else:
                    plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_"+det+"_rel_V3_diffs.jpg")
                    auxfunc.plt_two_2Dimgandhist(reldiffpv3_img, notnan_reldiffpv3, info_img, info_hist,
                                                 plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)


        else:
            msg = "NO plots were made because show_figs and save_figs were both set to False. \n"
            print(msg)
            log_msgs.append(msg)


    # If all tests passed then pytest will be marked as PASSED, else it will be FAILED
    FINAL_TEST_RESULT = "FAILED"
    for sl, testdir in total_test_result.items():
        for t, tr in testdir.items():
            if tr == "FAILED":
                FINAL_TEST_RESULT = "FAILED"
                msg = "\n * The test of "+t+" for slit "+sl+" FAILED."
                print(msg)
                log_msgs.append(msg)
            else:
                FINAL_TEST_RESULT = "PASSED"
                msg = "\n * The test of "+t+" for slit "+sl+" PASSED."
                print(msg)
                log_msgs.append(msg)

    if FINAL_TEST_RESULT == "PASSED":
        msg = "\n *** Final result for assign_wcs test will be reported as PASSED *** \n"
    else:
        msg = "\n *** Final result for assign_wcs test will be reported as FAILED *** \n"
    print(msg)
    log_msgs.append(msg)


    return FINAL_TEST_RESULT, log_msgs





if __name__ == '__main__':

    # This is a simple test of the code
    pipeline_path = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline"

    # input parameters that the script expects
    #data_dir = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline/testing_data/FS_FULL_FRAME/G235H_opaque/491_results"
    #data_dir = pipeline_path+"/testing_data/FS_FULL_FRAME/G140H_opaque"
    ##infile_name = data_dir+"/jwdata0010010_11010_0001_NRS1_assign_wcs.fits"   # for G140H
    #infile_name = pipeline_path+"/testing_data/FS_ALLSLITS/G140H_opaque/gain_scale_NRS1_assign_wcs.fits"
    #infile_name = pipeline_path+"/testing_data/FS_ALLSLITS/G140M_F070LP/final_output_caldet1_NRS1_assign_wcs.fits"
    #infile_name = pipeline_path+"/testing_data/FS_ALLSLITS/G235H_F170LP/ALLSLITS_g235H_gain_scale_NRS1_assignwcsstep.fits"
    infile_name = pipeline_path+"/testing_data/FS_ALLSLITS/G235H_F170LP/final_output_caldet1_NRS1_assign_wcs.fits"
    #infile_name = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline/testing_data/FS_FULL_FRAME/G140H_opaque/jwdata0010010_11010_0001_NRS1_assign_wcs.fits"
    #esa_files_path=pipeline_path+"/build7/test_data/ESA_intermediary_products/RegressionTestData_CV3_March2017_FixedSlit/"
    #esa_files_path = "/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/FS_CV3_cutouts/ESA_Int_products"
    esa_files_path = "/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/FS_CV3/ESA_Int_products"

    # print pipeline version
    import jwst
    print("\n  ** using pipeline version: ", jwst.__version__, "** \n")

    # Run the principal function of the script
    result, log_msgs = compare_wcs(infile_name, esa_files_path=esa_files_path, show_figs=False, save_figs=True,
                                    threshold_diff=1.0e-7, debug=False)




