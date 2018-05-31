import numpy as np
import os
from astropy.io import fits
from astropy import wcs

from jwst.assign_wcs import nirspec
from jwst import datamodels
from . import auxiliary_functions as auxfunc


"""
This script compares pipeline WCS info with ESA results for FIXED SLIT.

"""


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "2.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Completely changed script to use the datamodel instead of the compute_world_coordinates
#                         script, and added new routines for plot making and statistics calculations.


def compare_wcs(infile_name, esa_files_path=None, show_figs=True, save_figs=False, threshold_diff=1.0e-7, debug=False):
    """
    This function does the WCS comparison from the world coordinates calculated using the
    compute_world_coordinates.py script with the ESA files. The function calls that script.

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

    """

    # get grating and filter info from the rate file header
    det = fits.getval(infile_name, "DETECTOR", 0)
    print('infile_name=', infile_name)
    lamp = fits.getval(infile_name, "LAMP", 0)
    grat = fits.getval(infile_name, "GRATING", 0)
    filt = fits.getval(infile_name, "FILTER", 0)
    print ("from assign_wcs file  -->     Detector:", det, "   Grating:", grat, "   Filter:", filt, "   Lamp:", lamp)

    # get the datamodel from the assign_wcs output file
    img = datamodels.ImageModel(infile_name)

    # loop over the slits
    sltname_list = ["S200A1", "S200A2", "S400A1", "S1600A1", "S200B1"]
    sci_ext_list = auxfunc.get_sci_extensions(infile_name)
    print ('sci_ext_list=', sci_ext_list, '\n')

    # mapping the ESA slit names to pipeline names
    map_slit_names = {'SLIT_A_1600' : 'S1600A1',
                      'SLIT_A_200_1': 'S200A1',
                      'SLIT_A_200_2': 'S200A2',
                      'SLIT_A_400':   'S400A1',
                      'SLIT_B_200':   'S200B1',
                      }

    # list to determine if pytest is passed or not
    total_test_result = []

    # loop over the slits
    if det != "NRS2":
        sltname_list.pop(len(sltname_list)-1)

    # check if data is BOTS
    if fits.getval(infile_name, "EXP_TYPE", 0) == "NRS_BRIGHTOBJ":
        sltname_list = ["S1600A1"]

    for pipeslit in sltname_list:
        print ("\nWorking with slit: ", pipeslit)

        # Get the ESA trace
        #raw_data_root_file = "NRSV84600010001P0000000002101_4_491_SE_2016-01-17T17h34m08.fits"  # for testing with G140M FULLFRAME
        #raw_data_root_file = "NRSSMOS-MOD-G1H-02-5344031756_1_491_SE_2015-12-10T03h25m56.fits"  # for testing with G140H FULLFRAME
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
            print("Using NID = ", nid)
        else:
            nid = None
        esafile = auxfunc.get_esafile(esa_files_path, raw_data_root_file, "FS", specifics, nid=nid)

        # skip the test if the esafile was not found
        if esafile == "ESA file not found":
            print(" * compare_wcs_fs.py is exiting because the corresponding ESA file was not found.")
            print("   -> The WCS test is now set to skip and no plots will be generated. ")
            FINAL_TEST_RESULT = "skip"
            return FINAL_TEST_RESULT

        # Open the trace in the esafile
        print ("Using this ESA file: \n", esafile)
        with fits.open(esafile) as esahdulist:
            #print ("* ESA file contents ")
            #esahdulist.info()
            esa_slit_id = map_slit_names[esahdulist[0].header['SLITID']]
            # first check is esa_slit == to pipe_slit?
            if pipeslit == esa_slit_id:
                print("\n -> Same slit found for pipeline and ESA data: ", pipeslit, "\n")
            else:
                print("\n -> Missmatch of slits for pipeline and ESA data: ", pipeslit, esa_slit_id, "\n")

            # Assign variables according to detector
            if det == "NRS1":
                esa_flux = fits.getdata(esafile, "DATA1")
                esa_wave = fits.getdata(esafile, "LAMBDA1")
                esa_slity = fits.getdata(esafile, "SLITY1")
                esa_msax = fits.getdata(esafile, "MSAX1")
                esa_msay = fits.getdata(esafile, "MSAY1")
                pyw = wcs.WCS(esahdulist['LAMBDA1'].header)
            if det == "NRS2":
                try:
                    esa_flux = fits.getdata(esafile, "DATA2")
                    esa_wave = fits.getdata(esafile, "LAMBDA2")
                    esa_slity = fits.getdata(esafile, "SLITY2")
                    esa_msax = fits.getdata(esafile, "MSAX2")
                    esa_msay = fits.getdata(esafile, "MSAY2")
                    pyw = wcs.WCS(esahdulist['LAMBDA2'].header)
                except:
                    IndexError
                    print("\n * compare_wcs_fs.py is exiting because there are no extensions that match detector NRS2 in the ESA file.")
                    print("   -> The WCS test is now set to skip and no plots will be generated. \n")
                    FINAL_TEST_RESULT = "skip"
                    return FINAL_TEST_RESULT


        # get the WCS object for this particular slit
        wcs_slit = nirspec.nrs_wcs_set_input(img, pipeslit)

        # if we want to print all available transforms, uncomment line below
        #print(wcs_slit)

        # The WCS object attribute bounding_box shows all valid inputs, i.e. the actual area of the data according
        # to the slit. Inputs outside of the bounding_box return NaN values.
        #bbox = wcs_slit.bounding_box
        #print('wcs_slit.bounding_box: ', wcs_slit.bounding_box)

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

        # Compute pipeline RA, DEC, and lambda
        pra, pdec, pwave = wcs_slit(esax-1, esay-1)   # => RETURNS: RA, DEC, LAMBDA (lam *= 10**-6 to convert to microns)
        pwave *= 10**-6
        # calculate and print statistics for slit-y and x relative differences
        tested_quantity = "Wavelength Difference"
        rel_diff_pwave_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_wave, pwave, tested_quantity)
        rel_diff_pwave_img, notnan_rel_diff_pwave, notnan_rel_diff_pwave_stats = rel_diff_pwave_data
        test_result = auxfunc.does_median_pass_tes(tested_quantity, notnan_rel_diff_pwave_stats[1], threshold_diff)
        total_test_result.append(test_result)

        # get the transforms for pipeline slit-y
        det2slit = wcs_slit.get_transform('detector', 'slit_frame')
        slitx, slity, _ = det2slit(esax-1, esay-1)
        # calculate and print statistics for slit-y and x relative differences
        tested_quantity = "Slit-Y Difference"
        rel_diff_pslity_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_slity, slity, tested_quantity)
        rel_diff_pslity_img, notnan_rel_diff_pslity, notnan_rel_diff_pslity_stats = rel_diff_pslity_data
        test_result = auxfunc.does_median_pass_tes(tested_quantity, notnan_rel_diff_pslity_stats[1], threshold_diff)
        total_test_result.append(test_result)

        # do the same for MSA x, y and V2, V3
        detector2msa = wcs_slit.get_transform("detector", "msa_frame")
        pmsax, pmsay, _ = detector2msa(esax-1, esay-1)   # => RETURNS: msaX, msaY, LAMBDA (lam *= 10**-6 to convert to microns)
        # MSA-x
        tested_quantity = "MSA_X Difference"
        reldiffpmsax_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_msax, pmsax, tested_quantity)
        reldiffpmsax_img, notnan_reldiffpmsax, notnan_reldiffpmsax_stats = reldiffpmsax_data
        test_result = auxfunc.does_median_pass_tes(tested_quantity, notnan_reldiffpmsax_stats[1], threshold_diff)
        total_test_result.append(test_result)
        # MSA-y
        tested_quantity = "MSA_Y Difference"
        reldiffpmsay_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, esa_msay, pmsay, tested_quantity)
        reldiffpmsay_img, notnan_reldiffpmsay, notnan_reldiffpmsay_stats = reldiffpmsay_data
        test_result = auxfunc.does_median_pass_tes(tested_quantity, notnan_reldiffpmsay_stats[1], threshold_diff)
        total_test_result.append(test_result)

        # V2 and V3
        detector2v2v3 = wcs_slit.get_transform("detector", "v2v3")
        pv2, pv3, _ = detector2v2v3(esax-1, esay-1)   # => RETURNS: v2, v3, LAMBDA (lam *= 10**-6 to convert to microns)
        # get the corresponding ESA arr of interest
        #esav2, esav3 = ct.coords_transf("forward", det, filter_input, avgx3, avgy3, tilt, arcsecs, debug)
        #restricted_ev2, restricted_ev3 = [], []
        tested_quantity = "V2 difference"
        #reldiffpv2_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, restricted_esav2, pv2, tested_quantity)
        #reldiffpv2_img, notnan_reldiffpv2, notnan_reldiffpv2_stats = reldiffpv2_data
        #test_result = auxfunc.does_median_pass_tes(tested_quantity, notnan_reldiffpv2_stats[1], threshold_diff)
        #total_test_result.append(test_result)
        tested_quantity = "V3 difference"
        #reldiffpv3_data = auxfunc.get_reldiffarr_and_stats(threshold_diff, esa_slity, restricted_esav3, pv3, tested_quantity)
        #reldiffpv3_img, notnan_reldiffpv3, notnan_reldiffpv3_stats = reldiffpv3_data
        #test_result = auxfunc.does_median_pass_tes(tested_quantity, notnan_reldiffpv3_stats[1], threshold_diff)
        #total_test_result.append(test_result)

        # PLOTS
        if show_figs or save_figs:
            # set the common variables
            basenameinfile_name = os.path.basename(infile_name)
            main_title = filt+"   "+grat+"   SLIT="+pipeslit+"\n"
            bins = 15   # binning for the histograms
            #             lolim_x, uplim_x, lolim_y, uplim_y
            plt_origin = None

            # Wavelength
            title = main_title+r"Relative wavelength difference = $\Delta \lambda$"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta \lambda$ = ($\lambda_{pipe} - \lambda_{ESA}) / \lambda_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_rel_diff_pwave_stats]
            plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_rel_wave_diffs.jpg")
            auxfunc.plt_two_2Dimgandhist(rel_diff_pwave_img, notnan_rel_diff_pwave, info_img, info_hist,
                                         plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

            # Slit-y
            title = main_title+r"Relative slit position = $\Delta$slit_y"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$slit_y = (slit_y$_{pipe}$ - slit_y$_{ESA}$)/slit_y$_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_rel_diff_pslity_stats]
            plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_rel_slitY_diffs.jpg")
            auxfunc.plt_two_2Dimgandhist(rel_diff_pslity_img, notnan_rel_diff_pslity, info_img, info_hist,
                                         plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

            # MSA-x
            title = main_title+r"Relative MSA-x Difference = $\Delta$MSA_x"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$MSA_x = (MSA_x$_{pipe}$ - MSA_x$_{ESA}$)/MSA_x$_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_reldiffpmsax_stats]
            plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_rel_MSAx_diffs.jpg")
            auxfunc.plt_two_2Dimgandhist(reldiffpmsax_img, notnan_reldiffpmsax, info_img, info_hist,
                                         plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

            # MSA-y
            title = main_title+r"Relative MSA-y Difference = $\Delta$MSA_y"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$MSA_y = (MSA_y$_{pipe}$ - MSA_y$_{ESA}$)/MSA_y$_{ESA}$", "N"
            info_hist = [xlabel, ylabel, bins, notnan_reldiffpmsay_stats]
            plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_rel_MSAy_diffs.jpg")
            auxfunc.plt_two_2Dimgandhist(reldiffpmsay_img, notnan_reldiffpmsay, info_img, info_hist,
                                         plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)
            """
            # V2
            title = main_title+r"Relative V2 Difference = $\Delta$V2"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$V2 = (V2$_{pipe}$ - V2$_{ESA}$)/V2$_{ESA}$", "N"
            hist_data = notnan_reldiffpv2
            info_hist = [xlabel, ylabel, bins, notnan_reldiffpv2_stats]
            plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_rel_V2_diffs.jpg")
            auxfunc.plt_two_2Dimgandhist(reldiffpv2_img, hist_data, info_img, info_hist,
                                         plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

            # V3
            title = main_title+r"Relative V3 Difference = $\Delta$V3"+"\n"
            info_img = [title, "x (pixels)", "y (pixels)"]
            xlabel, ylabel = r"Relative $\Delta$V3 = (V3$_{pipe}$ - V3$_{ESA}$)/V3$_{ESA}$", "N"
            hist_data = notnan_reldiffpv3
            info_hist = [xlabel, ylabel, bins, notnan_reldiffpmsay_stats]
            plt_name = infile_name.replace(basenameinfile_name, pipeslit+"_rel_V3_diffs.jpg")
            auxfunc.plt_two_2Dimgandhist(reldiffpv3_img, hist_data, info_img, info_hist,
                                         plt_name=plt_name, plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)
            """

        else:
            print ("NO plots were made because show_figs and save_figs were both set to False. \n")


    # If all tests passed then pytest will be marked as PASSED, else it will be FAILED
    FINAL_TEST_RESULT = "PASSED"
    slittests = ["Wavelength", "Slit-y", "MSA-X", "MSA-Y"]#, "V2", "V3"]
    for i, tr in enumerate(total_test_result):
        if tr == "FAILED":
            FINAL_TEST_RESULT = "FAILED"
            # 4 tests per slit, and there are 4 slits, total of 16 results for det=NRS1
            if i < 4:#6:
                failedslit = sltname_list[0]
                failedtest = slittests[i]
            elif (i>3) and (i<8):#(i>5) and (i<12):
                failedslit = sltname_list[1]
                j = i - 4
                failedtest = slittests[j]
            elif (i>7) and (i<12):#(i>11) and (i<18):
                failedslit = sltname_list[2]
                j = i - 8
                failedtest = slittests[j]
            elif (i>11) and (i<16):#(i>17) and (i<24):
                failedslit = sltname_list[3]
                j = i - 12
                failedtest = slittests[j]
            if det == "NRS2":
                # 4 tests per slit, and there are 5 slits, total of 20 results
                if i > 15:#i > 23:
                    failedslit = sltname_list[4]
                    j = i - 16
                    failedtest = slittests[j]
            print("\n * Test of", failedtest, "FAILED for slit", failedslit)

    if FINAL_TEST_RESULT == "PASSED":
        print("\n *** Final result for assign_wcs test will be reported as PASSED *** \n")
    else:
        print("\n *** Final result for assign_wcs test will be reported as FAILED *** \n")


    return FINAL_TEST_RESULT







if __name__ == '__main__':

    # This is a simple test of the code
    pipeline_path = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline"

    # input parameters that the script expects
    #data_dir = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline/build7.1/part2/FS_FULL_FRAME/G235H_opaque/491_results"
    #data_dir = pipeline_path+"/build7.1/part2/FS_FULL_FRAME/G140M_opaque/491_results"
    data_dir = pipeline_path+"/build7.1/part2/FS_ALLSLITS/G140M_F070LP/491_results"
    infile_name = data_dir+"/gain_scale_assign_wcs.fits"
    #esa_files_path=pipeline_path+"/build7/test_data/ESA_intermediary_products/RegressionTestData_CV3_March2017_FixedSlit/"
    #esa_files_path = "/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/FS_CV3_cutouts/ESA_Int_products"
    esa_files_path = "/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/FS_CV3/ESA_Int_products"

    # print pipeline version
    import jwst
    print("\n  ** using pipeline version: ", jwst.__version__, "** \n")

    # Run the principal function of the script
    result = compare_wcs(infile_name, esa_files_path=esa_files_path, show_figs=False, save_figs=True,
                         threshold_diff=1.0e-7, debug=False)




