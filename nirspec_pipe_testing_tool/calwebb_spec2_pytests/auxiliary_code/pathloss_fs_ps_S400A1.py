import time
import os
from astropy import wcs
import numpy as np
from astropy.io import fits
import scipy
import pdb

import jwst
from gwcs import wcstools
from jwst import datamodels
from . import auxiliary_functions as auxfunc

from astropy.visualization import (ImageNormalize, AsinhStretch)
import matplotlib
import matplotlib.pyplot as plt

"""
This script tests the FS pipeline pathloss step output for a Point Source.
"""

# HEADER
__author__ = "T King"
__version__ = "1.2"

# HISTORY
# October 19, 2019 - Version 1.0: initial version started
# February 5, 2020 - Version 1.1: Tests pass
# February 26, 2020 - Version 1.2: Mostly pep8 compliant


def get_ps_uni_extensions(fits_file_name, is_point_source):
    """
    This functions obtains all the point source or uniform source
    extensions in the given file

    Args:
        fits_file_name: str, name of the fits file of interest
        is_point_source: boolean, true if point source; false if uniform source
    Returns:
        sci_list: list of the numbers of the science extensions
    """
    hdulist = fits.open(fits_file_name)
    ps_dict = {}
    uni_dict = {}
    s = 0
    for ext, hdu in enumerate(hdulist):
        if hdu.name == "PS":
            try:
                sltname = hdu.header["APERTURE"]
                ps_dict[sltname] = ext
            except KeyError:
                sltname = "Slit_"+repr(s+1)
                ps_dict[sltname] = ext
        if hdu.name == "UNI":
            try:
                sltname = hdu.header["APERTURE"]
                uni_dict[sltname] = ext
            except KeyError:
                sltname = "Slit_"+repr(s+1)
                uni_dict[sltname] = ext
    return ps_dict, uni_dict


def pathtest(step_input_filename, reffile, comparison_filename, writefile=True, show_figs=True,
             save_figs=False, plot_name=None, threshold_diff=1e-7,
             debug=False):
    """
    This function calculates the difference between the pipeline and the
    calculated pathloss values. The functions use the output of the
    compute_world_coordinates.py script.
    Args:
        step_input_filename: str, name of the output fits file from the
        source type step (with full path) 
        reffile: str, path to the pathloss FS reference fits file
        comparison_filename: str, path to comparison pipeline pathloss file
        writefile: boolean, if True writes the fits files of the calculated
        pathloss and difference image.
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots (the 3 plots can be saved or not
        independently with the function call)
        plot_name: string, desired name (if name is not given, the plot
        function will name the plot by default)
        threshold_diff: float, threshold difference between pipeline output
        and ESA file
        debug: boolean, if true a series of print statements will show
        on-screen
    Returns:
        - 1 plot, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to threshold.
        - log_msgs: list, all print statements are captured in this variable
    """

    log_msgs = []

    # start the timer
    pathtest_start_time = time.time()

    # get info from the rate file header
    det = fits.getval(step_input_filename, "DETECTOR", 0)
    msg = 'step_input_filename='+step_input_filename
    print(msg)
    log_msgs.append(msg)
    exptype = fits.getval(step_input_filename, "EXP_TYPE", 0)
    grat = fits.getval(step_input_filename, "GRATING", 0)
    filt = fits.getval(step_input_filename, "FILTER", 0)

    msg = "path_loss_file  -->  Grating:"+grat+"   Filter:"+filt+"   EXP_TYPE:"+exptype
    print(msg)
    log_msgs.append(msg)

    msg = "Using reference file: "+reffile
    print(msg)
    log_msgs.append(msg)

    is_point_source = True

    # get the datamodel from the assign_wcs output file
    extract2d_wcs_file = '/Users/tking/documents/king_FS_NRS1_V2/final_output_caldet1_NRS1_extract_2d.fits'
    model = datamodels.MultiSlitModel(extract2d_wcs_file)

    if writefile:
        # create the fits list to hold the calculated pathloss values for
        # each slit
        hdu0 = fits.PrimaryHDU()
        outfile = fits.HDUList()
        outfile.append(hdu0)

        # create the fits list to hold the image of pipeline-calculated
        # difference values
        hdu0 = fits.PrimaryHDU()
        compfile = fits.HDUList()
        compfile.append(hdu0)

    # list to determine if pytest is passed or not
    total_test_result = []

    # loop over the slits
    sltname_list = ["S200A1", "S200A2", "S400A1", "S1600A1"]
    msg = "Now looping through the slits. This may take a while... "
    print(msg)
    log_msgs.append(msg)
    if det == "NRS2":
        sltname_list.append("S200B1")

    # but check if data is BOTS
    if fits.getval(step_input_filename, "EXP_TYPE", 0) == "NRS_BRIGHTOBJ":
        sltname_list = ["S1600A1"]

    # get all the science extensions
    ps_uni_ext_list = get_ps_uni_extensions(reffile, is_point_source)

    # get files
    print("""Checking if files exist & obtaining datamodels.
          This takes a few minutes...""")
    if os.path.isfile(comparison_filename):
        if debug:
            print('Comparison file does exist.')
    else:
        result_msg = 'Comparison file does NOT exist. Skipping pathloss test.'
        print(result_msg)
        log_msgs.append(result_msg)
        result = 'skip'
        return result, result_msg, log_msgs

    # get the comparison data model
    pathloss_pipe = datamodels.open(comparison_filename)
    if debug:
        print('got comparison datamodel!')

    if os.path.isfile(step_input_filename):
        if debug:
            print('Input file does exist.')
    else:
        result_msg = 'Input file does NOT exist. Skipping pathloss test.'
        log_msgs.append(result_msg)
        result = 'skip'
        return result, result_msg, log_msgs
    # get the input data model
    pl = datamodels.open(step_input_filename)
    if debug:
        print('got input datamodel!')

    # loop through the wavelengths
    msg = "Looping through the wavelengths... "  # SLITS?
    print(msg)
    log_msgs.append(msg)

    slit_val = 0
    for slit, pipe_slit in zip(pl.slits, pathloss_pipe.slits):
        slit_val = slit_val+1

        # only testing S400A1 in this script
        if slit_val == 1:
            continue
        if slit_val == 2:
            continue
        if slit_val == 4:
            continue

        slit_id = pipe_slit.name
        print('Working with slitlet ', slit_id)
        if slit.name == slit_id:
            msg = """Slitlet name in fits file previous to pathloss
            and in pathloss output file are the same.
            """
            log_msgs.append(msg)
            print(msg)
        else:
            msg = """* Missmatch of slitlet names in fits file previous to
            pathloss and in pathloss output file. Skipping test.
            """
            result = 'skip'
            log_msgs.append(msg)
            return result, msg, log_msgs

        mode = "FS"

        if debug:
            print("grat = ", grat)

        continue_pl_test = False
        if fits.getval(step_input_filename, "EXP_TYPE", 0) == "NRS_BRIGHTOBJ":
            slit = model
            continue_pl_test = True
        else:
            for slit_in_MultiSlitModel in pl.slits:
                if slit_in_MultiSlitModel.name == slit_id:
                    slit = slit_in_MultiSlitModel
                    continue_pl_test = True
                    break

        if not continue_pl_test:
            continue
        else:
            try:
                if is_point_source is True:
                    ext = ps_uni_ext_list[0][slit_id]
                    print("Retrieved point source extension")
                elif is_point_source is False:
                    ext = ps_uni_ext_list[1][slit_id]
                    print("WARNING: Retrieved uniform source extension")
            except KeyError:
                ext = sltname_list.index(slit_id)
                print("Unable to retrieve extension.")

        wcs_obj = slit.meta.wcs

        # get the wavelength
        x, y = wcstools.grid_from_bounding_box(wcs_obj.bounding_box, step=(1, 1), center=True)
        ra, dec, wave = wcs_obj(x, y)   # wave is in microns
        wave_sci = wave * 10**(-6)  # microns --> meters

        # get positions of source in file:
        slit_x = slit.source_xpos
        slit_y = slit.source_ypos
        if debug:
            print("slit_x, slit_y", slit_x, slit_y)
        
        reffile = '/Users/tking/jwst_nirspec_pathloss_0001.fits'
        if slit_val == 3:
            if is_point_source:
                ext = 1
            else:
                ext = 3
            reffile = "/Users/tking/Desktop/jwst-nirspec-a400.plrf/jwst-nirspec-a400.plrf.fits"


        plcor_ref_ext = fits.getdata(reffile, ext)
        if debug:
            print("ext:", ext)
        hdul = fits.open(reffile)
        plcor_ref = hdul[1].data
        w = wcs.WCS(hdul[1].header)
        
        # make cube
        w1, y1, x1 = np.mgrid[:plcor_ref.shape[0], : plcor_ref.shape[1],
                              :plcor_ref.shape[2]]
        slitx_ref, slity_ref, wave_ref = w.all_pix2world(x1, y1, w1, 0)

        previous_sci = slit.data

        # print("Test for srctype == extract 2d:", previous_sci == extract_sci)
        # print("any differences:", np.where(previous_sci != extract_sci))
        pipe_correction = pipe_slit.pathloss

        # SET UP MANUALLY TO TEST CORRECTION AT NONZERO POINT
        # slit_x = 0.2
        # slit_y = 0.2
        # if debug:
        #     print("""WARNING: Using manually set slit_x and slit_y!
        # The pipeline correction will not use manually set values and
        # thus the residuals will change
        # """)

        correction_array = np.array([])
        lambda_array = np.array([])

        wave_sci_flat = wave_sci.reshape(wave_sci.size)
        wave_ref_flat = wave_ref.reshape(wave_ref.size)

        ref_xy = np.column_stack((slitx_ref.reshape(slitx_ref.size), slity_ref.reshape(slitx_ref.size)))

        
        # loop through slices in lambda from reference file
        shape = 0
        for lambda_val in wave_ref_flat:
            # loop through every lambda value
            # flattened so that looping works smoothly
            shape = shape + 1
            index = np.where(wave_ref[:, 0, 0] == lambda_val)
            # index of closest lambda value in reffile to given sci lambda
            #   took index of only the first slice of wave_ref because
            #   the others were repetitive & we got extra indices
            # take slice where lambda=index:
            plcor_slice = plcor_ref_ext[index[0][0]].reshape(plcor_ref_ext[index[0][0]].size)
            # do 2d interpolation to get a single correction factor for each slice
            corr_val = scipy.interpolate.griddata(ref_xy[:plcor_slice.size], plcor_slice,
                                                  np.asarray([slit_x, slit_y]),
                                                  method='linear')
            # append values from loop to create a vector of correction factors
            correction_array = np.append(correction_array, corr_val[0])
            # map to array with corresponding lambda
            lambda_array = np.append(lambda_array, lambda_val)

        # get correction value for each pixel
        corr_vals = np.interp(wave_sci_flat, lambda_array, correction_array)
        corr_vals = corr_vals.reshape(wave_sci.shape)
        corrected_array = previous_sci/corr_vals

        # set up generals for all the plots
        font = {'weight': 'normal',
                'size': 7}
        matplotlib.rc('font', **font)

        # Plots:
        # my correction values
        fig = plt.figure()
        plt.subplot(221)
        norm = ImageNormalize(corr_vals)
        plt.imshow(corr_vals, norm  = norm, aspect = 10.0, origin='lower', cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Calculated Correction')
        plt.colorbar()
        # pipe correction
        plt.subplot(222)
        norm = ImageNormalize(pipe_correction)
        plt.imshow(pipe_correction, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.title("Pipeline Correction")
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.colorbar()
        # residuals
        corr_residuals = pipe_correction - corr_vals
        plt.subplot(223)
        norm = ImageNormalize(corr_residuals)
        plt.imshow(corr_residuals, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Correction residuals')
        plt.colorbar()
        # my science data after
        plt.subplot(224)
        norm = ImageNormalize(corrected_array)
        plt.imshow(corrected_array, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.title('Calculated after pathloss data')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.colorbar()
        fig.suptitle("FS PS Pathloss Correction Test for slit " + str(slit_id))

        if save_figs:
            plt_name = step_input_filename + "_Pathloss_test_slitlet_"+str(mode) + "_" + str(slit_id) + "_" + str(slit_val) + ".jpg"
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_figs:
            plt.show()
        plt.clf()

        # create fits file to hold the calculated pathloss for each slit
        if writefile:
            msg = "Saving the fits files with the calculated pathloss for each slit..."
            print(msg)
            log_msgs.append(msg)

            # this is the file to hold the image of pipeline-calculated difference values
            outfile_ext = fits.ImageHDU(corr_vals, name=slit_id)
            outfile.append(outfile_ext)

            # this is the file to hold the image of pipeline-calculated difference values
            compfile_ext = fits.ImageHDU(corr_residuals, name=slit_id)
            compfile.append(compfile_ext)

        if debug:
            if corr_residuals[~np.isnan(corr_residuals)].size == 0:
                msg1 = " * Unable to calculate statistics because difference array has all values as NaN. Test will be set to FAILED."
                print(msg1)
                log_msgs.append(msg1)
                test_result = "FAILED"
            else:
                msg = "Calculating statistics... "
                print(msg)
                log_msgs.append(msg)
                corr_residuals = corr_residuals[np.where((corr_residuals != 999.0) & (corr_residuals < 0.1) & (corr_residuals > -0.1))]   # ignore outliers
                if corr_residuals.size == 0:
                    msg1 = " * Unable to calculate statistics because difference array has all outlier values. Test will be set to FAILED."
                    print(msg1)
                    log_msgs.append(msg1)
                    test_result = "FAILED"
                else:
                    stats_and_strings = auxfunc.print_stats(corr_residuals, "Difference", float(threshold_diff), abs=True)
                    stats, stats_print_strings = stats_and_strings
                    corr_residuals_mean, corr_residuals_median, corr_residuals_std = stats
                    for msg in stats_print_strings:
                        log_msgs.append(msg)

                    # This is the key argument for the assert pytest function
                    median_diff = False
                    if abs(corr_residuals_median) <= float(threshold_diff):
                        median_diff = True
                    if median_diff:
                        test_result = "PASSED"
                    else:
                        test_result = "FAILED"

                    msg = " *** Result of the test: "+test_result+"\n"
                    print(msg)
                    log_msgs.append(msg)
                    total_test_result.append(test_result)

    if writefile:
        outfile_name = step_input_filename.replace("srctype", "_calcuated_FS_PS_pathloss")
        compfile_name = step_input_filename.replace("srctype", "_comparison_FS_PS_pathloss")

        # create the fits list to hold the calculated flat values for each slit
        outfile.writeto(outfile_name, overwrite=True)

        # this is the file to hold the image of pipeline-calculated difference values
        compfile.writeto(compfile_name, overwrite=True)

        msg = "\nFits file with calculated pathloss values of each slit saved as: "
        print(msg)
        log_msgs.append(msg)
        print(outfile_name)
        log_msgs.append(outfile_name)

        msg = "Fits file with comparison (pipeline pathloss - calculated pathloss) saved as: "
        print(msg)
        log_msgs.append(msg)
        print(compfile_name)
        log_msgs.append(compfile_name)
   
    # If all tests passed then pytest is marked as PASSED, else it is FAILED
    FINAL_TEST_RESULT = False
    for t in total_test_result:
        if t == "FAILED":
            FINAL_TEST_RESULT = False
            break
        else:
            FINAL_TEST_RESULT = True

    if FINAL_TEST_RESULT:
        msg = "\n *** Final pathloss test result reported as PASSED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "All slits PASSED path_loss test."
    else:
        msg = "\n *** Final pathloss test result reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "One or more slits FAILED path_loss test."

    # end the timer
    pathloss_end_time = time.time() - pathtest_start_time
    if pathloss_end_time > 60.0:
        pathloss_end_time = pathloss_end_time/60.0  # in minutes
        pathloss_tot_time = "* Script FS_PS.py took ", repr(pathloss_end_time) + " minutes to finish."
        if pathloss_end_time > 60.0:
            pathloss_end_time = pathloss_end_time/60.  # in hours
            pathloss_tot_time = "* Script FS_PS.py took ", repr(pathloss_end_time) + " hours to finish."
    else:
        pathloss_tot_time = "* Script FS_PS.py took ", repr(pathloss_end_time) + " seconds to finish."
    print(pathloss_tot_time)
    log_msgs.append(pathloss_tot_time)

    return FINAL_TEST_RESULT, result_msg, log_msgs


if __name__ == '__main__':

    print("\n  ** using pipeline version: ", jwst.__version__, "** \n")

    pipeline_path = '/Users/tking/nirspec_pipe_testing_tool'

    # input parameters that the script expects
    step_input_filename = '/Users/tking/Documents/FS_NRS1_V2_test_run_pipeline/final_output_caldet1_NRS1_srctype.fits'
    reffile = '/Users/tking/jwst_nirspec_pathloss_0001.fits'
    comparison_filename = '/Users/tking/Documents/FS_NRS1_V2_test_run_pipeline/final_output_caldet1_NRS1_pathloss.fits'

    # name of the output images
    writefile = True

    # set the names of the resulting plots
    plot_name = None  # "FS_pathtest_histogram.pdf"

    # Run the principal function of the script
    median_diff = pathtest(step_input_filename, reffile, comparison_filename, writefile=writefile,
                           show_figs=False, save_figs=False,
                           plot_name=plot_name, threshold_diff=1.0e-7,
                           debug=True)
