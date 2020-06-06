import time
import os
from astropy import wcs
import numpy as np
from astropy.io import fits
import scipy
import pdb

import jwst
from jwst.pathloss import pathloss
from gwcs import wcstools
from jwst import datamodels
from jwst.assign_wcs import util
from . import auxiliary_functions as auxfunc
from astropy.visualization import ImageNormalize
import matplotlib.pyplot as plt

"""
This script tests the MOS pipeline pathloss step output for a Point Source.
"""

# HEADER
__author__ = "T King"
__version__ = "1.1"

# HISTORY
# October 19, 2019 - Version 1.0: initial version started
# February 25, 2020 - Version 1.1: Updted for Pep8 compliancy
# June 5, 2020 - Version 1.2: Included options to test additional interpolation methods

def get_corr_val(lambda_val, wave_ref, ref_ext, ref_xy, slit_x, slit_y):
    """Perform interpolation to get the correction value."""
    index = np.where(wave_ref[:, 0, 0] == lambda_val)
    plcor_slice = ref_ext[index[0][0]].reshape(ref_ext[index[0][0]].size)
    corr_val = scipy.interpolate.griddata(ref_xy[:plcor_slice.size], plcor_slice,
                                          np.asarray([slit_x, slit_y]),
                                          method='linear')
    return corr_val[0]

def get_corr_val_cubic(lambda_val, wave_ref, ref_ext, ref_xy, slit_x, slit_y, first_interp_method):
    """Perofrm an alternative interpolation in order to compare interpolation methods."""
    index = np.where(wave_ref[:, 0, 0] == lambda_val)
    plcor_slice = ref_ext[index[0][0]].reshape(ref_ext[index[0][0]].size)
    f = scipy.interpolate.interp2d(ref_xy[:plcor_slice.size][:,0], ref_xy[:plcor_slice.size][:,1], plcor_slice, kind='linear')
    corr_val = f(slit_x, slit_y)
    return corr_val[0]

def get_mos_ps_uni_extensions(fits_file_name, is_point_source):
    """
    This functions obtains all point or extended source extensions in ref file
    Args:
        fits_file_name: str, name of the fits file of interest
        is_point_source: boolean, true if point source; false if extended source
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


def pathtest(step_input_filename, reffile, comparison_filename, writefile=True,
             show_figs=True, save_figs=False, plot_name=None,
             threshold_diff=1.0e-7, debug=False):
    """
    This function calculates the difference between the pipeline and
    the calculated pathloss values. The functions use the output of
    the compute_world_coordinates.py script.
    Args:
        step_input_filename: str, name of the output fits file from
                             the sourcetype step (with full path)
        reffile: str, path to the pathloss MOS reference fits file
        comparison_filename: str, path to comparison pipeline pathloss file
        writefile: boolean, if True writes the fits files of the calculated
                   pathloss and difference images
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots (the 3 plots can be saved or
                   not independently with the function call)
        plot_name: string, desired name (if name is not given,
                   the plot function will name the plot by default)
        threshold_diff: float, threshold difference between pipeline output
                        and ESA file
        debug: boolean, if true, print statements will show on-screen
    Returns:
        - 1 plot, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to threshold.
        - log_msgs: list, all print statements are captured in this variable
    """

    log_msgs = []

    # start the timer
    pathtest_start_time = time.time()

    # get info from the input sourcetype file header
    msg = 'step_input_filename='+step_input_filename
    print(msg)
    log_msgs.append(msg)
    exptype = fits.getval(step_input_filename, "EXP_TYPE", 0)
    grat = fits.getval(step_input_filename, "GRATING", 0)
    filt = fits.getval(step_input_filename, "FILTER", 0)

    msg = "pathloss file:  Grating:"+grat+"  Filter:"+filt+" EXP_TYPE:"+exptype
    print(msg)
    log_msgs.append(msg)

    msg = "Using reference file: "+reffile
    print(msg)
    log_msgs.append(msg)

    if writefile:
        # create fits list to hold the calculated pathloss values for each slit
        hdu0 = fits.PrimaryHDU()
        outfile = fits.HDUList()
        outfile.append(hdu0)

        # create fits list to hold pipeline-calculated difference values
        hdu0 = fits.PrimaryHDU()
        compfile = fits.HDUList()
        compfile.append(hdu0)

    # list to determine if pytest is passed or not
    total_test_result = []

    # get all the science extensions
    is_point_source = True
    print("Retrieving exensions")
    ps_uni_ext_list = get_mos_ps_uni_extensions(reffile, is_point_source)

    # get files
    print("""Checking if files exist & obtaining datamodels.
          This takes a few minutes...""")
    if os.path.isfile(comparison_filename):
        if debug:
            msg = 'Comparison file does exist.'
            print(msg)
    else:
        result_msg = """Comparison file does NOT exist.
                     Pathloss test will be skipped."""
        print(result_msg)
        log_msgs.append(result_msg)
        result = 'skip'
        return result, result_msg, log_msgs

    # get the comparison data model
    pathloss_pipe = datamodels.open(comparison_filename)
    if debug:
        print('Retrieved comparison datamodel.')

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

    # loop through the slits
    msg = "Looping through the slits... "
    print(msg)
    log_msgs.append(msg)

    slit_val = 1
    for slit, pipe_slit in zip(pl.slits, pathloss_pipe.slits):
        try:
            nshutters = util.get_num_msa_open_shutters(slit.shutter_state)
            if is_point_source:
                if nshutters == 3:
                    shutter_key = "MOS1x3"
                elif nshutters == 1:
                    shutter_key = "MOS1x1"
                ext = ps_uni_ext_list[0][shutter_key]
                print("Retrieved point source extension {}".format(ext))
            if is_point_source is False:
                if nshutters == 1:
                    shutter_key = "MOS1x1"
                elif nshutters == 3:
                    shutter_key = "MOS1x3"
                ext = ps_uni_ext_list[1][shutter_key]
                print("Retrieved extended source extension {}".format(ext))
        except KeyError:
            print("Unable to retrieve extension. Using 1, but may be 5")
            ext = 1  # or 5 (1: one slit open. 5: three adjacent slits open)

        mode = "MOS"

        slit_id = pipe_slit.name
        print('Working with slitlet ', slit_id)

        if slit.name == slit_id:
            msg = """Slitlet name in fits file previous to pathloss
                   and in pathloss output file are the same."""
            log_msgs.append(msg)
            print(msg)
        else:
            msg = """* Missmatch of slitlet names in fits file previous
                  to pathloss and in pathloss output file. Skipping test."""
            result = 'skip'
            log_msgs.append(msg)
            return result, msg, log_msgs

        wcs_obj = slit.meta.wcs

        # get the wavelength
        x, y = wcstools.grid_from_bounding_box(wcs_obj.bounding_box,
                                               step=(1, 1), center=True)
        ra, dec, wave = wcs_obj(x, y)   # wave is in microns

        # get positions of source in file:
        slit_x = slit.source_xpos
        slit_y = slit.source_ypos*(-1)  # Scaling introduced for time being because error in assign_wcs transformation
        if debug:
            print("slit_x, slit_y (" + str(slit_x) + ", " + str(slit_y) + ")")

        ref_ext = fits.getdata(reffile, ext)
        hdul = fits.open(reffile)

        # plcor_ref = hdul[1].data
        if debug:
            print("ref_ext.shape", ref_ext.shape)
        w = wcs.WCS(hdul[1].header)

        # make cube
        w1, y1, x1 = np.mgrid[:ref_ext.shape[0], :ref_ext.shape[1],
                              :ref_ext.shape[2]]
        slitx_ref, slity_ref, wave_ref = w.all_pix2world(x1, y1, w1, 0)

        comp_sci = pipe_slit.data
        previous_sci = slit.data

        pipe_correction = pipe_slit.pathloss

        # Set up source position manually to test correction at nonzero point:
        # pipe_x = -0.2
        # pipe_y = -0.2

        # slit_x = pipe_x
        # slit_y = pipe_y
        # print("""WARNING: Using manually set slit_x and slit_y: ({}, {})!
        #       The pipeline correction will not use manually set values
        #       and thus the residuals will change""".format(slit_x, slit_y))

        wave_sci = wave * 10**(-6)  # microns --> meters
        wave_sci_flat = wave_sci.reshape(wave_sci.size)
        wave_ref_flat = wave_ref.reshape(wave_ref.size)

        ref_xy = np.column_stack((slitx_ref.reshape(slitx_ref.size),
                                 slity_ref.reshape(slitx_ref.size)))

        correction_list = [(get_corr_val(lambda_val, wave_ref, ref_ext,
                           ref_xy, slit_x, slit_y)) for
                           lambda_val in wave_ref_flat]
        correction_array = np.asarray(correction_list)

        # Option to test alternative interpolation method:
        # first_interp_method='linear'
        # second_interp_method = 'linear'
        
        # if first_interp_method == 'linear':
        #     correction_array_cubic = correction_array
        # else:
        #     correction_list_cubic = [(get_corr_val_cubic(lambda_val, wave_ref, ref_ext,
        #                              ref_xy, slit_x, slit_y, first_interp_method)) for lambda_val in wave_ref_flat]
        #     correction_array_cubic = np.asarray(correction_list_cubic)
        
        lambda_array = wave_ref_flat

        # get correction value for each pixel
        corr_vals = np.interp(wave_sci_flat, lambda_array, correction_array)
        # corr_vals_cubic = np.interp(wave_sci_flat, lambda_array, correction_array_cubic)
        corr_vals = corr_vals.reshape(wave_sci.shape)
        # corr_vals_cubic = corr_vals_cubic.reshape(wave_sci.shape)
        corrected_array = previous_sci/corr_vals

        # Plots:
        # my correction values
        fig = plt.figure(figsize=(12, 10))
        plt.subplot(221)
        norm = ImageNormalize(corr_vals)
        plt.imshow(corr_vals, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
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
        # residuals (pipeline correction-my correction)
        corr_residuals = pipe_correction-corr_vals
        plt.subplot(223)
        norm = ImageNormalize(corr_residuals)
        plt.imshow(corr_residuals, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Correction residuals')
        plt.colorbar()
        # Calculated Corrected Array
        plt.subplot(224)
        norm = ImageNormalize(corrected_array)
        plt.imshow(corrected_array, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Calculated Corrected Array')
        plt.colorbar()
        # #cubic v linear
        # plt.subplot(324)
        # norm = ImageNormalize(corr_vals_cubic-corr_vals)
        # plt.imshow(wave_sci, vmin = -0.00025, vmax=0.00025, aspect=10.0, origin='lower',
        #            cmap='viridis')
        # plt.xlabel('dispersion in pixels')
        # plt.ylabel('y in pixels')
        # plt.title('corr_vals_1+'+first_interp_method+'2'+second_interp_method+'-corr_vals_linear')
        # plt.colorbar()
        
        fig.suptitle("MOS PS at ({}, {}) Pathloss Correction Testing Slit ".format(slit_x, slit_y)
                     + str(slit_id))

        if save_figs:
            plt_name = step_input_filename + "_Pathloss_test_slitlet_" + str(mode) + "_" + str(slit_val) + ".png"
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_figs:
            plt.show()
        elif not save_figs and not show_figs:
            msg = "Not making plots because both show_figs and save_figs were set to False."
            if debug:
                print(msg)
            log_msgs.append(msg)
        elif not save_figs:
            msg = "Not saving plots because save_figs was set to False."
            if debug:  
                print(msg)
            log_msgs.append(msg)
        plt.close()

        # do not overwrite plots
        slit_val = slit_val + 1

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
                msg1 = """Unable to calculate statistics because difference
                       array has all values as NaN. Test will be set to FAILED.
                       """
                print(msg1)
                log_msgs.append(msg1)
                test_result = "FAILED"
                # delfg_mean, delfg_median, delfg_std = np.nan, np.nan, np.nan
                # stats = [delfg_mean, delfg_median, delfg_std]
            else:
                msg = "Calculating statistics... "
                print(msg)
                log_msgs.append(msg)
                # ignore outliers:
                corr_residuals = corr_residuals[np.where((corr_residuals != 999.0)
                                                & (corr_residuals < 0.1)
                                                & (corr_residuals > -0.1)
                                                & ~np.isnan(corr_residuals))]
                if corr_residuals.size == 0:
                    msg1 = """Unable to calculate statistics because
                              difference array has all outlier values.
                              Test will be set to FAILED."""
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
        outfile_name = step_input_filename.replace("srctype", "MOS_PS_calcuated_FS_UNI_pathloss")
        compfile_name = step_input_filename.replace("srctype", "MOS_PS_comparison_FS_UNI_pathloss")

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

    # If all tests passed pytest will be marked PASSED, else FAILED
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
        pathloss_tot_time = "* Script MSA.py took ", repr(pathloss_end_time)+" minutes to finish."
        if pathloss_end_time > 60.0:
            pathloss_end_time = pathloss_end_time/60.  # in hours
            pathloss_tot_time = "* Script MSA.py took ", repr(pathloss_end_time)+" hours to finish."
    else:
        pathloss_tot_time = "* Script MSA.py took ", repr(pathloss_end_time)+" seconds to finish."
    print(pathloss_tot_time)
    log_msgs.append(pathloss_tot_time)

    return FINAL_TEST_RESULT, result_msg, log_msgs


if __name__ == '__main__':

    # print pipeline version
    print("\n  ** using pipeline version: ", jwst.__version__, "** \n")
    print("Script takes a few minutes to run per slit")

    # input parameters that script expects
    working_dir = '/Users/tking/nirspec_pipe_testing_tool/calwebb_spec2_pytests'
    comparison_filename = '/Users/tking/Documents/king_MOS_grating/final_output_caldet1_NRS1_new_twofew_pathloss.fits'
    step_input_filename = '/Users/tking/Documents/king_MOS_grating/final_output_caldet1_NRS1_new_twofew_srctype.fits'
    reffile = '/Users/tking/jwst_nirspec_pathloss_0002.fits'

    # name of the output images
    writefile = True

    # set the names of the resulting plots
    plot_name = None  # eg, "MSA_pltest_histogram.pdf"

    # Run the principal function of the script
    median_diff = pathtest(step_input_filename, reffile, comparison_filename,
                           writefile=writefile, show_figs=False,
                           save_figs=True, plot_name=plot_name,
                           threshold_diff=1.0e-7, debug=True)
