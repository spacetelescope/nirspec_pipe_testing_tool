import time
import os
from astropy import wcs
import numpy as np
from astropy.io import fits
import scipy
import math
import argparse
import sys

import jwst
from gwcs import wcstools
from jwst import datamodels
from jwst.assign_wcs import nirspec
from . import auxiliary_functions as auxfunc

from astropy.visualization import (ImageNormalize, AsinhStretch)
from scipy.interpolate import interp1d

import matplotlib
import matplotlib.pyplot as plt

"""
This script tests the IFU pipeline pathloss step output for a Point Source.
"""

# HEADER
__author__ = "T King & M Pena-Guerrero"
__version__ = "1.4"

# HISTORY
# Oct 19, 2019 - Version 1.0: initial version started
# Feb 12, 2020 - Version 1.1: Tests pass for every slit
# Feb 26, 2020 - Mostly pep8 compliant
# June 11, 2020 - Version 1.2: Added changes to be able to run within NPTT
# September 25, 2020 - Version 1.3: Added option to use either data model or fits file as input for the test
# January 2021 - Version 1.4: Implemented option to use datamodels instead of fits files as input


def pathtest(step_input_filename, reffile, comparison_filename,
             writefile=True, show_figs=False, save_figs=True,
             threshold_diff=1.0e-7, debug=False):
    """
    This function calculates the difference between the pipeline and
    calculated pathloss values.
    Args:
        step_input_filename: str, full path name of sourcetype output fits file
        reffile: str, path to the pathloss IFU reference fits file
        comparison_filename: str, path to comparison pipeline pathloss file
        writefile: boolean, if True writes calculated flat and
                   difference image fits files
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots
        threshold_diff: float, threshold diff between pipeline and ESA file
        debug: boolean, if true print statements will show on-screen
    Returns:
        - 1 plot, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to threshold.
        - log_msgs: list, all print statements are captured in this variable
    """

    log_msgs = []

    # start the timer
    pathtest_start_time = time.time()

    # get info from the rate file header
    msg = 'step_input_filename='+step_input_filename
    print(msg)
    log_msgs.append(msg)
    exptype = fits.getval(step_input_filename, "EXP_TYPE", 0)
    grat = fits.getval(step_input_filename, "GRATING", 0)
    filt = fits.getval(step_input_filename, "FILTER", 0)

    msg = "pathloss file: Grating:"+grat+" Filter:"+filt+" EXP_TYPE:"+exptype
    print(msg)
    log_msgs.append(msg)

    msg = "Using reference file: "+reffile
    print(msg)
    log_msgs.append(msg)

    is_point_source = True

    if writefile:
        # create the fits list to hold the calculated pathloss values for each slit
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
    if is_point_source:
        ext = 1  # for all PS IFU

    # get files
    print('Checking files exist & obtaining datamodels. Takes a few mins...')
    if isinstance(comparison_filename, str):
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
        ifu_pipe_model = datamodels.open(comparison_filename)
        if debug:
            print('got comparison datamodel!')
    else:
        ifu_pipe_model = comparison_filename

    if os.path.isfile(step_input_filename):
        if debug:
            print('Input file does exist.')
    else:
        result_msg = 'Input file does NOT exist. Skipping pathloss test.'
        log_msgs.append(result_msg)
        result = 'skip'
        return result, result_msg, log_msgs

    # get the input data model
    ifu_input_model = datamodels.open(step_input_filename)
    if debug:
        print('got input datamodel!')

    # get slices (instead of using .slit)
    pl_ifu_slits = nirspec.nrs_ifu_wcs(ifu_input_model)
    print("got input slices")

    plcor_ref_ext = fits.getdata(reffile, ext)
    hdul = fits.open(reffile)
    plcor_ref = hdul[1].data
    w = wcs.WCS(hdul[1].header)

    w1, y1, x1 = np.mgrid[:plcor_ref.shape[0], :plcor_ref.shape[1], :plcor_ref.shape[2]]
    slitx_ref, slity_ref, wave_ref = w.all_pix2world(x1, y1, w1, 0)

    # these are full 2048 * 2048 files:
    previous_sci = fits.getdata(step_input_filename, "SCI")
    comp_sci = ifu_pipe_model.data
    pathloss_divided = comp_sci/previous_sci

    # Can manually test correction at nonzero point
    # slit_x = -0.3
    # slit_y = 0.3
    # if debug:
    #     print("""WARNING: Using manually set slit_x and slit_y! 
    #           The pipeline correction will not use manually set values 
    #           and thus the residuals will change""")

    # set up generals for all plots
    font = {'weight': 'normal',
            'size': 10}
    matplotlib.rc('font', **font)

    # loop through the slices
    msg = " Looping through the slices... "
    print(msg)
    log_msgs.append(msg)

    slit_list = np.ndarray.tolist(np.arange(0, 30))
    for slit, slit_num, pipeslit in zip(pl_ifu_slits, slit_list):
        print("working with slice {}".format(slit_num))

        # NIRSpec implementation
        x, y = wcstools.grid_from_bounding_box(slit.bounding_box, step=(1, 1))
        ra, dec, wave = slit(x, y)

        slit_x = 0  # This assumption is made for IFU sources only
        slit_y = 0

        if debug:
            print("slit_x, slit_y", slit_x, slit_y)

        correction_array = np.array([])
        lambda_array = np.array([])

        wave_sci = wave * 10**(-6)  # microns --> meters
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

        box = slit.bounding_box
        small_y = box[0][0]
        big_y = box[0][1]
        small_x = box[1][0]
        big_x = box[1][1]

        left = int(math.trunc(small_x))
        right = int(math.ceil(big_x))
        bottom = int(math.trunc(small_y))
        top = int(math.ceil(big_y))

        full_cut2slice = previous_sci[left:right, bottom:top]
        print("shapes:", full_cut2slice.shape, corr_vals.shape)

        if full_cut2slice.shape != corr_vals.shape:
            value = 0
            while (((full_cut2slice.shape[0] - corr_vals.shape[0]) != 0) or
                   ((full_cut2slice.shape[1] - corr_vals.shape[1]) != 0)) and \
                    value < 7:  # can delete second criteria once all pass
                if value == 6:
                    print("WARNING: may be in infinite loop!")
                x_amount_off = full_cut2slice.shape[0]-corr_vals.shape[0]
                if x_amount_off >= 1:
                    if x_amount_off % 2 == 0:  # need to add two values, so do one on each side
                        right = right-1
                        left = left+1
                        print("ALTERED SHAPE OF SLICE: V1")
                        value = value + 1
                        print("left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                    else:  # just add one value
                        left = left+1
                        print("ALTERED SHAPE OF SLICE: V2")
                        value = value+1
                        print("left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                elif x_amount_off <= -1:
                    if x_amount_off % 2 == 0:
                        right = right+1
                        left = left-1
                        print("ALTERED SHAPE OF SLICE: V3")
                        value = value+1
                        print("left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                    else:
                        left = left-1
                        print("ALTERED SHAPE OF SLICE: V4")
                        value = value+1
                        print("left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                y_amount_off = full_cut2slice.shape[1]-corr_vals.shape[1]
                if y_amount_off >= 1:
                    if y_amount_off % 2 == 0:
                        bottom = bottom-1
                        top = top+1
                        print("ALTERED SHAPE OF SLICE: V5")
                        value = value+1
                        print("left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                    else:
                        bottom = bottom+1
                        print("ALTERED SHAPE OF SLICE: V6")
                        value = value+1
                        print("left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                elif y_amount_off <= -1:
                    if y_amount_off % 2 == 0:
                        top = top+1
                        bottom = bottom-1
                        print("ALTERED SHAPE OF SLICE: V7")
                        value = value+1
                        print("left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                    else:
                        bottom = bottom-1
                        print("ALTERED SHAPE OF SLICE: V8")
                        value = value+1
                        print("left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                full_cut2slice = previous_sci[left:right, bottom:top]
                print("final left {}, right {}, top {}, bottom {}".format(left, right, top, bottom))
                print("NEW SHAPE OF SLICE: {} and corr_vals.shape: {}".format(full_cut2slice.shape, corr_vals.shape))

        if full_cut2slice.shape != corr_vals.shape:
            print("shapes did not match! full_cut2slice: {}, corr_vals {}".format(full_cut2slice.shape,
                                                                                  corr_vals.shape))
            continue

        corrected_array = full_cut2slice/corr_vals

        pipe_correction = pathloss_divided[left:right, bottom:top]
        if pipe_correction.shape != corr_vals.shape:
            print("shapes did not match! pipe_correction: {}, corr_vals {}".format(pipe_correction.shape,
                                                                                   corr_vals.shape))
            continue

        prev_sci_slit = previous_sci[left:right, bottom:top]
        if prev_sci_slit.shape != corr_vals.shape:
            print("shapes did not match! prev_sci_slit: {}, corr_vals {}".format(prev_sci_slit.shape, corr_vals.shape))
            continue

        comp_sci_slit = comp_sci[left:right, bottom:top]
        if comp_sci_slit.shape != corr_vals.shape:
            print("shapes did not match! comp_sci_slit: {}, corr_vals {}".format(comp_sci_slit.shape, corr_vals.shape))
            continue

        # Plots:
        step_input_filepath = step_input_filename.replace(".fits", "")
        # my correction values
        fig = plt.figure(figsize=(15, 15))
        plt.subplot(221)
        norm = ImageNormalize(corr_vals)
        plt.imshow(corr_vals, vmin=0.999995, vmax=1.000005, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Calculated Correction')
        plt.colorbar()
        # pipe correction
        plt.subplot(222)
        norm = ImageNormalize(pipe_correction)
        plt.imshow(pipe_correction, vmin=0.999995, vmax=1.000005, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Pipe Correction')
        plt.colorbar()
        # residuals (pipe correction - my correction)
        if pipe_correction.shape == corr_vals.shape:
            corr_residuals = pipe_correction - corr_vals
            plt.subplot(223)
            norm = ImageNormalize(corr_residuals)
            plt.imshow(corr_residuals, vmin=-0.000000005, vmax=0.000000005,
                       aspect=10.0, origin='lower', cmap='viridis')
            plt.xlabel('dispersion in pixels')
            plt.ylabel('y in pixels')
            plt.title('Correction residuals')
            plt.colorbar()
        # my science data after pathloss
        plt.subplot(224)
        norm = ImageNormalize(corrected_array)
        plt.imshow(corrected_array, vmin=0, vmax=300, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.title('My slit science data after pathloss')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.colorbar()
        fig.suptitle("IFU PS Pathloss Correction Testing" + str(corrected_array.shape))
        fig.tight_layout(pad=3.0)

        if show_figs:
            plt.show()
        if save_figs:
            plt_name = step_input_filepath + "_Pathloss_test_slitlet_IFU_PS_" + str(slit_num) + ".png"
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)

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
        plt.clf()

        # create fits file to hold the calculated pathloss for each slit
        if writefile:
            msg = "Saving the fits files with the calculated pathloss for each slit..."
            print(msg)
            log_msgs.append(msg)

            # this is the file to hold the image of pipeline-calculated difference values
            outfile_ext = fits.ImageHDU(corr_vals, name=str(slit_num))
            outfile.append(outfile_ext)

            # this is the file to hold the image of pipeline-calculated difference values
            compfile_ext = fits.ImageHDU(corr_residuals, name=str(slit_num))
            compfile.append(compfile_ext)

        if corr_residuals[~np.isnan(corr_residuals)].size == 0:
            msg1 = " * Unable to calculate statistics because difference array has all values as NaN. " \
                   "Test will be set to FAILED."
            print(msg1)
            log_msgs.append(msg1)
            test_result = "FAILED"
        else:
            msg = "Calculating statistics... "
            print(msg)
            log_msgs.append(msg)
            corr_residuals = corr_residuals[np.where((corr_residuals != 999.0) & (corr_residuals < 0.1) &
                                                     (corr_residuals > -0.1))]   # ignore outliers
            if corr_residuals.size == 0:
                msg1 = """ * Unable to calculate statistics because
                           difference array has all outlier values.
                           Test will be set to FAILED."""
                print(msg1)
                log_msgs.append(msg1)
                test_result = "FAILED"
            else:
                stats_and_strings = auxfunc.print_stats(corr_residuals, "Difference", float(threshold_diff),
                                                        absolute=True)
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

    hdul.close()
    ifu_pipe_model.close()
    ifu_input_model.close()

    if writefile:
        outfile_name = step_input_filename.replace("srctype", "_calcuated_pathloss")
        compfile_name = step_input_filename.replace("srctype", "_comparison_pathloss")

        # create the fits list to hold the calculated pathloss values for each slit
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

    # If all tests passed then pytest will be marked as PASSED, else FAILED
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
        pathloss_tot_time = "* Script IFU_PS.py took ", repr(pathloss_end_time)+" minutes to finish."
        if pathloss_end_time > 60.0:
            pathloss_end_time = pathloss_end_time/60.  # in hours
            pathloss_tot_time = "* Script IFU_PS.py took ", repr(pathloss_end_time)+" hours to finish."
    else:
        pathloss_tot_time = "* Script IFU_PS.py took ", repr(pathloss_end_time)+" seconds to finish."
    print(pathloss_tot_time)
    log_msgs.append(pathloss_tot_time)

    return FINAL_TEST_RESULT, result_msg, log_msgs


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("step_input_filename",
                        action='store',
                        default=None,
                        help='Name of input fits file prior to assign_wcs step, i.e. blah_rate.fits')
    parser.add_argument("reffile",
                        action='store',
                        default=None,
                        help='Path and name of reference file to use for the test.')
    parser.add_argument("comparison_filename",
                        action='store',
                        default=None,
                        help='Path and name the comparison file, i.e. the pathloss output file')
    parser.add_argument("-w",
                        dest="writefile",
                        action='store_false',
                        default=True,
                        help='Use flag -w to NOT write files with calculated correction.')
    parser.add_argument("-f",
                        dest="save_figs",
                        action='store_false',
                        default=True,
                        help='Use flag -f to NOT save final figures.')
    parser.add_argument("-s",
                        dest="show_figs",
                        action='store_true',
                        default=False,
                        help='Use flag -s to show final figures.')
    parser.add_argument("-t",
                        dest="threshold_diff",
                        action='store',
                        default=1.0e-07,
                        type=float,
                        help='Use flag -t to change the default threshold (currently set to 1.0e-07).')
    parser.add_argument("-d",
                        dest="debug",
                        action='store_true',
                        default=False,
                        help='Use flag -d to turn on debug mode.')
    args = parser.parse_args()

    # Set variables
    step_input_filename = args.step_input_filename
    reffile = args.reffile
    comparison_filename = args.comparison_filename
    writefile = args.writefile
    save_figs = args.save_figs
    show_figs = args.show_figs
    threshold_diff = args.threshold_diff
    debug = args.debug

    # print pipeline version
    print("\n  ** using pipeline version: ", jwst.__version__, "** \n")

    # Run the principal function of the script
    pathtest(step_input_filename, reffile, comparison_filename, writefile=writefile,
             show_figs=show_figs, save_figs=save_figs, threshold_diff=threshold_diff,
             debug=debug)


if __name__ == '__main__':
    sys.exit(main())

