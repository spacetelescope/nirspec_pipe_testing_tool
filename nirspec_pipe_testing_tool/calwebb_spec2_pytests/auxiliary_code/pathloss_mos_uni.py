import time
import os
from astropy import wcs
import numpy as np
from astropy.io import fits
import scipy
import pdb
import argparse
import sys

import jwst
from gwcs import wcstools
from jwst import datamodels
from jwst.assign_wcs import util

from . import auxiliary_functions as auxfunc

from astropy.visualization import (ImageNormalize, AsinhStretch)
from scipy.interpolate import interp1d

import matplotlib
import matplotlib.pyplot as plt

"""
This script tests the MSA pipeline pathloss step output for a Uniform Source.
"""

# HEADER
__author__ = "T King & M Pena-Guerrero"
__version__ = "1.3"

# HISTORY
# October 19, 2019 - Version 1.0: initial version started
# January 14, 2020 - Version 1.1: passes tests
# February 26, 2020 - Version 1.2: mostly pep8 compliant
# September 25, 2020 - Version 1.3: Added option to use either data model or fits file as input for the test


def get_mos_ps_uni_extensions(fits_file_name, is_point_source):
    """
    This functions obtains all the point source or extended source extensions
    in the given file
    Args:
        fits_file_name: str, name of the fits file of interest
        is_point_source: boolean, true if point source; false if extended source
    Returns:
        ps_dict, uni_dict: list of the numbers of the science extensions
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


def pathtest(step_input_filename, reffile, comparison_filename,
             writefile=True, show_figs=False, save_figs=True, threshold_diff=1.0e-7, debug=False):
    """
    This function calculates the difference between the pipeline and
    calculated pathloss values. The functions use the output of sourcetype.
    Args:
        step_input_filename: str, name of previous pipeline step output fits file
        reffile: str, path to the pathloss MSA UNI reference fits files
        comparison_filename: str, path to comparison pipeline pathloss file
        writefile: boolean, if True writes fits files of calculated pathloss
                   and difference images
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots
        threshold_diff: float, threshold difference between pipeline output
                        & ESA file
        debug: boolean, if true print statements will show on-screen
    Returns:
        - 1 plot, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to threshold.
        - log_msgs: list, all print statements are captured in this variable
    """

    log_msgs = []

    # start the timer
    pathtest_start_time = time.time()

    # get info from the input previous pipeline step file/datamodel
    print("Checking if files exist and obtaining datamodels. This takes a few minutes...")
    if isinstance(step_input_filename, str):
        if os.path.isfile(step_input_filename):
            if debug:
                print('Input file does exist.')
            msg = 'step_input_filename='+step_input_filename
            print(msg)
            log_msgs.append(msg)

            # get the input data model
            pl = datamodels.open(step_input_filename)
            if debug:
                print('got input datamodel!')
        else:
            result_msg = 'Input file does NOT exist. Skipping pathloss test.'
            log_msgs.append(result_msg)
            result = 'skip'
            return result, result_msg, log_msgs
    else:
        pl = step_input_filename

    # get comparison data
    if isinstance(comparison_filename, str):
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

    else:
        pathloss_pipe = comparison_filename

    # get info from data model
    det = pl.meta.instrument.detector
    lamp = pl.meta.instrument.lamp_state
    grat = pl.meta.instrument.grating
    filt = pl.meta.instrument.filter
    exptype = pl.meta.exposure.type

    msg = "from datamodel  -->     Detector: " + det + "   Grating: " + grat + "   Filter: " + \
          filt + "   Lamp: " + lamp + "   EXP_TYPE: " + exptype
    print(msg)
    log_msgs.append(msg)

    # get the reference files
    msg = "Using reference file: "+reffile
    print(msg)
    log_msgs.append(msg)

    if writefile:
        # create the fits list to hold the calculated flat values for each slit
        hdu0 = fits.PrimaryHDU()
        outfile = fits.HDUList()
        outfile.append(hdu0)

        # create fits list to hold image of pipeline-calculated diff values
        hdu0 = fits.PrimaryHDU()
        compfile = fits.HDUList()
        compfile.append(hdu0)

    # list to determine if pytest is passed or not
    total_test_result = []

    # loop through the slits
    msg = " Looping through the slits... "
    print(msg)
    log_msgs.append(msg)

    slit_val = 0
    for slit, pipe_slit in zip(pl.slits, pathloss_pipe.slits):
        slit_val = slit_val+1

        mode = "MOS"

        is_point_source = False

        print("Retrieving extensions")
        ps_uni_ext_list = get_mos_ps_uni_extensions(reffile, is_point_source)

        slit_id = slit.name
        try:
            nshutters = util.get_num_msa_open_shutters(slit.shutter_state)
            if is_point_source:
                if nshutters == 3:
                    shutter_key = "MOS1x3"
                elif nshutters == 1:
                    shutter_key = "MOS1x1"
                ext = ps_uni_ext_list[0][shutter_key]
                print("Retrieved point source extension")
            if is_point_source is False:
                if nshutters == 1:
                    shutter_key = "MOS1x1"
                elif nshutters == 3:
                    shutter_key = "MOS1x3"
                ext = ps_uni_ext_list[1][shutter_key]
                print("Retrieved extended source extension {}".format(ext))
        except KeyError:
            print("Unable to retrieve extension. Using ext 3, but may be 7")
            ext = 3

        wcs_obj = slit.meta.wcs

        x, y = wcstools.grid_from_bounding_box(wcs_obj.bounding_box, step=(1, 1), center=True)
        ra, dec, wave = slit.meta.wcs(x, y)
        wave_sci = wave * 10**(-6)   # microns --> meters

        plcor_ref_ext = fits.getdata(reffile, ext)
        print("plcor_ref_ext.shape", plcor_ref_ext.shape)

        hdul = fits.open(reffile)

        plcor_ref = hdul[1].data
        w = wcs.WCS(hdul[1].header)

        w1, y1, x1 = np.mgrid[:plcor_ref.shape[0], : plcor_ref.shape[1],
                              :plcor_ref.shape[2]]
        slitx_ref, slity_ref, wave_ref = w.all_pix2world(x1, y1, w1, 0)

        comp_sci = pipe_slit.data
        previous_sci = slit.data

        pipe_correction = pipe_slit.pathloss_uniform

        # set up generals for all the plots
        font = {'weight': 'normal',
                'size': 10}
        matplotlib.rc('font', **font)

        corr_vals = np.interp(wave_sci, wave_ref[:, 0, 0], plcor_ref_ext)
        corrected_array = previous_sci/corr_vals

        # plots:
        # my correction values
        fig = plt.figure()
        ax = plt.gca()
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_xaxis().get_major_formatter().set_scientific(False)
        # calculated correction values
        plt.subplot(221)
        norm = ImageNormalize(corr_vals)
        plt.imshow(corr_vals, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.xlabel('x in pixels')
        plt.ylabel('y in pixels')
        plt.title('Calculated Correction')
        plt.colorbar()
        # pipeline correction values
        plt.subplot(222)
        norm = ImageNormalize(pipe_correction)
        plt.imshow(pipe_correction, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.xlabel('x in pixels')
        plt.ylabel('y in pixels')
        plt.title('Pipeline Correction')
        plt.colorbar()
        # residuals (pipe correction - my correction)
        corr_residuals = pipe_correction - corr_vals
        plt.subplot(223)
        norm = ImageNormalize(corr_residuals)
        plt.ticklabel_format(useOffset=False)
        plt.imshow(corr_residuals, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.xlabel('x in pixels')
        plt.ylabel('y in pixels')
        plt.title('Correction residuals')
        plt.colorbar()
        # my science data after pathloss
        plt.subplot(224)
        norm = ImageNormalize(corrected_array)
        plt.imshow(corrected_array, norm=norm, aspect=10.0, origin='lower',
                   cmap='viridis')
        plt.title('Corrected Data After Pathloss')
        plt.xlabel('x in pixels')
        plt.ylabel('y in pixels')
        plt.colorbar()
        fig.suptitle("MOS UNI Pathloss Calibration Testing")

        if show_figs:
            plt.show()
        if save_figs:
            step_input_filepath = step_input_filename.replace(".fits", "")
            plt_name = step_input_filepath+"Pathloss_test_slitlet_" + str(mode) + "_UNI_" + str(slit_id) + ".png"
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        plt.close()

        ax = plt.subplot(212)
        plt.hist(corr_residuals[~np.isnan(corr_residuals)], bins=100,
                 range=(-0.00000013, 0.00000013))
        plt.title('Residuals Histogram')
        plt.xlabel("Correction Value")
        plt.ylabel("Number of Occurences")
        nanind = np.isnan(corr_residuals)  # get all the nan indexes
        notnan = ~nanind  # get all the not-nan indexes
        arr_mean = np.mean(corr_residuals[notnan])
        arr_median = np.median(corr_residuals[notnan])
        arr_stddev = np.std(corr_residuals[notnan])
        plt.axvline(arr_mean, label="mean = %0.3e" % (arr_mean), color="g")
        plt.axvline(arr_median, label="median = %0.3e" % (arr_median),
                    linestyle="-.", color="b")
        str_arr_stddev = "stddev = {:0.3e}".format(arr_stddev)
        ax.text(0.73, 0.67, str_arr_stddev, transform=ax.transAxes,
                fontsize=16)
        plt.legend()
        plt.minorticks_on()

        # Show and/or save figures
        if save_figs:
            plt_name = step_input_filepath + "Pathlosstest_MOS_UNI_slitlet_" + slit_id + ".png"
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

        if corr_residuals[~np.isnan(corr_residuals)].size == 0:
            msg1 = """Unable to calculate statistics because difference
                    array has all values as NaN.
                    Test will be set to FAILED."""
            print(msg1)
            log_msgs.append(msg1)
            test_result = "FAILED"
        else:
            msg = "Calculating statistics... "
            print(msg)
            log_msgs.append(msg)
            # ignore outliers:
            corr_residuals = corr_residuals[np.where((corr_residuals != 999.0)
                                                     & (corr_residuals < 0.1)
                                                     & (corr_residuals > -0.1))]
            if corr_residuals.size == 0:
                msg1 = """ * Unable to calculate statistics because
                       difference array has all outlier values.
                       Test will be set to FAILED."""
                print(msg1)
                log_msgs.append(msg1)
                test_result = "FAILED"
            else:
                stats_and_strings = auxfunc.print_stats(corr_residuals, "Difference", float(threshold_diff), absolute=True)
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
        outfile_name = step_input_filename.replace("srctype", det+"_calcuated_pathloss")
        compfile_name = step_input_filename.replace("srctype", det+"_comparison_pathloss")

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

    # If all tests passed then pytest will be marked as PASSED, else it will be FAILED
    FINAL_TEST_RESULT = False
    for t in total_test_result:
        if t == "FAILED":
            FINAL_TEST_RESULT = False
            break
        else:
            FINAL_TEST_RESULT = True

    if FINAL_TEST_RESULT:
        msg = "\n *** Final result for path_loss test will be reported as PASSED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "All slits PASSED path_loss test."
    else:
        msg = "\n *** Final result for path_loss test will be reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "One or more slits FAILED path_loss test."

    # end the timer
    pathloss_end_time = time.time() - pathtest_start_time
    if pathloss_end_time > 60.0:
        pathloss_end_time = pathloss_end_time/60.0  # in minutes
        pathloss_tot_time = "* Script msa_uni.py took ", repr(pathloss_end_time)+" minutes to finish."
        if pathloss_end_time > 60.0:
            pathloss_end_time = pathloss_end_time/60.  # in hours
            pathloss_tot_time = "* Script msa_uni.py took ", repr(pathloss_end_time)+" hours to finish."
    else:
        pathloss_tot_time = "* Script msa_uni.py took ", repr(pathloss_end_time)+" seconds to finish."
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

