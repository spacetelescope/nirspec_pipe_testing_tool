import time
import os
from astropy import wcs
import numpy as np
from astropy.io import fits
import argparse
import sys
import io

import jwst
from gwcs import wcstools
from jwst import datamodels

from . import auxiliary_functions as auxfunc
from astropy.visualization import (ImageNormalize, AsinhStretch)
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use("macOSX")

"""
This script tests the FS pipeline pathloss step output for a Uniform Source.
"""

# HEADER
__author__ = "T King & M Pena-Guerrero"
__version__ = "1.5"

# HISTORY
# Oct 19, 2019 - Version 1.0: initial version started
# Dec 2019, - Version 1.1: initial version passes test
# Feb 26 2020 - Version 1.2: mostly pep8 compliant
# June 8, 2020 - Version 1.3: Added changes to be able to run within NPTT
# September 25, 2020 - Version 1.4: Added option to use either data model or fits file as input for the test, and
#                      the option to provide an extract_2d file to the function
# January 2021 - Version 1.5: Implemented option to use datamodels instead of fits files as input


# put in auxiliary fxns
def get_ps_uni_extensions(fits_file_name, is_point_source):
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


def pathtest(step_input_filename, reffile, comparison_filename,
             writefile=True, show_figs=False, save_figs=True,
             threshold_diff=1.0e-7, debug=False):
    """
    This function calculates the difference between the pipeline and
    calculated pathloss values.
    Args:
        step_input_filename: str, full path name of sourcetype output fits file
        reffile: str, path to the pathloss FS reference fits file
        comparison_filename: str, path to pipeline-generated pathloss fits file
        writefile: boolean, if True writes the fits files of
                   calculated flat and difference images
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, whether to save the plots or not
        plot_name: string, desired name. If not given, plot has default name
        threshold_diff: float, threshold difference between
                        pipeline output and comparison file
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
    # For the moment, the pipeline is using the wrong reference file for slit 400A1, so read file that
    # re-processed with the right reference file and open corresponding data model
    # BUT we are skipping this since this is not in the released candidate of the pipeline
    """
    pathloss_400a1 = step_input_filename.replace("srctype.fits", "pathloss_400A1.fits")
    pathloss_pipe_400a1 = datamodels.open(pathloss_400a1)
    """
    if debug:
        print('got comparison datamodel!')
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

    print('Checking files exist & obtaining datamodels, takes a few mins...')
    # get the comparison data model
    if isinstance(comparison_filename, str):
        if os.path.isfile(comparison_filename):
            if debug:
                print('Comparison file does exist.')
        else:
            result_msg = 'Comparison file does NOT exist. Skipping pathloss test.'
            log_msgs.append(result_msg)
            result = 'skip'
            return result, result_msg, log_msgs

        pathloss_pipe = datamodels.open(comparison_filename)
    else:
        pathloss_pipe = comparison_filename
    # For the moment, the pipeline is using the wrong reference file for slit 400A1, so read file that
    # re-processed with the right reference file and open corresponding data model
    if os.path.isfile(step_input_filename.replace("srctype.fits", "pathloss_400A1.fits")):
        pathloss_400a1 = step_input_filename.replace("srctype.fits", "pathloss_400A1.fits")
        pathloss_pipe_400a1 = datamodels.open(pathloss_400a1)
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

    sltname_list = ["S200A1", "S200A2", "S400A1", "S1600A1"]
    if det == "NRS2":
        sltname_list.append("S200B1")

    # but check if data is BOTS
    if exptype == "NRS_BRIGHTOBJ":
        sltname_list = ["S1600A1"]

    # get all the science extensions
    is_point_source = False
    ps_uni_ext_list = get_ps_uni_extensions(reffile, is_point_source)

    slit_val = 0
    for slit, pipe_slit in zip(pl.slits, pathloss_pipe.slits):
        slit_val = slit_val+1
        slit_id = slit.name
        #if slit_id == 'S400A1':
        #    continue
        continue_pathloss_test = False
        if exptype == "NRS_BRIGHTOBJ":
            if isinstance(step_input_filename, str):
                extract2d_wcs_file = step_input_filename.replace("srctype.fits", "extract_2d.fits")
                model = datamodels.MultiSlitModel(extract2d_wcs_file)
            else:
                model = pl
            slit = model
            continue_pathloss_test = True
        else:
            for slit_in_MultiSlitModel in pl.slits:
                if slit_in_MultiSlitModel.name == slit_id:
                    slit = slit_in_MultiSlitModel
                    continue_pathloss_test = True
                    break

        if not continue_pathloss_test:
            continue
        else:
            try:
                if is_point_source is True:
                    ext = ps_uni_ext_list[0][slit_id]
                    print("Retrieved point source extension")
                elif is_point_source is False:
                    ext = ps_uni_ext_list[1][slit_id]
                    print("Retrieved extended source extension for {}".format(slit_val))
            except KeyError:
                # gets index associted with slit if issue above
                ext = sltname_list.index(slit_id)
                print("Unable to retrieve extension.")

        wcs_obj = slit.meta.wcs

        # get the wavelength
        x, y = wcstools.grid_from_bounding_box(wcs_obj.bounding_box,
                                               step=(1, 1), center=True)
        ra, dec, wave = wcs_obj(x, y)
        wave_sci = wave * 10**(-6)   # microns --> meters

        # adjustments for S400A1
        if slit_id == "S400A1":
            if is_point_source:
                ext = 1
            else:
                ext = 3
                print("Got uniform source extension frome extra reference file")
            reffile2use = "jwst-nirspec-a400.plrf.fits"
        else:
            reffile2use = reffile

        print("Using reference file {}".format(reffile2use))
        plcor_ref_ext = fits.getdata(reffile2use, ext)
        hdul = fits.open(reffile2use)

        plcor_ref = hdul[1].data
        w = wcs.WCS(hdul[1].header)

        w1, y1, x1 = np.mgrid[:plcor_ref.shape[0], : plcor_ref.shape[1],
                              :plcor_ref.shape[2]]
        slitx_ref, slity_ref, wave_ref = w.all_pix2world(x1, y1, w1, 0)

        previous_sci = slit.data
        if slit_id == 'S400A1':
            if pathloss_pipe_400a1 is not None:
                for pipe_slit_400a1 in pathloss_pipe_400a1.slits:
                    if pipe_slit_400a1.name == "S400A1":
                        comp_sci = pipe_slit_400a1.data
                        pipe_correction = pipe_slit_400a1.pathloss
                        break
                    else:
                        continue
        else:
            comp_sci = pipe_slit.data
            pipe_correction = pipe_slit.pathloss_uniform
        if len(pipe_correction) == 0:
            print("Pipeline pathloss correction in datamodel is empty. Skipping testing this slit.")
            continue

        # set up generals for all the plots
        font = {'weight': 'normal',
                'size': 7}
        matplotlib.rc('font', **font)

        corr_vals = np.interp(wave_sci, wave_ref[:, 0, 0], plcor_ref_ext)
        corrected_array = previous_sci/corr_vals

        # Plots:
        if save_figs:
            step_input_filepath = step_input_filename.replace(".fits", "")
        # my correction values
        fig = plt.figure()
        plt.subplot(321)
        norm = ImageNormalize(corr_vals)
        plt.imshow(corr_vals, norm=norm, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Calculated Correction')
        plt.colorbar()
        # pipe corerction
        plt.subplot(322)
        norm = ImageNormalize(pipe_correction)
        plt.imshow(pipe_correction, norm=norm, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Pathloss Correction Comparison')
        plt.colorbar()
        # residuals (pipe correction - my correction)
        corr_residuals = pipe_correction - corr_vals
        plt.subplot(323)
        norm = ImageNormalize(corr_residuals)
        plt.imshow(corr_residuals, norm=norm, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Correction residuals')
        plt.colorbar()
        # pipe science data before
        plt.subplot(324)
        norm = ImageNormalize(previous_sci)
        plt.imshow(previous_sci, norm=norm, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Normalized pipeline science data before pathloss')
        plt.colorbar()
        # pipe science data after
        plt.subplot(325)
        norm = ImageNormalize(comp_sci)
        plt.imshow(comp_sci, norm=norm, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.title('Normalized pipeline science data after pathloss')
        plt.colorbar()
        # pipe science data after pathloss
        plt.subplot(326)
        norm = ImageNormalize(corrected_array)
        plt.imshow(corrected_array, norm=norm, aspect=10.0,
                   origin='lower', cmap='viridis')
        plt.title('My science data after pathloss')
        plt.xlabel('dispersion in pixels')
        plt.ylabel('y in pixels')
        plt.colorbar()
        fig.suptitle("FS UNI Pathloss Correction Testing for {}".format(slit_id))

        # add space between the subplots
        fig.subplots_adjust(wspace=0.9)

        # Show and/or save figures
        if show_figs:
            plt.show()
        if save_figs:
            plt_name = step_input_filepath + "_Pathloss_test_slit_" + str(slit_id) + "_FS_extended.png"
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
            outfile_ext = fits.ImageHDU(corr_vals, name=slit_id)
            outfile.append(outfile_ext)

            # this is the file to hold the image of pipeline-calculated difference values
            compfile_ext = fits.ImageHDU(corr_residuals, name=slit_id)
            compfile.append(compfile_ext)

        # Histogram
        ax = plt.subplot(212)
        plt.hist(corr_residuals[~np.isnan(corr_residuals)],
                 bins=100, range=(-0.00000013, 0.00000013))
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
        ax.text(0.73, 0.67, str_arr_stddev, transform=ax.transAxes, fontsize=7)
        plt.legend()
        plt.minorticks_on()

        # Show and/or save figures
        if save_figs:
            plt_name = step_input_filename.replace(".fits", "") + "_Pathlosstest_slitlet_" + slit_id + ".png"
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
                msg1 = " * Unable to calculate statistics because difference array has all outlier values. Test " \
                       "will be set to FAILED."
                if debug:
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
                if debug:
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
        pathloss_tot_time = "* Script FS_UNI.py took ", repr(pathloss_end_time)+" minutes to finish."
        if pathloss_end_time > 60.0:
            pathloss_end_time = pathloss_end_time/60.  # in hours
            pathloss_tot_time = "* Script FS_UNI.py took ", repr(pathloss_end_time)+" hours to finish."
    else:
        pathloss_tot_time = "* Script FS_UNI.py took ", repr(pathloss_end_time)+" seconds to finish."
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
    parser.add_argument("-e",
                        dest="extract2d_file",
                        action='store',
                        default=None,
                        help='Path and name the extract_2d file, i.e. the output of the extract_2d step')
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
    extract2d_file = args.extract2d_file
    writefile = args.writefile
    save_figs = args.save_figs
    show_figs = args.show_figs
    threshold_diff = args.threshold_diff
    debug = args.debug

    # print pipeline version
    print("\n  ** using pipeline version: ", jwst.__version__, "** \n")

    # Run the principal function of the script
    pathtest(step_input_filename, reffile, comparison_filename, extract2d_file=extract2d_file, writefile=writefile,
             show_figs=show_figs, save_figs=save_figs, threshold_diff=threshold_diff,
             debug=debug)


if __name__ == '__main__':
    sys.exit(main())

