import os
import time
import numpy as np
from astropy.io import fits
from astropy import wcs
from collections import OrderedDict

from jwst import datamodels

from astropy.visualization import (ImageNormalize, AsinhStretch)
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from gwcs import wcstools

from . import auxiliary_functions as auxfunc


"""
This script tests the pipeline bar shadow step output for MOS data. It is the scripted and generalized version of the 
jupyter notebook barshadow.ipynb written by James Muzerolle in August of 2019.
"""


# HEADER
__author__ = "M. A. Pena-Guerrero & J. Muzerolle"
__version__ = "1.0"

# HISTORY
# Nov 2019 - Version 1.0: initial version completed



def run_barshadow_tests(plfile, bsfile, barshadow_threshold_diff=0.05, save_final_figs=False, show_final_figs=False,
                        save_intermediary_figs=False, show_intermediary_figs=False, write_barshadow_files=False,
                        debug=False):
    """

    Args:
        plfile: string, 2D spectra output prior to the bar shadow step (e.g., extract_2d or pathloss product)
        bsfile: string, read in 2D spectra output from the bar shadow step
        barshadow_threshold_diff: float, this value comes from the document ESA-JWST-SCI-NRS-TN-2016-016.pdf, it is
                                    an arbitrary error of the reference file of 0.0025 absolute error or 5% relative
                                    error (no justification provided)
        save_final_figs: boolean, if True the final figures with corresponding histograms will be saved
        show_final_figs: boolean, if True the final figures with corresponding histograms will be shown
        save_intermediary_figs: boolean, if True the intermediary figures with corresponding histograms will be saved
        show_intermediary_figs: boolean, if True the intermediary figures with corresponding histograms will be shown
        debug: boolean

    Returns:

    """

    # start the list of messages that will be added to the log file
    log_msgs = []

    # start the timer
    barshadow_test_start_time = time.time()

    # read in 2D spectra output prior to the bar shadow step (e.g., extract_2d or pathloss product)
    print('Checking if files exist and obtaining datamodels, this takes a few minutes...')
    if os.path.isfile(plfile):
        if debug:
            print('Extract_2d file does exist.')
    else:
        result_msg = 'Extract_2d file does NOT exist. Barshadow test will be skipped.'
        log_msgs.append(result_msg)
        result = 'skip'
        return result, result_msg, log_msgs

    # get the data model
    pl = datamodels.open(plfile)
    if debug:
        print('got extract_2d datamodel!')

    # read in 2D spectra output from the bar shadow step
    if os.path.isfile(bsfile):
        if debug:
            print('Bar shadow file does exist.')
    else:
        result_msg = 'Barshadow file does NOT exist. Barshadow test will be skipped.'
        log_msgs.append(result_msg)
        result = 'skip'
        return result, result_msg, log_msgs

    bs = datamodels.open(bsfile)
    if debug:
        print('got barshadow datamodel!')

    # list to determine if pytest is passed or not
    total_test_result = OrderedDict()

    if write_barshadow_files:
        # create the fits list to hold the image of the correction values
        hdu0 = fits.PrimaryHDU()
        outfile = fits.HDUList()
        outfile.append(hdu0)

        # create the fits list to hold the image of the comparison values
        hdu0 = fits.PrimaryHDU()
        complfile = fits.HDUList()
        complfile.append(hdu0)

    # loop over the slitlets in both files
    print('Looping over open slitlets...')
    for plslit, bsslit in zip(pl.slits, bs.slits):
        # check that slitlet name of the data from the pathloss or extract_2d and the barshadow datamodels are the same
        slit_id = bsslit.name
        print('Working with slitlet ', slit_id)
        if plslit.name == bsslit.name:
            msg = 'Slitlet name in fits file previous to barshadow and in barshadow output file are the same.'
            log_msgs.append(msg)
            print(msg)
        else:
            msg = '* Missmatch of slitlet names in fits file previous to barshadow and in barshadow output file. Skipping test.'
            result = 'skip'
            log_msgs.append(msg)
            return result, msg, log_msgs


        # obtain the data from the pathloss or extract_2d and the barshadow datamodels
        plsci = plslit.data
        bssci = bsslit.data

        if debug:
            print('plotting the data for both input files...')

        # set up generals for all the plots
        font = {  # 'family' : 'normal',
                'weight': 'normal',
                'size': 16}
        matplotlib.rc('font', **font)

        plt.figure(figsize=(12, 10))
        # Top figure
        plt.subplot(211)
        norm = ImageNormalize(plsci,vmin=0.,vmax=500.,stretch=AsinhStretch())
        plt.imshow(plsci, norm=norm, aspect=10.0, origin='lower', cmap='viridis')
        plt.title('Normalized science data before barshadow step for slitlet '+slit_id)
        # Bottom figure
        plt.subplot(212)
        norm=ImageNormalize(bssci,vmin=0.,vmax=500.,stretch=AsinhStretch())
        plt.imshow(bssci,norm=norm,aspect=10.0,origin='lower',cmap='viridis')
        plt.title('Normalized barshadow science data for slitlet '+slit_id)
        # Show and/or save figures
        file_path = bsfile.replace(os.path.basename(bsfile), "")
        file_basename = os.path.basename(bsfile.replace("_barshadow.fits", ""))
        if save_intermediary_figs:
            t = (file_basename, "Barshadowtest_NormSciData_slitlet" + slit_id + ".jpg")
            plt_name = "_".join(t)
            plt_name = os.path.join(file_path, plt_name)
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_intermediary_figs:
            plt.show()
        plt.close()

        # calculate spatial profiles for both products
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9,9))
        plt.subplots_adjust(hspace=0.5)
        fig.subplots_adjust(wspace=0.6)
        point1 = [355, 375]
        plprof1 = np.median(plsci[:,point1[0]:point1[1]],1)
        bsprof1 = np.median(bssci[:,point1[0]:point1[1]],1)
        # only use pixels that are not NaN
        x1 = np.squeeze(np.nonzero(~np.isnan(bsprof1)))
        ax1.plot(x1,plprof1[x1])
        ax1.set_title('Before barshadow array slice 1')
        ax1.set_xlabel('x (pixels)')
        ax1.set_ylabel('y (pixels)')
        if debug:
            print('ax1 std_dev/mean = ', np.nanstd(plprof1[x1])/np.nanmean(plprof1[x1]))

        point2 = [1190, 1210]
        plprof2 = np.median(plsci[:,point2[0]:point2[1]],1)
        bsprof2 = np.median(bssci[:,point2[0]:point2[1]],1)
        x2 = np.squeeze(np.nonzero(~np.isnan(bsprof2)))
        ax2.plot(x2,plprof2[x2])
        ax2.set_title('Before barshadow array slice 2')
        ax2.set_xlabel('x (pixels)')
        ax2.set_ylabel('y (pixels)')
        if debug:
            print('ax2 std_dev/mean = ', np.nanstd(plprof2[x2])/np.nanmean(plprof2[x2]))
        
        ax3.plot(x1,bsprof1[x1])
        ax3.set_title('Barshadow array slice 1')
        ax3.set_xlabel('x (pixels)')
        ax3.set_ylabel('y (pixels)')
        if debug:
            print('ax3 std_dev/mean = ', np.nanstd(bsprof1)/np.nanmean(bsprof1[x1]))
        
        ax4.plot(x2,bsprof2[x2])
        if debug:
            print('ax4 std_dev/mean = ', np.nanstd(bsprof2)/np.nanmean(bsprof2[x2]))
        ax4.set_title('Barshadow array slice 2')
        ax4.set_xlabel('x (pixels)')
        ax4.set_ylabel('y (pixels)')

        fig.suptitle('Spatial profiles before correction for slitlet '+slit_id, fontsize=20)

        # Show and/or save figures
        if save_intermediary_figs:
            t = (file_basename, "Barshadowtest_SpatialProfilesBe4correction_slitlet" + slit_id + ".jpg")
            plt_name = "_".join(t)
            plt_name = os.path.join(file_path, plt_name)
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_intermediary_figs:
            plt.show()
        plt.close()


        ### compare pipeline correction values with independent calculation

        # get the bar shadow corrections from the step product
        bscor_pipe = bsslit.barshadow

        # get correction from independent calculation
        msg = 'Calculating barshadow correction...'
        log_msgs.append(msg)
        print(msg)
        # Create x, y indices using the Trace WCS
        x, y = wcstools.grid_from_bounding_box(bsslit.meta.wcs.bounding_box, step=(1, 1))
        if debug:
            print('x = ', x)

        # derive the slity_y values per pixel
        wcsobj = bsslit.meta.wcs
        det2slit = wcsobj.get_transform('detector','slit_frame')
        bsslitx, bsslity, bswave = det2slit(x,y)
        # scale the slit_y values by 1.15 to take into account the shutter pitch
        bsslity = bsslity/1.15

        # compute bar shadow corrections independently, given the wavelength and slit_y from the data model
        # get the reference file (need the mos1x1 for this internal lamp case, where each shutter was extracted separately)
        #if bsslit.shutter_state == 'x':
        ref_file = '/grp/jwst/wit4/nirspec/CDP3/05_Other_Calibrations/5.3_BarShadow/referenceFilesBS-20160401/jwst-nirspec-mos1x1.bsrf.fits'
        if bsslit.shutter_state == '1':
            ref_file = '/grp/jwst/wit4/nirspec/CDP3/05_Other_Calibrations/5.3_BarShadow/referenceFilesBS-20160401/jwst-nirspec-mos1x3.bsrf.fits'
        if debug:
            '''    shutter_state : str ----- ``Slit.shutter_state`` attribute - a combination of
                                    possible values: ``1`` - open shutter, ``0`` - closed shutter, ``x`` - main shutter
            '''
            print('slit.shutter_state = ', bsslit.shutter_state)
        msg = 'Reference file used for barshadow calculation: '+ref_file
        log_msgs.append(msg)
        print(msg)
        hdul = fits.open(ref_file)
        bscor_ref = hdul[1].data
        w = wcs.WCS(hdul[1].header)
        y1, x1 = np.mgrid[:bscor_ref.shape[0], : bscor_ref.shape[1]]
        lam_ref, slity_ref = w.all_pix2world(x1,y1,0)

        # for slit wcs, interpolate over the reference file values
        lam_ref = lam_ref.reshape(bscor_ref.size)
        slity_ref = slity_ref.reshape(bscor_ref.size)
        pixels_ref = np.column_stack((lam_ref,slity_ref))
        bscor_ref = bscor_ref.reshape(bscor_ref.size)
        bswave_ex = bswave.reshape(bswave.size)
        indxs = ~np.isnan(bswave_ex)
        bsslity_ex = bsslity.reshape(bsslity.size)
        xyints = np.column_stack((bswave_ex[indxs],bsslity_ex[indxs]))
        bscor = np.empty(bswave_ex.size)
        bscor[:] = np.nan
        bscor[indxs] = griddata(pixels_ref, bscor_ref, xyints, method='linear')
        bscor = bscor.reshape(bswave.shape[0],bswave.shape[1])
        if debug:
            print('bscor.shape = ', bscor.shape)
        msg = 'Calculation of barshadow correction done.'
        log_msgs.append(msg)
        print(msg)

        shutter_status = bsslit.shutter_state
        fi = shutter_status.find('x')
        if debug:
            print('fi = ', fi)
        nax2 = hdul[1].header['NAXIS2']
        cv1 = hdul[1].header['CRVAL1']
        cd1 = hdul[1].header['CDELT1']
        cd2 = hdul[1].header['CDELT2']
        shutter_height = 1./cd2
        fi2 = nax2-shutter_height*(1+fi)
        if debug:
            print('nax2, fi2, shutter_height:', nax2, fi2, shutter_height)
        yrow = fi2 + bsslity*shutter_height
        wcol = (bswave-cv1)/cd1
        #print(yrow[9,1037],wcol[9,1037])
        if debug:
            print('np.shape(yrow)=', np.shape(yrow))
        point3 = [10, np.shape(yrow)[1]-50]
        print(yrow[point3[0], point3[1]],wcol[point3[0], point3[1]])

        fig = plt.figure(figsize=(12, 10))
        # Top figure
        plt.subplot(211)
        plt.imshow(bscor,vmin=0.,vmax=1.,aspect=10.0,origin='lower',cmap='viridis')
        plt.title('Calculated Correction')
        plt.colorbar()
        # Bottom figure
        plt.subplot(212)
        plt.imshow(bscor_pipe,vmin=0.,vmax=1.,aspect=10.0,origin='lower',cmap='viridis')
        plt.title('Pipeline Correction')
        plt.colorbar()

        fig.suptitle('Barshadow correction comparison for slitlet '+slit_id, fontsize=20)

        # Show and/or save figures
        if save_intermediary_figs:
            t = (file_basename, "Barshadowtest_CorrectionComparison_slitlet" + slit_id + ".jpg")
            plt_name = "_".join(t)
            plt_name = os.path.join(file_path, plt_name)
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_intermediary_figs:
            plt.show()
        plt.close()

        if debug:
            #print('bscor_pipe[9,1037],bswave[9,1037],bsslity[9,1037],bscor[9,1037]: ',
            #      bscor_pipe[9,1037],bswave[9,1037],bsslity[9,1037],bscor[9,1037])
            print('bscor_pipe[point3[0], point3[1]],bswave[point3[0], point3[1]],bsslity[point3[0], point3[1]],bscor[point3[0], point3[1]]: ',
                  bscor_pipe[point3[0], point3[1]],bswave[point3[0], point3[1]],bsslity[point3[0], point3[1]],bscor[point3[0], point3[1]])


        print('Creating final barshadow test plot...')
        reldiff = (bscor_pipe-bscor)/bscor
        if debug:
            print('np.nanmean(reldiff),np.nanstd(reldiff) : ', np.nanmean(reldiff),np.nanstd(reldiff))
        fig = plt.figure(figsize=(12, 10))
        # Top figure - 2D plot
        plt.subplot(211)
        plt.imshow(reldiff,vmin=-0.01,vmax=0.01,aspect=10.0,origin='lower',cmap='viridis')
        plt.colorbar()
        plt.title('Relative differences')
        plt.xlabel('x (pixels)')
        plt.ylabel('y (pixels)')
        # Bottom figure - histogram
        ax = plt.subplot(212)
        plt.hist(reldiff[~np.isnan(reldiff)],bins=100,range=(-0.1,0.1))
        plt.xlabel('(Pipeline_correction - Calculated_correction) / Calculated_correction')
        plt.ylabel('N')
        # add vertical line at mean and median
        nanind = np.isnan(reldiff)  # get all the nan indexes
        notnan = ~nanind  # get all the not-nan indexes
        arr_mean = np.mean(reldiff[notnan])
        arr_median = np.median(reldiff[notnan])
        arr_stddev = np.std(reldiff[notnan])
        plt.axvline(arr_mean, label="mean = %0.3e" % (arr_mean), color="g")
        plt.axvline(arr_median, label="median = %0.3e" % (arr_median), linestyle="-.", color="b")
        str_arr_stddev = "stddev = {:0.3e}".format(arr_stddev)
        ax.text(0.73, 0.67, str_arr_stddev, transform=ax.transAxes, fontsize=16)
        plt.legend()
        plt.minorticks_on()

        fig.suptitle('Barshadow correction relative differences for slitlet '+slit_id, fontsize=20)

        # Show and/or save figures
        if save_final_figs:
            t = (file_basename, "Barshadowtest_RelDifferences_slitlet" + slit_id + ".jpg")
            plt_name = "_".join(t)
            plt_name = os.path.join(file_path, plt_name)
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_final_figs:
            plt.show()
        plt.close()

        # Determine if median test is passed
        slitlet_test_result_list = []
        tested_quantity = 'barshadow_correction'
        stats = auxfunc.print_stats(reldiff[notnan], tested_quantity, barshadow_threshold_diff, abs=False, return_percentages=True)
        _, stats_print_strings, percentages = stats
        result = auxfunc.does_median_pass_tes(arr_median, barshadow_threshold_diff)
        slitlet_test_result_list.append({tested_quantity: result})
        for line in stats_print_strings:
            log_msgs.append(line)
        msg = " * Result of median test for slit "+slit_id+": "+result+"\n"
        print(msg)
        log_msgs.append(msg)

        tested_quantity = "percentage_greater_3threshold"
        result = auxfunc.does_median_pass_tes(percentages[1], 10)
        slitlet_test_result_list.append({tested_quantity: result})
        msg = " * Result of number of points greater than 3*threshold greater than 10%: "+result+"\n"
        print(msg)
        log_msgs.append(msg)

        tested_quantity = "percentage_greater_5threshold"
        result = auxfunc.does_median_pass_tes(percentages[2], 10)
        slitlet_test_result_list.append({tested_quantity: result})
        msg = " * Result of number of points greater than 5*threshold greater than 10%: "+result+"\n"
        print(msg)
        log_msgs.append(msg)

        # Make plots of normalized corrected data
        corrected = plsci/bscor
        plt.figure(figsize=(12, 10))
        norm=ImageNormalize(corrected,vmin=0.,vmax=500.,stretch=AsinhStretch())
        plt.imshow(corrected, norm=norm, aspect=10.0, origin='lower', cmap='viridis')
        plt.title('Normalized data before barshadow step with correction applied')
        plt.xlabel('Sci_data_before_barshadow / barshadow_calculated_correction')
        plt.ylabel('Normalized data')
        # Show and/or save figures
        if save_intermediary_figs:
            t = (file_basename, "Barshadowtest_CorrectedData_slitlet" + slit_id + ".jpg")
            plt_name = "_".join(t)
            plt_name = os.path.join(file_path, plt_name)
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_intermediary_figs:
            plt.show()
        plt.close()


        # calculate spatial profiles for both products
        fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(19,9))
        prof = np.median(corrected[:,point1[0]:point1[1]],1)
        x = np.arange(corrected.shape[0])
        ax1.plot(x,prof)
        ax1.set_title('Before barshadow array slice 1')
        ax1.set_xlabel('x (pixels)')
        ax1.set_ylabel('y (pixels)')
        if debug:
            print('np.nanstd(prof)/np.nanmean(prof) = ', np.nanstd(prof)/np.nanmean(prof))
        prof = np.median(corrected[:,point2[0]:point2[1]],1)
        x = np.arange(corrected.shape[0])
        ax2.plot(x,prof)
        ax2.set_title('Before barshadow array slice 2')
        ax2.set_xlabel('x (pixels)')
        ax2.set_ylabel('y (pixels)')
        if debug:
            print('np.nanstd(prof)/np.nanmean(prof) = ', np.nanstd(prof)/np.nanmean(prof))
        fig.suptitle('Corrected spatial profiles for slitlet '+slit_id, fontsize=20)
        # Show and/or save figures
        if save_intermediary_figs:
            t = (file_basename, "Barshadowtest_CorrectedSpatialProfiles_slitlet" + slit_id + ".jpg")
            plt_name = "_".join(t)
            plt_name = os.path.join(file_path, plt_name)
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_intermediary_figs:
            plt.show()
        plt.close()

        # store tests results in the total dictionary
        total_test_result[slit_id] = slitlet_test_result_list

        # create fits file to hold the calculated correction for each slitlet
        if write_barshadow_files:
            # this is the file to hold the image of the correction values
            outfile_ext = fits.ImageHDU(corrected, name=slit_id)
            outfile.append(outfile_ext)

            # this is the file to hold the image of pipeline-calculated difference values, the comparison
            complfile_ext = fits.ImageHDU(reldiff, name=slit_id)
            complfile.append(complfile_ext)

            # the file is not yet written, indicate that this slit was appended to list to be written
            msg = "Extension corresponing to slitlet "+slit_id+" appended to list to be written into calculated and comparison fits files."
            print(msg)
            log_msgs.append(msg)


    if debug:
        print('total_test_result = ', total_test_result)

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
        result_msg = "\n *** Final result for assign_wcs test will be reported as PASSED *** \n"
        print(result_msg)
        log_msgs.append(result_msg)
    else:
        result_msg = "\n *** Final result for assign_wcs test will be reported as FAILED *** \n"
        print(result_msg)
        log_msgs.append(result_msg)

    # end the timer
    barshadow_test_end_time = time.time() - barshadow_test_start_time
    if barshadow_test_end_time > 60.0:
        barshadow_test_end_time = barshadow_test_end_time/60.0  # in minutes
        barshadow_test_tot_time = "* Script flattest_mos.py took ", repr(barshadow_test_end_time)+" minutes to finish."
        if barshadow_test_end_time > 60.0:
            barshadow_test_end_time = barshadow_test_end_time/60.  # in hours
            barshadow_test_tot_time = "* Script flattest_mos.py took ", repr(barshadow_test_end_time)+" hours to finish."
    else:
        barshadow_test_tot_time = "* Script flattest_mos.py took ", repr(barshadow_test_end_time)+" seconds to finish."
    print(barshadow_test_tot_time)
    log_msgs.append(barshadow_test_tot_time)

    return FINAL_TEST_RESULT, result_msg, log_msgs




if __name__ == '__main__':

    # This is a simple test of the code
    pipeline_path = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline"

    # input parameters that the script expects
    plfile = pipeline_path+'/testing_data/MOS/G140M_LINE1/NRSV84600010001_G140M_LINE1_NRS2_extract_2d.fits'
    bsfile = pipeline_path+'/testing_data/MOS/G140M_LINE1/NRSV84600010001_G140M_LINE1_NRS2_barshadow.fits'

    # name of the output images
    writefile = True

    # This value comes from the document ESA-JWST-SCI-NRS-TN-2016-016.pdf, it is an arbitrary error of the
    # reference file of 0.0025 absolute error or 5% relative error (no justification)
    barshadow_threshold_diff = 0.0025

    # Run the principal function of the script
    run_barshadow_tests(plfile, bsfile, barshadow_threshold_diff=barshadow_threshold_diff, save_final_figs=True,
                        show_final_figs=False, save_intermediary_figs=True, show_intermediary_figs=False,
                        write_barshadow_files=False, debug=True)



