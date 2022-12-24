import time
import os
import argparse
import sys
import numpy as np
from astropy.io import fits
from glob import glob
from copy import deepcopy

from gwcs import wcstools
from jwst import datamodels

from . import auxiliary_functions as auxfunc

"""
This script tests the pipeline flat field step output for MOS data. It is the python version of the IDL script
(with the same name) written by James Muzerolle.
"""

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "2.8"


# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Completely changed script to use the datamodel instead of the compute_world_coordinates
#                         script, and added new routines for statistics calculations.
# Jun 2018 - Version 2.1: Changed extension numbers for the name of the extension in the D-, F-, and S-flats.
# Jun 2018 - Version 2.2: Removed function reverse_cols because it was not behaving as expected.
# Apr 2019 - Version 2.3: Implemented logging capability.
# May 2019 - Version 2.4: Implemented images of the residuals.
# Jun 2019 - Version 2.5: Updated name of interpolated flat to be the default pipeline name for this file.
# Sep 2019 - Version 2.6: Updated line to call model for SlitModel to work correctly with pipeline changes.
# Jan 2021 - Version 2.7: Implemented option to run with object instead of input fits file.
# Dec 2022 - Version 2.8: Fixed code to read new post-commissioning reference files in CRDS format.


def get_slit_wavdat(fastvartable, slit_id):
    """
    Extract the data for the appropriate slit from the fits table.
    Args:
        fastvartable: astropy table
        slit_id: string
    Returns:
        fastvar_wav: array, wavelengths for fast variation
        fastvar_dat: array
    """
    # obtain the table index for corresponding to the given slit
    idx = np.where(fastvartable['slit_name'] == slit_id)
    fastvar_wav = fastvartable['wavelength'][idx]
    fastvar_dat = fastvartable['data'][idx]
    return fastvar_wav[0], fastvar_dat[0]


def flattest(step_input_filename, dflat_path, sflat_path, fflat_path, writefile=True,
             show_figs=True, save_figs=False, interpolated_flat=None, threshold_diff=1.0e-7,
             output_directory=None, debug=False):
    """
    This function calculates the difference between the pipeline and the calculated flat field values.
    The functions uses the output of the compute_world_coordinates.py script.

    Args:
        step_input_filename: str, name of the output fits file from the 2d_extract step (with full path)
        dflat_path: str, path of where the D-flat reference fits files
        sflat_path: str, path of where the S-flat reference fits files
        fflat_path: str, path of where the F-flat reference fits files
        writefile: boolean, if True writes the fits files of the calculated flat and difference images
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots (the 3 plots can be saved or not independently with the function call)
        interpolated_flat: string, name of the on-the-fly interpolated pipeline flat
        threshold_diff: float, threshold difference between pipeline output and ESA file
        output_directory: None or string, path to the output_directory where to save the plots and output files
        debug: boolean, if true a series of print statements will show on-screen

    Returns:
        - 1 plot, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to threshold.
        - log_msgs: list, all print statements are captured in this variable

    """
    log_msgs = []

    # start the timer
    flattest_start_time = time.time()

    # get basic info from the input file
    if isinstance(step_input_filename, str):
        msg = 'step_input_filename=' + step_input_filename
        print(msg)
        log_msgs.append(msg)
        if "wavecorr" not in step_input_filename:
            wcs_file = step_input_filename.replace("_flat_field.fits", "_wavecorr.fits")
        else:
            wcs_file = step_input_filename
        model = datamodels.MultiSlitModel(wcs_file)

        # paths to save plots and files
        file_basename = os.path.basename(step_input_filename.replace(".fits", ""))
        if output_directory is not None:
            file_path = output_directory
        else:
            file_path = step_input_filename.replace(file_basename, "")

    else:
        file_path = output_directory
        model = step_input_filename
        file_basename = ''

    # read in the on-the-fly flat image
    if interpolated_flat is None:
        flatfile = step_input_filename.replace("flat_field.fits", "interpolatedflat.fits")
    else:
        flatfile = interpolated_flat
    # get all the science extensions in the flatfile
    sci_ext_list = auxfunc.get_sci_extensions(flatfile)

    # get basic info from model
    if isinstance(model, list):  # this was added for the validation notebooks to work
        model = model[0]
    det = model.meta.instrument.detector
    grat = model.meta.instrument.grating
    filt = model.meta.instrument.filter
    lamp = model.meta.instrument.lamp_state
    exp_type = model.meta.exposure.type.upper()

    msg = "flat_field_file  -->     Grating:" + grat + "   Filter:" + filt + "   LAMP:" + lamp
    print(msg)
    log_msgs.append(msg)

    # print the info about the reference files used:
    with datamodels.open(step_input_filename) as pipe_flat_field_mdl:
        msg0 = "\n * FOR COMPARISON PURPOSES, for file " + step_input_filename
        msg1 = "    DATE-OBS = " + str(pipe_flat_field_mdl.meta.observation.date)
        msg2 = "    Pipeline CRDS context: " + str(pipe_flat_field_mdl.meta.ref_file.crds.context_used)
        msg3 = "    Pipeline ref d-flat used:  " + str(pipe_flat_field_mdl.meta.ref_file.dflat.name)
        msg4 = "    Pipeline ref s-flat used:  " + str(pipe_flat_field_mdl.meta.ref_file.sflat.name)
        msg5 = "    Pipeline ref f-flat used:  " + str(pipe_flat_field_mdl.meta.ref_file.fflat.name) + "\n"
        print(msg0)
        print(msg1)
        print(msg2)
        print(msg3)
        print(msg4)
        print(msg5)
        log_msgs.append(msg0)
        log_msgs.append(msg1)
        log_msgs.append(msg2)
        log_msgs.append(msg3)
        log_msgs.append(msg4)
        log_msgs.append(msg5)

    # get the reference files
    # D-Flat
    if ".fits" not in dflat_path:
        # get all the fits files and find the appropriate detector
        fflatfiles = glob(os.path.join(dflat_path, "*.fits"))
        for ff in fflatfiles:
            if det.lower() in ff.lower():
                if debug:
                    print('Is detector ', det, ' in file ', ff)
                dfile = ff
                break
            else:
                result_msg = "Test skiped because there is no flat correspondence for the detector {} in directory {}".format(det,
                             dflat_path)
                print(msg)
                median_diff = "skip"
                return median_diff, result_msg, log_msgs
    else:
        dfile = dflat_path
    msg0 = " * This flat test is using the following reference files "
    msg1 = "    D-flat: " + dfile
    print(msg0)
    print(msg1)
    log_msgs.append(msg0)
    log_msgs.append(msg1)
    with fits.open(dfile) as dfile_hdu:
        dfim = dfile_hdu["SCI"].data
        dfimdq = dfile_hdu["DQ"].data
        dfrqe = dfile_hdu["FAST_VARIATION"].data
        dfhdr_sci = dfile_hdu["SCI"].header
    ns = np.shape(dfim)
    naxis3 = dfhdr_sci["NAXIS3"]
    if debug:
        print('np.shape(dfim) =', np.shape(dfim))
        print('np.shape(dfimdq) =', np.shape(dfimdq))

    # get the wavelength values
    dfwave = np.array([])
    for i in range(naxis3):
        t = ("PFLAT", str(i + 1))
        keyword = "_".join(t)
        dfwave = np.append(dfwave, fits.getval(dfile, keyword, 1))

    # S-flat
    mode = "FS"
    if filt == "F070LP":
        flat = "FLAT4"
    elif filt == "F100LP":
        flat = "FLAT1"
    elif filt == "F170LP":
        flat = "FLAT2"
    elif filt == "F290LP":
        flat = "FLAT3"
    elif filt == "CLEAR":
        flat = "FLAT5"
    else:
        msg = "No filter correspondence. Exiting the program."
        print(msg)
        log_msgs.append(msg)
        result_msg = "Test skiped because there is no flat correspondence for the filter in the data: {}".format(filt)
        median_diff = "skip"
        return median_diff, result_msg, log_msgs

    if ".fits" not in sflat_path:
        # get all the fits files and find the appropriate detector
        all_flatfiles = glob(os.path.join(sflat_path, "*.fits"))
        flatfiles = []
        for sflname in all_flatfiles:
            sflname_basename = os.path.basename(sflname)
            if det in sflname_basename and mode in sflname_basename:
                flatfiles.append(sflname)
        if debug:
            print('Files found in sflat path: ', len(flatfiles))
        for ff in flatfiles:
            if debug:
                print('looking for grating ', grat, ' in file ', ff)
            # look for the right configuration flat
            if grat in ff:
                if debug:
                    print('found grating, now looking for filter ', filt, ' in file ', ff)
                try:
                    if filt in ff:
                        sfile = ff
                        break
                    else:
                        result_msg = "Test skiped because there is no S-flat correspondence for configuration {} {} in directory {}".format(grat, filt, sflat_path)
                        print(msg)
                        median_diff = "skip"
                        return median_diff, result_msg, log_msgs
                except Exception as err:
                    print(err)
                    print('Trying OPAQUE as filter')
                    if 'OPAQUE' in ff:
                        sfile = ff
                        break
                    else:
                        result_msg = "Test skiped because there is no S-flat correspondence for configuration {} {} in directory {}".format(grat, 'OPAQUE', sflat_path)
                        print(msg)
                        median_diff = "skip"
                        return median_diff, result_msg, log_msgs
    else:
        sfile = sflat_path

    msg = "    S-flat: " + sfile
    print(msg)
    log_msgs.append(msg)
    with fits.open(sfile) as sfile_hdu:
        sfim = sfile_hdu["SCI"].data
        sfimdq = sfile_hdu["DQ"].data
        sffastvar = sfile_hdu["FAST_VARIATION"].data

    if debug:
        print("np.shape(sfim) = ", np.shape(sfim))
        print("np.shape(sfimdq) = ", np.shape(sfimdq))
        sf = fits.open(sfile)
        print(sf.info())
        sf.close()

    # F-Flat
    if ".fits" not in fflat_path:
        # get all the fits files and find the appropriate detector
        all_flatfiles = glob(os.path.join(fflat_path, "*.fits"))
        flatfiles = []
        for flname in all_flatfiles:
            flname_basename = os.path.basename(flname)
            if mode in flname_basename:
                flatfiles.append(flname)
        if debug:
            print('Files found in fflat path: ', len(flatfiles))
        for ff in flatfiles:
            # look for the right configuration flat
            if grat in ff:
                if debug:
                    print('found grating ', grat, ' in flat file')
                try:
                    if filt in ff:
                        if debug:
                            print('found filter ', filt, ' in flat file')
                        ffile = ff
                        break
                    else:
                        result_msg = "Test skiped because there is no F-flat correspondence for configuration {} {} in directory {}".format(grat, filt, fflat_path)
                        print(msg)
                        median_diff = "skip"
                        return median_diff, result_msg, log_msgs
                except Exception as err:
                    print(err)
                    print('Now trying OPAQUE as filter')
                    if 'OPAQUE' in ff:
                        ffile = ff
                        break
                    else:
                        result_msg = "Test skiped because there is no F-flat correspondence for configuration {} {} in directory {}".format(grat, 'OPAQUE', fflat_path)
                        print(msg)
                        median_diff = "skip"
                        return median_diff, result_msg, log_msgs
    else:
        ffile = fflat_path

    msg = "    F-flat: " + ffile
    print(msg)
    log_msgs.append(msg)
    with fits.open(ffile) as ffile_hdu:
        fferr = ffile_hdu["ERR"].data
        fffastvar = ffile_hdu["FAST_VARIATION"].data

    # now go through each pixel in the test data

    if writefile:
        # create the fits list to hold the calculated flat values for each slit
        hdu0 = fits.PrimaryHDU()
        outfile = fits.HDUList()
        outfile.append(hdu0)

        # create the fits list to hold the image of pipeline-calculated difference values
        hdu0 = fits.PrimaryHDU()
        complfile = fits.HDUList()
        complfile.append(hdu0)

    # list to determine if pytest is passed or not
    total_test_result = []

    # open the flatfile
    flatfile_hdu = fits.open(flatfile)

    # loop over the slits
    sltname_list = ["S200A1", "S200A2", "S400A1", "S1600A1"]
    msg = "\nNow looping through the slits. This may take a while... "
    print(msg)
    log_msgs.append(msg)
    if det == "NRS2":
        sltname_list.append("S200B1")

    # but check if data is BOTS
    if exp_type == "NRS_BRIGHTOBJ":
        sltname_list = ["S1600A1"]

    # do the loop over the slits
    for slit_id in sltname_list:
        continue_flat_field_test = False
        if exp_type == "NRS_BRIGHTOBJ":
            slit = model
            continue_flat_field_test = True
        else:
            for slit_in_MultiSlitModel in model.slits:
                if slit_in_MultiSlitModel.name == slit_id:
                    slit = slit_in_MultiSlitModel
                    continue_flat_field_test = True
                    break

        if not continue_flat_field_test:
            continue
        else:
            # select the appropriate S- and F-flat fast vector
            sfv_wav, sfv_dat = get_slit_wavdat(sffastvar, slit_id)
            ffv_wav, ffv_dat = get_slit_wavdat(fffastvar, slit_id)

            msg = "\nWorking with slit: " + slit_id
            print(msg)
            log_msgs.append(msg)

            # obtain corresponding occurrence of the SCI extension, i.e. for second occurrence of SCI then ext=2
            ext = 1
            if len(sci_ext_list) > 1:
                ei = sltname_list.index(slit_id)
                ext += ei

            print("exp_type = ", exp_type)
            print("SCI ext = ", ext)
            if isinstance(step_input_filename, str):
                ff = fits.open(step_input_filename)
                print(ff.info())
                ff.close()

            # get the wavelength
            x, y = wcstools.grid_from_bounding_box(slit.meta.wcs.bounding_box, step=(1, 1), center=True)
            ra, dec, wave = slit.meta.wcs(x, y)  # wave is in microns

            # get the subwindow origin
            px0 = slit.xstart - 1 + model.meta.subarray.xstart
            py0 = slit.ystart - 1 + model.meta.subarray.ystart
            msg = " Subwindow origin:   px0=" + repr(px0) + "   py0=" + repr(py0)
            print(msg)
            log_msgs.append(msg)
            n_p = np.shape(wave)
            nw = n_p[0] * n_p[1]
            nw1, nw2 = n_p[1], n_p[0]  # remember that x=nw1 and y=nw2 are reversed  in Python
            if debug:
                print(" nw1, nw2, nw = ", nw1, nw2, nw)

            delf = np.zeros([nw2, nw1]) + 999.0
            flatcor = np.zeros([nw2, nw1]) + 999.0
            flat_err = np.zeros([nw2, nw1]) + 999.0
            delflaterr = np.zeros([nw2, nw1])

            # read the pipeline-calculated flat image, using the corresponding SCI extension number
            pipeflat = flatfile_hdu["SCI", ext].data
            pipeflat_err = flatfile_hdu["ERR", ext].data

            # make sure the two arrays are the same shape
            if np.shape(flatcor) != np.shape(pipeflat):
                msg1 = 'WARNING -> Something went wrong, arrays are not the same shape:'
                msg2 = 'np.shape(flatcor) = ' + repr(np.shape(flatcor)) + '   np.shape(pipeflat) = ' + repr(
                    np.shape(pipeflat))
                msg3 = 'Mathematical operations will fail. Exiting the loop here and setting test as FAILED.'
                print(msg1)
                print(msg2)
                print(msg3)
                log_msgs.append(msg1)
                log_msgs.append(msg2)
                log_msgs.append(msg3)
                msg = " *** Result of the test: " + test_result + "\n"
                print(msg)
                log_msgs.append(msg)
                test_result = "FAILED"
                total_test_result.append(test_result)
                continue

            # loop through the wavelengths
            msg = " Looping through the wavelengths... "
            print(msg)
            log_msgs.append(msg)
            for j in range(nw1):  # in x
                for k in range(nw2):  # in y
                    if np.isfinite(wave[k, j]):  # skip if wavelength is NaN
                        # get thr full-frame pixel indeces for D- and S-flat image components
                        pind = [k + py0 - 1, j + px0 - 1]

                        # get the pixel bandwidth
                        if (j != 0) and (j < nw1 - 1):
                            if np.isfinite(wave[k, j + 1]) and np.isfinite(wave[k, j - 1]):
                                delw = 0.5 * (wave[k, j + 1] - wave[k, j - 1])
                            if np.isfinite(wave[k, j + 1]) and not np.isfinite(wave[k, j - 1]):
                                delw = wave[k, j + 1] - wave[k, j]
                            if not np.isfinite(wave[k, j + 1]) and np.isfinite(wave[k, j - 1]):
                                delw = wave[k, j] - wave[k, j - 1]
                        if j == 0:
                            delw = wave[k, j + 1] - wave[k, j]
                        if j == nw - 1:
                            delw = wave[k, j] - wave[k, j - 1]

                        # integrate over D-flat fast vector
                        dfrqe_wav = dfrqe["wavelength"][0]
                        dfrqe_rqe = dfrqe["data"][0]
                        iw = np.where((dfrqe_wav >= wave[k, j] - delw / 2.) & (dfrqe_wav <= wave[k, j] + delw / 2.))
                        int_tab = auxfunc.idl_tabulate(dfrqe_wav[iw], dfrqe_rqe[iw])
                        first_dfrqe_wav, last_dfrqe_wav = dfrqe_wav[iw[0]][0], dfrqe_wav[iw[0]][-1]
                        dff = int_tab / (last_dfrqe_wav - first_dfrqe_wav)

                        if debug:
                            print("k, py0, j, px0 = ", k, py0, j,px0)
                            print("pind = [k + py0 - 1, j + px0 - 1] = ", pind)
                            print("np.shape(dfrqe_wav) : ", np.shape(dfrqe_wav))
                            print("np.shape(dfrqe_rqe) : ", np.shape(dfrqe_rqe))
                            print("dfimdq[pind[0],[pind[1]] : ", dfimdq[pind[0], pind[1]])
                            print("for dfrqe_wav, np.shape(iw) = ", np.shape(iw))
                            print("np.shape(dfrqe_wav) = ", np.shape(dfrqe_wav[iw]))
                            print("np.shape(dfrqe_rqe) = ", np.shape(dfrqe_rqe[iw]))
                            print("int_tab = ", int_tab)
                            print("np.shape(dfim) = ", np.shape(dfim))
                            print("dff = ", dff)

                        # interpolate over D-flat cube
                        iloc = auxfunc.idl_valuelocate(dfwave, wave[k, j])[0]
                        if dfwave[iloc] > wave[k, j]:
                            iloc -= 1
                        ibr = [iloc]
                        if iloc != len(dfwave) - 1:
                            ibr.append(iloc + 1)
                        # get the values in the z-array at indeces ibr, and x=pind[1] and y=pind[0]
                        zz = dfim[:, pind[0], pind[1]][ibr]
                        # now determine the length of the array with only the finite numbers
                        zzwherenonan = np.where(np.isfinite(zz))
                        kk = np.size(zzwherenonan)
                        dfs = 1.0
                        if (wave[k, j] <= max(dfwave)) and (wave[k, j] >= min(dfwave)) and (kk == 2):
                            dfs = np.interp(wave[k, j], dfwave[ibr], zz[zzwherenonan])

                        # check DQ flags
                        if dfimdq[pind[0]][pind[1]] != 0:
                            dfs = 1.0

                        if debug:
                            print("wave[k, j] = ", wave[k, j])
                            print("iloc = ", iloc)
                            print("ibr = ", ibr)
                            print("np.interp(wave[k, j], dfwave[ibr], zz[zzwherenonan]) = ",
                                  np.interp(wave[k, j], dfwave[ibr], zz[zzwherenonan]))
                            print("dfs = ", dfs)

                        # integrate over S-flat fast vector
                        up_wav = wave[k, j] + delw / 2.0
                        lo_wav = wave[k, j] - delw / 2.0
                        iw = np.where((sfv_wav >= lo_wav) & (sfv_wav <= up_wav))
                        sff = 1.0
                        if np.size(iw) > 2:
                            int_tab = auxfunc.idl_tabulate(sfv_wav[iw], sfv_dat[iw])
                            first_sfv_wav, last_sfv_wav = sfv_wav[iw[0]][0], sfv_wav[iw[0]][-1]
                            sff = int_tab / (last_sfv_wav - first_sfv_wav)
                        else:
                            # find the nearest point in the reference table and use that as center plus the one previous
                            # and following it to create array to interpolate.
                            nearest_wav, nearest_wav_idx = auxfunc.find_nearest(sfv_wav, wave[k, j])
                            nearest_sfv = sfv_dat[nearest_wav_idx]
                            prev_nearest_wav, prev_nearest_sfv = sfv_wav[nearest_wav_idx-1], sfv_dat[nearest_wav_idx-1]
                            foll_nearest_wav, foll_nearest_sfv = sfv_wav[nearest_wav_idx+1], sfv_dat[nearest_wav_idx+1]
                            nearest_wav_arr = np.array([prev_nearest_wav, nearest_wav, foll_nearest_wav])
                            nearest_sfv_arr = np.array([prev_nearest_sfv, nearest_sfv, foll_nearest_sfv])

                            sff = np.interp(wave[k, j], nearest_wav_arr, nearest_sfv_arr)
                            if debug:
                                print("No wavelengths found to integrate over, trying an interpolation to find sff: ")
                                print("  prev_nearest_wav, prev_nearest_sfv: ", prev_nearest_wav, prev_nearest_sfv)
                                print("  nearest_wav, nearest_sfv: ", nearest_wav, nearest_sfv)
                                print("  foll_nearest_wav, foll_nearest_sfv: ", foll_nearest_wav, foll_nearest_sfv)
                        # get s-flat pixel-dependent correction
                        sfs = 1.0
                        if sfimdq[pind[0], pind[1]] == 0:
                            sfs = sfim[pind[0], pind[1]]

                        if debug:
                            print("wave[k, j] - delw / 2.0, wave[k, j] + delw / 2.0 = ", wave[k, j] - delw / 2.0, wave[k, j] + delw / 2.0)
                            print("for sfv_wav, iw = ", iw)
                            print("np.size(iw), np.shape(iw) = ", np.size(iw), np.shape(iw))
                            print("np.shape(sfv_wav) = ", np.shape(sfv_wav))
                            print("np.shape(sfv_dat) = ", np.shape(sfv_dat))
                            print("int_tab = ", int_tab)
                            print("sff = ", sff)
                            print("sfimdq[pind[0], pind[1]] = ", sfimdq[pind[0], pind[1]])
                            print("sfim[pind[0], pind[1]] = ", sfim[pind[0], pind[1]])
                            print("sfs = ", sfs)

                        # integrate over F-flat fast vector
                        # reference file blue cutoff is 1 micron, so need to force solution for shorter wavs
                        fff = 1.0
                        fff_err = 0.0
                        if wave[k, j] - delw / 2.0 >= 1.0:
                            iw = np.where((ffv_wav >= wave[k, j] - delw / 2.0) & (ffv_wav <= wave[k, j] + delw / 2.0))
                            if np.size(iw) > 1:
                                int_tab = auxfunc.idl_tabulate(ffv_wav[iw], ffv_dat[iw])
                                first_ffv_wav, last_ffv_wav = ffv_wav[iw[0]][0], ffv_wav[iw[0]][-1]
                                fff = int_tab / (last_ffv_wav - first_ffv_wav)
                                # f-flat corresponding error estimation by following same logic with average of
                                # value+err and value-err
                                if fferr.size != 0:
                                    int_tab_pluserr = auxfunc.idl_tabulate(ffv_wav[iw], ffv_dat[iw]+fferr[k, j])
                                    fff_pluserr = int_tab_pluserr / (last_ffv_wav - first_ffv_wav)
                                    int_tab_minuserr = auxfunc.idl_tabulate(ffv_wav[iw], ffv_dat[iw]-fferr[k, j])
                                    fff_minuserr = int_tab_minuserr / (last_ffv_wav - first_ffv_wav)
                                    fff_err = (fff_pluserr - fff_minuserr)/2
                            else:
                                # find the nearest point in the reference table and use that as center plus the one previous
                                # and following it to create array to interpolate.
                                nearest_wav, nearest_wav_idx = auxfunc.find_nearest(ffv_wav, wave[k, j])
                                nearest_ffv = ffv_dat[nearest_wav_idx]
                                prev_nearest_wav, prev_nearest_ffv = ffv_wav[nearest_wav_idx-1], ffv_dat[nearest_wav_idx-1]
                                foll_nearest_wav, foll_nearest_ffv = ffv_wav[nearest_wav_idx+1], ffv_dat[nearest_wav_idx+1]
                                nearest_wav_arr = np.array([prev_nearest_wav, nearest_wav, foll_nearest_wav])
                                nearest_ffv_arr = np.array([prev_nearest_ffv, nearest_ffv, foll_nearest_ffv])
                                fff = np.interp(wave[k, j], nearest_wav_arr, nearest_ffv_arr)
                                if debug:
                                    print("No wavelengths found to integrate over, trying an interpolation to find fff: ")
                                    print("  prev_nearest_wav, prev_nearest_ffv: ", prev_nearest_wav, prev_nearest_ffv)
                                    print("  nearest_wav, nearest_ffv: ", nearest_wav, nearest_ffv)
                                    print("  foll_nearest_wav, foll_nearest_ffv: ", foll_nearest_wav, foll_nearest_ffv)
                                # f-flat corresponding error estimation by following same logic with average of
                                # value+err and value-err, where err is value at wavelength[k, j]
                                if fferr.size != 0:
                                    nearest_ffv_arr = np.array([ffv_dat[nearest_wav_idx-1] + fferr[k, j],
                                                                ffv_dat[nearest_wav_idx] + fferr[k, j],
                                                                ffv_dat[nearest_wav_idx+1] + fferr[k, j]])
                                    fff_pluserr = np.interp(wave[k, j], nearest_wav_arr, nearest_ffv_arr)
                                    nearest_ffv_arr = np.array([ffv_dat[nearest_wav_idx-1] - fferr[k, j],
                                                                ffv_dat[nearest_wav_idx] - fferr[k, j],
                                                                ffv_dat[nearest_wav_idx+1] - fferr[k, j]])
                                    fff_minuserr = np.interp(wave[k, j], nearest_wav_arr, nearest_ffv_arr)
                                    fff_err = (fff_pluserr - fff_minuserr)/2

                        # calculate the total value for the flat correction
                        flatcor[k, j] = dff * dfs * sff * sfs * fff

                        # calculate the corresponding error propagation
                        # If the errors were independent from each other
                        # flat_err[k, j] = np.sqrt( dff_err**2/dff**2 + dfs_err**2/dfs**2 + sff_err**2/sff**2 + fff_err**2/fff**2 ) * flatcor[k, j]
                        # However, the errors of the D-flat and S-flat are correlated to the F-flat because it corrects for
                        # any field- and wavelength-dependent effects of the OTE and the FORE optics
                        # (also, the S-flat uncertainy is currently being overestimated as an arbitrary 10%), hence, the
                        # uncertainty from the S- and D-flats will be ignored and we assume that the F-flat uncertainties dominate
                        flat_err[k, j] = np.sqrt( fff_err**2/fff**2 ) * flatcor[k, j]

                        if debug:
                            print("np.shape(iw) =", np.shape(iw))
                            print("np.shape(ffv_wav) = ", np.shape(ffv_wav))
                            print("np.shape(ffv_dat) = ", np.shape(ffv_dat))
                            print("wave[k, j] - delw / 2.0 = ", wave[k, j] - delw / 2.0 )
                            print("fff = ", fff)
                            print("flatcor[k, j] = dff * dfs * sff * sfs * fff = ", flatcor[k, j])
                            print("dff, dfs, sff, sfs, fff:", dff, dfs, sff, sfs, fff)

                        try:
                            # Difference between pipeline and calculated values
                            delf[k, j] = pipeflat[k, j] - flatcor[k, j]
                            # difference between pipeline errors array and calculated values
                            delflaterr[k, j] = pipeflat_err[k, j] - flat_err[k, j]

                            # Remove all pixels with values=1 (outside slit boundaries) for statistics
                            if pipeflat[k, j] == 1:
                                delf[k, j] = 999.0
                                delflaterr[k, j] = 999.0
                            if np.isnan(wave[k, j]):
                                flatcor[k, j] = 1.0  # no correction if no wavelength
                                flat_err[k, j] = 999.0

                            if debug:
                                print("Difference between pipeline and calculated values:  delf[k, j] = ", delf[k, j])
                                print("k, j = ", k, j)
                                print()
                        except IndexError:
                            print('* ERROR: Index error for  [j, k] ', j, k)

            if debug:
                no_999 = delf[np.where(delf != 999.0)]
                print("np.shape(no_999) = ", np.shape(no_999))
                alldelf = delf.flatten()
                print("median of the whole array: ", np.median(alldelf))
                print("median, stdev in delf: ", np.median(no_999), np.std(no_999))
                neg_vals = no_999[np.where(no_999 < 0.0)]
                print("neg_vals = ", np.shape(neg_vals))
                print("np.shape(delf) = ", np.shape(delf))
                print()

            nanind = np.isnan(delf)  # get all the nan indexes
            notnan = ~nanind  # get all the not-nan indexes
            delf = delf[notnan]  # get rid of NaNs
            delflaterr = delflaterr[notnan]
            if debug:
                print('nanind = ', np.shape(nanind), nanind)
                print('notnan = ', np.shape(notnan), notnan)
                print('delf = ', np.shape(delf),  delf)

            if delf.size == 0:
                msg1 = " * Unable to calculate statistics because difference array has all values as NaN. Test will be set to FAILED."
                print(msg1)
                log_msgs.append(msg1)
                test_result = "FAILED"
                delfg_mean, delfg_median, delfg_std = np.nan, np.nan, np.nan
                stats = [delfg_mean, delfg_median, delfg_std]
            else:
                msg = "Calculating statistics... "
                print(msg)
                log_msgs.append(msg)
                # find out which outliers work better
                outlier_threshold = 0.1
                if len(delf[np.where((delf != 999.0) & (delf < outlier_threshold) & (delf > -1*outlier_threshold))]) <= 1:
                    # the differences are much larger than 0.1, need a larger outlier threshold
                    # sort the differences and return the 5th largest value
                    outlier_threshold = sorted(list(np.abs(delf[np.where(delf != 999.0)])), reverse=True)[5]
                # ignore outliers
                delfg = delf[np.where((delf != 999.0) & (delf < outlier_threshold) & (delf > -1*outlier_threshold))]
                delflaterr = delflaterr[np.where((delf != 999.0) & (delf < outlier_threshold) & (delf > -1*outlier_threshold))]
                if debug:
                    print('outlier_threshold = ', outlier_threshold)
                    print('delfg = ', np.shape(delfg), delfg)
                if delfg.size == 0:
                    msg1 = " * Unable to calculate statistics because difference array has all outlier values. Test will be set to FAILED."
                    print(msg1)
                    log_msgs.append(msg1)
                    test_result = "FAILED"
                    delfg_mean, delfg_median, delfg_std = np.nan, np.nan, np.nan
                    stats = [delfg_mean, delfg_median, delfg_std]
                else:
                    stats_and_strings = auxfunc.print_stats(delfg, "Flat Difference", float(threshold_diff), absolute=True)
                    stats, stats_print_strings = stats_and_strings
                    delfg_mean, delfg_median, delfg_std = stats
                    for msg in stats_print_strings:
                        log_msgs.append(msg)
                    _ = auxfunc.print_stats(delflaterr, "Flat Error Difference", float(threshold_diff), absolute=True)

                    # This is the key argument for the assert pytest function
                    median_diff = False
                    if abs(delfg_median) <= float(threshold_diff):
                        median_diff = True
                    if median_diff:
                        test_result = "PASSED"
                    else:
                        test_result = "FAILED"

            msg = " *** Result of the test: " + test_result + "\n"
            print(msg)
            log_msgs.append(msg)
            total_test_result.append(test_result)

            # make histogram
            flatcor_copy = deepcopy(flatcor)
            flatcor_copy[np.where(flatcor_copy == 999.0)] = np.nan
            flat_err_copy = deepcopy(flat_err)
            flat_err_copy[np.where(flat_err_copy == 999.0)] = np.nan
            if show_figs or save_figs:

                # set plot variables
                main_title = filt + "   " + grat + "   SLIT=" + slit_id + "\n"
                bins = None  # binning for the histograms, if None the function will select them automatically
                #             lolim_x, uplim_x, lolim_y, uplim_y
                plt_origin = None

                # Residuals img and histogram
                title = main_title + "Residuals"
                info_img = [title, "x (pixels)", "y (pixels)"]
                xlabel, ylabel = "flat$_{pipe}$ - flat$_{calc}$", "N"
                info_hist = [xlabel, ylabel, bins, stats]
                if delfg.size != 0 and delfg[1] is np.nan:
                    msg = "Unable to create plot of pipeline - calulated flat values."
                    print(msg)
                    log_msgs.append(msg)
                else:
                    if output_directory is not None:
                        t = (file_basename, "FS_flattest_" + slit_id + "_histogram.png")
                        plt_name = os.path.join(file_path, "_".join(t))
                    else:
                        plt_name = os.path.join(os.getcwd(), "FS_flattest_" + det + "_" + slit_id + "_histogram.png")
                        print("No output_directory was provided. Figures will be saved in current working directory:")
                        print(plt_name + "\n")
                    difference_img = (pipeflat - flatcor_copy)  # /flatcor
                    if debug:
                        print('np.shape(pipeflat), pipeflat: ', np.shape(pipeflat), pipeflat[np.where(pipeflat != 1.0)])
                        print('np.shape(flatcor), flatcor: ', np.shape(flatcor), flatcor[np.where(flatcor != 999.0)])
                        print('np.shape(difference_img), difference_img: ', np.shape(difference_img), ~np.isnan(difference_img))

                    # set the range of values to be shown in the image, will affect color scale
                    if delfg_std <= 0.0:
                        vminmax = [-5 * delfg_std, 5 * delfg_std]
                    else:
                        vminmax = None
                    auxfunc.plt_two_2Dimgandhist(difference_img, delfg, info_img, info_hist, plt_name=plt_name,
                                                 vminmax=vminmax,
                                                 plt_origin=plt_origin, show_figs=show_figs, save_figs=save_figs)

                # now make the plot for the errors comparison
                main_title = "ERRORS Comparison for " + filt + "   " + grat + "   SLIT=" + slit_id + "\n"
                title = main_title + "Residuals"
                info_img = [title, "x (pixels)", "y (pixels)"]
                xlabel, ylabel = "flat_err$_{pipe}$ - flat_err$_{calc}$", "N"
                info_hist = [xlabel, ylabel, bins, stats]
                if delflaterr.size != 0 and delflaterr[1] is np.nan:
                    msg = "Unable to create plot of pipeline - calulated flat error values."
                    print(msg)
                    log_msgs.append(msg)
                else:
                    if output_directory is not None:
                        t = (file_basename, "FS_flattest_error_" + slit_id + "_histogram.png")
                        plt_name = os.path.join(file_path, "_".join(t))
                    else:
                        plt_name = os.path.join(os.getcwd(), "FS_flattest_error_" + det + "_" + slit_id + "_histogram.png")
                        print("No output_directory was provided. Figures will be saved in current working directory:")
                        print(plt_name + "\n")
                    difference_img = pipeflat_err - flat_err_copy
                    auxfunc.plt_two_2Dimgandhist(difference_img, delflaterr, info_img, info_hist, plt_name=plt_name,
                                                 vminmax=None, plt_origin=plt_origin,
                                                 show_figs=show_figs, save_figs=save_figs)

            elif not save_figs and not show_figs:
                msg = "Not making plots because both show_figs and save_figs were set to False."
                print(msg)
                log_msgs.append(msg)
            elif not save_figs:
                msg = "Not saving plots because save_figs was set to False."
                print(msg)
                log_msgs.append(msg)

            # create fits file to hold the calculated flat for each slit
            if writefile:
                msg = "Saving the fits files with the calculated flat for each slit..."
                print(msg)
                log_msgs.append(msg)

                # this is the file to hold the image of calculated flat correctoion values
                #flatcor_copy[np.where(flatcor == 999.0)] = 1.0
                flatcor_copy[np.isnan(flatcor_copy)] = 1.0
                outfile_ext = fits.ImageHDU(flatcor_copy, name=slit_id)
                outfile.append(outfile_ext)
                # this is the file to hold the corresponding error extensions
                outfile_ext = fits.ImageHDU(flat_err_copy, name=slit_id+'_ERR')
                outfile.append(outfile_ext)

                # this is the file to hold the image of pipeline-calculated difference values
                complfile_ext = fits.ImageHDU(delf, name=slit_id)
                complfile.append(complfile_ext)

                # the file is not yet written, indicate that this slit was appended to list to be written
                msg = "Extension " + repr(
                    i) + " appended to list to be written into calculated and comparison fits files."
                print(msg)
                log_msgs.append(msg)

    if writefile:
        outfile_name = flatfile.replace("interpolatedflat.fits", "flat_calc.fits")
        complfile_name = flatfile.replace("interpolatedflat.fits", "flat_comp.fits")

        # create the fits list to hold the calculated flat values for each slit
        outfile.writeto(outfile_name, overwrite=True)

        # this is the file to hold the image of pipeline-calculated difference values
        complfile.writeto(complfile_name, overwrite=True)

        msg = "\nFits file with calculated flat values of each slit saved as: "
        print(msg)
        log_msgs.append(msg)
        print(outfile_name)
        log_msgs.append(outfile_name)

        msg = "Fits file with comparison (pipeline flat - calculated flat) saved as: "
        print(msg)
        log_msgs.append(msg)
        print(complfile_name)
        log_msgs.append(complfile_name)

    # If all tests passed then pytest will be marked as PASSED, else it will be FAILED
    FINAL_TEST_RESULT = False
    for t in total_test_result:
        if t == "FAILED":
            FINAL_TEST_RESULT = False
            break
        else:
            FINAL_TEST_RESULT = True

    if FINAL_TEST_RESULT:
        msg = "\n *** Final result for flat_field test will be reported as PASSED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "All slits PASSED flat_field test."
    else:
        msg = "\n *** Final result for flat_field test will be reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "One or more slits FAILED flat_field test."

    # end the timer
    flattest_end_time = time.time() - flattest_start_time
    if flattest_end_time > 60.0:
        flattest_end_time = flattest_end_time / 60.0  # in minutes
        flattest_tot_time = "* Script flattest_fs.py took ", repr(flattest_end_time) + " minutes to finish."
        if flattest_end_time > 60.0:
            flattest_end_time = flattest_end_time / 60.  # in hours
            flattest_tot_time = "* Script flattest_fs.py took ", repr(flattest_end_time) + " hours to finish."
    else:
        flattest_tot_time = "* Script flattest_fs.py took ", repr(flattest_end_time) + " seconds to finish."
    print(flattest_tot_time)
    log_msgs.append(flattest_tot_time)

    return FINAL_TEST_RESULT, result_msg, log_msgs


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("step_input_filename",
                        action='store',
                        default=None,
                        help='Name of input fits file prior to assign_wcs step, i.e. blah_rate.fits')
    parser.add_argument("dflat_path",
                        action='store',
                        default=None,
                        help='Path and name of D-flat file.')
    parser.add_argument("sflat_path",
                        action='store',
                        default=None,
                        help='Path and name of S-flat file.')
    parser.add_argument("fflat_path",
                        action='store',
                        default=None,
                        help='Path and name of F-flat file.')
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
                        default=9.999e-05,
                        type=float,
                        help='Use flag -t to change the default threshold (currently set to 9.999e-05).')
    parser.add_argument("-o",
                        dest="output_directory",
                        action='store',
                        default=None,
                        help='Use flag -o to provide the output_directory to save the figures and files.')
    parser.add_argument("-d",
                        dest="debug",
                        action='store_true',
                        default=False,
                        help='Use flag -d to turn on debug mode.')
    args = parser.parse_args()

    # Set variables
    step_input_filename = args.step_input_filename
    dflat_path = args.dflat_path
    sflat_path = args.sflat_path
    fflat_path = args.fflat_path
    writefile = args.writefile
    save_figs = args.save_figs
    show_figs = args.show_figs
    threshold_diff = args.threshold_diff
    output_directory = args.output_directory
    debug = args.debug

    # print pipeline version
    import jwst

    print("\n  ** using pipeline version: ", jwst.__version__, "** \n")

    # Run the principal function of the script
    flattest(step_input_filename, dflat_path=dflat_path, sflat_path=sflat_path,
             fflat_path=fflat_path, writefile=writefile, show_figs=show_figs, save_figs=save_figs,
             threshold_diff=threshold_diff, output_directory=output_directory, debug=debug)


if __name__ == '__main__':
    sys.exit(main())

