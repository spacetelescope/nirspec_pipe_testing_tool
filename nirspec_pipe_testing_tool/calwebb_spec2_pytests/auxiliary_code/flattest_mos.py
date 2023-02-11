import time
import os
import argparse
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from copy import deepcopy

from jwst import datamodels

from . import auxiliary_functions as auxfunc


"""
This script tests the pipeline flat field step output for MOS data. It is the python version of the IDL script
(with the same name) written by James Muzerolle.
"""


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "3.9"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Completely changed script to use the datamodel instead of the compute_world_coordinates
#                         script, and added new routines for statistics calculations.
# Jun 2018 - Version 2.1: Changed extension numbers for the name of the extension in the D-, F-, and S-flats.
# Jun 2018 - Version 3.0: Change the loop over the pixels to go over both indices instead of flattening arrays.
# Jun 2018 - Version 3.1: Removed function reverse_cols because it was not behaving as expected.
# Aug 2018 - Version 3.2: Fixed bugs per Phil Hodge recommendations.
# Apr 2019 - Version 3.3: Implemented capability to return logging messages.
# May 2019 - Version 3.4: Implemented images of the residuals.
# Jun 2019 - Version 3.5: Updated name of interpolated flat to be the default pipeline name for this file.
# Sep 2020 - Version 3.6: Fixed code to match latest pipeline changes to fix MOS flat field: set flat-d, -s, -f to 1
#                         if the data is Nan or 0, then DQ flag is set to DO_NOT_USE + NO_FLAT_FIELD. The code now
#                         only looks at the flat DQ flag set to DO_NOT_USE in determining if the data should be used.
# Jan 2021 - Version 3.7: Implemented option to run with object instead of input fits file.
# Sep 2021 - Version 3.8: Changing wavelength array to be read from model.slit instead of wcs
#                         (as recommended in Jira issue https://jira.stsci.edu/browse/JP-2225)
# Oct 2022 - Version 3.9: Fixed code to read new post-commissioning reference files in CRDS format.


def flattest(step_input_filename, dflat_path, sflat_path, fflat_path, msa_shutter_conf,
             writefile=False, show_figs=True, save_figs=False, interpolated_flat=None,
             threshold_diff=1.0e-14, debug=False):
    """
    This function does the WCS comparison from the world coordinates calculated using the
    compute_world_coordinates.py script with the ESA files. The function calls that script.

    Args:
        step_input_filename: str, name of the output fits file from the flat field step (with full path)
        dflatref_path: str, path of where the D-flat reference fits files
        sflat_path: str, path of where the S-flat reference fits files
        fflat_path: str, path of where the F-flat reference fits files
        msa_shutter_conf: str, full path and name of the MSA configuration fits file
        writefile: boolean, if True writes the fits files of the calculated flat and difference images
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots (the 3 plots can be saved or not independently with the function call)
        interpolated_flat: string, name of the on-the-fly interpolated pipeline flat
        threshold_diff: float, threshold difference between pipeline output and ESA file
        debug: boolean, if true a series of print statements will show on-screen

    Returns:
        - 1 plot, if told to save and/or show.
        - median_diff: Boolean, True if smaller or equal to 1e-14
        - log_msgs: list, all print statements are captured in this variable

    """

    log_msgs = []

    # start the timer
    flattest_start_time = time.time()

    # get info from the rate file header
    if isinstance(step_input_filename, str):
        msg = 'step_input_filename=' + step_input_filename
        print(msg)
        log_msgs.append(msg)
        if "wavecorr" not in step_input_filename:
            wcs_file = step_input_filename.replace("_flat_field", "_wavecorr")
        else:
            wcs_file = step_input_filename
        model = datamodels.MultiSlitModel(wcs_file)

    else:
        model = step_input_filename
    
    # read in the on-the-fly flat image
    if interpolated_flat is None:
        flatfile = step_input_filename.replace("_flat_field", "_interpolatedflat")
    else:
        flatfile = interpolated_flat
    # get all the science extensions in the flatfile
    sci_ext_name_list, hdu_idx_list = auxfunc.get_sci_extensions(flatfile, lists=True)

    # get basic info from model
    if isinstance(model, list):  # this was added for the validation notebooks to work
        model = model[0]

    # get basic info from model
    det = model.meta.instrument.detector
    grat = model.meta.instrument.grating
    filt = model.meta.instrument.filter
    lamp = model.meta.instrument.lamp_state
    exptype = model.meta.exposure.type.upper()

    msg = "flat_field_file  -->     Grating:" + grat + "   Filter:" + filt + "   LAMP:" + lamp
    print(msg)
    log_msgs.append(msg)

    if isinstance(step_input_filename, str):
        file_path = step_input_filename.replace(os.path.basename(step_input_filename), "")
        file_basename = os.path.basename(step_input_filename.replace(".fits", ""))
    else:
        file_path = os.path.dirname(os.path.realpath(__file__))
        file_basename = grat + "_" + filt + "_" + det

    # define the mode
    if "msa" in exptype.lower():
        mode = "MOS"
    else:
        mode = "not MOS data"

    if isinstance(step_input_filename, str):
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
    if not os.path.isfile(dflat_path):
        result_msg = "Test skiped because the D-flat provided does not exist: {}".format(dflat_path)
        print(msg)
        median_diff = "skip"
        return median_diff, result_msg, log_msgs
    dfile = dflat_path
    msg0 = " * This flat test is using the following reference files "
    msg1 = "    D-flat: " + dfile
    print(msg0)
    print(msg1)
    log_msgs.append(msg0)
    log_msgs.append(msg1)
    with fits.open(dfile) as dfhdu:
        dfhdr_sci = dfhdu["SCI"].header
        dfim = dfhdu["SCI"].data
        dfimdq = dfhdu["DQ"].data
        dfimerr = dfhdu["ERR"].data
        dfrqe = dfhdu["FAST_VARIATION"].data
    ns = np.shape(dfim)
    naxis3 = dfhdr_sci["NAXIS3"]
    # get the wavelength values
    dfwave = np.array([])
    for i in range(naxis3):
        keyword = "_".join(("PFLAT", str(i+1)))
        dfwave = np.append(dfwave, dfhdr_sci[keyword])

    # S-flat
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
        # This is the key argument for the assert pytest function
        result_msg = "Test skiped because there is no flat correspondence for the filter in the data: {}".format(filt)
        median_diff = "skip"
        return median_diff, result_msg, log_msgs

    if not os.path.isfile(sflat_path):
        result_msg = "Test skiped because the S-flat provided does not exist: {}".format(sflat_path)
        print(msg)
        median_diff = "skip"
        return median_diff, result_msg, log_msgs
    sfile = sflat_path
    msg = "    S-flat: " + sfile
    print(msg)
    log_msgs.append(msg)
    with fits.open(sfile) as sfhdu:
        sfim = sfhdu["SCI"].data
        sfimdq = sfhdu["DQ"].data
        sfimerr = sfhdu["ERR"].data
        sffastvar = sfhdu["FAST_VARIATION"].data
        sfhdu_sci = sfhdu["SCI"].header
    sfv_wav, sfv_dat = auxfunc.get_slit_wavdat(sffastvar, 'ANY')

    # get the wavelength values for sflat cube
    sfimwave = np.array([])
    naxis3 = sfhdu_sci["NAXIS3"]
    for i in range(0, naxis3):
        if i+1 < 10:
            keyword = "".join(("FLAT_0", str(i+1)))
        else:
            keyword = "".join(("FLAT_", str(i+1)))
        if debug:
            print("S-flat -> using ", keyword)
        try:
            sfimwave = np.append(sfimwave, sfhdu_sci[keyword])
        except KeyError:
            continue

    # F-Flat
    if not os.path.isfile(fflat_path):
        result_msg = "Test skiped because the F-flat provided does not exist: {}".format(fflat_path)
        print(msg)
        median_diff = "skip"
        return median_diff, result_msg, log_msgs
    ffile = fflat_path
    msg = "    F-flat: " + ffile
    print(msg)
    log_msgs.append(msg)

    # assign the correct quadrant variables accordingly
    # Quadrant 1
    ffhdu = fits.open(ffile)
    ffsq1 = ffhdu['SCI', 1].data   # first occurrence of the SCI extension
    ffhdr_sci = ffhdu['SCI', 1].header
    naxis3 = ffhdr_sci["NAXIS3"]
    ffsallwave = np.array([])
    for i in range(0, naxis3):
        if i <= 9:
            suff = "".join(("0", str(i)))
        else:
            suff = str(i)
        t = ("FLAT", suff)
        keyword = "_".join(t)
        if debug:
            print(" F-flat -> ", keyword)
        ffsallwave = np.append(ffsallwave, ffhdr_sci[keyword])
        # this array is independant from the quadrant
    ffserrq1 = ffhdu["ERR", 1].data
    ffsdqq1 = ffhdu["DQ", 1].data
    ffvq1 = ffhdu["FAST_VARIATION", 1].data
    # Quadrant 2
    ffsq2 = ffhdu["SCI", 2].data
    ffserrq2 = ffhdu["ERR", 2].data
    ffsdqq2 = ffhdu["DQ", 2].data
    ffvq2 = ffhdu["FAST_VARIATION", 2].data
    # Quadrant 3
    ffsq3 = ffhdu["SCI", 3].data
    ffserrq3 = ffhdu["ERR", 3].data
    ffsdqq3 = ffhdu["DQ", 3].data
    ffvq3 = ffhdu["FAST_VARIATION", 3].data
    # Quadrant 4
    ffsq4 = ffhdu["SCI", 4].data
    ffserrq4 = ffhdu["ERR", 4].data
    ffsdqq4 = ffhdu["DQ", 4].data
    ffvq4 = ffhdu["FAST_VARIATION", 4].data
    ffhdu.close()

    # now prepare the output files and structures to go through each pixel in the test data

    if writefile:
        # create the fits list to hold the image of the correction values
        hdu0 = fits.PrimaryHDU()
        outfile = fits.HDUList()
        outfile.append(hdu0)

        # create the fits list to hold the image of the comparison values
        hdu0 = fits.PrimaryHDU()
        complfile = fits.HDUList()
        complfile.append(hdu0)

    # list to determine if pytest is passed or not
    total_test_result = []

    # get the slitlet info, needed for the F-Flat
    ext_shutter_info = "SHUTTER_INFO"  # this is extension 2 of the msa file, that has the shutter info
    slitlet_info = fits.getdata(msa_shutter_conf, ext_shutter_info)
    sltid = slitlet_info.field("SLITLET_ID")

    # open the flatfile
    flatfile_hdu = fits.open(flatfile)

    #debug = True
    # loop over the 2D subwindows and read in the WCS values
    total_slits = len(model.slits)
    for si, slit in enumerate(model.slits):
        slit_id = slit.name
        msg = "\nWorking with slit ID: "+slit_id+" - which is "+repr(si+1)+" out of "+repr(total_slits)
        print(msg)
        log_msgs.append(msg)
        # get the science extension occurrence in the pipeline calculated flat
        ext = sci_ext_name_list.index(slit_id)+1

        # get the wavelength
        # y, x = np.mgrid[:slit.data.shape[0], :slit.data.shape[1]]
        # ra, dec, wave = slit.meta.wcs(x, y)   # wave is in microns, commented because it's not wavecor product
        wave = slit.wavelength  # uses the object from step wavecor

        # get the subwindow origin
        px0 = slit.xstart - 1 + model.meta.subarray.xstart
        py0 = slit.ystart - 1 + model.meta.subarray.ystart
        msg = " Subwindow origin:   px0="+repr(px0)+"   py0="+repr(py0)
        print(msg)
        log_msgs.append(msg)
        n_p = np.shape(wave)
        nw = n_p[0]*n_p[1]
        nw1, nw2 = n_p[1], n_p[0]   # remember that x=nw1 and y=nw2 are reversed in Python
        if debug:
            print("nw = ", nw)

        # Define the arrays to hold calculations, cor=correction calculation,
        # err=uncdertanties, del=pipeline-calculation  -> 999.0 to know what values to ignore
        delf = np.zeros([nw2, nw1]) + 999.0
        flatcor = np.zeros([nw2, nw1]) + 999.0
        flat_err = np.zeros([nw2, nw1])
        delflaterr = np.zeros([nw2, nw1])

        for j, s in enumerate(sltid):
            if s == int(slit_id):
                im = j
                # get the shutter with the source in it
                if slitlet_info.field("BACKGROUND")[im] == "N":
                    isrc = j
        # changes suggested by Phil Hodge
        quad = slit.quadrant  # slitlet_info.field("SHUTTER_QUADRANT")[isrc]
        row = slit.xcen  # slitlet_info.field("SHUTTER_ROW")[isrc]
        col = slit.ycen  # slitlet_info.field("SHUTTER_COLUMN")[isrc]
        slitlet_id = repr(row)+"_"+repr(col)
        msg = 'silt_id='+repr(slit_id)+"   quad="+repr(quad)+"   row="+repr(row)+"   col="+repr(col)+\
              "   slitlet_id="+repr(slitlet_id)
        print(msg)
        log_msgs.append(msg)

        # get the relevant F-flat reference data
        if quad == 1:
            ffsall = ffsq1
            ffsalldq = ffsdqq1
            ffv = ffvq1
            fferr = ffserrq1
        if quad == 2:
            ffsall = ffsq2
            ffsalldq = ffsdqq2
            ffv = ffvq2
            fferr = ffserrq2
        if quad == 3:
            ffsall = ffsq3
            ffsalldq = ffsdqq3
            ffv = ffvq3
            fferr = ffserrq3
        if quad == 4:
            ffsall = ffsq4
            ffsalldq = ffsdqq4
            ffv = ffvq4
            fferr = ffserrq4
        # Extract the data for the appropriate slit from the fits table
        ffv_wav, ffv_dat = auxfunc.get_slit_wavdat(ffv, 'ANY')

        # loop through the pixels
        msg = "Now looping through the pixels, this will take a while ... "
        print(msg)
        log_msgs.append(msg)
        wave_shape = np.shape(wave)
        for j in range(nw1):   # in x
            for k in range(nw2):   # in y
                if np.isfinite(wave[k, j]):   # skip if wavelength is NaN
                    # get the pixel indexes
                    jwav = wave[k, j]
                    pind = [k+py0-1, j+px0-1]
                    #if debug:
                    #    print('j, k, jwav, px0, py0 : ', j, k, jwav, px0, py0)
                    #    print('pind = ', pind)

                    # get the pixel bandwidth
                    if (j != 0) and (j < nw1-1):
                        if np.isfinite(wave[k, j+1]) and np.isfinite(wave[k, j-1]):
                            delw = 0.5 * (wave[k, j+1] - wave[k, j-1])
                        if np.isfinite(wave[k, j+1]) and not np.isfinite(wave[k, j-1]):
                            delw = wave[k, j+1] - wave[k, j]
                        if not np.isfinite(wave[k, j+1]) and np.isfinite(wave[k, j-1]):
                            delw = wave[k, j] - wave[k, j-1]
                    if j == 0:
                        delw = wave[k, j+1] - wave[k, j]
                    if j == nw-1:
                        delw = wave[k, j] - wave[k, j-1]

                    #if debug:
                    #    print("wave[k, j+1], wave[k, j-1] : ", np.isfinite(wave[k, j+1]), wave[k, j+1], wave[k, j-1])
                    #    print("delw = ", delw)

                    # read the pipeline-calculated flat image
                    pipeflat = flatfile_hdu['SCI', ext].data
                    pipeflat_err = flatfile_hdu['ERR', ext].data

                    # define temporary variables and set wavelength range to operate in
                    dff, dfs, sff, sfs, fff, ffs = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
                    dff_err, dfs_err, sff_err, sfs_err, fff_err, ffs_err = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                    sflat_dqflags_ok, dflat_dqflags_ok, fflat_dqflags_ok = True, True, True
                    if (wave[k, j] >= 0.6) and (wave[k, j] <= 5.3):
                        # integrate over dflat fast vector
                        dfrqe_wav = dfrqe["wavelength"][0]
                        dfrqe_rqe = dfrqe["data"][0]
                        # find the nearest point in the reference table and use that as center plus the one previous
                        # and following it to create array to interpolate.
                        dff = auxfunc.interp_close_pts(wave[k, j], dfrqe_wav, dfrqe_rqe, debug)
                        # the corresponding error is 0.0 because we currently have no information on this

                        # interpolate over dflat cube
                        iloc = auxfunc.idl_valuelocate(dfwave, wave[k, j])[0]
                        if dfwave[iloc] > wave[k, j]:
                            iloc -= 1
                        ibr = [iloc]
                        if iloc != len(dfwave)-1:
                            ibr.append(iloc+1)
                        # get the values in the z-array at indices ibr, and x=pind[1] and y=pind[0]
                        zz = dfim[:, pind[0], pind[1]][ibr]
                        # now determine the length of the array with only the finite numbers
                        zzwherenonan = np.where(np.isfinite(zz))
                        kk = np.size(zzwherenonan)
                        if (wave[k, j] <= max(dfwave)) and (wave[k, j] >= min(dfwave)) and (kk == 2):
                            dfs = np.interp(wave[k, j], dfwave[ibr], zz[zzwherenonan])
                            # calculate corresponding error
                            try:
                                dfs_err = np.interp(wave[k, j], dfwave[ibr], dfimerr[:, pind[0], pind[1]][ibr])
                            except IndexError:   # meaning that the ERR extension does not have the same size array
                                dfs_err = dfimerr[pind[0], pind[1]]

                        # check DQ flags for d-flat
                        if dfimdq[pind[0], pind[1]] != 0:
                            # print("d-flat: DQ flag 1 = DO NOT USE forcing -> dfs=1.0   or  ",
                            # "DQ flag 4 = NO_FLAT_FIELD, also forcing -> dfs=1.0")
                            dfs, dfs_err = 1.0, 0.0
                            dflat_dqflags_ok = False

                        # integrate over s-flat fast vector
                        # find the nearest point in the reference table and use that as center plus the one previous
                        # and following it to create array to interpolate.
                        sff = auxfunc.interp_close_pts(wave[k, j], sfv_wav, sfv_dat, debug)
                        # the corresponding error is 0.0 because we currently have no information on this

                        # interpolate s-flat cube
                        # first, find the index and value of the closest element to the wav on interest
                        iloc = auxfunc.idl_valuelocate(sfimwave, wave[k, j])[0]
                        if sfimwave[iloc] > wave[k, j]:
                            iloc -= 1
                        ibr = [iloc]
                        if iloc != len(sfimwave)-1:
                            ibr.append(iloc+1)
                        # get the values in the z-array at indices ibr, and x=pind[1] and y=pind[0]
                        zz = sfim[:, pind[0], pind[1]][ibr]
                        # now determine the length of the array with only the finite numbers
                        zzwherenonan = np.where(np.isfinite(zz))
                        kk = np.size(zzwherenonan)
                        if (wave[k, j] <= max(sfimwave)) and (wave[k, j] >= min(sfimwave)) and (kk == 2):
                            sfs = np.interp(wave[k, j], sfimwave[ibr], zz[zzwherenonan])
                            # calculate corresponding error
                            try:
                                sfs_err = np.interp(wave[k, j], dfwave[ibr], sfimerr[:, pind[0], pind[1]][ibr])
                            except IndexError:   # meaning that the ERR extension does not have the same size array
                                sfs_err = sfimerr[pind[0], pind[1]]

                        # check DQ flags for s-flat
                        kk = np.where(sfimdq[:, col-1, row-1] != 0)
                        if np.size(kk) >= 1:
                            # print("s-flat: DQ flag DO NOT USE  or  NO_FLAT_FIELD force sfs=1.0")
                            sfs, sfs_err = 1.0, 0.0
                            sflat_dqflags_ok = False

                        # integrate over f-flat fast vector
                        # find the nearest point in the reference table and use that as center plus the one previous
                        # and following it to create array to interpolate.
                        fff = auxfunc.interp_close_pts(wave[k, j], ffv_wav, ffv_dat, debug)
                        # the corresponding error is 0.0 because we currently have no information on this

                        # interpolate over f-flat cube
                        ffs = np.interp(wave[k, j], ffsallwave, ffsall[:, col-1, row-1])
                        # get corresponding error estimation
                        try:
                            ffs_err = fferr[pind[0], pind[1]]
                        except IndexError:   # this means that the ERR extension is not an array
                            ffs_err = 0.0

                        # check DQ flags for f-flat
                        kk = np.where(ffsalldq[:, col-1, row-1] != 0)
                        if np.size(kk) >= 1:
                            ffs = 1.0
                            ffs_err = 0.0
                            fflat_dqflags_ok = False

                        # add correction
                        flatcor[k, j] = 1.0
                        if sflat_dqflags_ok and dflat_dqflags_ok and fflat_dqflags_ok:
                            flatcor[k, j] = dff * dfs * sff * sfs * fff * ffs
                            if np.isnan(flatcor[k, j]):
                                flatcor[k, j] = 1.0
                        if flatcor[k, j] <= 0.0 or pipeflat[k, j] == 1:  # follow the pipeline:
                            flatcor[k, j] = 1.0

                        # calculate the corresponding error propagation
                        try:
                            dff_err2 = dff_err**2/dff**2
                        except ZeroDivisionError:
                            dff_err2 = 0.0
                        try:
                            dfs_err2 = dfs_err**2/dfs**2
                        except ZeroDivisionError:
                            dfs_err2 = 0.0
                        try:
                            sff_err2 = sff_err**2/sff**2
                        except ZeroDivisionError:
                            sff_err2 = 0.0
                        try:
                            sfs_err2 = sfs_err**2/sfs**2
                        except ZeroDivisionError:
                            sfs_err2 = 0.0
                        try:
                            fff_err2 = fff_err**2/fff**2
                        except ZeroDivisionError:
                            fff_err2 = 0.0
                        try:
                            ffs_err2 = ffs_err**2/ffs**2
                        except ZeroDivisionError:
                            ffs_err2 = 0.0
                        if np.isnan(dff_err2):
                            dff_err2 = 0.0
                        if np.isnan(dfs_err2):
                            dfs_err2 = 0.0
                        if np.isnan(sff_err2):
                            sff_err2 = 0.0
                        if np.isnan(sfs_err2):
                            sfs_err2 = 0.0
                        if np.isnan(fff_err2):
                            fff_err2 = 0.0
                        if np.isnan(ffs_err2):
                            ffs_err2 = 0.0
                        error_sq_sum = dff_err2 + dfs_err2 + sff_err2 + sfs_err2 + fff_err2 + ffs_err2
                        flat_err[k, j] = np.sqrt( error_sq_sum ) * flatcor[k, j]

                        if (pind[1]-px0+1 == 9999) and (pind[0]-py0+1 == 9999):
                            if debug:
                                print("pind = ", pind)
                                print("wave[k, j] = ", wave[k, j])
                                print("dfs, dff = ", dfs, dff)
                                print("sfs, sff = ", sfs, sff)
                                print("fff, ffs = ", fff, ffs)
                                print("flatcor[k, j] = ", flatcor[k, j])

                            msg = "Making the plot fot this slitlet..."
                            print(msg)
                            log_msgs.append(msg)
                            # make plot
                            font = {'weight': 'normal',
                                    'size': 16}
                            matplotlib.rc('font', **font)
                            fig = plt.figure(1, figsize=(12, 10))
                            plt.subplots_adjust(hspace=.4)
                            ax = plt.subplot(111)
                            xmin = wave[k, j]-0.01
                            xmax = wave[k, j]+0.01
                            plt.xlim(xmin, xmax)
                            plt.plot(dfwave, dfim[:, pind[0], pind[1]], linewidth=7, marker='D', color='k',
                                     label="dflat_im")
                            plt.plot(wave[k, j], dfs, linewidth=7, marker='D', color='r')
                            plt.plot(dfrqe_wav, dfrqe_rqe, linewidth=7, marker='D', c='k', label="dflat_vec")
                            plt.plot(wave[k, j], dff, linewidth=7, marker='D', color='r')
                            plt.plot(sfimwave, sfim[:, pind[0], pind[1]], linewidth=7, marker='D', color='k',
                                     label="sflat_im")
                            plt.plot(wave[k, j], sfs, linewidth=7, marker='D', color='r')
                            plt.plot(sfv_wav, sfv_dat, linewidth=7, marker='D', color='k', label="sflat_vec")
                            plt.plot(wave[k, j], sff, linewidth=7, marker='D', color='r')
                            # add legend
                            box = ax.get_position()
                            ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])
                            ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
                            plt.minorticks_on()
                            plt.tick_params(axis='both', which='both', bottom=True, top=True, right=True, direction='in',
                                            labelbottom=True)
                            plt.show()
                            msg = "Exiting the program. Unable to calculate statistics. Test set to be SKIPPED."
                            print(msg)
                            log_msgs.append(msg)
                            plt.close()
                            result_msg = "Unable to calculate statistics. Test set be SKIP."
                            median_diff = "skip"
                            return median_diff, result_msg, log_msgs

                        # Difference between pipeline and calculated values
                        delf[k, j] = pipeflat[k, j] - flatcor[k, j]
                        # difference between pipeline errors array and calculated values
                        delflaterr[k, j] = pipeflat_err[k, j] - flat_err[k, j]

                        if debug:
                            # Uncomment to print where the differeces are too big, and see where we differ with the pipeline
                            if abs(delf[k, j]) >= 1.0 or abs(delflaterr[k, j]) >= 1.0:
                                print('using SCI extension of interpolatedflat number: ', ext)
                                print("wave[k, j] = ", wave[k, j])
                                print("dfs, dff = ", dfs, dff)
                                print("sfs, sff = ", sfs, sff)
                                print("fff, ffs = ", fff, ffs)
                                print('dflat dq flag = ', dfimdq[pind[0], pind[1]])
                                print('sflat dq flag = ', sfimdq[:, pind[0], pind[1]][ibr])
                                print('pind[0], pind[1] = ', pind[0], pind[1])
                                print('x, y, pipeflat, calcflat, diff: ')
                                print(j+1, k+1, pipeflat[k, j], flatcor[k, j], delf[k, j])
                                print('is the calc point somwhere in the pipeflat?', np.where(pipeflat == flatcor[k, j]))
                                print('pipeflat_err, calcflat_err: ')
                                print(pipeflat_err[k, j], flat_err[k, j])
                                print('dff_err, dfs_err, sff_err, sfs_err, fff_err, ffs_err :')
                                print(dff_err, dfs_err, sff_err, sfs_err, fff_err, ffs_err)
                                print()
                                input()

                        # Remove all pixels with values=1 (outside slit boundaries) for statistics
                        if pipeflat[k, j] == 1:
                            delf[k, j] = 999.0
                            delflaterr[k, j] = 999.0

        # only keep points in slitlet
        delfg = delf[np.where(delf != 999.0)]
        delflaterr = delflaterr[np.where(delf != 999.0)]
        # attempt remove outliers, for better statistics, only use points where pipe-calc <= 1.0
        outliers_idx = np.where(np.absolute(delfg) <= 1.0)
        if len(outliers_idx) >= len(delfg)/2.0:
            # the remaining points more than half the original number, remove outliers
            delfg = delfg[outliers_idx]
        # same with the error differences
        outliers_idx = np.where(np.absolute(delflaterr) <= 1.0)
        if len(outliers_idx) >= len(delflaterr)/2.0:
            delflaterr = delflaterr[outliers_idx]
        test_result = "FAILED"
        if delf.size == 0:
            msg1 = " * Unable to calculate statistics because difference array has all values as NaN."
            msg2 = "   Test will be set to FAILED and NO plots will be made."
            print(msg1)
            print(msg2)
            log_msgs.append(msg1)
            log_msgs.append(msg2)
        else:
            msg = "Calculating statistics... "
            print(msg)
            log_msgs.append(msg)
            if delfg.size == 0:
                msg1 = " * Unable to calculate statistics because difference array has all outlier values."
                msg2 = "   Test will be set to FAILED and NO plots will be made."
                print(msg1)
                print(msg2)
                log_msgs.append(msg1)
                log_msgs.append(msg2)
            else:
                stats_and_strings = auxfunc.print_stats(delfg, "Flat Difference", float(threshold_diff), absolute=True)
                stats, stats_print_strings = stats_and_strings
                delfg_mean, delfg_median, delfg_std = stats
                _ = auxfunc.print_stats(delflaterr, "Flat Error Difference", float(threshold_diff), absolute=True)

                # This is the key argument for the assert pytest function
                median_diff = False
                if abs(delfg_median) <= float(threshold_diff):
                    median_diff = True
                if median_diff:
                    test_result = "PASSED"
                else:
                    test_result = "FAILED"

                # make histogram
                flatcor_copy = deepcopy(flatcor)
                flatcor_copy[np.where(flatcor_copy == 999.0)] = np.nan
                flat_err_copy = deepcopy(flat_err)
                flat_err_copy[np.where(flatcor == 999.0)] = np.nan
                if save_figs or show_figs:
                    msg = "Making histogram plot for this slitlet..."
                    print(msg)
                    log_msgs.append(msg)
                    # set the plot variables
                    main_title = filt+"   "+grat+"   SLIT="+slit_id+"\n"
                    bins = None   # binning for the histograms, if None the function will select them automatically
                    #             lolim_x, uplim_x, lolim_y, uplim_y
                    plt_origin = None

                    # Residuals img and histogram
                    title = main_title+"Residuals"
                    info_img = [title, "x (pixels)", "y (pixels)"]
                    xlabel, ylabel = "flat$_{pipe}$ - flat$_{calc}$", "N"
                    info_hist = [xlabel, ylabel, bins, stats]
                    if delfg[1] is np.nan:
                        msg = "Unable to create plot of relative wavelength difference."
                        print(msg)
                        log_msgs.append(msg)
                    else:
                        t = (file_basename, "MOS_flattest_"+slitlet_id+"_histogram.png")
                        plt_name = "_".join(t)
                        plt_name = os.path.join(file_path, plt_name)
                        difference_img = pipeflat - flatcor_copy
                        # set the range of values to be shown in the image, will affect color scale
                        vminmax = [-5*delfg_std, 5*delfg_std]
                        auxfunc.plt_two_2Dimgandhist(difference_img, delfg, info_img, info_hist, plt_name=plt_name,
                                                     vminmax=vminmax, plt_origin=plt_origin, show_figs=show_figs,
                                                     save_figs=save_figs)

                elif not save_figs and not show_figs:
                    msg = "Not making plots because both show_figs and save_figs were set to False."
                    print(msg)
                    log_msgs.append(msg)
                elif not save_figs:
                    msg = "Not saving plots because save_figs was set to False."
                    print(msg)
                    log_msgs.append(msg)

        msg = " *** Result of the test: "+test_result+"\n"
        print(msg)
        log_msgs.append(msg)
        total_test_result.append(test_result)

        # create fits file to hold the calculated flat for each slit
        if writefile:
            flatcor_copy[np.isnan(flatcor_copy)] = 1.0
            # this is the file to hold the image of the correction values
            outfile_ext = fits.ImageHDU(flatcor_copy, name=slitlet_id)
            outfile.append(outfile_ext)
            # this is the file to hold the corresponding error extensions
            outfile_ext = fits.ImageHDU(flat_err_copy, name=slitlet_id+'_ERR')
            outfile.append(outfile_ext)

            # this is the file to hold the image of the comparison values
            complfile_ext = fits.ImageHDU(difference_img, name=slitlet_id)
            complfile.append(complfile_ext)

            # the file is not yet written, indicate that this slit was appended to list to be written
            msg = "Extension corresponding to slitlet "+slitlet_id+" appended to list to be written into calculated " \
                                                                   "and comparison fits files."
            print(msg)
            log_msgs.append(msg)

    flatfile_hdu.close()

    if writefile:
        outfile_name = flatfile.replace("interpolatedflat.fits", "flat_calc.fits")
        complfile_name = flatfile.replace("interpolatedflat.fits", "flat_comp.fits")

        # this is the file to hold the image of pipeline-calculated difference values
        outfile.writeto(outfile_name, overwrite=True)

        # this is the file to hold the image of pipeline-calculated difference values
        complfile.writeto(complfile_name, overwrite=True)

        msg = "\nFits file with calculated flat values of each slit saved as: "
        print(msg)
        print(outfile_name)
        log_msgs.append(msg)
        log_msgs.append(outfile_name)

        msg = "Fits file with comparison (pipeline flat - calculated flat) saved as: "
        print(msg)
        print(complfile_name)
        log_msgs.append(msg)
        log_msgs.append(complfile_name)

    # If all tests passed then pytest will be marked as PASSED, else it will be FAILED
    FINAL_TEST_RESULT = True
    for t in total_test_result:
        if t == "FAILED":
            FINAL_TEST_RESULT = False
            break
    if FINAL_TEST_RESULT:
        msg = "\n *** Final result for flat_field test will be reported as PASSED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "All slitlets PASSED flat_field test."
    else:
        msg = "\n *** Final result for flat_field test will be reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "One or more slitlets FAILED flat_field test."

    # Total error calculation according to equation taken from:
    # https://jwst-pipeline.readthedocs.io/en/latest/jwst/flatfield/main.html
    auxfunc.calc_flat_total_errs(step_input_filename, show_plts=show_figs, save_plts=save_figs)

    # end the timer
    flattest_end_time = time.time() - flattest_start_time
    if flattest_end_time > 60.0:
        flattest_end_time = flattest_end_time/60.0  # in minutes
        flattest_tot_time = "* Script flattest_mos.py took ", repr(flattest_end_time)+" minutes to finish."
        if flattest_end_time > 60.0:
            flattest_end_time = flattest_end_time/60.  # in hours
            flattest_tot_time = "* Script flattest_mos.py took ", repr(flattest_end_time)+" hours to finish."
    else:
        flattest_tot_time = "* Script flattest_mos.py took ", repr(flattest_end_time)+" seconds to finish."
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
    parser.add_argument("msa_shutter_conf",
                        action='store',
                        default=None,
                        help='Name of the MSA shutter configuration file in pipeline format, e.g. blah_msa.fits')
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
    msa_shutter_conf = args.msa_shutter_conf
    writefile = args.writefile
    save_figs = args.save_figs
    show_figs = args.show_figs
    threshold_diff = args.threshold_diff
    debug = args.debug

    # Run the principal function of the script
    flattest(step_input_filename, dflat_path=dflat_path, sflat_path=sflat_path, fflat_path=fflat_path,
             msa_shutter_conf=msa_shutter_conf, writefile=writefile, show_figs=show_figs, save_figs=save_figs,
             threshold_diff=threshold_diff, debug=debug)


if __name__ == '__main__':
    sys.exit(main())
