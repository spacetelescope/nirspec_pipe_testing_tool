import time
import os
import argparse
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

from jwst import datamodels

from . import auxiliary_functions as auxfunc


"""
This script tests the pipeline flat field step output for MOS data. It is the python version of the IDL script
(with the same name) written by James Muzerolle.
"""


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "3.7"

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


def flattest(step_input_filename, dflat_path, sflat_path, fflat_path, msa_shutter_conf,
             writefile=False, show_figs=True, save_figs=False, interpolated_flat=None,
             threshold_diff=1.0e-14, debug=False):
    """
    This function does the WCS comparison from the world coordinates calculated using the
    compute_world_coordinates.py script with the ESA files. The function calls that script.

    Args:
        step_input_filename: str, name of the output fits file from the 2d_extract step (with full path)
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
        if "extract_2d" not in step_input_filename:
            extract2d_wcs_file = step_input_filename.replace("_flat_field", "_extract_2d")
        else:
            extract2d_wcs_file = step_input_filename
        model = datamodels.MultiSlitModel(extract2d_wcs_file)

    else:
        model = step_input_filename

    # read in the on-the-fly flat image
    if interpolated_flat is None:
        flatfile = step_input_filename.replace("_flat_field", "_interpolatedflat")
    else:
        flatfile = interpolated_flat
    # get all the science extensions in the flatfile
    sci_ext_list = auxfunc.get_sci_extensions(flatfile)

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
        model = step_input_filename
        file_path = os.path.dirname(os.path.realpath(__file__))
        file_basename = grat + "_" + filt + "_" + det

    # define the mode
    if "msa" in exptype.lower():
        mode = "MOS"
    else:
        mode = "not MOS data"

    # get the reference files
    # D-Flat
    if ".fits" not in dflat_path:
        dflat_ending = "f_01.03.fits"
        t = (dflat_path, "nrs1", dflat_ending)
        dfile = "_".join(t)
        if det == "NRS2":
            dfile = dfile.replace("nrs1", "nrs2")
    else:
        dfile = dflat_path
    msg = "".join(["Using D-flat: ", dfile])
    print(msg)
    log_msgs.append(msg)
    dfim = fits.getdata(dfile, "SCI")
    dfimdq = fits.getdata(dfile, "DQ")
    # need to flip/rotate the image into science orientation
    ns = np.shape(dfim)
    dfim = np.transpose(dfim, (0, 2, 1))   # keep in mind that 0,1,2 = z,y,x in Python, whereas =x,y,z in IDL
    dfimdq = np.transpose(dfimdq)
    if det == "NRS2":
        # rotate science data by 180 degrees for NRS2
        dfim = dfim[..., ::-1, ::-1]
        dfimdq = dfimdq[..., ::-1, ::-1]
    naxis3 = fits.getval(dfile, "NAXIS3", "SCI")

    # get the wavelength values
    dfwave = np.array([])
    for i in range(naxis3):
        keyword = "_".join(("PFLAT", str(i+1)))
        dfwave = np.append(dfwave, fits.getval(dfile, keyword, "SCI"))
    dfrqe = fits.getdata(dfile, 2)

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

    if ".fits" not in sflat_path:
        sflat_ending = "f_01.01.fits"
        t = (sflat_path, grat, "OPAQUE", flat, "nrs1", sflat_ending)
        sfile = "_".join(t)
        if det == "NRS2":
            sfile = sfile.replace("nrs1", "nrs2")
    else:
        sfile = sflat_path

    msg = "Using S-flat: "+sfile
    print(msg)
    log_msgs.append(msg)

    if mode not in sflat_path:
        msg = "Wrong path in for mode S-flat. This script handles mode " + mode + "only."
        print(msg)
        log_msgs.append(msg)
        # This is the key argument for the assert pytest function
        result_msg = "Wrong path in for mode S-flat. Test skiped because mode is not FS."
        median_diff = "skip"
        return median_diff, result_msg, log_msgs

    sfim = fits.getdata(sfile, "SCI")#1)
    sfimdq = fits.getdata(sfile, "DQ")#3)

    # need to flip/rotate image into science orientation
    sfim = np.transpose(sfim, (0, 2, 1))
    sfimdq = np.transpose(sfimdq, (0, 2, 1))
    if det == "NRS2":
        # rotate science data by 180 degrees for NRS2
        sfim = sfim[..., ::-1, ::-1]
        sfimdq = sfimdq[..., ::-1, ::-1]

    # get the wavelength values for sflat cube
    sfimwave = np.array([])
    naxis3 = fits.getval(sfile, "NAXIS3", "SCI")
    for i in range(0, naxis3):
        if i+1 < 10:
            keyword = "".join(("FLAT_0", str(i+1)))
        else:
            keyword = "".join(("FLAT_", str(i+1)))
        if debug:
            print("S-flat -> using ", keyword)
        try:
            sfimwave = np.append(sfimwave, fits.getval(sfile, keyword, "SCI"))
        except KeyError:
            continue
    sfv = fits.getdata(sfile, 5)

    # F-Flat
    if ".fits" not in fflat_path:
        fflat_ending = "01.01.fits"
        ffile = "_".join((fflat_path, filt, fflat_ending))
    else:
        ffile = fflat_path
    sci_ext = "SCI"
    dq_ext = "DQ"
    err_ext = "ERR"
    fast_var_ext = "FAST_VARIATION"

    if mode not in fflat_path:
        msg = "Wrong path in for mode F-flat. This script handles mode " + mode + "only."
        print(msg)
        log_msgs.append(msg)
        # This is the key argument for the assert pytest function
        median_diff = "skip"
        return median_diff, msg, log_msgs

    msg = "Using F-flat: "+ffile
    print(msg)
    log_msgs.append(msg)
    Q_in_ext = False
    try:
        ffsq1 = fits.getdata(ffile, sci_ext)
    except KeyError:
        sci_ext = "SCI_Q1"
        ffsq1 = fits.getdata(ffile, sci_ext)
        Q_in_ext = True
    naxis3 = fits.getval(ffile, "NAXIS3", sci_ext)
    ffswaveq1 = np.array([])
    for i in range(0, naxis3):
        if i <= 9:
            suff = "".join(("0", str(i)))
        else:
            suff = str(i)
        t = ("FLAT", suff)
        keyword = "_".join(t)
        if debug:
            print("1. F-flat -> ", keyword)
        ffswaveq1 = np.append(ffswaveq1, fits.getval(ffile, keyword, sci_ext))
    if Q_in_ext:
        ffserrq1 = fits.getdata(ffile, "ERR_Q1")
        ffsdqq1 = fits.getdata(ffile, "DQ_Q1")
        ffvq1 = fits.getdata(ffile, "Q1")
        ffsq2 = fits.getdata(ffile, "SCI_Q2")
    else:
        ffserrq1 = fits.getdata(ffile, err_ext)
        ffsdqq1 = fits.getdata(ffile, dq_ext)
        ffvq1 = fits.getdata(ffile, fast_var_ext)
        ffsq2 = fits.getdata(ffile, sci_ext, 2)
    ffswaveq2 = np.array([])
    for i in range(0, naxis3):
        if i <= 9:
            suff = "".join(("0", str(i)))
        else:
            suff = str(i)
        t = ("FLAT", suff)
        keyword = "_".join(t)
        if debug:
            print("2. F-flat -> using ", keyword)
        if Q_in_ext:
            ffswaveq2 = np.append(ffswaveq2, fits.getval(ffile, keyword, "SCI_Q2"))
        else:
            ffswaveq2 = np.append(ffswaveq2, fits.getval(ffile, keyword, sci_ext, 2))
    if Q_in_ext:
        ffserrq2 = fits.getdata(ffile, "ERR_Q2")
        ffsdqq2 = fits.getdata(ffile, "DQ_Q2")
        ffvq2 = fits.getdata(ffile, "Q2")
        ffsq3 = fits.getdata(ffile, "SCI_Q3")
    else:
        ffserrq2 = fits.getdata(ffile, err_ext, 2)
        ffsdqq2 = fits.getdata(ffile, dq_ext, 2)
        ffvq2 = fits.getdata(ffile, fast_var_ext, 2)
        ffsq3 = fits.getdata(ffile, sci_ext, 3)
    ffswaveq3 = np.array([])
    for i in range(0, naxis3):
        if i <= 9:
            suff = "".join(("0", str(i)))
        else:
            suff = str(i)
        t = ("FLAT", suff)
        keyword = "_".join(t)
        if debug:
            print("3. F-flat -> using ", keyword)
        if Q_in_ext:
            ffswaveq3 = np.append(ffswaveq3, fits.getval(ffile, keyword, "SCI_Q3"))
        else:
            ffswaveq3 = np.append(ffswaveq3, fits.getval(ffile, keyword, sci_ext, 3))
    if Q_in_ext:
        ffserrq3 = fits.getdata(ffile, "ERR_Q3")
        ffsdqq3 = fits.getdata(ffile, "DQ_Q3")
        ffvq3 = fits.getdata(ffile, "Q3")
        ffsq4 = fits.getdata(ffile, "SCI_Q4")
    else:
        ffserrq3 = fits.getdata(ffile, err_ext, 3)
        ffsdqq3 = fits.getdata(ffile, dq_ext, 3)
        ffvq3 = fits.getdata(ffile, fast_var_ext, 3)
        ffsq4 = fits.getdata(ffile, sci_ext, 4)
    ffswaveq4 = np.array([])
    for i in range(0, naxis3):
        if i <= 9:
            suff = "0"+str(i)
        else:
            suff = str(i)
        keyword = "FLAT_"+suff
        if debug:
            print("4. F-flat -> using ", keyword)
        if Q_in_ext:
            ffswaveq4 = np.append(ffswaveq4, fits.getval(ffile, keyword, "SCI_Q4"))
        else:
            ffswaveq4 = np.append(ffswaveq4, fits.getval(ffile, keyword, sci_ext, 4))
    if Q_in_ext:
        ffserrq4 = fits.getdata(ffile, "ERR_Q4")
        ffsdqq4 = fits.getdata(ffile, "DQ_Q4")
        ffvq4 = fits.getdata(ffile, "Q4")
    else:
        ffserrq4 = fits.getdata(ffile, err_ext, 4)
        ffsdqq4 = fits.getdata(ffile, dq_ext, 4)
        ffvq4 = fits.getdata(ffile, fast_var_ext, 4)

    # now go through each pixel in the test data

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

    # loop over the 2D subwindows and read in the WCS values
    for slit in model.slits:
        slit_id = slit.name
        msg = "\nWorking with slit: "+slit_id
        print(msg)
        log_msgs.append(msg)
        ext = sci_ext_list[slit_id]   # this is for getting the science extension in the pipeline calculated flat

        # get the wavelength
        y, x = np.mgrid[:slit.data.shape[0], :slit.data.shape[1]]
        ra, dec, wave = slit.meta.wcs(x, y)   # wave is in microns

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

        delf = np.zeros([nw2, nw1]) + 999.0
        flatcor = np.zeros([nw2, nw1]) + 999.0

        # get the slitlet info, needed for the F-Flat
        ext_shutter_info = "SHUTTER_INFO"   # this is extension 2 of the msa file, that has the shutter info
        slitlet_info = fits.getdata(msa_shutter_conf, ext_shutter_info)
        sltid = slitlet_info.field("SLITLET_ID")
        for j, s in enumerate(sltid):
            if s == int(slit_id):
                im = j
                # get the shutter with the source in it
                if slitlet_info.field("BACKGROUND")[im] == "N":
                    isrc = j
        # changes suggested by Phil Hodge
        quad = slit.quadrant  #slitlet_info.field("SHUTTER_QUADRANT")[isrc]
        row = slit.xcen  #slitlet_info.field("SHUTTER_ROW")[isrc]
        col = slit.ycen  #slitlet_info.field("SHUTTER_COLUMN")[isrc]
        slitlet_id = repr(row)+"_"+repr(col)
        msg = 'silt_id='+repr(slit_id)+"   quad="+repr(quad)+"   row="+repr(row)+"   col="+repr(col)+\
              "   slitlet_id="+repr(slitlet_id)
        print(msg)
        log_msgs.append(msg)

        # get the relevant F-flat reference data
        if quad == 1:
            ffsall = ffsq1
            ffsallwave = ffswaveq1
            ffsalldq = ffsdqq1
            ffv = ffvq1
        if quad == 2:
            ffsall = ffsq2
            ffsallwave = ffswaveq2
            ffsalldq = ffsdqq2
            ffv = ffvq2
        if quad == 3:
            ffsall = ffsq3
            ffsallwave = ffswaveq3
            ffsalldq = ffsdqq3
            ffv = ffvq3
        if quad == 4:
            ffsall = ffsq4
            ffsallwave = ffswaveq4
            ffsalldq = ffsdqq4
            ffv = ffvq4

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
                    if debug:
                        print('j, k, jwav, px0, py0 : ', j, k, jwav, px0, py0)
                        print('pind = ', pind)

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

                    if debug:
                        print("wave[k, j+1], wave[k, j-1] : ", np.isfinite(wave[k, j+1]), wave[k, j+1], wave[k, j-1])
                        print("delw = ", delw)

                    # integrate over dflat fast vector
                    dfrqe_wav = dfrqe.field("WAVELENGTH")
                    dfrqe_rqe = dfrqe.field("RQE")
                    iw = np.where((dfrqe_wav >= wave[k, j]-delw/2.0) & (dfrqe_wav <= wave[k, j]+delw/2.0))
                    dff = 1.0
                    if np.size(iw) >= 1:
                        int_tab = auxfunc.idl_tabulate(dfrqe_wav[iw[0]], dfrqe_rqe[iw[0]])
                        first_dfrqe_wav, last_dfrqe_wav = dfrqe_wav[iw[0]][0], dfrqe_wav[iw[0]][-1]
                        dff = int_tab/(last_dfrqe_wav - first_dfrqe_wav)

                    if debug:
                        print("np.shape(dfrqe_wav) : ", np.shape(dfrqe_wav))
                        print("np.shape(dfrqe_rqe) : ", np.shape(dfrqe_rqe))
                        print("dfimdq[pind[0]][pind[1]] : ", dfimdq[pind[0]][pind[1]])
                        print("np.shape(iw) =", np.shape(iw))
                        print("np.shape(dfrqe_wav[iw[0]]) = ", np.shape(dfrqe_wav[iw[0]]))
                        print("np.shape(dfrqe_rqe[iw[0]]) = ", np.shape(dfrqe_rqe[iw[0]]))
                        print("int_tab=", int_tab)
                        print("dff = ", dff)

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
                    dfs = 1.0

                    if (wave[k, j] <= max(dfwave)) and (wave[k, j] >= min(dfwave)) and (kk == 2):
                        dfs = np.interp(wave[k, j], dfwave[ibr], zz[zzwherenonan])

                    # check DQ flags for d-flat
                    if dfimdq[pind[0], pind[1]] == 1:
                        # print("d-flat: DQ flag DO NOT USE forcing dfs=1.0")
                        dfs = 1.0

                    # force d-flat to 1.0 for fixing flat field issue with pipeline
                    #dff, dfs = 1.0, 1.0

                    # integrate over s-flat fast vector
                    sfv_wav = sfv.field("WAVELENGTH")
                    sfv_dat = sfv.field("DATA")
                    iw = np.where((sfv_wav >= wave[k, j]-delw/2.0) & (sfv_wav <= wave[k, j]+delw/2.0))
                    sff = 1.0
                    if np.size(iw) >= 1:
                        int_tab = auxfunc.idl_tabulate(sfv_wav[iw], sfv_dat[iw])
                        first_sfv_wav, last_sfv_wav = sfv_wav[iw[0]][0], sfv_wav[iw[0]][-1]
                        sff = int_tab/(last_sfv_wav - first_sfv_wav)

                    # interpolate s-flat cube
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
                    sfs = 1.0
                    if (wave[k, j] <= max(sfimwave)) and (wave[k, j] >= min(sfimwave)) and (kk == 2):
                        sfs = np.interp(wave[k, j], sfimwave[ibr], zz[zzwherenonan])

                    # check DQ flags for s-flat
                    kk = np.where(sfimdq[:, pind[0], pind[1]][ibr] == 1)
                    if np.size(kk) >= 1:
                        # print("s-flat: DQ flag DO NOT USE forcing sfs=1.0")
                        sfs = 1.0

                    # force s-flat to 1.0 for fixing flat field issue with pipeline
                    #sff, sfs = 1.0, 1.0

                    # integrate over f-flat fast vector
                    # reference file wavelength range is from 0.6 to 5.206 microns, so need to force
                    # solution to 1 for wavelengths outside that range
                    ffv_wav = ffv.field("WAVELENGTH")
                    ffv_dat = ffv.field("DATA")
                    fff = 1.0
                    if (wave[k, j]-delw/2.0 >= 0.6) and (wave[k, j]+delw/2.0 <= 5.206):
                        iw = np.where((ffv_wav >= wave[k, j]-delw/2.0) & (ffv_wav <= wave[k, j]+delw/2.0))
                        if np.size(iw) > 1:
                            int_tab = auxfunc.idl_tabulate(ffv_wav[iw], ffv_dat[iw])
                            first_ffv_wav, last_ffv_wav = ffv_wav[iw[0]][0], ffv_wav[iw[0]][-1]
                            fff = int_tab/(last_ffv_wav - first_ffv_wav)

                    # interpolate over f-flat cube
                    ffs = np.interp(wave[k, j], ffsallwave, ffsall[:, col-1, row-1])

                    # check DQ flags for f-flat
                    kk = np.where(ffsalldq[:, col-1, row-1] == 1)
                    if np.size(kk) >= 1:
                        ffs = 1.0

                    # force f-flat to 1.0 for fixing flat field issue with pipeline
                    #fff, ffs = 1.0, 1.0

                    # add correction
                    flatcor[k, j] = dff * dfs * sff * sfs * fff * ffs

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

                    if debug:
                        print("wave[k, j] = ", wave[k, j])
                        print("dfs, dff = ", dfs, dff)
                        print("sfs, sff = ", sfs, sff)
                        print("fff, ffs = ", fff, ffs)
                        print("flatcor[k, j] = ", flatcor[k, j])

                    # read the pipeline-calculated flat image
                    # there are four extensions in the flatfile: SCI, DQ, ERR, WAVELENGTH
                    pipeflat = fits.getdata(flatfile, ext)

                    try:
                        # Difference between pipeline and calculated values
                        delf[k, j] = pipeflat[k, j] - flatcor[k, j]

                        # Remove all pixels with values=1 (outside slit boundaries) for statistics
                        if pipeflat[k, j] == 1:
                            delf[k, j] = 999.0
                        if np.isnan(wave[k, j]):
                            flatcor[k, j] = 1.0   # no correction if no wavelength

                        if debug:
                            print("flatcor[k, j] = ", flatcor[k, j])
                            print("delf[k, j] = ", delf[k, j])
                    except:
                        IndexError

        nanind = np.isnan(delf)   # get all the nan indexes
        notnan = ~nanind   # get all the not-nan indexes
        delf = delf[notnan]   # get rid of NaNs
        delf_shape = np.shape(delf)
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
            delfg = delf[np.where((delf != 999.0) & (delf < 0.1) & (delf > -0.1))]   # ignore outliers
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

                # This is the key argument for the assert pytest function
                median_diff = False
                if abs(delfg_median) <= float(threshold_diff):
                    median_diff = True
                if median_diff:
                    test_result = "PASSED"
                else:
                    test_result = "FAILED"

                if save_figs or show_figs:
                    # make histogram
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
                        difference_img = (pipeflat - flatcor)#/flatcor
                        in_slit = np.logical_and(difference_img < 900.0, difference_img > -900.0)  # ignore out of slit
                        difference_img[~in_slit] = np.nan   # Set values outside the slit to NaN
                        nanind = np.isnan(difference_img)   # get all the nan indexes
                        difference_img[nanind] = np.nan   # set all nan indexes to have a value of nan
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
            # this is the file to hold the image of the correction values
            outfile_ext = fits.ImageHDU(flatcor.reshape(wave_shape), name=slitlet_id)
            outfile.append(outfile_ext)

            # this is the file to hold the image of the comparison values
            complfile_ext = fits.ImageHDU(delf.reshape(delf_shape), name=slitlet_id)
            complfile.append(complfile_ext)

            # the file is not yet written, indicate that this slit was appended to list to be written
            msg = "Extension corresponding to slitlet "+slitlet_id+" appended to list to be written into calculated " \
                                                                   "and comparison fits files."
            print(msg)
            log_msgs.append(msg)

    if writefile:
        outfile_name = flatfile.replace("interpolatedflat.fits", "_flat_calc.fits")
        complfile_name = flatfile.replace("interpolatedflat.fits", "_flat_comp.fits")

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
