from __future__ import print_function, division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

import auxiliary_functions as auxfunc


"""
This script tests the pipeline flat field step output for MOS data. It is the python version of the IDL script
(with the same name) written by James Muzerolle.
"""


def reverse_cols(arr):
    """
    This function permutates the last column of the array with the first, e.g. a = [4,5,6]
    b = reverse_cols(a) = [6,5,4].
    Args:
        arr: numpy array

    Returns:
        rev_arr: numpy array with first and last columns reversed
    """
    last_idx = np.shape(arr)[-1]-1
    permutation = [last_idx]
    for i, a in enumerate(arr):
        if (i != 0) and (i != last_idx):
            permutation.append(i)
        if i == last_idx:
            permutation.append(0)
    p = np.argsort(permutation)
    rev_arr = arr[:, p]
    return rev_arr


def flattest(step_input_filename, dflatref_path=None, sfile_path=None, fflat_path=None, writefile=False,
             show_figs=True, save_figs=False, plot_name=None, threshold_diff=1.0e-14, debug=False):
    """
    This function calculates the difference between the pipeline and the calculated flat field values.
    The functions uses the output of the compute_world_coordinates.py script.

    Args:
        step_input_filename: str, name of the output fits file from the 2d_extract step (with full path)
        dflatref_path: str, path of where the D-flat reference fits files
        sfile_path: str, path of where the S-flat reference fits files
        fflat_path: str, path of where the F-flat reference fits files
        writefile: boolean, if True writes the fits files of the calculated flat and difference images
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots (the 3 plots can be saved or not independently with the function call)
        plot_name: string, desired name (if name is not given, the plot function will name the plot by
                    default)
        threshold_diff: float, threshold difference between pipeline output and ESA file
        debug: boolean, if true a series of print statements will show on-screen

    Returns:
        - 1 plot, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to threshold.

    """

    # get info from the rate file header
    det = fits.getval(step_input_filename, "DETECTOR", 0)
    print('step_input_filename=', step_input_filename)
    exptype = fits.getval(step_input_filename, "EXP_TYPE", 0)
    grat = fits.getval(step_input_filename, "GRATING", 0)
    filt = fits.getval(step_input_filename, "FILTER", 0)
    print ("rate_file  -->     Grating:", grat, "   Filter:", filt, "   EXP_TYPE:", exptype)

    # read in the on-the-fly flat image
    flatfile = step_input_filename.replace("2d_flat_field.fits", "intflat.fits")
    #flatfile = step_input_filename    # this is for testing purposes only

    # get the reference files
    # D-Flat
    dflat_ending = "f_01.03.fits"
    t = (dflatref_path, "nrs1", dflat_ending)
    dfile = "_".join(t)
    if det == "NRS2":
        dfile = dfile.replace("nrs1", "nrs2")
    dfim = fits.getdata(dfile, 1)
    dfimdq = fits.getdata(dfile, 4)
    # need to flip/rotate the image into science orientation
    ns = np.shape(dfim)
    dfim = np.transpose(dfim, (0, 2, 1))   # keep in mind that 0,1,2 = z,y,x in Python, whereas =x,y,z in IDL
    dfimdq = np.transpose(dfimdq)
    if det == "NRS2":
        dfim = reverse_cols(dfim)
        dfim = dfim[::-1]
        dfimdq = reverse_cols(dfimdq)
        dfimdq = dfimdq[::-1]
    naxis3 = fits.getval(dfile, "NAXIS3", 1)

    # get the wavelength values
    dfwave = np.array([])
    for i in range(naxis3):
        t = ("PFLAT", str(i+1))
        keyword = "_".join(t)
        dfwave = np.append(dfwave, fits.getval(dfile, keyword, 1))
    dfrqe = fits.getdata(dfile, 2)

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
        print ("No filter correspondence. Exiting the program.")
        exit()
    sflat_ending = "f_01.01.fits"
    if mode in sfile_path:
        t = (sfile_path, grat, "OPAQUE", flat, "nrs1", sflat_ending)
        sfile = "_".join(t)
    else:
        print ("Wrong path in for mode S-flat. This script handles mode ", mode, "only.")
        # This is the key argument for the assert pytest function
        median_diff = None
        return median_diff

    if debug:
        print ("grat = ", grat)
        print ("flat = ", flat)
        print ("sfile used = ", sfile)

    if det == "NRS2":
        sfile = sfile.replace("nrs1", "nrs2")
    sfim = fits.getdata(sfile, 1)
    sfimdq = fits.getdata(sfile, 3)

    # need to flip/rotate image into science orientation
    sfim = np.transpose(sfim, (0, 2, 1))
    sfimdq = np.transpose(sfimdq, (0, 2, 1))
    if det == "NRS2":
        sfim = reverse_cols(sfim)
        sfim = sfim[::-1]
        sfimdq = reverse_cols(sfimdq)
        sfimdq = sfimdq[::-1]

    if det == "NRS1":
        sfv_a2001 = fits.getdata(sfile, 5)
        sfv_a2002 = fits.getdata(sfile, 6)
        sfv_a400 = fits.getdata(sfile, 7)
        sfv_a1600 = fits.getdata(sfile, 8)
    elif det == "NRS2":
        sfv_b200 = fits.getdata(sfile, 5)

    # F-Flat
    fflat_ending = "01.01.fits"
    if mode in fflat_path:
        ffile = "_".join((fflat_path, filt, fflat_ending))
    else:
        print ("Wrong path in for mode F-flat. This script handles mode ", mode, "only.")
        # This is the key argument for the assert pytest function
        median_diff = None
        return median_diff
    ffv = fits.getdata(ffile, 1)

    # now go through each pixel in the test data
    wc_file_name = step_input_filename.replace("_flat_field.fits", "_world_coordinates.fits")
    wc_hdulist = fits.open(wc_file_name)
    # loop over the slits
    for i, _ in enumerate(wc_hdulist):
        ext = i+1
        if ext >= len(wc_hdulist):
            break
        slit_id = fits.getval(wc_file_name, "SLIT", ext)

        # select the appropriate S-flat fast vector
        if slit_id == "S200A1":
            sfv = sfv_a2001
        if slit_id == "S200A2":
            sfv = sfv_a2002
        if slit_id == "S400A1":
            sfv = sfv_a400
        if slit_id == "S1600A1":
            sfv = sfv_a1600
        if slit_id == "S200B1":
            sfv = sfv_b200

        # get the wavelength
        wc_data = fits.getdata(wc_file_name, ext)
        wave = wc_data[0, :, :]

        # get the subwindow origin
        px0 = int(fits.getval(wc_file_name, "CRVAL1", ext))-1
        py0 = int(fits.getval(wc_file_name, "CRVAL2", ext))-1
        n_p = np.shape(wave)
        nw = n_p[0]*n_p[1]
        nw1, nw2 = n_p[1], n_p[0]   # remember that x and y are reversed  in Python
        if debug:
            print ("subwindow origin:   px0=",px0, "   py0=", py0)
            #print ("nw = ", nw)
        delf = np.zeros([n_p[1], n_p[0]]) + 999.0
        flatcor = np.zeros([n_p[1], n_p[0]]) + 999.0

        # loop through the wavelengths
        for j in range(nw1):   # in x
            for k in range(nw2):   # in y
                if np.isfinite(wave[k, j]):   # skip if wavelength is NaN
                    # get thr full-frame pixel indeces for D- and S-flat image components
                    # **** this does not account for subarrays!
                    pind = [k+py0-1, j+px0-1]
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

                    # integrate over D-flat fast vector
                    dfrqe_wav = dfrqe.field("WAVELENGTH")
                    dfrqe_rqe = dfrqe.field("RQE")
                    iw = np.where((dfrqe_wav >= wave[k, j]-delw/2.) & (dfrqe_wav <= wave[k, j]+delw/2.))
                    int_tab = auxfunc.idl_tabulate(dfrqe_wav[iw], dfrqe_rqe[iw])
                    first_dfrqe_wav, last_dfrqe_wav = dfrqe_wav[iw[0]][0], dfrqe_wav[iw[0]][-1]
                    dff = int_tab/(last_dfrqe_wav - first_dfrqe_wav)

                    if debug:
                        #print ("np.shape(dfrqe_wav) : ", np.shape(dfrqe_wav))
                        #print ("np.shape(dfrqe_rqe) : ", np.shape(dfrqe_rqe))
                        #print ("dfimdq[pind[0]][pind[1]] : ", dfimdq[pind[0]][pind[1]])
                        #print ("np.shape(iw) =", np.shape(iw))
                        #print ("np.shape(dfrqe_wav[iw[0]]) = ", np.shape(dfrqe_wav[iw[0]]))
                        #print ("np.shape(dfrqe_rqe[iw[0]]) = ", np.shape(dfrqe_rqe[iw[0]]))
                        #print ("int_tab=", int_tab)
                        print ("dff = ", dff)

                    # interpolate over D-flat cube
                    iloc = auxfunc.idl_valuelocate(dfwave, wave[k, j])[0]
                    ibr = [iloc]
                    if iloc != len(dfwave)-1:
                        ibr.append(iloc+1)
                    # get the values in the z-array at indeces ibr, and x=pind[1] and y=pind[0]
                    zz = dfim[:, pind[0], pind[1]][[ibr]]
                    # now determine the length of the array with only the finite numbers
                    zzwherenonan = np.where(np.isfinite(zz))
                    kk = np.size(zzwherenonan)
                    dfs = 1.0
                    if (wave[k, j] <= max(dfwave)) and (wave[k, j] >= min(dfwave)) and (kk == 2):
                        dfs = np.interp(wave[k, j], dfwave[ibr], zz[zzwherenonan])
                    # check DQ flags
                    if dfimdq[pind[0]][pind[1]] != 0:
                        dfs = 1.0

                    # integrate over S-flat fast vector
                    sfv_wav = sfv.field("WAVELENGTH")
                    sfv_dat = sfv.field("DATA")
                    iw = np.where((sfv_wav >= wave[k, j]-delw/2.0) & (sfv_wav <= wave[k, j]+delw/2.0))
                    sff = 1.0
                    if np.size(iw) >= 2:
                        int_tab = auxfunc.idl_tabulate(sfv_wav[iw], sfv_dat[iw])
                        first_sfv_wav, last_sfv_wav = sfv_wav[iw[0]][0], sfv_wav[iw[0]][-1]
                        sff = int_tab/(last_sfv_wav - first_sfv_wav)
                    # get s-flat pixel-dependent correction
                    sfs = 1.0
                    if sfimdq[pind[0], pind[1]] == 0:
                        sfs = sfim[pind[0], pind[1]]

                    # integrate over F-flat fast vector
                    # reference file blue cutoff is 1 micron, so need to force solution for shorter wavs
                    ffv_wav = ffv.field("WAVELENGTH")
                    ffv_dat = ffv.field("DATA")
                    fff = 1.0
                    if wave[k, j]-delw/2.0 >= 1.0:
                        iw = np.where((ffv_wav >= wave[k, j]-delw/2.0) & (ffv_wav <= wave[k, j]+delw/2.0))
                        if np.size(iw) > 1:
                            int_tab = auxfunc.idl_tabulate(ffv_wav[iw], ffv_dat[iw])
                            first_ffv_wav, last_ffv_wav = ffv_wav[iw[0]][0], ffv_wav[iw[0]][-1]
                            fff = int_tab/(last_ffv_wav - first_ffv_wav)

                    flatcor[k, j] = dff * dfs * sff * sfs * fff

                # read the pipeline-calculated flat image
                pipeflat = fits.getdata(flatfile, ext+(ext-1)*2)

                # Difference between pipeline and calculated values
                delf[k, j] = pipeflat[k, j] - flatcor[k, j]

                # Remove all pixels with values=1 (outside slit boundaries) for statistics
                if pipeflat[k, j] == 1:
                    delf[j] = 999.0
                else:
                    flatcor[k, j] = 1.0   # no correction if no wavelength


        delfg = delf[np.where((delf != 999.0) & (delf < 0.1) & (delf > -0.1))]   # ignore outliers
        delfg_median, delfg_std = np.median(delfg), np.std(delfg)
        print ("median, stdev in flat value differences: ", delfg_median, delfg_std)

        # This is the key argument for the assert pytest function
        median_diff = False
        if delfg_median <= float(threshold_diff):
            median_diff = True

        # make histogram
        font = {#'family' : 'normal',
                'weight' : 'normal',
                'size'   : 16}
        matplotlib.rc('font', **font)
        alpha = 0.2
        fontsize = 15
        fig = plt.figure(1, figsize=(8, 6))
        plt.subplots_adjust(hspace=.4)
        ax = plt.subplot(111)
        t = (filt, grat, "  SLIT"+slit_id)
        plt.title("_".join(t))
        plt.xlabel("flat$_{pipe}$ - flat$_{calc}$")
        plt.ylabel("N")
        xmin = delfg_median - delfg_std*5
        xmax = delfg_median + delfg_std*5
        plt.xlim(xmin, xmax)
        x_median = "median = {:0.3}".format(delfg_median)
        x_stddev = "stddev = {:0.3}".format(delfg_std)
        ax.text(0.7, 0.9, x_median, transform=ax.transAxes, fontsize=fontsize)
        ax.text(0.7, 0.83, x_stddev, transform=ax.transAxes, fontsize=fontsize)
        #_, _, _ = ax.hist(delfg, bins=int((xmax-xmin)/40.), histtype='bar', ec='k', facecolor="red", alpha=alpha)
        _, _, _ = ax.hist(delfg, bins=30, histtype='bar', ec='k', facecolor="red", alpha=alpha)

        if save_figs:
            if plot_name is None:
                file_basename = step_input_filename.replace(".fits", "")
                t = (file_basename, "FS_flattest_histogram.pdf")
                plot_name = "_".join(t)
            plt.savefig(plot_name)
            print ('\n Plot saved: ', plot_name)
        if show_figs:
            plt.show()
        plt.close()


        # create fits file to hold the calculated flat for each slit
        if writefile:
            outfile = step_input_filename.replace("2d_flat_field.fits", det+"_flat_calc.fits")
            main_hdr = fits.getheader(step_input_filename, 0)
            #ext_hdr = fits.getheader(msafile, 2)
            hdu = fits.PrimaryHDU(main_hdr)
            hdulist = fits.HDUList([hdu, flatcor])
            hdulist.writeto(outfile)

            # this is the file to hold the image of pipeline-calculated difference values
            complfile = step_input_filename.replace("2d_flat_field.fits", det+"_flat_comp.fits")
            hdu = fits.PrimaryHDU(main_hdr)
            hdulist = fits.HDUList([hdu, delf])
            hdulist.writeto(complfile)


    return median_diff



if __name__ == '__main__':

    # This is a simple test of the code
    pipeline_path = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline"

    # input parameters that the script expects
    auxiliary_code_path = pipeline_path+"/src/pytests/calwebb_spec2_pytests/auxiliary_code"
    step_input_filename = "jwtest1010001_01101_00001_NRS1_rate_short_assign_wcs_extract_2d_flat_field.fits"
    #dflatref_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.2_D_Flat/nirspec_dflat"
    sfile_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.3_S_Flat/FS/nirspec_FS_sflat"
    fflat_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.1_F_Flat/FS/nirspec_FS_fflat"
    dflatref_path = "nirspec_dflat"
    #sfile_path = "nirspec_FS_sflat"
    #fflat_path = "nirspec_FS_fflat"

    # name of the output images
    writefile = True

    # set the names of the resulting plots
    plot_name = "FS_flattest_histogram.pdf"

    # Run the principal function of the script
    median_diff = flattest(step_input_filename, dflatref_path=dflatref_path, sfile_path=sfile_path,
                           fflat_path=fflat_path, writefile=writefile, show_figs=True, save_figs=False,
                           plot_name=plot_name, threshold_diff=1.0e-14, debug=False)

