from __future__ import print_function, division
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

from .. auxiliary_code import auxiliary_functions as auxfunc


"""
This script tests the pipeline flat field step output for IFU data. It is the python version of the IDL script
(with the same name) written by James Muzerolle, and changes on it made by Ben Sargent.
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


def mk_hist(title, delfg, delfg_median, delfg_std, save_figs, show_figs, plot_name):
    # create histogram
    font = {#'family' : 'normal',
            'weight' : 'normal',
            'size'   : 16}
    matplotlib.rc('font', **font)
    alpha = 0.2
    fontsize = 15
    fig = plt.figure(1, figsize=(12, 10))
    plt.subplots_adjust(hspace=.4)
    ax = plt.subplot(111)
    plt.title(title)
    plt.xlabel("flat$_{pipe}$ - flat$_{calc}$")
    plt.ylabel("N")
    xmin = min(delfg) - (max(delfg) - min(delfg))*0.1
    xmax = max(delfg) + (max(delfg) - min(delfg))*0.1
    plt.xlim(xmin, xmax)
    x_median = "median = {:0.3}".format(delfg_median)
    x_stddev = "stddev = {:0.3}".format(delfg_std)
    ax.text(0.7, 0.9, x_median, transform=ax.transAxes, fontsize=fontsize)
    ax.text(0.7, 0.83, x_stddev, transform=ax.transAxes, fontsize=fontsize)
    #_, _, _ = ax.hist(delfg, bins=int((xmax-xmin)/40.), histtype='bar', ec='k', facecolor="red", alpha=alpha)
    _, _, _ = ax.hist(delfg, bins=30, histtype='bar', ec='k', facecolor="red", alpha=alpha)

    if save_figs:
        if plot_name is None:
            t = (title, ".pdf")
            plot_name = "".join(t)
        plt.savefig(plot_name)
        print ('\n Plot saved: ', plot_name)
    if show_figs:
        plt.show()
    plt.close()



def flattest(step_input_filename, dflatref_path=None, sfile_path=None, fflat_path=None, writefile=False,
             mk_all_slices_plt=False, show_figs=True, save_figs=False, plot_name=None,
             threshold_diff=1.0e-14, debug=False):
    """
    This function calculates the difference between the pipeline and the calculated flat field values.
    The functions uses the output of the compute_world_coordinates.py script.

    Args:
        step_input_filename: str, name of the output fits file from the 2d_extract step (with full path)
        dflatref_path: str, path of where the D-flat reference fits files
        sfile_path: str, path of where the S-flat reference fits files
        fflat_path: str, path of where the F-flat reference fits files
        msa_conf_root: str, path to where the MSA configuration fits file lives
        writefile: boolean, if True writes the fits files of the calculated flat and difference images
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots (the 3 plots can be saved or not independently with the function call)
        plot_name: string, desired name (if name is not given, the plot function will name the plot by
                    default)
        threshold_diff: float, threshold difference between pipeline output and ESA file
        debug: boolean, if true a series of print statements will show on-screen

    Returns:
        - 1 plot, if told to save and/or show.
        - median_diff: Boolean, True if smaller or equal to 1e-14

    """

    # get info from the rate file header
    det = fits.getval(step_input_filename, "DETECTOR", 0)
    exptype = fits.getval(step_input_filename, "EXP_TYPE", 0)
    grat = fits.getval(step_input_filename, "GRATING", 0)
    filt = fits.getval(step_input_filename, "FILTER", 0)
    file_basename = step_input_filename.replace(".fits", "")
    print('step_input_filename=', step_input_filename)
    print ("rate_file  -->     Grating:", grat, "   Filter:", filt, "   EXP_TYPE:", exptype)

    # read in the on-the-fly flat image
    flatfile = step_input_filename.replace("2d_flat_field.fits", "intflat.fits")
    #flatfile = step_input_filename    # this is for testing purposes only
    pipeflat = fits.getdata(flatfile, 1)

    # get the reference files
    # D-Flat
    dflat_ending = "f_01.03.fits"
    dfile = dflatref_path+"_nrs1_"+dflat_ending
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
        keyword = "PFLAT_"+str(i+1)
        dfwave = np.append(dfwave, fits.getval(dfile, keyword, 1))
    dfrqe = fits.getdata(dfile, 2)

    # S-flat
    tsp = exptype.split("_")
    mode = tsp[1]
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
        # This is the key argument for the assert pytest function
        median_diff = None
        return median_diff

    sflat_ending = "f_01.01.fits"
    sfile = sfile_path+"_"+grat+"_OPAQUE_"+flat+"_nrs1_"+sflat_ending

    if debug:
        print ("grat = ", grat)
        print ("flat = ", flat)
        print ("sfile used = ", sfile)

    if det == "NRS2":
        sfile = sfile.replace("nrs1", "nrs2")
    sfim = fits.getdata(sfile, 1)
    sfimdq = fits.getdata(sfile, 3)

    # need to flip/rotate image into science orientation
    sfim = np.transpose(sfim)
    sfimdq = np.transpose(sfimdq)
    if det == "NRS2":
        sfim = reverse_cols(sfim)
        sfim = sfim[::-1]
        sfimdq = reverse_cols(sfimdq)
        sfimdq = sfimdq[::-1]

    # get the wavelength values for sflat cube
    naxis3 = fits.getval(sfile, "NAXIS3", 1)
    sfimwave = np.array([])
    for i in range(0, naxis3):
        if i+1 < 10:
            keyword = "FLAT_0"+str(i+1)
        else:
            keyword = "FLAT_"+str(i+1)
        sfimwave = np.append(sfimwave, fits.getval(sfile, keyword, 1))
    sfv = fits.getdata(sfile, 5)

    # F-Flat
    fflat_ending = "_01.01.fits"
    if mode in fflat_path:
        ffile = fflat_path+"_"+filt+fflat_ending
    else:
        print ("Wrong path in for mode F-flat. This script handles mode ", mode, "only.")
        # This is the key argument for the assert pytest function
        median_diff = None
        return median_diff
    ffv = fits.getdata(ffile, 1)

    # go through each pixel in the test data
    wc_file_name = step_input_filename.replace("_flat_field.fits", "_world_coordinates.fits")
    wc_hdulist = fits.open(wc_file_name)
    # loop over the slices and read in the WCS values
    for ext, _ in enumerate(wc_hdulist):
        slice_id = fits.getval(wc_file_name, "SLIT", 2)
        wc_data = fits.getdata(wc_file_name, ext)
        # get the wavelength
        wave = wc_data[0, :, :]

        # get the subwindow origin (technically no subwindows for IFU, but need this for comparing to the
        # full frame on-the-fly flat image).
        px0 = int(fits.getval(wc_file_name, "CRVAL1", ext))
        py0 = int(fits.getval(wc_file_name, "CRVAL2", ext))
        n_p = np.shape(wave)
        nx, ny = n_p[1], n_p[0]
        nw = nx * ny
        if debug:
            print ("subwindow origin:   px0=",px0, "   py0=", py0)
            print ("nw = ", nw)
        delf = np.zeros([nw]) + 999.0
        flatcor = np.zeros([nw]) + 999.0
        sffarr = np.zeros([nw])

        # loop through the wavelengths
        print ("looping through the wavelngth, this may take a little time ... ")
        flat_wave = wave.flatten()
        for j in range(nw-1):
            if np.isfinite(flat_wave[j]):   # skip if wavelength is NaN
                # get the pixel indeces
                jwav = flat_wave[j]
                t=np.where(wave == jwav)
                pind = [t[0][0]+py0, t[1][0]+px0]   # pind =[pixel_y, pixe_x] in python, [x, y] in IDL
                if debug:
                    print ('j, jwav, px0, py0 : ', j, jwav, px0, py0)
                    print ('pind = ', pind)

                # get the pixel bandwidth **this needs to be modified for prism, since the dispersion is not linear!**
                delw = 0.0
                if (j!=0) and (int((j-1)/nx)==int(j/nx)) and (int((j+1)/nx)==int(j/nx)) and np.isfinite(flat_wave[j+1]) and np.isfinite(flat_wave[j-1]):
                    delw = 0.5 * (flat_wave[j+1] - flat_wave[j-1])
                if (j==0) or not np.isfinite(flat_wave[j-1]) or (int((j-1)/nx) != int(j/nx)):
                    delw = 0.5 * (flat_wave[j+1] - flat_wave[j])
                if (j==nw-1) or not np.isfinite(flat_wave[j+1]) or (int((j+1)/nx) != int(j/nx)):
                    delw = 0.5 (flat_wave[j] - flat_wave[j-1])

                if debug:
                    #print ("(j, (j-1), nx, (j-1)/nx, (j+1), (j+1)/nx)", j, (j-1), nx, int((j-1)/nx), (j+1), int((j+1)/nx))
                    #print ("np.isfinite(flat_wave[j+1]), np.isfinite(flat_wave[j-1])", np.isfinite(flat_wave[j+1]), np.isfinite(flat_wave[j-1]))
                    #print ("flat_wave[j+1], flat_wave[j-1] : ", np.isfinite(flat_wave[j+1]), flat_wave[j+1], flat_wave[j-1])
                    print ("delw = ", delw)

                # integrate over D-flat fast vector
                dfrqe_wav = dfrqe.field("WAVELENGTH")
                dfrqe_rqe = dfrqe.field("RQE")
                iw = np.where((dfrqe_wav >= jwav-delw/2.0) & (dfrqe_wav <= jwav+delw/2.0))
                int_tab = auxfunc.idl_tabulate(dfrqe_wav[iw[0]], dfrqe_rqe[iw[0]])
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
                dfs = 1.0
                if dfimdq[pind[0], pind[1]] == 0:
                    dfs = np.interp(jwav, dfwave, dfim[:, pind[0], pind[1]])

                # integrate over S-flat fast vector
                sfv_wav = sfv.field("WAVELENGTH")
                sfv_dat = sfv.field("DATA")
                if (jwav < 5.3) and (jwav > 0.6):
                    iw = np.where((sfv_wav >= jwav-delw/2.0) & (sfv_wav <= jwav+delw/2.0))
                    sff = sfv_dat
                    if np.size(iw) != 0:
                        int_tab = auxfunc.idl_tabulate(sfv_wav[iw], sfv_dat[iw])
                        first_sfv_wav, last_sfv_wav = sfv_wav[iw[0]][0], sfv_wav[iw[0]][-1]
                        sff = int_tab/(last_sfv_wav - first_sfv_wav)

                # get s-flat pixel-dependent correction
                sfs = 1.0
                if sfimdq[pind[0], pind[1]] == 0:
                    sfs = sfim[pind[0], pind[1]]

                # integrate over f-flat fast vector
                # reference file blue cutoff is 1 micron, so need to force solution for shorter wavs
                ffv_wav = ffv.field("WAVELENGTH")
                ffv_dat = ffv.field("DATA")
                fff = 1.0
                if jwav-delw/2.0 >= 1.0:
                    iw = np.where((ffv_wav >= jwav-delw/2.0) & (ffv_wav <= jwav+delw/2.0))
                    if np.size(iw) > 1:
                        int_tab = auxfunc.idl_tabulate(ffv_wav[iw], ffv_dat[iw])
                        first_ffv_wav, last_ffv_wav = ffv_wav[iw[0]][0], ffv_wav[iw[0]][-1]
                        fff = int_tab/(last_ffv_wav - first_ffv_wav)
                    else:
                        fff = ffv_dat

                flatcor[j] = dff * dfs * sff * sfs * fff
                sffarr[j] = sff

                # Difference between pipeline and calculated values
                delf[j] = pipeflat[pind[0]-py0, pind[1]-px0] - flatcor[j]

                # Remove all pixels with values=1 (mainly inter-slit pixels) for statistics
                if pipeflat[pind[0]-py0, pind[1]-px0] == 1:
                    delf[j] = 999.0
                else:
                    flatcor[j] = 1.0   # no correction if no wavelength

                if debug:
                    print ("pind = ", pind)
                    print ("jwav = ", jwav)
                    print ("dfs, dff = ", dfs, dff)
                    print ("sfs, sff = ", sfs, sff)

                delfg = delf[np.where((delf != 999.0) & (delf >= -0.1))]   # ignore outliers
                delfg_median, delfg_std = np.median(delfg), np.std(delfg)
                print ("median and stdev in flat value differences for slice #: ", i)
                print ("median = ", delfg_median, "    stdev =", delfg_std)

                # make the slice plot
                if (show_figs) or (save_figs):
                    # create histogram
                    t = (file_basename, det, "IFUflatcomp_hist", repr(i))
                    title = ("_".join(t))
                    mk_hist(title, delfg, delfg_median, delfg_std, save_figs, show_figs, plot_name)

                if debug:
                    print ("dfs = ", dfs)
                    print ("sff = ", sff)
                    print ("sfs = ", sfs)

        wc_hdulist.close()

        delfg = delf[np.where((delf != 999.0) & (delf >= -0.1))]   # ignore outliers
        delfg_median, delfg_std = np.median(delfg), np.std(delfg)
        print ("median, stdev in flat value differences: ", delfg_median, delfg_std)

        if mk_all_slices_plt:
            # create histogram
            t = (file_basename, det, "IFU_flatcomp_hist")
            title = ("_".join(t))
            mk_hist(title, delfg, delfg_median, delfg_std, save_figs, show_figs, plot_name)


        # This is the key argument for the assert pytest function
        median_diff = False
        if delfg_median <= float(threshold_diff):
            median_diff = True


        # create fits file to hold the calculated flat for each slit
        if writefile:
            outfile = step_input_filename.replace("2d_flat_field.fits", det+"_flat_calc.fits")
            main_hdr = fits.getheader(step_input_filename, 0)
            #ext_hdr = fits.getheader(msafile, 2)
            hdu = fits.PrimaryHDU(main_hdr)
            hdulist = fits.HDUList([hdu, flatcor.reshape(n_p)])
            hdulist.writeto(outfile)
            # this is the file to hold the image of pipeline-calculated difference values
            complfile = step_input_filename.replace("2d_flat_field.fits", det+"_flat_comp.fits")
            hdu = fits.PrimaryHDU(main_hdr)
            hdulist = fits.HDUList([hdu, delf.reshape(n_p)])
            hdulist.writeto(complfile)


    return median_diff



if __name__ == '__main__':

    # This is a simple test of the code
    pipeline_path = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline"

    # input parameters that the script expects
    auxiliary_code_path = pipeline_path+"/src/pytests/calwebb_spec2_pytests/auxiliary_code"
    step_input_filename = "jwtest1010001_01101_00001_NRS1_rate_short_assign_wcs_extract_2d_flat_field.fits"
    #dflatref_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.2_D_Flat/nirspec_dflat"
    #sfile_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.3_S_Flat/IFU/nirspec_IFU_sflat"
    #fflat_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.1_F_Flat/IFU/nirspec_IFU_fflat"
    dflatref_path = "nirspec_dflat"
    sfile_path = "nirspec_IFU_sflat"
    fflat_path = "nirspec_IFU_fflat"

    # name of the output images
    writefile = True

    # set the names of the resulting plots
    plot_name = "IFU_flattest_histogram.pdf"

    # Run the principal function of the script
    median_diff = flattest(step_input_filename, dflatref_path=dflatref_path, sfile_path=sfile_path,
                           fflat_path=fflat_path, writefile=writefile, mk_all_slices_plt=False,
                           show_figs=True, save_figs=False, plot_name=plot_name, threshold_diff=1.0e-14, debug=False)

