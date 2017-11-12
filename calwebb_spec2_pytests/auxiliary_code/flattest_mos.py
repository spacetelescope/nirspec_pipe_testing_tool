from __future__ import print_function, division
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import fits

import auxiliary_functions as auxfunc


"""
This script tests the pipeline flat field step output for MOS data.
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


def flattest(step_input_filename, dflatref_path=None, sfile_path=None, fflat_path=None, msa_conf_root=None,
             writefile=False, show_figs=True, save_figs=False, plot_name=None, threshold_diff=1.0e-14, debug=False):
    """
    This function does the WCS comparison from the world coordinates calculated using the
    compute_world_coordinates.py script with the ESA files. The function calls that script.

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
        - 3 plots, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to 1e-14

    """

    # get info from the rate file header
    det = fits.getval(step_input_filename, "DETECTOR", 0)
    print('step_input_filename=', step_input_filename)
    msametfl = fits.getval(step_input_filename, "MSAMETFL", 0)
    lamp = fits.getval(step_input_filename, "LAMP", 0)
    exptype = fits.getval(step_input_filename, "EXP_TYPE", 0)
    grat = fits.getval(step_input_filename, "GRATING", 0)
    filt = fits.getval(step_input_filename, "FILTER", 0)
    print ("rate_file  -->     Grating:", grat, "   Filter:", filt, "   Lamp:", lamp)

    # read in the on-the-fly flat image
    #flatfile = step_input_filename.replace("2d_flat_field.fits", "intflat.fits")
    flatfile = step_input_filename    # this is for testing purposes only
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
    for i in xrange(0, naxis3):
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
        exit()
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
    sfim = np.transpose(sfim, (0, 2, 1))
    sfimdq = np.transpose(sfimdq, (0, 2, 1))
    if det == "NRS2":
        sfim = reverse_cols(sfim)
        sfim = sfim[::-1]
        sfimdq = reverse_cols(sfimdq)
        sfimdq = sfimdq[::-1]

    # get the wavelength values for sflat cube
    naxis3 = fits.getval(sfile, "NAXIS3", 1)
    sfimwave = np.array([])
    for i in xrange(0, naxis3):
        if i+1 < 10:
            keyword = "FLAT_0"+str(i+1)
        else:
            keyword = "FLAT_"+str(i+1)
        sfimwave = np.append(sfimwave, fits.getval(sfile, keyword, 1))
    sfv = fits.getdata(sfile, 5)

    # F-Flat
    fflat_ending = "_01.01.fits"
    ffile = fflat_path+"_"+filt+fflat_ending
    ffsq1 = fits.getdata(ffile, 1)
    naxis3 = fits.getval(ffile, "NAXIS3", 1)
    ffswaveq1 = np.array([])
    for i in xrange(0, naxis3):
        if i <= 9 :
            suff = "0"+str(i)
        else:
            suff = str(i)
        keyword = "FLAT_"+suff
        ffswaveq1 = np.append(ffswaveq1, fits.getval(ffile, keyword, 1))
    ffserrq1 = fits.getdata(ffile, 2)
    ffsdqq1 = fits.getdata(ffile, 3)
    ffvq1 = fits.getdata(ffile, 4)
    ffsq2 = fits.getdata(ffile, 1)
    ffswaveq2 = np.array([])
    for i in xrange(0, naxis3):
        if i <= 9:
            suff = "0"+str(i)
        else:
            suff = str(i)
        keyword = "FLAT_"+suff
        ffswaveq2 = np.append(ffswaveq2, fits.getval(ffile, keyword, 1))
    ffserrq2 = fits.getdata(ffile, 2)
    ffsdqq2 = fits.getdata(ffile, 3)
    ffvq2 = fits.getdata(ffile, 4)
    ffsq3 = fits.getdata(ffile, 1)
    ffswaveq3 = np.array([])
    for i in xrange(0, naxis3):
        if i <= 9 :
            suff = "0"+str(i)
        else:
            suff = str(i)
        keyword = "FLAT_"+suff
        ffswaveq3 = np.append(ffswaveq3, fits.getval(ffile, keyword, 1))
    ffserrq3 = fits.getdata(ffile, 2)
    ffsdqq3 = fits.getdata(ffile, 3)
    ffvq3 = fits.getdata(ffile, 4)
    ffsq4 = fits.getdata(ffile, 1)
    ffswaveq4 = np.array([])
    for i in xrange(0, naxis3):
        if i <= 9:
            suff = "0"+str(i)
        else:
            suff = str(i)
        keyword = "FLAT_"+suff
        ffswaveq4 = np.append(ffswaveq4, fits.getval(ffile, keyword, 1))
    ffserrq4 = fits.getdata(ffile, 2)
    ffsdqq4 = fits.getdata(ffile, 3)
    ffvq4 = fits.getdata(ffile, 4)

    # go through each pixel in the test data
    wc_file_name = step_input_filename.replace("_flat_field.fits", "_world_coordinates.fits")
    wc_hdulist = fits.open(wc_file_name)
    # loop over the 2D subwindows and read in the WCS values
    for i, _ in enumerate(wc_hdulist):
        ext = i+1
        if ext >= len(wc_hdulist):
            break
        slit_id = fits.getval(wc_file_name, "SLIT", ext)
        wc_data = fits.getdata(wc_file_name, ext)
        # get the wavelength
        wave = wc_data[0, :, :]

        # get the subwindow origin
        px0 = int(fits.getval(wc_file_name, "CRVAL1", ext))-1
        py0 = int(fits.getval(wc_file_name, "CRVAL2", ext))-1
        n_p = np.shape(wave)
        nw = n_p[0]*n_p[1]
        if debug:
            print ("subwindow origin:   px0=",px0, "   py0=", py0)
            #print ("nw = ", nw)
        delf = np.zeros([nw]) + 999.0
        flatcor = np.zeros([nw]) + 999.0

        # get the slitlet info, needed for the F-Flat
        ext_shutter_info = "SHUTTER_INFO"   # this is extension 2 of the msa file, that has the shutter info
        msafile = os.path.join(msa_conf_root, msametfl)
        slitlet_info = fits.getdata(msafile, ext_shutter_info)
        sltid = slitlet_info.field("SLITLET_ID")
        for j, s in enumerate(sltid):
            if s == int(slit_id):
                im = j
                # get the shutter with the source in it
                if slitlet_info.field("BACKGROUND")[im] == "N":
                    isrc = j
        quad = slitlet_info.field("SHUTTER_QUADRANT")[im]
        row = slitlet_info.field("SHUTTER_ROW")[im]
        col = slitlet_info.field("SHUTTER_COLUMN")[im]

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
        print ("looping through the pixels, this may take a little time ... ")
        flat_wave = wave.flatten()
        for j in xrange(nw-1):
            if np.isfinite(flat_wave[j]):   # skip if wavelength is NaN
                # get the pixel indeces
                jwav = flat_wave[j]
                t=np.where(wave == jwav)
                pind = [t[0][0]+py0, t[1][0]+px0]
                if debug:
                    print ('j, jwav, px0, py0 : ', j, jwav, px0, py0)
                    print ('pind = ', pind)

                # get the pixel bandwidth **this needs to be modified for prism, since the dispersion is not linear!**
                delw = 0.0
                if (j!=0) and (int((j-1)/n_p[1])==int(j/n_p[1])) and (int((j+1)/n_p[1])==int(j/n_p[1])) and np.isfinite(flat_wave[j+1]) and np.isfinite(flat_wave[j-1]):
                    delw = 0.5 * (flat_wave[j+1] - flat_wave[j-1])
                if (j==0) or not np.isfinite(flat_wave[j-1]) or (int((j-1)/n_p[1]) != int(j/n_p[1])):
                    delw = flat_wave[j+1] - flat_wave[j]
                if (j==nw-1) or not np.isfinite(flat_wave[j+1]) or (int((j+1)/n_p[1]) != int(j/n_p[1])):
                    delw = flat_wave[j] - flat_wave[j-1]

                if debug:
                    #print ("(j, (j-1), n_p[1], (j-1)/n_p[1], (j+1), n_p[1], (j+1)/n_p[1])", j, (j-1), n_p[1], int((j-1)/n_p[1]), (j+1), n_p[1], int((j+1)/n_p[1]))
                    #print ("np.isfinite(flat_wave[j+1]), np.isfinite(flat_wave[j-1])", np.isfinite(flat_wave[j+1]), np.isfinite(flat_wave[j-1]))
                    #print ("flat_wave[j+1], flat_wave[j-1] : ", np.isfinite(flat_wave[j+1]), flat_wave[j+1], flat_wave[j-1])
                    print ("delw = ", delw)

                # integrate over dflat fast vector
                dfrqe_wav = dfrqe.field("WAVELENGTH")
                dfrqe_rqe = dfrqe.field("RQE")
                iw = np.where((dfrqe_wav >= flat_wave[j]-delw/2.0) & (dfrqe_wav <= flat_wave[j]+delw/2.0))
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

                # interpolate over dflat cube
                iloc = auxfunc.idl_valuelocate(dfwave, flat_wave[j])[0]
                ibr = [iloc]
                if iloc != len(dfwave)-1:
                    ibr.append(iloc+1)
                # get the values in the z-array at indeces ibr, and x=pind[1] and y=pind[0]
                zz = dfim[:, pind[0], pind[1]][[ibr]]
                # now determine the length of the array with only the finite numbers
                zzwherenonan = np.where(np.isfinite(zz))
                kk = np.size(zzwherenonan)
                dfs = 1.0
                if (flat_wave[j] <= max(dfwave)) and (flat_wave[j] >= min(dfwave)) and (kk == 2):
                    dfs = np.interp(flat_wave[j], dfwave[ibr], zz[zzwherenonan])
                # check DQ flags
                if dfimdq[pind[0]][pind[1]] != 0:
                    dfs = 1.0

                # integrate over S-flat fast vector
                sfv_wav = sfv.field("WAVELENGTH")
                sfv_dat = sfv.field("DATA")
                iw = np.where((sfv_wav >= flat_wave[j]-delw/2.0) & (sfv_wav <= flat_wave[j]+delw/2.0))
                sff = 1.0
                if np.size(iw) != 0:
                    int_tab = auxfunc.idl_tabulate(sfv_wav[iw], sfv_dat[iw])
                    first_sfv_wav, last_sfv_wav = sfv_wav[iw[0]][0], sfv_wav[iw[0]][-1]
                    sff = int_tab/(last_sfv_wav - first_sfv_wav)

                # interpolate s-flat cube
                iloc = auxfunc.idl_valuelocate(sfimwave, flat_wave[j])[0]
                ibr = [iloc]
                if iloc != len(sfimwave)-1:
                    ibr.append(iloc+1)
                # get the values in the z-array at indeces ibr, and x=pind[1] and y=pind[0]
                zz = sfim[:, pind[0], pind[1]][[ibr]]
                # now determine the length of the array with only the finite numbers
                zzwherenonan = np.where(np.isfinite(zz))
                kk = np.size(zzwherenonan)
                sfs = 1.0
                if (flat_wave[j] <= max(sfimwave)) and (flat_wave[j] >= min(sfimwave)) and (kk == 2):
                    sfs = np.interp(flat_wave[j], sfimwave[ibr], zz[zzwherenonan])
                # check DQ flags
                kk = np.where(sfimdq[:, pind[0], pind[1]][[ibr]] == 0)
                if np.size(kk) != 2:
                    sfs = 1.0

                # integrate over f-flat fast vector
                # reference file wavelength range is from 0.6 to 5.206 microns, so need to force
                # solution to 1 for wavelengths outside that range
                ffv_wav = ffv.field("WAVELENGTH")
                ffv_dat = ffv.field("DATA")
                fff = 1.0
                if (flat_wave[j]-delw/2.0 >= 0.6) and (flat_wave[j]+delw/2.0 <= 5.206):
                    iw = np.where((ffv_wav >= flat_wave[j]-delw/2.0) & (ffv_wav <= flat_wave[j]+delw/2.0))
                    if np.size(iw) > 1:
                        int_tab = auxfunc.idl_tabulate(ffv_wav[iw], ffv_dat[iw])
                        first_ffv_wav, last_ffv_wav = ffv_wav[iw[0]][0], ffv_wav[iw[0]][-1]
                        fff = int_tab/(last_ffv_wav - first_ffv_wav)
                # interpolate over f-flat cube
                ffs = np.interp(flat_wave[j], ffsallwave, ffsall[:, col-1, row-1])
                flatcor[j] = dff*dfs*sff*sfs*fff*ffs
                if (pind[1]-px0+1 == 9999) and (pind[0]-py0+1 == 9999):
                    print ("pind = ", pind)
                    print ("flat_wave[j] = ", flat_wave[j])
                    print ("dfs, dff = ", dfs, dff)
                    print ("sfs, sff = ", sfs, sff)
                    # make plot
                    font = {#'family' : 'normal',
                            'weight' : 'normal',
                            'size'   : 16}
                    matplotlib.rc('font', **font)
                    fig = plt.figure(1, figsize=(12, 10))
                    plt.subplots_adjust(hspace=.4)
                    ax = plt.subplot(111)
                    xmin = flat_wave[j]-0.01
                    xmax = flat_wave[j]+0.01
                    plt.xlim(xmin, xmax)
                    #plt.ylim(min(dfim[:, pind[0], pind[1]])*0.9, max(dfim[:, pind[0], pind[1]])*1.1)
                    plt.plot(dfwave, dfim[:, pind[0], pind[1]], linewidth=7, marker='D', color='k', label="dflat_im")
                    plt.plot(flat_wave[j], dfs, linewidth=7, marker='D', color='r')
                    plt.plot(dfrqe_wav, dfrqe_rqe, linewidth=7, marker='D', c='k', label="dflat_vec")
                    plt.plot(flat_wave[j], dff, linewidth=7, marker='D', color='r')
                    plt.plot(sfimwave, sfim[:, pind[0], pind[1]], linewidth=7, marker='D', color='k', label="sflat_im")
                    plt.plot(flat_wave[j], sfs, linewidth=7, marker='D', color='r')
                    plt.plot(sfv_wav, sfv_dat, linewidth=7, marker='D', color='k', label="sflat_vec")
                    plt.plot(flat_wave[j], sff, linewidth=7, marker='D', color='r')
                    # add legend
                    box = ax.get_position()
                    ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])
                    ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
                    plt.minorticks_on()
                    plt.tick_params(axis='both', which='both', bottom='on', top='on', right='on', direction='in', labelbottom='on')
                    plt.show()
                    print ("Exiting the program.")
                    plt.close()
                    exit()

                if debug:
                    print ("dfs = ", dfs)
                    print ("sff = ", sff)
                    print ("sfs = ", sfs)
                    print ("ffs = ", ffs)

                # Difference between pipeline and calculated values
                delf[j] = pipeflat[pind[0]-py0, pind[1]-px0] - flatcor[j]
                #print("pind[0], px0, pind[1], flatcor[j], delf[j] : ", pind[0], px0, pind[1], flatcor[j], delf[j])
                #print("np.shape(pipeflat) = ", np.shape(pipeflat))

                # Remove all pixels with values=1 (mainly inter-slit pixels) for statistics
                if pipeflat[pind[0]-py0, pind[1]-px0] == 1:
                    delf[j] = 999.0
                else:
                    flatcor[j] = 1.0   # no correction if no wavelength


        delfg = delf[np.where((delf != 999.0) & (delf >= -0.1))]   # ignore outliers
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
        plt.title(filt+" "+grat+" SLIT_"+slit_id)
        plt.xlabel("flat$_{pipe}$ - flat$_{calc}$")
        plt.ylabel("N")
        xmin = delfg_median - delfg_std*5
        xmax = delfg_median + delfg_std*5
        plt.xlim(xmin, xmax)
        x_median = "median = {:0.3}".format(delfg_median)
        x_stddev = "stddev = {:0.3}".format(delfg_std)
        ax.text(0.7, 0.9, x_median, transform=ax.transAxes, fontsize=fontsize)
        ax.text(0.7, 0.83, x_stddev, transform=ax.transAxes, fontsize=fontsize)
        #print("np.shape(delfg) = ", np.shape(delfg))
        #_, _, _ = ax.hist(delfg, bins=int((xmax-xmin)/40.), histtype='bar', ec='k', facecolor="red", alpha=alpha)
        _, _, _ = ax.hist(delfg, bins=30, histtype='bar', ec='k', facecolor="red", alpha=alpha)

        if save_figs:
            if fig_name is None:
                file_basename = step_input_filename.replace(".fits", "")
                fig_name = file_basename+"_MOSflattest_histogram.jpg"
            plt.savefig(fig_name)
            print ('\n Plot saved: ', fig_name)
        if show_figs:
            plt.show()
        plt.close()


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
    msa_conf_root = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline/build7/test_data/MOS_CV3/complete_pipeline_testset"
    #dflatref_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.2_D_Flat/nirspec_dflat"
    #sfile_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.3_S_Flat/MOS/nirspec_MOS_sflat"
    #fflat_path = "/grp/jwst/wit4/nirspec/CDP3/04_Flat_field/4.1_F_Flat/MOS/nirspec_MOS_fflat"
    dflatref_path = "nirspec_dflat"
    sfile_path = "nirspec_MOS_sflat"
    fflat_path = "nirspec_MOS_fflat"

    # name of the output images
    writefile = True

    # set the names of the resulting plots
    plot_name = "jwtest1010001_01101_00001_MOS_flattest_histogram.jpg"

    # Run the principal function of the script
    median_diff = flattest(step_input_filename, dflatref_path=dflatref_path, sfile_path=sfile_path,
                           fflat_path=fflat_path, msa_conf_root=msa_conf_root, writefile=writefile,
                           show_figs=True, save_figs=False, plot_name=plot_name, threshold_diff=1.0e-14, debug=False)

