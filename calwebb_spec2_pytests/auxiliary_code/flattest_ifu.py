from __future__ import print_function, division
import numpy as np
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import fits


"""
This script tests the pipeline flat field step output.
"""


def get_esafile(auxiliary_code_path, det, grat, filt, sltname_list, esa_files_path):
    """
    This function gets the ESA file corresponding to the input given.
    Args:
        auxiliary_code_path: str, path where to find the auxiliary code. If not set the code will assume
                            it is in the the auxiliary code directory
        det: str, e.g "NRS1"
        grat: str, grating
        filt: str, filter
        sltname_list: list, slit from data extension
        esa_files_path: str, full path of where to find all ESA intermediary products to make comparisons for the tests

    Returns:
        esafile: str, full path of the ESA file corresponding to input given
    """

    # check if a specific file needs to be used
    if ".fits" in esa_files_path:
        return esa_files_path

    # get the corresponding ESA file to the input file
    # to do this, the script needs the python dictionary of the CV3 data
    sys.path.append(auxiliary_code_path)
    import CV3_testdata_used4build7
    if det == "NRS1":
        file4detector = 0
    elif det == "NRS2":
        file4detector = 1
    for NID, nid_dict_key in CV3_testdata_used4build7.CV3_testdata_dict["FS"]["NID"].items():
        if nid_dict_key["grism"] == grat:
            if nid_dict_key["filter"] == filt:
                CV3filename = nid_dict_key["CV3filename"][file4detector]
                print ("NID of ESA file:", NID)
                print("CV3filename =", CV3filename)
    for sltname in sltname_list:
        # change the format of the string to match the ESA trace
        sltname = sltname.split("S")[1]
        if sltname[-1] == "A1":
            sltname = "A_"+sltname.split("A")[0]+"_1"
        elif sltname[-1] == "A2":
            sltname = "A_"+sltname.split("A")[0]+"_2"
        elif sltname[-1] == "A":
            sltname = "A_"+sltname.split("A")[0]+"_"
        elif sltname[-1] == "B":
            sltname = "B_"+sltname.split("B")[0]+"_"

    # the ESA direcoty names use/follow their name conventions
    ESA_dir_name = CV3filename.split("_")[0].replace("NRS", "")+"_"+NID+"_JLAB88"
    esafile_directory = esa_files_path+ESA_dir_name+"/"+ESA_dir_name+"_trace_SLIT"

    # to match current ESA intermediary files naming convention
    esafile_basename = "Trace_SLIT_"+sltname+ESA_dir_name+".fits"
    print ("Using this ESA file: \n", "Directory =", esafile_directory, "\n", "File =", esafile_basename)
    esafile = os.path.join(esafile_directory, esafile_basename)
    return esafile


def mk_plots(title, show_figs=True, save_figs=False, info_fig1=None, info_fig2=None,
             histogram=False, deltas_plt=False, fig_name=None):
    """
    This function makes all the plots of the script.
    Args:
        title: str, title of the plot
        show_figs: boolean, show figures on screen or not
        save_figs: boolean, save figures or not
        info_fig1: list, arrays, number of bins, and limits for the first figure in the plot
        info_fig2: list, arrays, number of bins, and limits for the second figure in the plot
        histogram: boolean, are the figures in the plot histograms
        deltas_plt: boolean, regular plot
        fig_name: str, name of plot

    Returns:
        It either shows the resulting figure on screen and saves it, or one of the two.
    """
    font = {#'family' : 'normal',
            'weight' : 'normal',
            'size'   : 16}
    matplotlib.rc('font', **font)
    fig = plt.figure(1, figsize=(12, 10))
    plt.subplots_adjust(hspace=.4)
    alpha = 0.2
    fontsize = 15

    # FIGURE 1
    # number in the parenthesis are nrows, ncols, and plot number, numbering in next row starts at left
    ax = plt.subplot(211)
    if histogram:
        xlabel1, ylabel1, xarr1, yarr1, xmin, xmax, bins, x_median, x_stddev = info_fig1
        x_median = "median = {:0.3}".format(x_median)
        x_stddev = "stddev = {:0.3}".format(x_stddev)
        plt.title(title)
        plt.xlabel(xlabel1)
        plt.ylabel(ylabel1)
        plt.xlim(xmin, xmax)
        ax.text(0.7, 0.9, x_median, transform=ax.transAxes, fontsize=fontsize)
        ax.text(0.7, 0.83, x_stddev, transform=ax.transAxes, fontsize=fontsize)
        n, bins, patches = ax.hist(xarr1, bins=bins, histtype='bar', ec='k', facecolor="red", alpha=alpha)
    if deltas_plt:
        title1, xlabel1, ylabel1, xarr1, yarr1, xdelta, x_median, x_stddev = info_fig1
        plt.title(title1)
        plt.xlabel(xlabel1)
        plt.ylabel(ylabel1)
        mean_minus_1half_std = x_median - 1.5*x_stddev
        mean_minus_half_std = x_median - 0.5*x_stddev
        mean_plus_half_std = x_median + 0.5*x_stddev
        mean_plus_1half_std = x_median + 1.5*x_stddev
        '''
        for xd, xi, yi in zip(xdelta, xarr1, yarr1):
            if xd > mean_plus_1half_std:
                plt.plot(xi, yi, linewidth=7, marker='D', color='red')#, label="")
            if xd < mean_minus_1half_std:
                plt.plot(xi, yi, linewidth=7, marker='D', color='fuchsia')#, label="")
            if (xd > mean_minus_1half_std) and (xd < mean_minus_half_std):
                plt.plot(xi, yi, linewidth=7, marker='D', color='blue')#, label="")
            if (xd > mean_minus_half_std) and (xd < mean_plus_half_std):
                plt.plot(xi, yi, linewidth=7, marker='D', color='lime')#, label="")
            if (xd > mean_plus_half_std) and (xd < mean_plus_1half_std):
                plt.plot(xi, yi, linewidth=7, marker='D', color='black')#, label="")
        '''
        idx_red = np.where(xdelta > mean_plus_1half_std)
        idx_fuchsia = np.where(xdelta < mean_minus_1half_std)
        idx_blue = np.where((xdelta > mean_minus_1half_std) & (xdelta < mean_minus_half_std))
        idx_lime = np.where((xdelta > mean_minus_half_std) & (xdelta < mean_plus_half_std))
        idx_black = np.where((xdelta > mean_plus_half_std) & (xdelta < mean_plus_1half_std))
        plt.plot(xarr1[idx_red], yarr1[idx_red], linewidth=7, marker='D', color='red')#, label="")
        plt.plot(xarr1[idx_fuchsia], yarr1[idx_fuchsia], linewidth=7, marker='D', color='fuchsia')#, label="")
        plt.plot(xarr1[idx_blue], yarr1[idx_blue], linewidth=7, marker='D', color='blue')#, label="")
        plt.plot(xarr1[idx_lime], yarr1[idx_lime], linewidth=7, marker='D', color='lime')#, label="")
        plt.plot(xarr1[idx_black], yarr1[idx_black], linewidth=7, marker='D', color='black')#, label="")
        # add legend
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])
        #ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
        #plt.plot(xarr1, yarr1, linewidth=7)
    plt.minorticks_on()
    plt.tick_params(axis='both', which='both', bottom='on', top='on', right='on', direction='in', labelbottom='on')

    # FIGURE 2
    # number in the parenthesis are nrows, ncols, and plot number, numbering in next row starts at left
    ax = plt.subplot(212)
    if histogram:
        xlabel2, ylabel2, xarr2, yarr2, xmin, xmax, bins, y_median, y_stddev = info_fig2
        y_median = "median = {:0.3}".format(y_median)
        y_stddev = "stddev = {:0.3}".format(y_stddev)
        plt.xlabel(xlabel2)
        plt.ylabel(ylabel2)
        plt.xlim(xmin, xmax)
        ax.text(0.7, 0.9, y_median, transform=ax.transAxes, fontsize=fontsize)
        ax.text(0.7, 0.83, y_stddev, transform=ax.transAxes, fontsize=fontsize)
        n, bins, patches = ax.hist(xarr2, bins=bins, histtype='bar', ec='k', facecolor="red", alpha=alpha)
    if deltas_plt:
        title2, xlabel2, ylabel2, xarr2, yarr2, ydelta, y_median, y_stddev = info_fig2
        plt.title(title2)
        plt.xlabel(xlabel2)
        plt.ylabel(ylabel2)
        mean_minus_1half_std = y_median - 1.5*y_stddev
        mean_minus_half_std = y_median - 0.5*y_stddev
        mean_plus_half_std = y_median + 0.5*y_stddev
        mean_plus_1half_std = y_median + 1.5*y_stddev
        '''
        for yd, xi, yi in zip(ydelta, xarr2, yarr2):
            if yd > mean_plus_1half_std:
                plt.plot(xi, yi, linewidth=7, marker='D', color='red')#, label="")
            if yd < mean_minus_1half_std:
                plt.plot(xi, yi, linewidth=7, marker='D', color='fuchsia')#, label="")
            if (yd > mean_minus_1half_std) and (yd < mean_minus_half_std):
                plt.plot(xi, yi, linewidth=7, marker='D', color='blue')#, label="")
            if (yd > mean_minus_half_std) and (yd < mean_plus_half_std):
                plt.plot(xi, yi, linewidth=7, marker='D', color='lime')#, label="")
            if (yd > mean_plus_half_std) and (yd < mean_plus_1half_std):
                plt.plot(xi, yi, linewidth=7, marker='D', color='black')#, label=r"$\mu+0.5*\sigma$ > $\mu+1.5*\sigma$")
        '''
        idx_red = np.where(ydelta > mean_plus_1half_std)
        idx_fuchsia = np.where(ydelta < mean_minus_1half_std)
        idx_blue = np.where((ydelta > mean_minus_1half_std) & (ydelta < mean_minus_half_std))
        idx_lime = np.where((ydelta > mean_minus_half_std) & (ydelta < mean_plus_half_std))
        idx_black = np.where((ydelta > mean_plus_half_std) & (ydelta < mean_plus_1half_std))
        plt.plot(xarr2[idx_red], yarr2[idx_red], linewidth=7, marker='D', color='red')#, label="")
        plt.plot(xarr2[idx_fuchsia], yarr2[idx_fuchsia], linewidth=7, marker='D', color='fuchsia')#, label="")
        plt.plot(xarr2[idx_blue], yarr2[idx_blue], linewidth=7, marker='D', color='blue')#, label="")
        plt.plot(xarr2[idx_lime], yarr2[idx_lime], linewidth=7, marker='D', color='lime')#, label="")
        plt.plot(xarr2[idx_black], yarr2[idx_black], linewidth=7, marker='D', color='black')#, label="")
        # add legend
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])
        #ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
        #plt.plot(xarr2, yarr2, linewidth=7)
    plt.tick_params(axis='both', which='both', bottom='on', top='on', right='on', direction='in', labelbottom='on')
    plt.minorticks_on()
    if save_figs:
        if histogram:
            if fig_name is None:
                fig_name = "FS_wcs_histogram.jpg"
        if deltas_plt:
            if fig_name is None:
                fig_name = "FS_wcs_Deltas.jpg"
        print ('\n Plot saved: ', fig_name)
    if show_figs:
        plt.show()
    plt.close()


def flattest(infile_name, esa_files_path=None, auxiliary_code_path=None,
                show_figs=True, save_figs=False, plot_names=None, threshold_diff=1.0e-14, debug=False):
    """
    This function does the comparison from the pipeline flat field output and
    the intermediary ESA products.

    Args:
        infile_name: str, name of the output fits file from the 2d_extract step (with full path)
        esa_files_path: str, full path of where to find all ESA intermediary products to make comparisons for the tests
        auxiliary_code_path: str, path where to find the auxiliary code. If not set the code will assume
                            it is in the the auxiliary code directory
        show_figs: boolean, whether to show plots or not
        save_figs: boolean, save the plots (the 3 plots can be saved or not independently with the function call)
        plot_names: list of 3 strings, desired names (if names are not given, the plot function will name the plots by
                    default)
        threshold_diff: float, threshold difference between pipeline output and ESA file
        debug: boolean, if true a series of print statements will show on-screen

    Returns:
        - 2 plots, if told to save and/or show them.
        - median_diff: Boolean, True if smaller or equal to 1e-14

    """

    # read in the output of the flat step for relevant info
    print('infile_name=', infile_name)
    det = fits.getval(infile_name, "DETECTOR", 0)
    lamp = fits.getval(infile_name, "LAMP", 0)
    grat = fits.getval(infile_name, "GRATING", 0)
    filter = fits.getval(infile_name, "FILTER", 0)
    exptype = fits.getval(infile_name, "EXP_TYPE", 0)
    print ("flat field fits file  -->     Detector:", det, "   Grating:", grat, "   Filter:", filt, "   Lamp:", lamp)


    # read in the on-the-fly flat image
    basenameinfile_name = os.path.basename(infile_name)
    fileID = basenameinfile_name.split("_")[0]   # obtain the id of the file
    working_dir_path = os.getcwd()
    flatfile = working_dir_path+"/"+fileID+"_"+det+"_uncal_rate_assign_intflat.fits"
    pipeflat =
    
    
    
    
    
    
    if auxiliary_code_path is None:
        auxiliary_code_path = "./"

    #compute_world_coordinates.compute_world_coordinates(infile_name)

    # The world coordinate file was created but it needs to be renamed
    basenameinfile_name = os.path.basename(infile_name)
    fileID = basenameinfile_name.split("_")[0]   # obtain the id of the file
    working_dir_path = os.getcwd()
    wcoordfile = working_dir_path+"/"+fileID+"_world_coordinates.fits"
    #print (wcoordfile)
    # to move file to location of infile
    #cwc_fname = infile_name.replace(".fits", "_world_coordinates.fits")
    # to rename file within the working directory
    cwc_fname = basenameinfile_name.replace(".fits", "_world_coordinates.fits")
    print (cwc_fname)
    #os.system("mv "+wcoordfile+" "+cwc_fname)

    # loop over the slits
    sltname_list = []
    wchdu = fits.open(cwc_fname)
    #n_ext = len(wchdu)
    sci_ext_list = wcsfunc.get_sci_extensions(infile_name)
    print ('sci_ext_list=', sci_ext_list, '\n')

    for i, s_ext in enumerate(sci_ext_list):
        print("-> opening science extension =", s_ext, "  in ", infile_name)
        print("   which corresponds to ext:", i+1, " of file:", cwc_fname)
        hdr = wchdu[i+1].header

        # what is the slit of this exposure
        pslit = hdr["SLIT"]
        print("SLIT = ", pslit)

        # for matched spectrum, get the wavelength and Delta_Y values
        fdata = fits.getdata(infile_name, ext=s_ext)
        pwave = fdata[0,:]
        pdy = fdata[3,:]
        pskyx = fdata[1,:]
        pskyy = fdata[2,:]

        # get the origin of the subwindow
        px0 = fits.getval(infile_name, "SLTSTRT1", ext=s_ext)+fits.getval(infile_name, "SUBSTRT1", ext=0)-1
        py0 = fits.getval(infile_name, "SLTSTRT2", ext=s_ext)+fits.getval(infile_name, "SUBSTRT2", ext=0)-1
        sltname = fits.getval(infile_name, "SLTNAME", ext=s_ext)
        sltname_list.append(sltname)
        n_p = np.shape(fdata)
        npx = n_p[0]
        npy = n_p[1]
        print("npx+1=", npx+1, "px0_list=", px0)
        px = np.arange(1, npx+1)+np.array(px0)
        py = np.arange(1, npy+1)+np.array(py0)
        print  ("Pipeline subwindow corner pixel ID: ", px0, py0)

        # read in ESA data
        esafile = get_esafile(auxiliary_code_path, det, grat, filt, sltname_list, esa_files_path)
        esahdulist = fits.open(esafile)
        #print ("* ESA file contents ")
        #esahdulist.info()
        esahdr1 = esahdulist[1].header
        enext = []
        for ext in esahdulist:
            enext.append(ext)
        if det == "NRS1":
            eflux = fits.getdata(esafile, 1)
            ewave = fits.getdata(esafile, 4)
            edy = fits.getdata(esafile, 5)
        if det == "NRS2":
            eflux = fits.getdata(esafile, 6)
            ewave = fits.getdata(esafile, 9)
            edy = fits.getdata(esafile, 10)
        esahdulist.close()
        n_p = np.shape(eflux)
        nex = n_p[1]
        ney = n_p[0]
        # get the origin of the subwindow
        if det == "NRS1":
            ex0 = esahdr1["CRVAL1"] - esahdr1["CRPIX1"] + 1
            ey0 = esahdr1["CRVAL2"] - esahdr1["CRPIX2"] + 1
        else:
            ex0 = 2048.0 - (esahdr1["CRPIX1"] - esahdr1["CRVAL1"] + 1)
            ey0 = 2048.0 - (esahdr1["CRPIX2"] - esahdr1["CRVAL2"] + 1)
        ex = np.arange(nex) + ex0
        ey = np.arange(ney) + ey0
        print("ESA subwindow corner pixel ID: ", ex0, ey0)
        if debug:
            print("From ESA file: ")
            print("   ex0 =", ex0)
            print("   ex0+nex-1 =", ex0+nex-1)
            print("   ey0 =", ey0)
            print("   ey0+ney-1 =", ey0+ney-1)
            print("   ex=", ex, "   ey=", ey)

        # match up the correct elements in each data set
        subpx, subex = wcsfunc.do_idl_match(px, ex)
        subpy, subey = wcsfunc.do_idl_match(py, ey)
        print("matched elements in the 2D spectra: ", len(subex), len(subey))
        imp, ime = [], []
        for spy in subpy:
            im0 = subpx + npx * spy
            imp.append(im0)
        for sey in subey:
            im0 = subex + nex * sey
            ime.append(im0)
        imp, ime = np.array(imp), np.array(ime)
        imp, ime  = imp.flatten(), ime.flatten()

        # get the difference between the two in units of resels
        # do not include pixels where one or the other solution is 0 or NaN
        flat_pwave, flat_ewave = pwave.flatten(), ewave.flatten()
        ig = []
        for ip, ie, ig_i in zip(imp, ime, range(len(imp))):
            if all( [flat_pwave[ip] != 0 and flat_ewave[ie].size != 0 and np.isfinite(flat_pwave[ip]) and np.isfinite(flat_ewave[ip])] ):
                ig.append(ig_i)

        delwave = []
        for ig_i in ig:
            delw = flat_pwave[imp[ig_i]] - flat_ewave[ime[ig_i]]
            delwave.append(delw)
        delwave = np.array(delwave)
        for dw in delwave:
            if not np.isfinite(dw):
                print("Got a NaN in delwave array!, median and standard deviation will fail.")

        pxr, pyr = np.array([]), np.array([])
        for _ in range(npy):
            pxr = np.concatenate((pxr, px))
        pxr = pxr.astype(int)
        reshaped_py = py.reshape(npy, 1)
        for rpy_i in reshaped_py:
            for _ in range(npx):
                pyr = np.concatenate((pyr, rpy_i))
        pyr = pyr.astype(int)

        pxrg, pyrg, deldy = [], [], []
        flat_pdy, flat_edy = pdy.flatten(), edy.flatten()
        for ig_i in ig:
            pxrg_i = pxr[imp[ig_i]]
            pxrg.append(pxrg_i)
            pyrg_i = pyr[imp[ig_i]]
            pyrg.append(pyrg_i)
            deldy_i = flat_pdy[imp[ig_i]] - flat_edy[ime[ig_i]]
            deldy.append(deldy_i)
        pxrg, pyrg, deldy = np.array(pxrg), np.array(pyrg), np.array(deldy)
        for d in deldy:
            if not np.isfinite(d):
                print("Got a NaN in deldy array!, median and standard deviation will fail.")

        # get the median and standard deviations
        median_diff = False
        if len(delwave) != 0:
            delwave_median, delwave_stddev = np.median(delwave), np.std(delwave)
            deldy_median, deldy_stddev = np.median(deldy), np.std(deldy)
            print("\n  delwave:   median =", delwave_median, "   stdev =", delwave_stddev)
            print("\n  deldy:   median =", deldy_median, "   stdev =", deldy_stddev)

            # This is the key argument for the assert pytest function
            if delwave_median <= threshold_diff:
                median_diff = True

            # PLOTS
            if plot_names is not None:
                hist_name, deltas_name = plot_names

            # HISTOGRAM
            title = filt+"   "+grat+"   SLIT="+sltname
            xmin1 = min(delwave) - (max(delwave)-min(delwave))*0.1
            xmax1 = max(delwave) + (max(delwave)-min(delwave))*0.1
            xlabel1, ylabel1 = r"$\lambda_{pipe}$ - $\lambda_{ESA}$ (10$^{-10}$m)", "N"
            yarr = None
            bins = 15
            info_fig1 = [xlabel1, ylabel1, delwave, yarr, xmin1, xmax1, bins, delwave_median, delwave_stddev]
            xmin2 = min(deldy) - (max(deldy)-min(deldy))*0.1
            xmax2 = max(deldy) + (max(deldy)-min(deldy))*0.1
            xlabel2, ylabel2 = r"$\Delta y_{pipe}$ - $\Delta y_{ESA}$ (relative slit position)", "N"
            info_fig2 = [xlabel2, ylabel2, deldy, yarr, xmin2, xmax2, bins, deldy_median, deldy_stddev]
            mk_plots(title, info_fig1=info_fig1, info_fig2=info_fig2, show_figs=show_figs, save_figs=save_figs,
                     histogram=True, fig_name=hist_name)

            # DELTAS PLOT
            title = ""
            title1, xlabel1, ylabel1 = r"$\Delta \lambda$", "x (pixels)", "y (pixels)"
            info_fig1 = [title1, xlabel1, ylabel1, pxrg, pyrg, delwave, delwave_median, delwave_stddev]
            title2, xlabel2, ylabel2 = r"$\Delta$ Flux", "x (pixels)", "y (pixels)"
            info_fig2 = [title2, xlabel2, ylabel2, pxrg, pyrg, deldy, deldy_median, deldy_stddev]
            mk_plots(title, info_fig1=info_fig1, info_fig2=info_fig2, show_figs=show_figs, save_figs=save_figs,
                     deltas_plt=True, fig_name=deltas_name)
        else:
            print(" * Delta_wavelength array is emtpy. No plots being made. \n")

    return median_diff



if __name__ == '__main__':

    # This is a simple test of the code
    pipeline_path = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline"

    # input parameters that the script expects
    auxiliary_code_path = pipeline_path+"/src/pytests/calwebb_spec2_pytests/auxiliary_code"
    infile_name = "jwtest1003001_01101_00001_NRS1_uncal_rate_assign_wcs_extract_2d_flat_field.fits"
    #esa_files_path=pipeline_path+"/build7/test_data/ESA_intermediary_products/RegressionTestData_CV3_March2017_IFU/"
    # if a specific file needs to be used
    esa_files_path = pipeline_path+"/build7/test_data/ESA_intermediary_products/RegressionTestData_CV3_March2017_IFU/SIMA-QUAL-04-B-6007022859_37668_JLAB88/SIMA-QUAL-04-B-6007022859_37668_JLAB88_trace_IFU/Trace_IFU_Slice_00_SIMA-QUAL-04-B-6007022859_37668_JLAB88.fits"

    # set the names of the resulting plots
    hist_name = "IFU_jwtest1003001_01101_00001_wcs_histogram.jpg"
    deltas_name = "IFU_jwtest1003001_01101_00001_wcs_deltas.jpg"
    plot_names = [hist_name, deltas_name]

    # Run the principal function of the script
    median_diff = compare_wcs(infile_name, esa_files_path=esa_files_path,
                              auxiliary_code_path=auxiliary_code_path,
                              plot_names=plot_names, show_figs=True,
                              save_figs=False, threshold_diff=1.0e-14)




