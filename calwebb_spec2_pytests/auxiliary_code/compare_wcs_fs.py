import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import fits

from jwst.assign_wcs.tools.nirspec import compute_world_coordinates
from . import auxiliary_functions as auxfunc


"""
This script compares pipeline WCS info with ESA results for FIXED SLIT.

slit = 'S200A1', 'S200A2', 'S400A1', 'S1600A1', 'S200B1'
det = 'NRS1' or 'NRS2'
subarray_origin = [SUBSTRT1, SUBSTRT2] from original image; needed since SUBSTRT in
                    world_coordinates file is not in full frame reference
"""


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
        info_fig3: list, arrays, number of bins, and limits for the second figure in the plot
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
    if histogram:
        ax = plt.subplot(211)
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
        ax = plt.subplot(211)
        title1, xlabel1, ylabel1, xarr1, yarr1, xdelta, x_median, x_stddev = info_fig1
        xarr1, yarr1 = xarr1[0], yarr1[0]
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
    if histogram:
        ax.xaxis.set_major_locator(MaxNLocator(6))
        from matplotlib.ticker import FuncFormatter
        def MyFormatter(x, lim):
            if x == 0:
                return 0
            return '{0}E{1}'.format(round(x/1e-10, 2), -10)
        majorFormatter = FuncFormatter(MyFormatter)
        ax.xaxis.set_major_formatter(majorFormatter)
    plt.tick_params(axis='both', which='both', bottom='on', top='on', right='on', direction='in', labelbottom='on')

    # FIGURE 2
    # number in the parenthesis are nrows, ncols, and plot number, numbering in next row starts at left
    if histogram:
        ax = plt.subplot(212)
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
        ax = plt.subplot(212)
        title2, xlabel2, ylabel2, xarr2, yarr2, ydelta, y_median, y_stddev = info_fig2
        xarr2, yarr2 = xarr2[0], yarr2[0]
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
        """
        ax = plt.subplot(313)
        title3, xlabel3, ylabel3, xarr3, yarr3, ydelta, y_median, y_stddev = info_fig3
        xarr3, yarr3 = xarr3[0], yarr3[0]
        plt.title(title3)
        plt.xlabel(xlabel3)
        plt.ylabel(ylabel3)
        mean_minus_1half_std = y_median - 1.5*y_stddev
        mean_minus_half_std = y_median - 0.5*y_stddev
        mean_plus_half_std = y_median + 0.5*y_stddev
        mean_plus_1half_std = y_median + 1.5*y_stddev
        idx_red = np.where(ydelta > mean_plus_1half_std)
        idx_fuchsia = np.where(ydelta < mean_minus_1half_std)
        idx_blue = np.where((ydelta > mean_minus_1half_std) & (ydelta < mean_minus_half_std))
        idx_lime = np.where((ydelta > mean_minus_half_std) & (ydelta < mean_plus_half_std))
        idx_black = np.where((ydelta > mean_plus_half_std) & (ydelta < mean_plus_1half_std))
        plt.plot(xarr3[idx_red], yarr3[idx_red], linewidth=7, marker='D', color='red')#, label="ydelta > median_plus_1half_std")
        plt.plot(xarr3[idx_fuchsia], yarr3[idx_fuchsia], linewidth=7, marker='D', color='fuchsia')#, label="ydelta < median_minus_1half_std")
        plt.plot(xarr3[idx_blue], yarr3[idx_blue], linewidth=7, marker='D', color='blue')#, label="median_plus_1half_std > ydelta > median_minus_1half_std")
        plt.plot(xarr3[idx_lime], yarr3[idx_lime], linewidth=7, marker='D', color='lime')#, label="median_plus_half_std > ydelta > median_minus_half_std")
        plt.plot(xarr3[idx_black], yarr3[idx_black], linewidth=7, marker='D', color='black')#, label="median_plus_1half_std > ydelta > median_minus_half_std")
        # add legend
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])
        #ax.legend(loc='upper right', bbox_to_anchor=(1, 1))
        #plt.plot(xarr3, yarr3, linewidth=7)
        """
    plt.tick_params(axis='both', which='both', bottom='on', top='on', right='on', direction='in', labelbottom='on')
    plt.minorticks_on()
    if save_figs:
        type_fig = "pdf"
        if histogram:
            if fig_name is None:
                fig_name = ".".join(("FS_wcs_histogram", type_fig))
        if deltas_plt:
            if fig_name is None:
                fig_name = ".".join(("FS_wcs_Deltas", type_fig))
        plt.savefig(fig_name)
        print ('\n Plot saved: ', fig_name)
    if show_figs:
        plt.show()
    plt.close()


def compare_wcs(infile_name, esa_files_path=None, auxiliary_code_path=None,
                show_figs=True, save_figs=False, plot_names=None, threshold_diff=1.0e-14, debug=False):
    """
    This function does the WCS comparison from the world coordinates calculated using the
    compute_world_coordinates.py script with the ESA files. The function calls that script.

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

    # get grating and filter info from the rate file header
    det = fits.getval(infile_name, "DETECTOR", 0)
    print('infile_name=', infile_name)
    lamp = fits.getval(infile_name, "LAMP", 0)
    grat = fits.getval(infile_name, "GRATING", 0)
    filt = fits.getval(infile_name, "FILTER", 0)
    print ("extract_2d  -->     Detector:", det, "   Grating:", grat, "   Filter:", filt, "   Lamp:", lamp)


    # Run compute_world_coordinates.py in order to produce the necessary file
    # !!! note that the code expects to be in the build environment !!!
    if auxiliary_code_path is None:
        auxiliary_code_path = "./"

    compute_world_coordinates.compute_world_coordinates(infile_name)

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
    cwc_fname = infile_name.replace(basenameinfile_name, cwc_fname)
    os.system("mv "+wcoordfile+" "+cwc_fname)

    # loop over the slits
    sltname_list = []
    wchdu = fits.open(cwc_fname)
    #n_ext = len(wchdu)
    sci_ext_list = auxfunc.get_sci_extensions(infile_name)
    print ('sci_ext_list=', sci_ext_list, '\n')

    for i, s_ext in enumerate(sci_ext_list):
        print("-> opening extension =", i+1, "  in ", cwc_fname)
        print("   which corresponds to science ext:", s_ext, " of file:", infile_name)
        hdr = wchdu[i+1].header

        # what is the slit of this exposure
        pslit = hdr["SLIT"]
        print("SLIT = ", pslit)

        # for matched spectrum, get the wavelength and Delta_Y values
        fdata = fits.getdata(cwc_fname, ext=i+1)
        pwave = fdata[0,:,:]
        pdy = fdata[3,:,:]
        pskyx = fdata[1,:,:]
        pskyy = fdata[2,:,:]
        #print('pskyx=', pskyx)
        #print('pskyy=', pskyy)
        #print('np.shape(fdata) = ', np.shape(fdata))

        # get the origin of the subwindow
        px0 = fits.getval(infile_name, "SLTSTRT1", ext=s_ext)+fits.getval(infile_name, "SUBSTRT1", ext=0)-1
        py0 = fits.getval(infile_name, "SLTSTRT2", ext=s_ext)+fits.getval(infile_name, "SUBSTRT2", ext=0)-1
        sltname = fits.getval(infile_name, "SLTNAME", ext=s_ext)
        sltname_list.append(sltname)
        n_p = np.shape(fdata)
        npx = n_p[2]
        npy = n_p[1]
        #print("npx=", npx, "npy=", npy)
        px = np.arange(npx)+np.array(px0)
        py = np.arange(npy)+np.array(py0)
        print  ("Pipeline subwindow corner pixel ID: ", px0, py0)

        # read in the ESA file using raw data root file name
        #raw_data_root_file = "NRSV00300060001P000000000210T_1_491_SE_2016-01-06T06h27m34.fits"  # for testing script
        _, raw_data_root_file = auxfunc.get_modeused_and_rawdatrt_PTT_cfg_file()
        specifics = sltname_list
        # check if ESA data is not in the regular directories
        NIDs = ["30055", "30205", "30133"]
        special_cutout_files = ["NRSSMOS-MOD-G1H-02-5344031756_1_491_SE_2015-12-10T03h25m56.fits",
                                "NRSSMOS-MOD-G2M-01-5344191938_1_491_SE_2015-12-10T19h29m26.fits",
                                "NRSSMOS-MOD-G3H-02-5344120942_1_492_SE_2015-12-10T12h18m25.fits"]
        if raw_data_root_file in special_cutout_files:
            nid = NIDs[special_cutout_files.index(raw_data_root_file)]
            print("Using NID = ", nid)
        else:
            nid = None
        esafile = auxfunc.get_esafile(esa_files_path, raw_data_root_file, "FS", specifics, nid=nid)

        # skip the test if the esafile was not found
        if esafile == "ESA file not found":
            print(" * compare_wcs_fs.py is exiting because the corresponding ESA file was not found.")
            print("   -> The WCS test is now set to skip and no plots will be generated. ")
            median_diff = "skip"
            return median_diff

        if not isinstance(esafile, list):
            esafile_list = [esafile]
        else:
            esafile_list = esafile

        slit = sltname.replace("S", "")
        if "A" in slit:
            slit = "_"+slit.split("A")[0]+"_"
        if "B" in slit:
            slit = "_"+slit.split("B")[0]+"_"

        # choose corresponding esa file
        for esafile in esafile_list:
            #print ("looking for slit: ", slit)
            #print (" in esafile: ", esafile)
            if slit in esafile:
                print ("Using this ESA file: \n", esafile)
                esahdulist = fits.open(esafile)
                #print ("* ESA file contents ")
                #esahdulist.info()
                break

        esahdr1 = esahdulist[1].header
        enext = []
        for ext in esahdulist:
            enext.append(ext)
        if det == "NRS1":
            eflux = fits.getdata(esafile, 1)
            ewave = fits.getdata(esafile, 4)
            edy = fits.getdata(esafile, 5)
        if det == "NRS2":
            try:
                eflux = fits.getdata(esafile, 6)
                ewave = fits.getdata(esafile, 9)
                edy = fits.getdata(esafile, 10)
            except:
                IndexError
                print(" * compare_wcs_fs.py is exiting because there are no extensions that match detector NRS2 in the ESA file.")
                print("   -> The WCS test is now set to skip and no plots will be generated. ")
                median_diff = "skip"
                return median_diff

        esahdulist.close()
        n_p = np.shape(eflux)
        #print("eflux=", eflux.flatten())
        nex = n_p[1]
        ney = n_p[0]
        # get the origin of the subwindow
        if det == "NRS1":
            ex0 = esahdr1["CRVAL1"] - esahdr1["CRPIX1"] + 1
            ey0 = esahdr1["CRVAL2"] - esahdr1["CRPIX2"] + 1
        else:
            ex0 = 2048.0 - (esahdr1["CRPIX1"] - esahdr1["CRVAL1"])
            ey0 = 2048.0 - (esahdr1["CRPIX2"] - esahdr1["CRVAL2"])
        ex = np.arange(nex) + ex0
        ey = np.arange(ney) + ey0
        print("ESA subwindow corner pixel ID: ", ex0, ey0)
        if debug:
            print("From ESA file: ")
            #print("   ex0 =", ex0)
            print("   ex0+nex-1 =", ex0+nex-1)
            #print("   ey0 =", ey0)
            print("   ey0+ney-1 =", ey0+ney-1)
            print("   ex=", len(ex), "   ey=", len(ey))
            print("   px=", len(px), "   py=", len(py))

        # match up the correct elements in each data set
        subpx, subex = auxfunc.do_idl_match(px, ex)
        subpy, subey = auxfunc.do_idl_match(py, ey)
        print("matched elements in the 2D spectra: ", len(subpx), len(subey))
        for sx in subpx:
            if px[sx] not in ex:
                print("aha! found it!  ", ex[sx])
        imp, ime = [], []
        for spy in subpy:
            im0 = subpx + npx * spy
            imp.append(im0)
        for sey in subey:
            im0 = subex + nex * sey
            ime.append(im0)
        imp, ime = np.array(imp), np.array(ime)
        imp, ime  = imp.flatten(), ime.flatten()

        #print('ime = ', np.shape(ime), '   imp =', np.shape(imp))
        #print('ime[0:9] = ', ime[0:10])

        # get the difference between the two in units of resels
        # do not include pixels where one or the other solution is 0 or NaN
        flat_pwave, flat_ewave = pwave.flatten(), ewave.flatten()
        # get the elements with index imp and ime
        #print('flat_pwave=', len(flat_pwave), '  flat_ewave=', len(flat_ewave))

        try:
            ig = np.where((flat_pwave[imp]!=0.0) & (flat_ewave[ime]!=0.0) & (np.isfinite(flat_pwave[imp])) & (np.isfinite(flat_ewave[imp])))
        except:
            IndexError
        #print('ig = ', np.shape(ig))

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
        try:
            for ig_i in ig:
                pxrg_i = pxr[imp[ig_i]]
                pxrg.append(pxrg_i)
                pyrg_i = pyr[imp[ig_i]]
                pyrg.append(pyrg_i)
                deldy_i = flat_pdy[imp[ig_i]] - flat_edy[ime[ig_i]]
                deldy.append(deldy_i)
            pxrg, pyrg, deldy = np.array(pxrg), np.array(pyrg), np.array(deldy)

            igy = np.where((flat_pdy[imp]!=0.0) & (flat_edy[ime]!=0.0) & (np.isfinite(flat_pdy[imp])) & (np.isfinite(flat_edy[imp])))
            #print('igy = ', np.shape(igy))

            delwave = (flat_pwave[imp[ig]]*1.0e-6 -flat_ewave[ime[ig]])
            deldy = flat_pdy[imp[igy]] - flat_edy[ime[igy]]
            #delflx = flat_pdy[imp[ig]] - flat_edy[ime[ig]]
            #print (flat_pdy[imp[ig]], flat_edy[ime[ig]])

            #print('np.shape(delwave) = ', np.shape(delwave))
            #print('np.shape(deldy) = ', np.shape(deldy))
        except:
            IndexError

        # get the median and standard deviations
        median_diff = False
        if (len(delwave) != 0) and (len(deldy) != 0):
            delwave_median, delwave_stddev = np.median(delwave), np.std(delwave)
            deldy_median, deldy_stddev = np.median(deldy), np.std(deldy)
            print("\n  delwave:   median =", delwave_median, "   stdev =", delwave_stddev)
            print("\n    deldy:   median =", deldy_median, "   stdev =", deldy_stddev)

            # This is the key argument for the assert pytest function
            if abs(delwave_median) <= threshold_diff:
                median_diff = True
            if median_diff:
                test_result = "PASSED"
            else:
                test_result = "FAILED"
            print (" *** Result of the test: ",test_result)

            # PLOTS
            if np.isfinite(delwave_median):
                print ("Making WCS plots...")
                #if plot_names is not None:
                #    hist_name, deltas_name = plot_names
                #else:
                #    #hist_name, deltas_name = None, None
                hist_name = infile_name.replace(".fits", "_"+sltname+"_wcs_histogram.pdf")
                deltas_name = infile_name.replace(".fits", "_"+sltname+"_wcs_Deltas.pdf")


                # HISTOGRAM
                if show_figs or save_figs:
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
                    title2, xlabel2, ylabel2 = "Relative slit position", "x (pixels)", "y (pixels)"
                    info_fig2 = [title2, xlabel2, ylabel2, pxrg, pyrg, deldy, deldy_median, deldy_stddev]
                    #title3, xlabel3, ylabel3 = r"$\Delta$ Flux", "x (pixels)", "y (pixels)"
                    #info_fig3 = [title3, xlabel3, ylabel3, pxrg, pyrg, delflx, deldy_median, deldy_stddev]
                    mk_plots(title, info_fig1=info_fig1, info_fig2=info_fig2,
                             show_figs=show_figs, save_figs=save_figs, deltas_plt=True, fig_name=deltas_name)

                    print("Done.")

                else:
                    print ("compare_wcs.py ran but NO plots were made because show_figs and save_figs were both set to False. \n")
            else:
                print("Not making plots because median is NaN.")


    return median_diff



if __name__ == '__main__':

    # This is a simple test of the code
    pipeline_path = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline"

    # input parameters that the script expects
    auxiliary_code_path = pipeline_path+"/src/pytests/calwebb_spec2_pytests/auxiliary_code"
    data_dir = "/Users/pena/Documents/PyCharmProjects/nirspec/pipeline/build7.1/part2/FS_FULL_FRAME/G140H_opaque/491_results/"
    infile_name = data_dir+"gain_scale_assign_wcs_extract_2d.fits"
    #esa_files_path=pipeline_path+"/build7/test_data/ESA_intermediary_products/RegressionTestData_CV3_March2017_FixedSlit/"
    esa_files_path = "/grp/jwst/wit4/nirspec_vault/prelaunch_data/testing_sets/b7.1_pipeline_testing/test_data_suite/FS_CV3_cutouts/ESA_Int_products"
    # if a specific file needs to be used
    #esa_files_path = pipeline_path+"/build7/test_data/ESA_intermediary_products/RegressionTestData_CV3_March2017_FixedSlit/V84600003001P0000000002104_39528_JLAB88/V84600003001P0000000002104_39528_JLAB88_trace_SLIT/Trace_SLIT_A_200_1_V84600003001P0000000002104_39528_JLAB88.fits"

    # set the names of the resulting plots
    hist_name = "FS_wcs_histogram.pdf"
    deltas_name = "FS_wcs_deltas.pdf"
    plot_names = None#[hist_name, deltas_name]

    # Run the principal function of the script
    median_diff = compare_wcs(infile_name, esa_files_path=esa_files_path, auxiliary_code_path=None, plot_names=None,
                              show_figs=False, save_figs=True, threshold_diff=9.9e-14, debug=False)




