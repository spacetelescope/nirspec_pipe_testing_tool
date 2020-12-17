import numpy as np
import os
from scipy import integrate
from scipy import interpolate
from astropy.io import fits
from glob import glob
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from decimal import Decimal

"""
This script contains the auxiliary functions that the wcs FS, MOS, and IFU WCS scripts use.
"""

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "2.2"


# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2018 - Version 2.0: Added routines for plot making and statistics calculations for assign_wcs with using
#                         the datamodel instead of the compute_world_coordinates script.
# Aug 2018 - Version 2.1: Added case to catch simulation rawdataroot names for the ESA files
# Sep 2019 - Version 2.2: Modified function to identify science extensions to work with build 7.3


def find_nearest(arr, value):
    """
    This function gives the content and the index in the array of the number that is closest to
    the value given.
    :param arr = 1-D numpy array
    :param value = float or integer
    :return: The array element closest to value and its index
    """
    idx = (np.abs(arr - value)).argmin()
    return arr[idx], idx


def get_sci_extensions(fits_file_name):
    """
    This functions obtains all the science extensions in the given file
    Args:
        fits_file_name: name of the fits file of interest

    Returns:
        sci_list: list of the numbers of the science extensions
    """
    hdulist = fits.open(fits_file_name)
    hdulist.info()
    sci_dict = {}
    s = 0
    for ext, hdu in enumerate(hdulist):
        if hdu.name == "SCI":
            try:
                sltname = hdu.header["SLTNAME"]
                sci_dict[sltname] = ext
            except KeyError:
                sltname = "Slit_" + repr(s + 1)
                sci_dict[sltname] = ext

    return sci_dict


def do_idl_match(arrA, arrB):
    """
    This function does the same that the IDL match function does. It finds the elements common
    in both arrays and it returns two arrays with the index of those elements in each array.
    (The arrays do not need to be the same length)
    Args:
        arrA: numpy array
        arrB: numpy array

    Returns:
        subA: numpy array of index of arrA from elements also present in arrB
        subB: numpy array of index of arrB from elements also present in arrA
    """
    # Find the index corresponding to the intersection elements and return them as arrays
    subA, subB = [], []
    for i, ai in enumerate(arrA):
        if ai in arrB:
            subA.append(i)
    for i, bi in enumerate(arrB):
        if bi in arrA:
            subB.append(i)
    return np.array(subA), np.array(subB)


def do_idl_rebin(a, *args):
    """
    * This function was copied from Example 2 in http://scipy-cookbook.readthedocs.io/items/Rebinning.html

    This acts identically to IDL's rebin command where all values in the original array are summed
    and divided amongst the entries in the new array. As in IDL, the new shape must be a factor of
    the old one. The ugly 'evList trick' builds and executes a python command of the form
    a.reshape(args[0],factor[0],).sum(1)/factor[0]
    a.reshape(args[0],factor[0],args[1],factor[1],).sum(1).sum(2)/factor[0]/factor[1]
    etc. This general form is extended to cover the number of required dimensions.

    This function rebins ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
     a=rand(6,4); b=rebin(a,3,2)
     a=rand(6); b=rebin(a,2)
    """
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape) / np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],' % (i, i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)' % (i + 1) for i in range(lenShape)] + \
             ['/factor[%d]' % i for i in range(lenShape)]
    # print (''.join(evList))
    return eval(''.join(evList))


def get_modeused_and_rawdatrt_PTT_cfg_file(infile_name):
    # get script directory and config name
    out_dir = os.path.dirname(infile_name)
    PPT_cfg_file = os.path.join(out_dir, "PTT_config.cfg")
    if not os.path.isfile(PPT_cfg_file):
        detector = fits.getval(infile_name, "DETECTOR")
        PPT_cfg_file = PPT_cfg_file.replace(".cfg", "_" + detector + ".cfg")
    with open(PPT_cfg_file, "r") as cfg:
        for i, line in enumerate(cfg.readlines()):
            if "#" not in line:
                if "mode_used" in line:
                    mode_used = line.split()[2]
                if "raw_data_root_file" in line:
                    raw_data_root_file = line.split()[2]
    return mode_used, raw_data_root_file


def get_esafile(esa_files_path, rawdatroot, mode, specifics, nid=None, debug=False):
    """
    This function gets the ESA file corresponding to the input given.
    Args:
        esa_files_path: str, top level of where the regression test data lives
        rawdatroot: str, name of the raw data file (the file ran in create_data)
        mode: string, either 'MOS', 'FS', or 'IFU'
        specifics: list, specific parameters needed for each mode
        nid: string, ESA NID of the raw data file used for the create_data script

    Returns:
        esafile: str, full path of the ESA file corresponding to input given
        log_msgs: list, all screen messages are captured here
    """

    def convert_sltname2esaformat(specifics, esafiles=None, debug=False):
        sltname_list = specifics
        if debug:
            print("specifics: ", specifics)
        for sltname in sltname_list:
            sltname = sltname.replace("S", "")
            if ("A1" in sltname) and ("200" in sltname):
                sltname = "A_" + sltname.split("A")[0] + "_1_"
                if esafiles is not None:
                    esafiles.append(sltname)
            if ("A2" in sltname) and ("200" in sltname):
                sltname = "A_" + sltname.split("A")[0] + "_2_"
                if esafiles is not None:
                    esafiles.append(sltname)
            if ("400" in sltname) or ("1600" in sltname):
                sltname = "A_" + sltname.split("A")[0] + "_"
                if esafiles is not None:
                    esafiles.append(sltname)
            if "B" in sltname:
                sltname = "B_" + sltname.split("B")[0] + "_"
                if esafiles is not None:
                    esafiles.append(sltname)
        if debug:
            print("the slit name used is: ", sltname)
            print(esafiles)
        if esafiles is not None:
            return esafiles
        else:
            return sltname

    def find_esafile_basename(specifics, jlab88_dir, debug=False):
        """
        This function simply avoids code repetition.
        """
        # get the right esa file according to the mode
        esafile_basename = "No match found for esafile"
        if "MOS" in mode:
            quad, row, col = specifics
            # add a 0 if necessary for convention purposes
            if col < 10:
                col = "00" + repr(col)
            elif col < 100:
                col = "0" + repr(col)
            else:
                col = repr(col)
            if row < 10:
                row = "00" + repr(row)
            elif row < 100:
                row = "0" + repr(row)
            else:
                row = repr(row)
            # to match current ESA intermediary files naming convention
            esafile_basename = "Trace_MOS_" + repr(quad) + "_" + row + "_" + col + "_" + jlab88_dir + ".fits"
            if debug:
                print("esafile_basename = ", esafile_basename)
        if "SLIT" in mode:
            esafiles = []
            esafiles = convert_sltname2esaformat(specifics, esafiles=esafiles)
            # to match current ESA intermediary files naming convention
            esafile_basename_list = []
            for sltname in esafiles:
                esafile_basename = "Trace_SLIT_" + sltname + jlab88_dir + ".fits"
                esafile_basename_list.append(esafile_basename)
            esafile_basename = esafile_basename_list[-1]
            if debug:
                print("esafile_basename = ", esafile_basename)
        if "IFU" in mode:
            IFUslice = specifics[0]
            # to match current ESA intermediary files naming convention
            esafile_basename = "Trace_IFU_Slice_" + IFUslice + "_" + jlab88_dir + ".fits"
        return esafile_basename

    # get the root name from rawdatroot keyword (e.g. NRSV96214001001P0000000002105_1_491_SE_2016-01-24T01h59m01.fits)
    esaroot = rawdatroot.split("_")[0]
    if "NRS" in esaroot:
        esaroot = esaroot.replace("NRS", "")  # captures all cases for ground observations
    else:
        esaroot = rawdatroot.replace(".fits", "")  # captures simulation data

    log_msgs = []

    # go into the esa_files_path directory and enter the the mode to get the right esafile
    # get all subdirectories within esa_files_path
    subdir_list = glob(esa_files_path + "/*")
    # check if there are subdirectories or not
    dirs_not_in_list = []
    if debug:
        print("List of all items and/or subdirectories in ", esa_files_path)
    for i, item in enumerate(subdir_list):
        if debug:
            print(i, item)
        if "List" in item:
            subdir_list.pop(i)
            continue
        if ".fits" in item:
            item_not_dir = True
        else:
            item_not_dir = False
        dirs_not_in_list.append(item_not_dir)

    same_nid_files = []
    jlab88_list = []
    esafile = "ESA file not found"
    nidrawfile = ""
    if all(item_not_dir == True for item_not_dir in dirs_not_in_list):
        raw_file_name = rawdatroot.split("_")[0].replace("NRS", "")
        for item in subdir_list:
            if raw_file_name in item:
                if "List" not in item:
                    nidrawfile = fits.getval(item, "GS_JOBID", 0).split("_")[1].replace("000", "")
                    if debug:
                        print("nidrawfile =", nidrawfile)
                    same_nid_files.append(item)
        nid = nidrawfile
    else:
        same_nid_files = []
        for item in subdir_list:
            # check if the file is at the first level using the NID
            if nid is not None:
                if ".fits" in item and "List" not in item:
                    nid2compare = fits.getval(item, "GS_JOBID", 0).split("_")[1].replace("000", "")
                    if debug:
                        print("NID_raw_data_file =", nid, "    nid2compare =", nid2compare)
                    if nid == nid2compare:
                        # collect all files with the same NID
                        same_nid_files.append(item)
            else:
                subdir = item
                if esaroot in subdir:
                    jlab88_list.append(subdir)

    # get the specific ESA file
    if mode == "FS":
        mode = "SLIT"

    # this is specific code for FS_CV3_cutout ESA direcotry, where some ESA files are at the top level dir
    if len(same_nid_files) != 0:
        esafiles = same_nid_files
        if debug:
            print("len(esafiles) = ", len(esafiles))
            print("esafiles = ", esafiles)
        if mode == "SLIT" and nid is not None:
            sltname = convert_sltname2esaformat(specifics, esafiles=None)
            for esaf in esafiles:
                if sltname in esaf:
                    esafile = esaf

    if len(jlab88_list) != 0 and nid is None:
        # If the file is not at the first level go into further directories
        if debug:
            print("jlab88_list=", jlab88_list)
        for jlab88_dir in jlab88_list:
            mode_dir = os.path.join(jlab88_dir, jlab88_dir.split("/")[-1] + "_trace_" + mode)
            if debug:
                print("using these specifics: ", specifics)
            esafile_basename = find_esafile_basename(specifics, jlab88_dir.split("/")[-1])
            if debug:
                print("Using this ESA file: \n", "Directory =", mode_dir, "\n", "File =", esafile_basename)
            if not isinstance(esafile_basename, list):
                esafile = os.path.join(mode_dir, esafile_basename)
                # check if we got the right esafile
                if not os.path.isfile(esafile):
                    esafile = "ESA file not found"
                    break
                try:
                    root_filename = fits.getval(esafile, "FILENAME", 0)
                    if rawdatroot.replace(".fits", "") in root_filename:
                        print(" * File name matches raw file used for create_data.")
                    else:
                        print(" * WARNING: Raw data file name used for create_data does not match esa root file name.")
                        print("            Expected: ", rawdatroot.replace(".fits", ""))
                        print("            Obtained: ", root_filename.split('.')[0])
                        break
                except KeyError:
                    print("* WARNING: PTT is unable to confirm the raw file name used matches the esa file root name.")
            else:
                esafile = []
                for esabase in esafile_basename:
                    esaf = os.path.join(mode_dir, esabase)
                    esafile.append(esaf)
                    # check if we got the right esafile
                    root_filename = fits.getval(esaf, "FILENAME", 0)
                    root_filename_msg = "root_filename = " + root_filename
                    rawdatroot_msg = "rawdatroot = " + rawdatroot
                    print(root_filename_msg)
                    print(rawdatroot_msg)
                    log_msgs.append(root_filename_msg)
                    log_msgs.append(rawdatroot_msg)
                    if rawdatroot.replace(".fits", "") in root_filename:
                        msg = " * File name matches raw file used for create_data."
                        print(msg)
                        log_msgs.append(msg)
                    else:
                        msg = " * WARNING: Raw data file name used for create_data does not match esa root file name."
                        print(msg)
                        log_msgs.append(msg)

    return esafile, log_msgs


def idl_tabulate(x, f, p=5):
    """
    This is a Python proxy to the IDL int_tabulate function taken from:
    https://stackoverflow.com/questions/14345001/idls-int-tabulate-scipy-equivalent

    Args:
        x: array
        f: array
        p: integer, integrator order

    Returns:
        ret: array, integrated values
    """

    def newton_cotes(x, f):
        if x.shape[0] < 2:
            return 0
        rn = (x.shape[0] - 1) * (x - x[0]) / (x[-1] - x[0])
        weights = integrate.newton_cotes(rn)[0]
        """
        # I added this part for the last remaining non 5 points, it will only use the available points
        lw, lf = len(weights), len(f)
        if lw != lf:
            last_weights = []
            for i, fi in enumerate(f):
                last_weights.append(weights[i])
            weights = np.array(last_weights)
        """
        dot_wf = np.dot(weights, f)
        return (x[-1] - x[0]) / (x.shape[0] - 1) * dot_wf

    ret = 0
    try:
        for idx in range(0, x.shape[0], p - 1):
            ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    except:
        IndexError
    return ret


def idl_valuelocate(arr, vals):
    """
    This function is equivalent to value_locate() in IDL.

    Args:
        arr: array where values will be located. ** This array MUST be in increasing order **
        vals: array, values that will be located in arr

    Returns:
        idx: list, indices of where the values would be located if placed in arr

    """

    if isinstance(vals, float):
        vals = [vals]

    idx_list = []
    for v in vals:
        # Find the index and value of the closest element to vals
        if (v > arr[0]) and (v < arr[-1]):
            arr_val, i = find_nearest(arr, v)
            """
            if arr[i] > v:
                i =-1
            """
            idx_list.append(i)
        else:
            idx_list.append(-1)
    # print ("v, i, arr[i] : ", v, i, arr[i])
    return idx_list


def interp_spline(x, y, atx):
    """
    Since spline is being deprecated in scipy, this function does a cubic spline interpolation.
    Args:
        x: array
        y: array
        atx: float

    Returns:
        float corresponding to atx
    """
    t = interpolate.splrep(x, y)
    return interpolate.splev(atx, t)


def construct_reldiff_img(x_size, y_size, forced_nan_idxs, edy_not_restricted, restricted_ig, relative_diff_arr):
    """
    This function constructs an image of the relative differences of the pipeline with respect to ESA.
    Args:
        x_size: integer, dimension of relative difference image in x-direction
        y_size: integer, dimension of relative difference image in y-direction
        forced_nan_idxs: numpy array, indices of all those values that have to be set to NaNs
        edy_not_restricted: numpy array, indices of values in the slit-y array that are outside of the restriction zone
        restricted_ig: numpy array, indices of values in the slit-y array that are inside of the restriction zone
        relative_diff_arr: numpy array, values of the relative differences (e.g. lambdas, msax1, etc.)

    Returns:
        reshaped_rel_diff = numpy array, image ready to be plotted with matplotlib.pyplot.plot.imshow
    """
    # Reconstructing the relative difference images
    ones_arr = np.ones((y_size, x_size))
    # convert to nan everything else but the comparison region
    ones_arr[ones_arr == 1.] = np.nan
    flatten_ones_arr = ones_arr.flatten()
    for fni in forced_nan_idxs:
        flatten_ones_arr[fni] = np.nan
    # add nans where the matched esa slit-y values are outside of the restriction
    for enr in edy_not_restricted:
        flatten_ones_arr[enr] = np.nan
    # add the relative difference values
    for idx, rig_i in enumerate(restricted_ig):
        flatten_ones_arr[rig_i] = relative_diff_arr[idx]
    # reshape to dimensions of matched wavelengths
    reshaped_rel_diff = np.reshape(flatten_ones_arr, (y_size, x_size))
    return reshaped_rel_diff


def compute_percentage(values, threshold):
    """
    This function computes the percentage of pixels above certain threshold.

    Args:
        values: numpy array of quantities to be evaluated
        threshold: float, threshold value to evaluate against

    Returns:
        res: list of float percentages of values above the threshold
    """
    values = values[~np.isnan(values)]
    n_total = values.size

    thresh = [threshold, 3 * threshold, 5 * threshold]
    res = []
    for i in thresh:
        n = np.logical_or(values > i, values < -i).nonzero()[0].size
        result = (n / n_total) * 100
        res.append(result)
    return res


def print_stats(arrX, xname, threshold_diff, absolute=False, return_percentages=False):
    """
    This function prints the statistics for the given arrays.
    Args:
        arrX: numpy array
        xname: string, name of arrX to be printed in legend of statistics
        threshold_diff: float, threshold value to determine if test failed or passed
        absolute: boolean, if True then legend will read 'absolute', otherwise it will read 'relative'

    Returns:
        x_stats: list, all quantities calculated for arrX
        stats_print_strings: list, contains all print statements
    """
    if absolute:
        type_of_calculations = "  Absolute"
        type_of_percentages = "absolute differences"
    else:
        type_of_calculations = "  Relative"
        type_of_percentages = "relative differences"
    # calculate statistics
    stats_print_strings = []
    arrX_mean, arrX_median, arrX_stdev = np.mean(arrX), np.median(arrX), np.std(arrX)
    print("\n", type_of_calculations, xname, " :   mean = %0.3e" % (arrX_mean), "   median = %0.3e" % (arrX_median),
          "   stdev = %0.3e" % (arrX_stdev))
    rel_max = np.max(arrX)
    rel_min = np.min(arrX)
    percentage_results = compute_percentage(arrX, threshold_diff)
    max_str = "    Maximum " + type_of_calculations + xname + " = %0.3e" % rel_max
    min_str = "    Minimum " + type_of_calculations + xname + " = %0.3e" % rel_min
    pix_percentage_greater_than_min = "    Percentage of pixels where median of " + type_of_percentages + \
                                      " is greater than: "
    threshold1 = "                            ->  1xtheshold = " + str(int(round(percentage_results[0], 0))) + "%"
    threshold3 = "                            ->  3xtheshold = " + str(int(round(percentage_results[1], 0))) + "%"
    threshold5 = "                            ->  5xtheshold = " + str(int(round(percentage_results[2], 0))) + "%"
    print(max_str)
    print(min_str)
    print(pix_percentage_greater_than_min)
    print(threshold1)
    print(threshold3)
    print(threshold5)
    stats_print_strings.append(max_str)
    stats_print_strings.append(min_str)
    stats_print_strings.append(pix_percentage_greater_than_min)
    stats_print_strings.append(threshold1)
    stats_print_strings.append(threshold3)
    stats_print_strings.append(threshold5)
    if int(round(percentage_results[1], 0)) > 10:
        msg = '\n *** WARNING: More than 10% of pixels have a median value greater than 3xthreshold!\n'
        print(msg)
        stats_print_strings.append(msg)
    if int(round(percentage_results[2], 0)) > 10:
        msg = '\n *** WARNING: More than 10% of pixels have a median value greater than 5xthreshold!\n'
        print(msg)
        stats_print_strings.append(msg)
    x_stats = [arrX_mean, arrX_median, arrX_stdev]
    if return_percentages:
        percentages = [int(round(percentage_results[0], 0)), int(round(percentage_results[1], 0)),
                       int(round(percentage_results[2], 0))]
        return x_stats, stats_print_strings, percentages
    else:
        return x_stats, stats_print_strings


def get_reldiffarr_and_stats(threshold_diff, edy, esa_arr, arr, arr_name, absolute=False):
    """
    This functions performs all the steps necessary to obtain the relative difference image and the corresponding
    statistics, considering only non-NaN elements.
    Args:
        threshold_diff: float, threshold value used to calculate statistics
        edy: numpy array, the array of the ESA slit-y positions
        esa_arr: numpy array, array of the ESA quantity to be investigated
        arr: numpy array, the array of the quantity to be investigated
        arr_name: string, name of the array to be printed in statistics
        absolute: boolean, default set to False, i.e. quantities are relative

    Returns:
        DATAMODEL_rel_diff: np array, 2D array ready for plotting (image)
        DATAMODEL_rel_diff[notnan]: np array, flat array of the relative differences for the restricted, not-nan values
        notnan_reldiffarr_stats: list, calculated statistics

    """
    # get rid of nans and restrict according to slit-y
    in_slit = np.logical_and(edy < .5, edy > -.5)
    arr[~in_slit] = np.nan  # Set lam values outside the slit to NaN
    nanind = np.isnan(arr)  # get all the nan indexes
    notnan = ~nanind  # get all the not-nan indexes
    arr[nanind] = np.nan  # set all nan indexes to have a value of nan
    # Set the values to NaN in the ESA lambda extension
    esa_arr[nanind] = np.nan

    # Compute the difference in wavelength
    stats_print_strings = []
    if absolute:
        DATAMODEL_diff = arr / esa_arr
    else:
        DATAMODEL_diff = (arr - esa_arr) / esa_arr
    if arr[notnan].size == 0:
        print(arr_name, " has all NaN values. Differences array will also be NaNs and statistics "
                        "calculations will fail.")
        notnan_reldiffarr_stats = np.nan, np.nan, np.nan
    else:
        # calculate and print stats
        if absolute:
            abs_difference = True
        else:
            abs_difference = False
        notnan_reldiffarr_stats, stats_print_strings = print_stats(DATAMODEL_diff[notnan], arr_name, threshold_diff,
                                                                   absolute=abs_difference)
    return DATAMODEL_diff, DATAMODEL_diff[notnan], notnan_reldiffarr_stats, stats_print_strings


def does_median_pass_tes(arr_median, threshold_diff):
    """
    This function determines if the given median is less than or equal to the given threshold
    Args:
        arr_median: float, calculated median
        threshold_diff: float, threshold to compare against

    Returns:
        test_result: string, result of the test, i.e. PASSED or FAILED
    """
    median_diff = False
    if arr_median == np.nan:
        return "FAILED"
    if abs(arr_median) <= threshold_diff:
        median_diff = True
    if median_diff:
        test_result = "PASSED"
    else:
        test_result = "FAILED"
    return test_result


def plt_two_2Dimgandhist(img, hist_data, info_img, info_hist, plt_name=None, plt_origin=None, limits=None, vminmax=None,
                         show_figs=False, save_figs=False):
    """
    This function creates and shows/saves one figure with the 2D plot for the given array and the
    corresponding histogram.
    Args:
        img: numpy array to be ploted in top figure
        hist_data: numpy array to be plotted in bottom figure
        info_img: list, contains title, label_x, and label_y information for top figure
        info_hist: list, contains xlabel, ylabel, bins, and stats information for bottom figure, where stats is
                    also a list containing the mean, median, and standard deviation for hist_data
        plt_name: None or string, if not None the string is the full path and name with which the figure will be saved
        plt_origin: None or list, if not None, code is expecting a list with the following content,
                    plt_origin = lolim_x, uplim_x, lolim_y, uplim_y, this re-sets the origin of the plot
        limits: list, limits of x- and y-axis
        vminmax: list, minimum and maximum values to show in the image (will affect the color scale)
        show_figs: boolean, if True the plot will be showed on-screen
        save_figs: boolean, if True the figure will be saved using the plt_name input

    Returns:
        Nothing.
    """
    # set up generals
    font = {'weight': 'normal',
            'size': 16}
    matplotlib.rc('font', **font)
    plt.figure(1, figsize=(12, 10))
    plt.subplots_adjust(hspace=0.4)
    alpha = 0.2
    fontsize = 15

    # Top figure - 2D plot
    ax = plt.subplot(211)
    if vminmax is None:
        vminmax = [None, None]
    if plt_origin is None:
        im = ax.imshow(img, aspect="auto", origin='lower', vmin=vminmax[0], vmax=vminmax[1])
    else:
        lolim_x, uplim_x, lolim_y, uplim_y = plt_origin
        im = ax.imshow(img, aspect="auto", origin='lower', extent=[lolim_x, uplim_x, lolim_y, uplim_y], vmin=vminmax[0],
                       vmax=vminmax[1])
    plt.tick_params(axis='both', which='both', bottom=True, top=True, right=True, direction='in', labelbottom=True)
    plt.minorticks_on()
    plt.colorbar(im, ax=ax)
    title, label_x, label_y = info_img
    plt.title(title)
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    if limits is not None:
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])

    # Bottom figure - histogram
    ax = plt.subplot(212)
    xlabel, ylabel, bins, stats = info_hist
    x_mean, x_median, x_stddev = stats
    str_x_stddev = "stddev = {:0.3e}".format(x_stddev)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    ax.text(0.73, 0.67, str_x_stddev, transform=ax.transAxes, fontsize=fontsize)
    if bins is None:
        try:
            xmin = x_median - x_stddev * 5
            xmax = x_median + x_stddev * 5
            binwidth = (xmax - xmin) / 40.
            bins = np.arange(xmin, xmax + binwidth, binwidth)
        except:
            ValueError
            bins = 15
    n, bins, patches = ax.hist(hist_data, bins=bins, histtype='bar', ec='k', facecolor="red", alpha=alpha)
    ax.xaxis.set_major_locator(MaxNLocator(8))
    from matplotlib.ticker import FuncFormatter

    def MyFormatter(x, lim):
        if x == 0:
            return 0
        return "%0.1E" % Decimal(x)

    majorFormatter = FuncFormatter(MyFormatter)
    ax.xaxis.set_major_formatter(majorFormatter)
    # add vertical line at mean and median
    plt.axvline(x_mean, label="mean = %0.3e" % (x_mean), color="g")
    plt.axvline(x_median, label="median = %0.3e" % (x_median), linestyle="-.", color="b")
    plt.legend()
    plt.minorticks_on()

    # Show and/or save figures
    if save_figs:
        if plt_name is not None:
            plt.savefig(plt_name)
        else:
            print(" * The figure was given no name. Please set the plt_name variable in the call to the function")
            print("   named plt_two_2Dimgs. FIGURE WILL NOT BE SAVED.")
    if show_figs:
        plt.show()

    plt.close()


def main():
    print("Functions in this module cannot be executed from the command line. To access functions, call from script.")
    pass


if __name__ == '__main__':
    import sys

    sys.exit(main())
