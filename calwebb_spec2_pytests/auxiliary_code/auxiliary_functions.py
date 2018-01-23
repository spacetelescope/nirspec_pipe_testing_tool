import numpy as np
import os
from scipy import integrate
from scipy import interpolate
from astropy.io import fits
from glob import glob


"""
This script contains the auxiliary functions that the wcs FS, MOS, and IFU WCS scripts use.
"""


def find_nearest(arr, value):
    '''
    This function gives the content and the index in the array of the number that is closest to
    the value given.
    :param arr = 1-D numpy array
    :param value = float or integer
    :return: The array element closest to value and its index
    '''
    idx=(np.abs(arr-value)).argmin()
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
    sci_list = []
    for ext, hdu in enumerate(hdulist):
        if hdu.name == "SCI":
            sci_list.append(ext)
    return sci_list


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
    '''
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
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    #print (''.join(evList))
    return eval(''.join(evList))


def get_esafile(esa_files_path, rawdatroot, mode, specifics):
    """
    This function gets the ESA file corresponding to the input given.
    Args:
        esa_files_path: str, top level of where the regression test data lives
        rawdatroot: str, name of the raw data file (the file ran in create_data)
        mode: string, either 'MOS', 'FS', or 'IFU'
        specifics: list, specific parameters needed for each mode

    Returns:
        esafile: str, full path of the ESA file corresponding to input given
    """

    def find_esafile_basename(specifics, jlab88_dir):
        """
        This function simply avoids code repetition.
        """
        # get the right esa file according to the mode
        esafile_basename = "No match found for esafile"
        if "MOS" in mode:
            quad,row, col = specifics
            # add a 0 if necessary for convention purposes
            if col < 99:
                col = "0"+str(col[0])
            else:
                col = str(col[0])
            if row < 99:
                row = "0"+str(row[0])
            else:
                row = str(row[0])
            # to match current ESA intermediary files naming convention
            esafile_basename = "Trace_MOS_"+str(quad[0])+"_"+row+"_"+col+"_"+jlab88_dir+".fits"
        if "SLIT" in mode:
            sltname_list = specifics
            esafiles = []
            for sltname in sltname_list:
                # change the format of the string to match the ESA trace
                sltname = sltname.replace("S", "")
                print("sltname after removing S: ", sltname)
                if ("A" in sltname) and ("200" in sltname):
                    sltname = "A_"+sltname.split("A")[0]+"_1_"
                    esafiles.append(sltname)
                if ("A2" in sltname) and ("200" in sltname):
                    sltname = "A_"+sltname.split("A")[0]+"_2_"
                    esafiles.append(sltname)
                if ("400" in sltname) or ("1600" in sltname):
                    sltname = "A_"+sltname.split("A")[0]+"_"
                    esafiles.append(sltname)
                if "B" in sltname:
                    sltname = "B_"+sltname.split("B")[0]+"_"
                    esafiles.append(sltname)
            # to match current ESA intermediary files naming convention
            esafile_basename_list = []
            for sltname in esafiles:
                esafile_basename = "Trace_SLIT_"+sltname+jlab88_dir+".fits"
                esafile_basename_list.append(esafile_basename)
            esafile_basename = esafile_basename_list
        if "IFU" in mode:
            IFUslice = specifics[0]
            # to match current ESA intermediary files naming convention
            esafile_basename = "Trace_IFU_Slice_"+IFUslice+"_"+jlab88_dir+".fits"
        return esafile_basename


    # get the root name from rawdatroot keyword (e.g. NRSV96214001001P0000000002105_1_491_SE_2016-01-24T01h59m01.fits)
    esaroot = rawdatroot.split("_")[0].replace("NRS", "")

    # go into the esa_files_path directory and enter the the mode to get the right esafile
    # get all subdirectories within esa_files_path
    subdir_list = glob(esa_files_path+"/*")
    jlab88_list = []
    esafile = "ESA file not found"
    for subdir in subdir_list:
        #subdir = subdir.split("/")[-1]
        if esaroot in subdir:
            jlab88_list.append(subdir)
    #print("jlab88_list=", jlab88_list)
    for jlab88_dir in jlab88_list:
        if mode == "FS":
            mode = "SLIT"
        mode_dir = os.path.join(jlab88_dir, jlab88_dir.split("/")[-1]+"_trace_"+mode)
        esafile_basename = find_esafile_basename(specifics, jlab88_dir.split("/")[-1])
        print ("Using this ESA file: \n", "Directory =", mode_dir, "\n", "File =", esafile_basename)
        if not isinstance(esafile_basename, list):
            esafile = os.path.join(mode_dir, esafile_basename)
            # check if we got the right esafile
            root_filename = fits.getval(esafile, "FILENAME", 0)
            print("root_filename = ", root_filename)
            print("rawdatroot = ", rawdatroot)
            if rawdatroot.replace(".fits", "") in root_filename:
                print (" * File name matches raw file used for create_data.")
                break
            else:
                print (" * WARNING: Raw data file name used for create_data does not match esa root file name.")
        else:
            esafile = []
            for esabase in esafile_basename:
                esaf = os.path.join(mode_dir, esabase)
                esafile.append(esaf)
                # check if we got the right esafile
                root_filename = fits.getval(esaf, "FILENAME", 0)
                print("root_filename = ", root_filename)
                print("rawdatroot = ", rawdatroot)
                if rawdatroot.replace(".fits", "") in root_filename:
                    print (" * File name matches raw file used for create_data.")
                else:
                    print (" * WARNING: Raw data file name used for create_data does not match esa root file name.")

    return esafile



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
    def newton_cotes(x, f) :
        if x.shape[0] < 2 :
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
    for idx in range(0, x.shape[0], p - 1) :
        ret += newton_cotes(x[idx:idx + p], f[idx:idx + p])
    return ret


def idl_valuelocate(arr, vals):
    """
    This function is equivalent to value_locate() in IDL.

    Args:
        arr: array where values will be located. ** This array MUST be in increasing order **
        vals: array, values that will be located in arr

    Returns:
        idx: list, indeces of where the values would be located if placed in arr

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
    #print ("v, i, arr[i] : ", v, i, arr[i])
    return idx_list


def interp_spline(x,y, atx):
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
