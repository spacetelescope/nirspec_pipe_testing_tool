from __future__ import print_function, division
from astropy.io import fits
import collections
import os

'''
This script contains functions frequently used in the test suite.
'''

def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default,
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]


def read_hdrfits(fits_file_name, info=False, show_hdr=False):
    '''
    This function reads the header fits file and returns a dictionary of the keywords with
    corresponding values. Keywords will be stored in the order they are read.
    Args:
        fits_file_name: full path with name of the header text file
        info: if True the function will show the contents and shapes of the file
        show_hdr: if True the function will print the header of the file

    Returns:
        hdrl: The header of the fits file
    '''
    #  Read the fits file
    hdulist = fits.open(fits_file_name)
    # print on screen what extensions are in the file
    if info:
        print ('\n FILE INFORMATION: \n')
        hdulist.info()
    # get and print header
    hdrl = hdulist[0].header
    if show_hdr:
        print ('\n FILE HEADER: \n')
        print (repr(hdrl))
    # close the fits file
    hdulist.close()
    return hdrl


def get_filedata(fits_file_name, extension=None):
    """
    This function gets the data from the science extension of the fits file provided.
    Args:
        fits_file_name: name of the fits file used as input for the pipeline
        extension: integer, if no number is given as input, the default is extension 1

    Returns:
        fdata: the data extension of the file
    """
    if extension is None:
        ext = 1
    else:
        ext = extension
    fdata = fits.getdata(fits_file_name, ext)
    return fdata


def get_keywd_val(fits_file_name, keywd, ext=0):
    """
    This function obtains the value corresponding to the given keyword.
    Args:
        fits_file_name: name of the fits file used as input for the pipeline
        keywd: keyword for which to obtain the value
        ext: extension in which the kwyword lives, by default it is set to the
             primary extension.

    Returns:
        keywd_val: the value corresponding to the inputed keyword
    """
    keywd_val = fits.getval(fits_file_name, keywd, ext)
    return keywd_val


def get_sci_extensions(fits_file_name):
    """
    This function obtains all the science extensions in the given file
    Args:
        fits_file_name: name of the fits file of interest

    Returns:
        sci_list: list of the numbers of the science extensions
    """
    hdulist = fits.open(fits_file_name)
    sci_list = []
    for ext, hdu in enumerate(hdulist):
        if hdu.name == "SCI":
            sci_list.append(ext)
    return sci_list


def get_step_inandout_filename(step, initial_input_file, steps_dict):
    """
    This function determines the corresponding input file name for the step (i.e. the pipeline expects a specific
    format and name for each step). This is according to which steps where set to True  in the configuration file.
    Args:
        step: string, pipeline step to be ran
        initial_input_file: the base name of the input file for calwebb_spec2
        steps_dict: dictionary, pipeline steps listed in the input configuration file

    Returns:
        step_input_filename: string, the base name of the input file for the specified step
        step_output_filename: string, the base name of the output file for the specified step
        in_file_suffix : string, the suffix added to the initial file name for the step input file
        out_file_suffix : string, the suffix added to the initial file name for the step output file
    """

    # dictionary of the steps and corresponding strings to be added to the file name after the step has ran
    step_string_dict = collections.OrderedDict()
    step_string_dict["bkg_subtract"] = "_subtract_images"
    step_string_dict["assign_wcs"] = "_assign_wcs"
    step_string_dict["imprint_subtract"] = "_imprint"
    step_string_dict["msa_flagging"] = "_msa_flag"
    step_string_dict["extract_2d"] = "_extract_2d"
    step_string_dict["flat_field"] = "_flat_field"
    step_string_dict["straylight"] = "_stray"
    step_string_dict["fringe"] = "_fringe"
    step_string_dict["pathloss"] = "_pathloss"
    step_string_dict["photom"] = "_photom"
    step_string_dict["resample_spec"] = "_resample"
    step_string_dict["cube_build"] = "_cube"
    step_string_dict["extract_1d"] = "_extract_1d"

    # get all the steps set to run
    steps_set_to_True = []
    for stp, val in steps_dict.items():
        if val:
            steps_set_to_True.append(stp)

    # order the steps set to True according to the ordered dictionary step_string_dict
    ordered_steps_set_to_True = []
    for key in step_string_dict:
        if key in steps_set_to_True:
            ordered_steps_set_to_True.append(key)

    # get the right input and output name according to the steps set to True in the configuration file
    step_input_filename, step_output_filename, in_file_suffix, out_file_suffix = initial_input_file, "", "", ""
    for i, stp in enumerate(ordered_steps_set_to_True):
        if stp == step:
            out_file_suffix = step_string_dict[stp]
            step_output_filename = initial_input_file.replace(".fits", out_file_suffix+".fits")
            break
        else:
            in_file_suffix = step_string_dict[ordered_steps_set_to_True[i]]
            step_input_filename = step_input_filename.replace(".fits", in_file_suffix+".fits")

    return in_file_suffix, out_file_suffix, step_input_filename, step_output_filename


def read_True_steps_suffix_map(txtfile_name_with_path):
    """
    This function reads the text file that contains all the steps in the configuration file, the
    corresponding suffix, and whether they were completed or not.

    Args:
        txtfile_name_with_path: string, full name and path of the text file

    Returns:
        steps_list: list, steps set to True in configuration file
        suffix_list: list, suffix for the output file corresponding to each step had it completed
        completion_list: list, strings of True or False depending on whether the step completed or not

    """
    steps_list, suffix_list, completion_list = [], [], []
    with open(txtfile_name_with_path, "r") as tf:
        for line in tf.readlines():
            if "#" not in line:
                info = line.split()
                steps_list.append(info[0])
                suffix_list.append(info[1])
                completion_list.append(info[2])
    return steps_list, suffix_list, completion_list


def add_completed_steps(True_steps_suffix_map, step, outstep_file_suffix, step_completed):
    """
    This function adds the completed steps along with the corresponding suffix of the output file name into a text file.
    Args:
        True_steps_suffix_map: string, full path of where the text file will be written into
        step: string, pipeline step just ran
        outstep_file_suffix: string, suffix added right before .fits to the input file
        step_completed: boolean, True if the step was completed and False if it was skiped

    Returns:
        nothing
    """
    print ("Map saved at: ", True_steps_suffix_map)
    line2write = "{:<20} {:<20} {:<20}".format(step, outstep_file_suffix, str(step_completed))
    print (line2write)
    with open(True_steps_suffix_map, "a") as tf:
        tf.write(line2write+"\n")


def get_correct_input_step_filename(initial_input_file, steps_list, suffix_list, completion_list):
    """
    This function gets the name of the input step file depending on the info stored in the steps mapping
    text file.

    Args:
        initial_input_file: string, name of the input file to initialize calwebb_spec2
        steps_list: list, steps set to True in configuration file
        suffix_list: list, suffix corresponding to the output file for the step
        completion_list: list, strings of booleans depending if the step ran to completion or not

    Returns:
        step_input_filename: string, base name of the input file for the step
    """
    step_input_filename = initial_input_file
    for i, stp in enumerate(steps_list):
        if str_to_bool(completion_list[i]):
            step_input_filename = step_input_filename.replace(".fits", suffix_list[i]+".fits")
    return step_input_filename


def str_to_bool(s):
    """
    This function converts a string of True or False into a boolean value.
    Args:
        s: string

    Returns:
        boolean value of s
    """
    if s == 'True':
         return True
    elif s == 'False':
         return False

def check_FS_true(output_hdul):
    """
    This function checks if the fits file is a Fixed Slit.
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        result: boolean, if true, the file is assumed to be Fixed Slit
    """
    result = False
    if "EXP_TYPE" in output_hdul:
        if "EXP_TYPE" == "NRS_FIXEDSLIT":
            result = True
    return result


def check_MOS_true(output_hdul):
    """
    This function checks if the fits file is Multi-Object Spectroscopy (MOS).
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        result: boolean, if true, the file is assumed to be MOS data
    """
    result = False
    if "EXP_TYPE" in output_hdul:
        if "EXP_TYPE" == "NRS_MSASPEC":
            result = True
    return result


def check_IFU_true(output_hdul):
    """
    This function checks if the fits file is IFU data.
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        result: boolean, if true, the file is assumed to be IFU data
    """
    result = False
    if "EXP_TYPE" in output_hdul:
        if "EXP_TYPE" == "NRS_IFU":
            result = True
    return result

def find_which_slit(output_hdul):
    """
    This function determines which Fixed Slit was used
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        s: string, slit used or is None if not in the list
    """
    # the order of this list corresponds to the
    slits = ["S200A1", "S200A2", "S200B1", "S400A1", "S1600A1"]
    if "FXD_SLIT" in output_hdul:
        for i, s in enumerate(slits):
            if "FXD_SLIT" == s:
                return i+1, s


def find_DETECTOR(output_hdul):
    """
    This function determines which detector was used
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        det: string, either NRS1 or NRS2
    """
    if "DETECTOR" in output_hdul:
        det = output_hdul["DETECTOR"]
        return det
