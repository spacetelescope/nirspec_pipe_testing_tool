import collections
import os
import configparser
import glob
import subprocess
from astropy.io import fits

'''
This script contains functions frequently used in the test suite.
'''


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed



def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default,
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]


def read_hdrfits(fits_file_name, info=False, show_hdr=False, ext=0):
    '''
    This function reads the header fits file and returns a dictionary of the keywords with
    corresponding values. Keywords will be stored in the order they are read.
    Args:
        fits_file_name: full path with name of the header text file
        info: if True the function will show the contents and shapes of the file
        show_hdr: if True the function will print the header of the file
        ext: integer, number of extension to be read

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
    hdrl = hdulist[ext].header
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


def get_step_inandout_filename(step, initial_input_file, steps_dict, working_directory=None):
    """
    This function determines the corresponding input file name for the step (i.e. the pipeline expects a specific
    format and name for each step). This is according to which steps where set to True  in the configuration file.
    Args:
        step: string, pipeline step to be ran
        initial_input_file: the base name of the input file for calwebb_spec2
        steps_dict: dictionary, pipeline steps listed in the input configuration file
        working_directory: string, path where the output files will be saved at

    Returns:
        step_input_filename: string, the base name of the input file for the specified step
        step_output_filename: string, the base name of the output file for the specified step
        in_file_suffix : string, the suffix added to the initial file name for the step input file
        out_file_suffix : string, the suffix added to the initial file name for the step output file
    """

    # dictionary of the steps and corresponding strings to be added to the file name after the step has ran
    step_string_dict = collections.OrderedDict()
    step_string_dict["assign_wcs"] = "_assign_wcs"
    step_string_dict["bkg_subtract"] = "_subtract_images"
    step_string_dict["imprint_subtract"] = "_imprint"
    #step_string_dict["msa_flagging"] = "_msa_flag"
    step_string_dict["extract_2d"] = "_extract_2d"
    step_string_dict["flat_field"] = "_flat_field"
    step_string_dict["straylight"] = "_stray"
    step_string_dict["fringe"] = "_fringe"
    step_string_dict["pathloss"] = "_pathloss"
    step_string_dict["barshadow"] = "_barshadow"
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
            if working_directory is not None  and  stp=="assign_wcs":
                initial_input_file_basename = os.path.basename(initial_input_file)
                step_output_filename = initial_input_file_basename.replace(".fits", out_file_suffix+".fits")
                step_output_filename = os.path.join(working_directory, step_output_filename)
            else:
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


def read_completion_to_full_run_map(full_run_map, step):
    """
    This function adds the completion boolean to the full run map file.
    Args:
        full_run_map: string, path and name of the text file to be updated
        step: string, name of the step to be updated

    Returns:
        step_product = string, name of the intermediary fits product of the step
    """
    step_product = ""
    with open(full_run_map, "r") as tf:
        for line in tf.readlines():
            if "#" not in line:
                if step in line:
                    step_product = line.split()[1]
    return step_product


def add_completed_steps(True_steps_suffix_map, step, outstep_file_suffix, step_completed, end_time):
    """
    This function adds the completed steps along with the corresponding suffix of the output file name into a text file.
    Args:
        True_steps_suffix_map: string, full path of where the text file will be written into
        step: string, pipeline step just ran
        outstep_file_suffix: string, suffix added right before .fits to the input file
        step_completed: boolean, True if the step was completed and False if it was skiped
        end_time: string, time it took for the step to run (in seconds)

    Returns:
        nothing
    """
    print ("Map saved at: ", True_steps_suffix_map)
    if (float(end_time)) > 60.0:
        end_time_min = float(end_time)/60.  # this is in minutes
        end_time = end_time+"  ="+repr(round(end_time_min, 2))+"min"
    line2write = "{:<20} {:<20} {:<20} {:<20}".format(step, outstep_file_suffix, str(step_completed), end_time)
    print (line2write)
    with open(True_steps_suffix_map, "a") as tf:
        tf.write(line2write+"\n")


def start_end_PTT_time(txt_name, start_time=None, end_time=None):
    """
    This function calculates and prints the starting/ending PTT running time in the True_steps_suffix_map.txt or
    full_run_map.txt file.
    Args:
        txt_name: string, path and name of the text file
        start_time: float, starting time
        end_time: float, ending time

    Returns:
        Nothing.
    """
    if start_time is not None:
        # start the timer to compute the step running time of PTT
        print("PTT starting time: ", repr(start_time), "\n")
        line2write = "{:<20} {:<20}".format('# Starting PTT running time: ', repr(start_time))
        with open(txt_name, "a") as tf:
            tf.write(line2write+"\n")

    if end_time is not None:
        # get the start time from the file
        with open(txt_name, "r") as tf:
            for line in tf.readlines():
                if "Starting PTT running time" in line:
                    PTT_start_time = float(line.split(":")[-1])
                    break
        # compute end the timer to compute PTT running time
        PTT_total_time = end_time - PTT_start_time   # this is in seconds
        if PTT_total_time > 60.0:
            PTT_total_time_not_in_sec = round(PTT_total_time / 60.0, 2)   # in minutes
            PTT_total_run_time = repr(PTT_total_time_not_in_sec)+"min"
        if PTT_total_time_not_in_sec > 60.0:
            PTT_total_time_not_in_sec = round(PTT_total_time / 60.0, 2)   # in hrs
            PTT_total_run_time = repr(PTT_total_time_not_in_sec)+"hr"
        print ("The total time for PTT to run (including pipeline) was "+repr(PTT_total_time)+" seconds.")
        if "full_run_map" not in txt_name:
            line2write = "{:<20} {:<20} {:<20} {:<20}".format('', '', 'PTT_total_run_time  ', repr(PTT_total_time)+'  ='+PTT_total_run_time)
        else:
            line2write = "{:<20} {:<20}".format('PTT_total_run_time  ', repr(PTT_total_time)+'  ='+PTT_total_run_time)
        print (line2write)
        with open(txt_name, "a") as tf:
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


def set_inandout_filenames(step, config):
    """
    Set the step input and output file names, and add the step suffix to the step map.
    Args:
        step: string, name of the step as in the pipeline
        config: object, this is the configuration file object

    Returns:
        step_input_filename = string
        step_output_filename = string
        in_file_suffix = string, suffix of the input file
        out_file_suffix = string, suffix of the output file
        True_steps_suffix_map = string, path to the suffix map
    """
    step_dict = dict(config.items("steps"))
    initial_input_file = config.get("calwebb_spec2_input_file", "input_file")
    True_steps_suffix_map = config.get("calwebb_spec2_input_file", "True_steps_suffix_map")
    pytests_directory = os.getcwd()
    True_steps_suffix_map = os.path.join(pytests_directory, True_steps_suffix_map)
    steps_list, suffix_list, completion_list = read_True_steps_suffix_map(True_steps_suffix_map)
    step_input_filename = get_correct_input_step_filename(initial_input_file, steps_list,
                                                                     suffix_list, completion_list)
    suffix_and_filenames = get_step_inandout_filename(step, initial_input_file, step_dict)
    in_file_suffix, out_file_suffix, _, _ = suffix_and_filenames
    step_output_filename = step_input_filename.replace(".fits", out_file_suffix+".fits")
    print ("step_input_filename = ", step_input_filename)
    print ("step_output_filename = ", step_output_filename)
    return step_input_filename, step_output_filename, in_file_suffix, out_file_suffix, True_steps_suffix_map


def read_info4outputhdul(config, step_info):
    """
    Unfold the variables from the return of the function set_inandout_filenames, and get variables from the
    configuration file.
    Args:
        config: object, this is the configuration file object
        step_info: list, this is the output from the function set_inandout_filenames

    Returns:
        set_inandout_filenames_info = list, variables from set_inandout_filenames function and configuration file
    """
    initiate_calwebb_spc2 = "calwebb_spec2_input_file"
    working_directory = config.get(initiate_calwebb_spc2, "working_directory")
    step, step_input_filename, output_file, in_file_suffix, outstep_file_suffix, True_steps_suffix_map = step_info
    txt_name = os.path.join(working_directory, True_steps_suffix_map)
    step_input_file = os.path.join(working_directory, step_input_filename)
    step_output_file = os.path.join(working_directory, output_file)
    run_calwebb_spec2 = config.getboolean("run_calwebb_spec2_in_full", "run_calwebb_spec2")
    set_inandout_filenames_info = [step, txt_name, step_input_file, step_output_file, run_calwebb_spec2, outstep_file_suffix]
    return set_inandout_filenames_info


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
        if output_hdul["EXP_TYPE"] == "NRS_FIXEDSLIT":
            result = True
    return result


def check_BOTS_true(output_hdul):
    """
    This function checks if the fits file is a Bright Object Time Series, BOTS.
    Args:
        output_hdul: the HDU list of the output header keywords

    Returns:
        result: boolean, if true, the file is assumed to be BOTS
    """
    result = False
    if "EXP_TYPE" in output_hdul:
        if output_hdul["EXP_TYPE"] == "NRS_BRIGHTOBJ":
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
        if "MSA" in output_hdul["EXP_TYPE"]:
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
        if "IFU" in output_hdul["EXP_TYPE"]:
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
    slits = ["S200A1", "S200A2", "S200B1", "S400A1", "S1600A1", ]
    if "FXD_SLIT" in output_hdul:
        for i, s in enumerate(slits):
            print("slit: ", s, i)
            if s in output_hdul["FXD_SLIT"]:
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


def get_time_to_run_pipeline(True_steps_suffix_map):
    """
    This function calculates the total running time by adding the individual time per step.
    Args:
        True_steps_suffix_map: string, full path of where the text file will be written into

    Returns:
        total_time: string, total calculate the total time by reading the time of each step from the map

    """
    #times_per_step = np.loadtxt(True_steps_suffix_map, comments="#", usecols=(3), unpack=True)
    times_per_step = []
    with open(True_steps_suffix_map, "r") as tf:
        for line in tf.readlines():
            if '#' in line:
                continue
            # this is for "msa_flagging" and "srctype"
            if "=" in line:
                line  = line.split("=")[0]
            if len(line.split()) < 4:
                t=line.split()[2]
            else:
                t=line.split()[3]
            t = float(t)
            times_per_step.append(t)
    total_time = sum(times_per_step)
    return total_time


def get_reffile_used(output_hdul):
    '''
    This function obtains the list of reference files used in the processing of step_outfile_fits. It only looks for
    the keywords that match the text file of expected reference files used, located at: /auxiliary_code/ref_files.txt
    Args:
        output_hdul: dictionary, main header of the step output fits file

    Returns:
        ref_files_used_so_far: dictionary, reference files used in processing of step_outfile_fits
    '''

    # list of the keywords to look for according to the step (NOT in order of processing)
    # these are only the reference keywords for which we have an expected correspnding file in the text file:
    # /auxiliary_code/ref_files.txt
    txt_file_ref_keys = ["R_MASK", "R_SATURA", "R_SUPERB", "R_LINEAR",  "R_DARK", "R_READNO", "R_GAIN", "R_FORE",
                         "R_DISPER", "R_CAMERA", "R_COLLIM", "R_FPA", "R_IFUFOR", "R_IFUPOS", "R_IFUSLI", "R_MSA",
                         "R_OTE", "R_WAVRAN", "R_PTHLOS", "R_PHOTOM", "R_EXTR1D", "R_DFLAT", "R_FFLAT", "R_SFLAT"]

    # get the list of reference file keywords set so far
    ref_files_used_so_far = {}
    for reffile_key in txt_file_ref_keys:
        if reffile_key in output_hdul:
            reffile_key_value = output_hdul[reffile_key]
            ref_files_used_so_far[reffile_key] = reffile_key_value
        else:
            print ("Keyword: ", reffile_key, "  not found in main header.")

    return ref_files_used_so_far


def get_latest_file(filetype, disregard_known_files=False):
    """
    This function gets the latest modified file of filetype. This function will look into the calwebb_spec2_pytests
    directory for the given file.
    Args:
        filetype: string, name/type of file type, e.g. *.txt, *.html, full_run_map.txt
        disregard_known_files: boolean, if True function will not look for True_steps_suffix_map.txt or full_run_map.txt

    Returns:
        latest_filetypefile: string, name of the latest file modified of type filetype
    """
    # get a list of all the filetype files in the calwebb_spec2_pytests dir
    list_of_filetypefiles = glob.glob(filetype)
    # find the latest of the filetype files but exclude known file names
    if disregard_known_files:
        if "True_steps_suffix_map.txt" in list_of_filetypefiles:
            idx = list_of_filetypefiles.index("True_steps_suffix_map.txt")
            list_of_filetypefiles.pop(idx)
        if "full_run_map.txt" in list_of_filetypefiles:
            idx = list_of_filetypefiles.index("full_run_map.txt")
            list_of_filetypefiles.pop(idx)
    latest_filetypefile = "File not found."
    if len(list_of_filetypefiles) > 0:
        latest_filetypefile = max(list_of_filetypefiles, key=os.path.getctime)
    #print("filetype, latest_filetypefile = ", filetype, latest_filetypefile)
    return latest_filetypefile


def convert_html2pdf():
    """
    This function converts the latest html file into a pdf.
    In order to work it needs this plugin (https://pypi.org/project/pdfkit/#description):
    pip install pdfkit
    and the OSX 64-bit download from: https://wkhtmltopdf.org/downloads.html
    """
    # get a list of all the html files in the calwebb_spec2_pytests dir
    latest_htmlfile = get_latest_file("*.html")
    # create the pdf output name
    pdf_file = latest_htmlfile.replace(".html", ".pdf")
    # convert the html report into a pdf file
    #import pdfkit
    #options = {'dpi': 96}
    #pdfkit.from_file(latest_htmlfile, pdf_file, options=options)
    print ("\n Converted ", latest_htmlfile, " to ", pdf_file, ". Both files are available in current directory. \n")


def move_latest_report_and_txt_2workdir():
    """
    This function moves the PTT output reporting files into the working directory.
    Args:
        Nothing

    Returns:
        Nothing
    """
    # get the working directory
    config = configparser.ConfigParser()
    config.read(['../calwebb_spec2_pytests/PTT_config.cfg'])
    working_dir = config.get("calwebb_spec2_input_file", "working_directory")
    # get a list of all the html and txt files in the calwebb_spec2_pytests dir
    latest_htmlfile = get_latest_file("*.html")   # this will pick up the report.html just created/modified
    latest_pdffile = latest_htmlfile.replace(".html", ".pdf")   # this will pick up the report.html converted to pdf
    latest_screenoutputtxtfile = get_latest_file("*.txt", disregard_known_files=True) # this should pick up the output_screen.txt
    latest_suffixmaptxtfile = get_latest_file("True_steps_suffix_map.txt")
    latest_fullrunmaptxtfile = get_latest_file("full_run_map.txt")
    # move these files into the working directory
    files2move = [latest_htmlfile, latest_pdffile, latest_screenoutputtxtfile,
                  latest_suffixmaptxtfile, latest_fullrunmaptxtfile]
    for f in files2move:
        if f != "File not found."  and  os.path.isfile(f):
            subprocess.run(["mv", f, working_dir])
            print(f, " moved to ", working_dir)


