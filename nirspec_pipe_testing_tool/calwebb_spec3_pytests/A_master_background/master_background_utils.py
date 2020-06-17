from nirspec_pipe_testing_tool.calwebb_spec2_pytests.auxiliary_code.reffile_test import create_rfile_test

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Jun 2020 - Version 1.0: initial version completed


def create_completed_steps_txtfile(txt_suffix_map, step_input_file):
    """
    This function creates the completed steps along with the corresponding suffix of the output file name into a text file.
    Args:
        txt_suffix_map: string, full path of where the text file will be written into
        step_input_file: string, name of the input file for the pipeline step

    Returns:
        Nothing. A text file will be created in the pytests directory where all steps will be added
    """
    # name of the text file to collect the step name and suffix
    line0 = "# {:<20}".format("Input file: "+step_input_file)
    line1 = "# {:<17} {:<20} {:<20} {:<20}".format("Step", "Added suffix", "Step complition", "Time to run [s]")
    print(line1)
    with open(txt_suffix_map, "w+") as tf:
        tf.write(line0+"\n")
        tf.write(line1+"\n")


def print_time2file(txt_name, end_time, string2print):
    """
    This function is only used in the case of running the pipeline in full. It changes the total running time to the
    appropriate units and returns a string of the right spacing.
    Args:
        txt_name: string, name of the suffix map file
        end_time: string, time it took to run
        string2print: string, phrase to be printed in the file before the time

    Returns:
        nothing
    """
    if float(end_time) > 60.0:
        total_time_min = repr(round(float(end_time)/60.0, 1))   # in minutes
        total_time = total_time_min+'min'
        if float(total_time_min) > 60.0:
            total_time_hr = repr(round(float(end_time)/60.0, 1))   # in hours
            total_time = total_time_hr+'hr'
        end_time = end_time+'  ='+total_time
    line2write = "{:<20} {:<20} {:<20} {:<20}".format('', '', string2print, end_time)
    print(line2write)
    with open(txt_name, "a") as tf:
        tf.write(line2write+"\n")


# VERIFICATION FUNCTIONS

def wavstart_exists(output_hdul):
    """
    This function checks that the keyword WAVSTART was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "WAVSTART" in output_hdul
    return result


# REFERENCE FILES checks

#camera_rfile_is_correct = create_rfile_test("R_CAMERA", "camera model")

