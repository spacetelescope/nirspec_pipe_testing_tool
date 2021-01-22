

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.1"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Jun 2020 - Version 1.1: implemented APT user input vs pipeline logic

"""
This file contains the functions which will be used to test the source_type step
of the JWST Calibration Pipeline.

"""


# VERIFICATION FUNCTIONS

def s_srctyp_exists(output_hdul):
    """
    This function checks that the keyword SRCTYPE was added.
    :param output_hdul:
    :return:
        result: boolean, true if the keyword was indeed added
        msg: string, output message
    """
    result = "S_SRCTYP" in output_hdul
    return result


def srctype_logic_correct(output_hdul):
    """
    Logic of the Source_type step keyword values according to the user input value recorded in the SRCTYAPT keyword.
    :param output_hdul:
    :return:
    """
    possible_apt_values = ['UNKNOWN', 'POINT', 'EXTENDED']
    # output_hdul = hdul, step_output_file, run_pytests, scihdul
    srctyapt = output_hdul[0]['SRCTYAPT']
    exptype = output_hdul[0]['EXP_TYPE']
    srctype = output_hdul[3]['SRCTYPE']
    result = False
    msg = "SRCTYPE keyword in SCI extension matches expected value."

    # unknown in APT, then defaults should be set in pipeline
    if srctyapt == possible_apt_values[0]:
        # default for FS and MOS should be point source, and for IFU extended
        if exptype == 'NRS_IFU':
            if srctype == 'EXTENDED':
                result = True
            else:
                msg = 'SRCTYPE keyword WRONG in SCI extension. With SRCTYAPT=UNKNOWN for IFU data ' \
                      'expected_value=EXTENDED, pipeline_value=' + srctype
        else:
            if srctype == 'POINT':
                result = True
            else:
                msg = 'SRCTYPE keyword WRONG in SCI extension. With SRCTYAPT=UNKNOWN for non-IFU data ' \
                      'expected_value=POINT, pipeline_value=' + srctype

    # point source in APT, then pipeline should be the same
    if srctyapt == possible_apt_values[1]:
        if srctype == 'POINT':
            result = True
        else:
            msg = 'SRCTYPE keyword WRONG in SCI extension. With SRCTYAPT=POINT expected_value=POINT, ' \
                  'pipeline_value=' + srctype

    # extended source in APT, then pipeline should be the same
    if srctyapt == possible_apt_values[2]:
        if srctype == 'EXTENDED':
            result = True
        else:
            msg = 'SRCTYPE keyword WRONG in SCI extension. With SRCTYAPT=EXTENDED expected_value=EXTENDED, ' \
                  'pipeline_value=' + srctype

    return result, msg

