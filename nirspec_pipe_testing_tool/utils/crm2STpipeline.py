import argparse
import numpy as np
import sys
from astropy.io import fits

# import the modules needed within our tool
from . import move_data2ext1
from . import level2b_hdr_keywd_check as lev2bcheck
from . import level2b_hdr_keywd_dict_map2sim as map2sim

"""

This script is a wrapper for the 2 scripts needed to convert the ESA simulations in CRM to ST pipeline to
ingest products.

Example usage:
    The code works from the terminal or called as a module.

    Terminal
        To create a NEW fits file with the updated header type:
        $ nptt_crm2STpipeline blah.fits MODE

        To add the reference pixels add at the end of the command -rfpx, e.g.
        $ nptt_path_to_this_script/crm2STpipeline blah.fits MODE -rfpx

    As a module
        # imports
        import nirspec_pipe_testing_tool as nptt
        
        # create the pipeline-ready count rate file
        stsci_pipe_ready_file = nptt.utils.crm2STpipeline.crm2pipe(input_fits_file, mode_used, add_ref_pix, only_update)
        
        # create dictionary of special arguments
        additional_args_dict = {'TITLE': proposal_title, 'TARGNAME': target_name, 'new_file': new_file}
        
        # modify the keyword values to match IPS information
        nptt.utils.level2b_hdr_keywd_dict_map2sim.match_IPS_keywords(stsci_pipe_ready_file, input_fits_file, 
                                                                     additional_args_dict=additional_args_dict)

    * NOTE: In all cases the MODE can be any of the following:
            FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf

"""

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.1"

# HISTORY
# May 2019 - Version 1.0: initial version completed
# Feb 2020 - Version 1.1: added part to match and replace keyword values from IPS file


def crm2pipe(input_fits_file, mode_used, add_ref_pix, only_update=False, subarray=None):
    """
    This function is the wraper for the scripts needed to convert a crm file to a pipeline ready product.
    Args:
        input_fits_file: string, name and path of the input file
        mode_used: string, possible values are FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf
        add_ref_pix: boolean, if True it will add the reference pixels for full frame
        only_update: boolean, if False the code will create a new file with updated header
        subarray: None or string, name of subarray to use

    Returns:
        out_fits: string, path and name of the final product is the file that is pipeline ready.
    """

    # set the subarray value to what the code expects
    if subarray is not None:
        if "full" in subarray.lower():
            subarray = "FULL-FRAME"
        elif "all" in subarray.lower():
            subarray = "ALLSLITS"
        elif "200a1" in subarray.lower():
            subarray = "S200A1"
        elif "200a2" in subarray.lower():
            subarray = "S200A2"
        elif "200b1" in subarray.lower():
            subarray = "S200B1"
        elif "400" in subarray.lower():
            subarray = "S400A1"
        elif "1024a" in subarray.lower():
            subarray = "SUB1024A"
        elif "1024b" in subarray.lower():
            subarray = "SUB1024B"
        elif "2048" in subarray.lower():
            subarray = "SUB2048"
        elif "32" in subarray.lower():
            subarray = "SUB32"
        elif "512" in subarray.lower():
            subarray = "SUB512"
        elif "512s" in subarray.lower():
            subarray = "SUB512S"
        else:
            print(" * Oh oh, subarray {} is not a recognized value.".format(subarray))
            print("   Recognized values: FULL, 200A1, 200A2, 400, 1600, 1024A, 1024B, 2048, 32, 512, 512S")
            print("   Exiting script.\n")
            exit()

    # Perform data move to the science extension if needed, i.e. if there is only 1 data extension
    hdulist = fits.open(input_fits_file)
    print("Contents of the input file")
    print(hdulist.info(), "\n")
    sci_ext, dq_ext, err_ext = False, False, False
    data_ext, var_ext, quality_ext, header_ext = True, True, True, True
    ext_name_list = []
    for hdu in hdulist:
        if hdu.name == "PRIMARY":
            continue
        if hdu.name == "SCI":
            sci_ext = True
        elif hdu.name == "DQ":
            dq_ext = True
        elif hdu.name == "ERR":
            err_ext = True
        if "DATA" not in hdu.name:
            data_ext = False
        elif "VAR" not in hdu.name:
            var_ext = False
        elif "QUALITY" not in hdu.name:
            quality_ext = False
        elif "HEADER" not in hdu.name:
            header_ext = False
        ext_name_list.append(hdu.name)

    # Determine which detector to use
    detectors = ["NRS1", "NRS2"]
    if all("1" in ext for ext in ext_name_list):
        #print("only 1 in file name")
        detectors = ["NRS1"]
    elif all("2" in ext for ext in ext_name_list):
        #print("only 2 in file name")
        detectors = ["NRS2"]
    elif "NRS1" in input_fits_file or "491" in input_fits_file:
        #print("NRS1 in file name")
        detectors = ["NRS1"]
    elif "NRS2" in input_fits_file or "492" in input_fits_file:
        #print("NRS2 in file name")
        detectors = ["NRS2"]
    if len(hdulist) == 5:
        detectors = [fits.getval(input_fits_file, "DET", 0)]

    if not data_ext and not var_ext and not quality_ext:
        print("The extension names do not match the expected: DATA, VAR, QUALITY or SCI, DQ, ERR.")
        print("Please make sure you are using the correct file. Exiting script.")
        exit()

    if not sci_ext and not dq_ext and not err_ext:
        print("Renaming extensions for ST pipeline...")
        for det in detectors:
            # rename extensions in the file to match expected names in the pipeline
            move_data2ext1.move_data(input_fits_file, det, add_ref_pix)
            st_pipe_file = input_fits_file.replace(".fits", "_" + det + "_modified.fits")

            # perform rotations expected in the pipeline
            print("Rotating data for ST pipeline ingestion.")
            st_pipe_file = rm_extra_exts_and_rotate(st_pipe_file, det)

            # Perform the keyword check on the file with the right number of extensions
            print("Fixing the header keywords for detector ", det)
            out_fits = lev2bcheck.check_lev2b_hdr_keywd(st_pipe_file, only_update, mode_used, det, subarray)

            print("Done with detector", det, "\n")

    else:
        print("The extension names are correct.")
        if len(hdulist) > 4:
            for det in detectors:
                print("Removing the extra extensions and rotating data for ST pipeline ingestion.")
                outfile_name = rm_extra_exts_and_rotate(input_fits_file, det)

                # Perform the keyword check on the file with the right number of extensions
                print("Fixing the header keywords")
                out_fits = lev2bcheck.check_lev2b_hdr_keywd(outfile_name, only_update, mode_used, det, subarray)

        else:
            print("No need to rename or modify file for ST pipeline ingestion. Exiting script.")
            exit()

    return out_fits


def transpose(arr, detector):
    transposed_array = np.transpose(arr)
    if "2" in detector:
        # rotate data by 180 degrees
        transposed_array = np.rot90(transposed_array)
        transposed_array = np.rot90(transposed_array)
    return transposed_array


def rm_extra_exts_and_rotate(input_fits_file, detector):
    """
    This function removes the extra extensions for the ones needed for ingestion in the STScI pipeline - these are
    SCI, DQ, and ERR. Rotate the data if necessary.
    Args:
        input_fits_file: string, name and path of the fits input file
        detector: list, expected to have one or both strings, NRS1, NRS2

    Returns:
        outfile_name: string, name and path of the fits file with the 3 extensions the pipeline expects

    """

    # create the fits list to hold the calculated flat values for each slit
    original_hdulist = fits.open(input_fits_file)

    # append the extensions with the right name and orientation
    outfile = fits.HDUList()
    outfile.append(original_hdulist[0])
    input_fits_file_data = fits.getdata(input_fits_file, "SCI")
    transposed_array = transpose(input_fits_file_data, detector)   # Correct the orientation for pipeline processing
    outfile.append(fits.ImageHDU(transposed_array, name="SCI"))

    input_fits_file_data = fits.getdata(input_fits_file, "ERR")
    transposed_array = transpose(input_fits_file_data, detector)   # Correct the orientation for pipeline processing
    outfile.append(fits.ImageHDU(transposed_array, name="ERR"))

    input_fits_file_data = fits.getdata(input_fits_file, "DQ")
    transposed_array = transpose(input_fits_file_data, detector)   # Correct the orientation for pipeline processing
    outfile.append(fits.ImageHDU(transposed_array, name="DQ"))

    # write the new output file
    if detector not in input_fits_file and 'modified' not in input_fits_file:
        outfile_name = input_fits_file.replace(".fits", "_" + detector + "_modified.fits")
    else:
        outfile_name = input_fits_file
    outfile.writeto(outfile_name, overwrite=True)

    return outfile_name


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("ips_file",
                        action='store',
                        default=None,
                        help='Name of IPS fits file, i.e. blah.fits')
    parser.add_argument("mode_used",
                        action='store',
                        default=None,
                        help='Observation mode options: FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, '
                             'msata, focus, mimf')
    parser.add_argument("-r",
                        dest="add_ref_pix",
                        action='store_true',
                        default=False,
                        help='Add the reference pixels.')
    parser.add_argument("-p",
                        dest="proposal_title",
                        action='store',
                        default=None,
                        help='Add the proposal title to the header keyword, i.e. -p=some_title')
    parser.add_argument("-t",
                        dest="target_name",
                        action='store',
                        default=None,
                        help='Add the target name to the header keyword, i.e. -t=some_target')
    parser.add_argument("-n",
                        dest="new_file",
                        action='store_true',
                        default=True,
                        help='Use -n if wanting to create a new file with updated header. Default is to update '
                             'header without creating a new file')
    parser.add_argument("-s",
                        dest="subarray",
                        action='store',
                        default=None,
                        help='Use -s to use a specific subarray.')
    args = parser.parse_args()

    # Set the variables
    ips_file = args.ips_file
    mode_used = args.mode_used
    add_ref_pix = args.add_ref_pix
    proposal_title = args.proposal_title
    target_name = args.target_name
    new_file = args.new_file
    subarray = args.subarray

    # Perform data move to the science extension and the keyword check on the file with the right number of extensions
    stsci_pipe_ready_file = crm2pipe(ips_file, mode_used, add_ref_pix, new_file, subarray)

    # create dictionary of command-line arguments
    additional_args_dict = {'TITLE': proposal_title,
                            'TARGNAME': target_name,
                            'new_file': new_file
                            }
    # modify the keyword values to match IPS information
    print('Matching IPS keyword values to corresponding STScI pipeline keywords...')
    map2sim.match_IPS_keywords(stsci_pipe_ready_file, ips_file, additional_args_dict=additional_args_dict)

    print('\n * Script  crm2STpipeline.py  finished * \n')


if __name__ == '__main__':
    sys.exit(main())
