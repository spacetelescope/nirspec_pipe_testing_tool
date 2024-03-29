import argparse
import shutil
import os
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

        To add the reference pixels add at the end of the command -rfpx, only old IFU simulations (<2020) needed this. 
        
        To change the name of the MSA metafile in the file header use the flag -m, e.g. -m=V962150_msa.fits

    As a module
        # imports
        import nirspec_pipe_testing_tool as nptt
        
        # create the pipeline-ready count rate file, create dictionary of special arguments, and
        # modify the keyword values to match IPS information
        
        # set the variables
        ips_file = '/path_to_crm_file/crm.fits'
        mode_used = 'MOS'  # One of FS, MOS, IFU, or BOTS
        add_ref_pix = False  # Add the reference pixels for IFU - old simulations (<2020) needed this
        proposal_title = 'some cool title'
        target_name = 'some target name'
        subarray = 'FULL-FRAME'  # name of the subarray to use
        new_file = False  # create a new file with the modified/fixed header
        msa_metafile = 'V962150_msa.fits'
          
        # run the script
        nptt.utils.crm2STpipeline.crm2STpipeline(ips_file, mode_used, add_ref_pix, proposal_title, target_name, 
                                                 subarray=subarray, new_file=new_file, msa_metafile=msa_metafile)

    * NOTE: In all cases the MODE can be any of the following: FS, MOS, IFU, BOTS

"""

# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.4"

# HISTORY
# May 2019 - Version 1.0: initial version completed
# Feb 2020 - Version 1.1: added part to match and replace keyword values from IPS file
# Feb 2021 - Version 1.2: implemented create a metafile for MOS data
# Apr 2021 - Version 1.3: implemented rotation of data depending on data array dimensions
# Jun 2022 - Version 1.4: corrected naming and implemented rotation_needed as a switch


def crm2pipe(input_fits_file, mode_used, add_ref_pix, new_file, subarray=None, msa_metafile='N/A',
             output_dir=None, force_rotation=False, rotation_needed=False, verbose=False):
    """
    This function is the wrapper for the scripts needed to convert a crm file to a pipeline ready product.
    Args:
        input_fits_file: string, name and path of the input file
        mode_used: string, possible values are FS, MOS, IFU, BOTS
        add_ref_pix: boolean, if True it will add the reference pixels for full frame
        new_file: boolean, if True the code will create a new file with updated header
        subarray: None or string, name of subarray to use
        msa_metafile: string, name of the MSA metafile
        output_dir: string, path to place the output file - if None output will be in same dir as input
        force_rotation: boolean, if True data will be rotated
        rotation_needed: boolean,  True value includes the case for 'RAW' data
        verbose: boolean

    Returns:
        out_fits: string, path and name of the final product is the file that is pipeline ready.
    """
    # determine if a new file is to be generated or only update
    only_update = True
    if new_file:
        only_update = False

    if 'mos' in mode_used.lower() or mode_used.lower() == 'msa':
        if msa_metafile != 'N/A':
            if verbose:
                print('(crm2STpipeline.crm2pipe:) MSA metafile will be set to ', msa_metafile)
        else:
            if verbose:
                print('(crm2STpipeline.crm2pipe:) WARNING! MSA metafile not specified, will be set to N/A')

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
        elif "1600" in subarray.lower():
            subarray = "S1600A1"
        else:
            print("(crm2STpipeline.crm2pipe:) WARNING! Uh oh, subarray {} is not a recognized value.".format(subarray))
            print("                           Recognized values: FULL, 200A1, 200A2, 400, 1600, 1024A, 1024B, 2048, "
                  "32, 512, 512S")
            print("                           Exiting script.\n")
            exit()

    # Perform data move to the science extension if needed, i.e. if there is only 1 data extension
    hdulist = fits.open(input_fits_file)
    if verbose:
        print("(crm2STpipeline.crm2pipe:) Contents of the input file ")
    hdulist.info()
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
        detectors = ["NRS1"]
    elif all("2" in ext for ext in ext_name_list):
        detectors = ["NRS2"]
    elif "NRS1" in input_fits_file.upper() or "491" in input_fits_file:
        detectors = ["NRS1"]
    elif "NRS2" in input_fits_file.upper() or "492" in input_fits_file:
        detectors = ["NRS2"]
    if len(hdulist) == 5:
        detectors = [fits.getval(input_fits_file, "DET", 0)]

    # determine if rotation is needed
    if not force_rotation:
        # rotation_needed = True  # this includes the case for value='raw'
        try:
            file_type = fits.getval(input_fits_file, "FILETYPE", 0)
            if 'cts' in file_type.lower() or 'crm' in file_type.lower():
                rotation_needed = False
            # ensure that in this case the data is rotated
            if 'uncal' in file_type.lower():
                rotation_needed = True
            print("(crm2STpipeline.crm2pipe:) keyword FILETYPE = ", file_type)
            print("(crm2STpipeline.crm2pipe:) rotation_needed = ", rotation_needed)
            if not rotation_needed:
                rotation_input = input("Force data to be rotated? [Y]/n ")
                if 'n' in rotation_input.lower():
                    rotation_needed = False
                else:
                    rotation_needed = True
        except KeyError:
            print("(crm2STpipeline.crm2pipe:) FILETYPE keyword not found. Data will be rotated.")
    else:
        rotation_needed = force_rotation
    if verbose:
        print("(crm2STpipeline.crm2pipe:) Will the data be rotated? ", rotation_needed)

    hdulist.close()

    if not data_ext and not var_ext and not quality_ext:
        print("(crm2STpipeline.crm2pipe:) ERROR! The extension names do not match the expected: DATA, VAR, QUALITY or "
              "SCI, DQ, ERR.")
        print("                           Please make sure you are using the correct file. Exiting script.")
        exit()

    if not sci_ext:   # and not dq_ext and not err_ext:
        print("(crm2STpipeline.crm2pipe:) Renaming extensions for ST pipeline...")
        for det in detectors:
            # rename extensions in the file to match expected names in the pipeline
            print("(crm2STpipeline.crm2pipe:) Renaming extensions and moving data for ST pipeline ingestion...")
            st_pipe_file = move_data2ext1.move_data(input_fits_file, det, add_ref_pix=add_ref_pix,
                                                    output_dir=output_dir)

            # perform rotations expected in the pipeline
            if rotation_needed:
                print("(crm2STpipeline.crm2pipe:) Rotating data for ST pipeline ingestion...")
                st_pipe_file = rm_extra_exts_and_rotate(st_pipe_file, det, output_dir=output_dir)

            # Perform the keyword check on the file with the right number of extensions
            print("(crm2STpipeline.crm2pipe:) Fixing the header keywords for detector ", det)
            out_fits = lev2bcheck.check_lev2b_hdr_keywd(st_pipe_file, only_update, mode_used, detector=det,
                                                        subarray=subarray, msa_metafile=msa_metafile,
                                                        mktxt=False, verbose=verbose)

            print("(crm2STpipeline.crm2pipe:) Done with detector", det, "\n")

    else:
        print("(crm2STpipeline.crm2pipe:) The extension names are correct.")
        if len(hdulist) > 4:
            for det in detectors:
                if rotation_needed:
                    print("(crm2STpipeline.crm2pipe:) Removing the extra extensions and rotating data for ST "
                          "pipeline ingestion.")
                    outfile_name = rm_extra_exts_and_rotate(input_fits_file, det, output_dir=output_dir)
                else:
                    outfile_name = input_fits_file

                # Perform the keyword check on the file with the right number of extensions
                print("(crm2STpipeline.crm2pipe:) Fixing the header keywords")
                out_fits = lev2bcheck.check_lev2b_hdr_keywd(outfile_name, only_update, mode_used, detector=det,
                                                            subarray=subarray, msa_metafile=msa_metafile,
                                                            mktxt=False, verbose=verbose)

        else:
            print("(crm2STpipeline.crm2pipe:) No need to rename or modify file for ST pipeline ingestion. "
                  "Exiting script.")
            exit()

    return out_fits


def transpose(arr, detector):
    transposed_array = np.transpose(arr)
    if "2" in detector:
        # rotate data by 180 degrees
        transposed_array = np.rot90(transposed_array)
        transposed_array = np.rot90(transposed_array)
    return transposed_array


def rm_extra_exts_and_rotate(input_fits_file, detector, output_dir=None):
    """
    This function removes the extra extensions for the ones needed for ingestion in the STScI pipeline - these are
    SCI, DQ, and ERR. Rotate the data if necessary.
    Args:
        input_fits_file: string, name and path of the fits input file
        detector: list, expected to have one or both strings, NRS1, NRS2
        output_dir: string, path to place the output file - if None output will be in same dir as input
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
    outfile_name = input_fits_file
    if detector.lower() not in outfile_name.lower():
        outfile_name = outfile_name.replace(".fits", "_" + detector + ".fits")
    if 'modified' not in outfile_name:
        outfile_name = outfile_name.replace(".fits", "_modified.fits")
    if output_dir is not None:
        outfile_name = os.path.join(output_dir, os.path.basename(outfile_name))
        if '.fits' in output_dir:
            outfile_name = output_dir
    outfile.writeto(outfile_name, overwrite=True)
    original_hdulist.close()
    return outfile_name


def crm2STpipeline(ips_file, mode_used, add_ref_pix, proposal_title, target_name, subarray=None, new_file=False,
                   msa_metafile='N/A', output_dir=None, force_rotation=False,
                   rotation_needed=False, verbose=False):
    """
    This function is the wrapper for the scripts needed to convert a crm file to a pipeline ready product.
    :param ips_file: string, full path to IPS crm fits file
    :param mode_used: string, possible values are FS, MOS, IFU, BOTS
    :param add_ref_pix: boolean, if True it will add the reference pixels for full frame
    :param proposal_title: string, proposal title to be added to the header
    :param target_name: string, name of the target to be added to the header
    :param subarray: None or string, name of subarray to use
    :param new_file: boolean, if True create a new file with updated header. Default is to update
                     header without creating a new file
    :param msa_metafile: string, name of the MSA metafile
    :param output_dir: string, path to place the output file - if None output will be in same dir as input
    :param force_rotation: boolean, if True data will be rotated to the science frame
    :param rotation_needed: boolean,  True value includes the case for 'RAW' data
    :param verbose: boolean
    :return:
    """
    # Perform data move to the science extension and the keyword check on the file with the right number of extensions
    stsci_pipe_ready_file = crm2pipe(ips_file, mode_used, add_ref_pix, new_file, subarray=subarray,
                                     msa_metafile=msa_metafile, output_dir=output_dir,
                                     force_rotation=force_rotation,
                                     rotation_needed=rotation_needed, verbose=verbose)

    # create dictionary of command-line arguments
    additional_args_dict = {'TITLE': proposal_title,
                            'TARGNAME': target_name,
                            'new_file': new_file
                            }

    # modify the keyword values to match IPS information
    map2sim.match_IPS_keywords(stsci_pipe_ready_file, ips_file, additional_args_dict=additional_args_dict,
                               verbose=verbose)
    return stsci_pipe_ready_file


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
                        help='Observation mode options: FS, MOS, IFU, BOTS')
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
                        action='store_false',
                        default=True,
                        help='Use -n if wanting to create a new file with updated header. Default is to update '
                             'header without creating a new file')
    parser.add_argument("-s",
                        dest="subarray",
                        action='store',
                        default=None,
                        help='Use -s to use a specific subarray.')
    parser.add_argument("-m",
                        dest="msa_metafile",
                        action='store',
                        default=None,
                        help='Use -m to use a specific MSA metafile. Only applies to MOS mode.')
    parser.add_argument("-o",
                        dest="output_dir",
                        action='store',
                        default=None,
                        help='Use -o to provide a path to place the output fits file.')
    parser.add_argument("-f",
                        dest="force_rotation",
                        action='store_true',
                        default=False,
                        help='Use -f to force the data to be rotated to the Science frame. Default '
                             'is to rotate only for raw-type data (which is in Detector frame).')
    parser.add_argument("-rn",
                        dest="rotation_needed",
                        action='store_true',
                        default=False,
                        help='Use -rn to indicate if rotaton is needed for this data. Default is False.')
    parser.add_argument("-v",
                        dest="verbose",
                        action='store_true',
                        default=False,
                        help='Use -v to print on-screen keywords and values.')
    args = parser.parse_args()

    # Set the variables
    ips_file = args.ips_file
    mode_used = args.mode_used
    add_ref_pix = args.add_ref_pix
    proposal_title = args.proposal_title
    target_name = args.target_name
    new_file = args.new_file
    subarray = args.subarray
    msa_metafile = args.msa_metafile
    output_dir = args.output_dir
    force_rotation = args.force_rotation
    rotation_needed = args.rotation_needed
    verbose = args.verbose

    # Run wrapper function
    crm2STpipeline(ips_file, mode_used, add_ref_pix, proposal_title, target_name, subarray=subarray,
                   new_file=new_file, msa_metafile=msa_metafile, output_dir=output_dir,
                   force_rotation=force_rotation, rotation_needed=rotation_needed, verbose=verbose)

    print('\n * Script  crm2STpipeline.py  finished * \n')


if __name__ == '__main__':
    sys.exit(main())
