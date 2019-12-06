import argparse
import numpy as np
from astropy.io import fits

# import the modules needed within our tool
import move_data2ext1
import level2b_hdr_keywd_check as lev2bcheck

"""

This script is a wraper for the 2 scripts needed to convert the ESA simulations in CRM to ST pipeline 
ingestible products. 

Example usage:
    The code works from the terminal or called as a module.
    
    Terminal 
        To create a NEW fits file with the updated header type:
        > python /path_to_this_script/crm2STpipeline.py blah.fits MODE

        To add the reference pixels add at the end of the command -rfpx, e.g.
        > python /path_to_this_script/crm2STpipeline.py blah.fits MODE -rfpx
    
    As a module 
        import crm2STpipeline 
        crm2STpipeline.crm2pipe(input_fits_file, mode_used, add_ref_pix, only_update)
    
    * NOTE: In all cases the MODE can be any of the following: 
            FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf

"""

def crm2pipe(input_fits_file, mode_used, add_ref_pix, only_update=False):
    """
    This function is the wraper for the scripts needed to convert a crm file to a pipeline ready product.
    Args:
        input_fits_file: string, name and path of the input file
        mode_used: string, possible values are FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf
        add_ref_pix: boolean, if True it will add the reference pixels for full frame
        only_update: boolean, if False the code will create a new file with updated header

    Returns:
        Nothing. The final product is the file that is pipeline ready.
    """

    # Perform data move to the science extension if needed, i.e. if there is only 1 data extension
    hdulist = fits.open(input_fits_file)
    print("Contents of the input file")
    print(hdulist.info(), "\n")
    sci_ext, dq_ext, err_ext = False, False, False
    data_ext, var_ext, quality_ext = True, True, True
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
            lev2bcheck.perform_check(st_pipe_file, only_update, mode_used, det)

            print("Done with detector", det, "\n")

    else:
        print("The extension names are correct.")
        if len(hdulist) > 4:
            for det in detectors:
                print("Removing the extra extensions and rotating data for ST pipeline ingestion.")
                outfile_name = rm_extra_exts_and_rotate(input_fits_file, det)

                # Perform the keyword check on the file with the right number of extensions
                print("Fixing the header keywords")
                lev2bcheck.perform_check(outfile_name, only_update, mode_used, det)

        else:
            print("No need to rename or modify file for ST pipeline ingestion. Exiting script.")
            exit()


def transpose(arr, detector):
    transposed_array = np.transpose(arr)
    if "2" in detector:
        # rotate dq data by 180 degrees
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



if __name__ == '__main__':

    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("fits_file",
                        action='store',
                        default=None,
                        help='Name of fits file, i.e. blah.fits')
    parser.add_argument("mode_used",
                        #dest="mode_used",
                        action='store',
                        default=None,
                        help='Observation mode used: FS, MOS, IFU, BOTS, dark, image, confirm, taconfirm, wata, msata, focus, mimf.')
    parser.add_argument("-rfpx",
                        dest="add_ref_pix",
                        action='store_true',
                        default=False,
                        help='Add the reference pixels.')
    args = parser.parse_args()

    # Set the variables
    input_fits_file = args.fits_file
    mode_used = args.mode_used
    add_ref_pix = args.add_ref_pix
    only_update = False

    # Perform data move to the science extension and the keyword check on the file with the right number of extensions
    crm2pipe(input_fits_file, mode_used, add_ref_pix, only_update)

    print ('\n * Script  crm2STpipeline.py  finished * \n')


