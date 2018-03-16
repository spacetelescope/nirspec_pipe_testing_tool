import os
import argparse
import subprocess
from glob import glob
from astropy.io import fits

# import the modules to prepare the data
import hdr_keywd_check as hkch


"""
This script will prepare the data to run through the pipeline via the create_data script and a header keyword check.

This script ** NEEDS ** to be run once PER fits file BEFORE running the pipeline. It will create a subdirectory in
order to run the create_data script, and then erase it (if asked to) once the job is done.

Example usage:
    The code works from the terminal within the pipeline conda environment.
    In the directory where you want to store the uncalibrated data, type:
        > python /path_to_this_script/prepare_data2run.py blah.fits MODE
    where MODE is either FS, IFU, or MOS
"""


def modify_PTT_cfg_file(fits_file, mode):
    """
    This function reads and modifies the config file with the mode used and the raw data file name.
    Args:
        fits_file: string, basename of the raw data file.fits
        mode: string, either FS, MOS, or IFU

    Returns:
        nothing
    """
    # get script directory and config name
    utils_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    PPT_cfg_file = utils_dir.replace("utils", "calwebb_spec2_pytests/cwspec2_config.cfg")
    with open(PPT_cfg_file, "r") as cfg:
        cfg_lines = []
        for line in cfg.readlines():
            if "#" not in line:
                if "mode_used" in line:
                    line = "mode_used = "+mode+"\n"
                if "raw_data_root_file" in line:
                    line = "raw_data_root_file = "+fits_file+"\n"
            cfg_lines.append(line)
    subprocess.run(["rm", PPT_cfg_file])
    with open(PPT_cfg_file, "w") as cfg:
        for line in cfg_lines:
            cfg.write(line)


# Get arguments to run script
parser = argparse.ArgumentParser(description='')
parser.add_argument("fits_file",
                    action='store',
                    default=None,
                    help='Name of fits file, i.e. blah.fits')
parser.add_argument("mode",
                    action='store',
                    default=None,
                    help='Mode in which the data was taken, i.e. FS, MOS, IFU')
parser.add_argument("-rm",
                    dest="rm_prep_data",
                    action='store_true',
                    default=False,
                    help='If -rm is used, the directory prep_data will be removed once the job is done.')
parser.add_argument("-u",
                    dest="only_update",
                    action='store_true',
                    default=False,
                    help='Use -u if NOT wanting to create a new file with updated header.')
args = parser.parse_args()

# Set the variables inputed from the command line
fits_file = args.fits_file
mode = args.mode
rm_prep_data = args.rm_prep_data
only_update = args.only_update

# Create the directory to run the create_data script
subprocess.call(["mkdir", "prep_data"])

# Put there the input file
subprocess.call(["cp", fits_file, "prep_data/"])

# Create the data.prop file
datapropfile = "prep_data/data.prop"
with open(datapropfile, "w") as dpf:
    dpf.write("<Proposal title='Test1'>\n")
    dpf.write("  <Observation>\n")
    dpf.write("    <Visit>\n")
    dpf.write("      <VisitGroup>\n")
    dpf.write("        <ParallelSequenceID>\n")
    dpf.write("          <Activity>\n")
    dpf.write("            <Exposure>\n")
    dpf.write("              <Detector>\n")
    dpf.write("                <base>"+fits_file+"</base>\n")
    dpf.write("                <subarray></subarray>\n")
    dpf.write("                <exp_type></exp_type>\n")
    dpf.write("              </Detector>\n")
    dpf.write("            </Exposure>\n")
    dpf.write("          </Activity>\n")
    dpf.write("        </ParallelSequenceID>\n")
    dpf.write("      </VisitGroup>\n")
    dpf.write("    </Visit>\n")
    dpf.write("  </Observation>\n")
    dpf.write("</Proposal>\n")

# run the create_data script
subprocess.run(["create_data", "prep_data"])

# if create_data was successfull the new uncalibrated file must have appeared in the prep_data directory
uncal_file = glob("prep_data/*_uncal.fits")[0]
if os.path.isfile(uncal_file):

    # move the fixed uncal file out of prep_data
    subprocess.run(["mv", uncal_file, "."])
    uncal_file = os.path.basename(uncal_file)
    fits_file = os.path.basename(fits_file)

    # add keywords with the name of the raw data fits file name , and mode used
    try:
        modify_PTT_cfg_file(fits_file, mode)
    except:
        IsADirectoryError

    # make sure the lamp, filter, and grating values are correctly propagated
    lamp = fits.getval(fits_file, 'CAA_LAMP')
    if "NO_LAMP" in lamp:
        try:
            filt = fits.getval(fits_file, "FWA_POS")
        except:
            print (" *** Unable to determine what was the FILTER used...  :(  ")
            print (" The FILTER keyword has to be set up manually in order for the pipeline to be able to process data.")
            filt = ""
            KeyError

    elif 'LINE1' in lamp:
        filt = 'F100LP'
    elif 'LINE2' in lamp:
        filt = 'F170LP'
    elif 'LINE3' in lamp:
        filt = 'F290LP'
    elif 'LINE4' in lamp:
        filt = 'CLEAR'
    elif 'FLAT4' in lamp:
        filt = 'F070LP'
    fits.setval(uncal_file, 'LAMP', 0, value=lamp)
    fits.setval(uncal_file, 'FILTER', 0, value=filt, after='DETECTOR')

    # add the missing header keywords and fix the format to the one the pipeline expects
    uncal_file = glob("*_uncal.fits")[0]
    hkch.perform_check(uncal_file, only_update, mode)

else:
    print("The *uncal.fits file does not exist. Exiting the code.")
    exit()

# if asked to, remove the prep_data directory
if rm_prep_data:
    subprocess.run(["rm", "-R", "prep_data"])

print ('\n * Script  prepare_data.py  finished * \n')
