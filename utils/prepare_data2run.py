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
        > python /path_to_this_script/prepare_data2run.py blah.fits
"""

# Get arguments to run script
parser = argparse.ArgumentParser(description='')
parser.add_argument("fits_file",
                    action='store',
                    default=None,
                    help='Name of fits file, i.e. blah.fits')
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

    # add the missing header keywords and fix the format to the one the pipeline expects
    uncal_file = glob("*_uncal.fits")[0]
    hkch.perform_check(uncal_file, only_update)

    # add a keyword with the name of the raw data fits file name
    fits.setval(uncal_file, 'rawdatroot', 0, value=uncal_file, after='DPSW_VER')

else:
    print("The *uncal.fits file does not exist. Exiting the code.")
    exit()

# if asked to, remove the prep_data directory
if rm_prep_data:
    subprocess.run(["rm", "-R", "prep_data"])

print ('\n * Script  prepare_data.py  finished * \n')