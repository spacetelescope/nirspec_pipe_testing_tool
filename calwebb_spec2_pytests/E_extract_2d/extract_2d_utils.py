import subprocess
import os
from collections import OrderedDict
from astropy.io import fits
from ..auxiliary_code import auxiliary_functions as auxfunc


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed


"""
This file contains the functions which will be used to test the extract_2d step
of the JWST Calibration Pipeline.

Selected keywords are checked to verify that the step ran through successfully.
"""


### VERIFICATION FUNCTIONS

def s_ext2d_exists(output_hdul):
    """
    This function checks that the keyword S_EXTR2D was added.
    Args:
        outout_hdul: the HDU list of the header keywords

    Returns:
        result: boolean, true if the keyword was indeed added
    """
    result = "S_EXTR2D" in output_hdul
    return result


def find_FSwindowcorners(infile_name, esa_files_path):
    """
	Find the slits corners of pipeline and ESA file and determine if they match.
	Args:
		infile_name - string, name of the input pipeline fits file
		esa_files_path - string, path to locate esa files

	Retrurs:
		result - string or dictionary
				if string, test is set to skip
				if dictionary, contains a boolean per slit test
	"""

    result = OrderedDict()

    #iterate over slits
    sci_ext_list = auxfunc.get_sci_extensions(infile_name)

    primary_header = fits.getheader(infile_name, ext=0)
    detector = primary_header["DETECTOR"]

    for i, s_ext in enumerate(sci_ext_list):
        sci_header = fits.getheader(infile_name, ext=s_ext)

        #grab corners of extracted subwindow
        px0 = sci_header["SLTSTRT1"] + primary_header["SUBSTRT1"] - 1
        py0 = sci_header["SLTSTRT2"] + primary_header["SUBSTRT2"] - 1
        pnx = sci_header["SLTSIZE1"]
        pny = sci_header["SLTSIZE2"]

        px1 = px0 + pnx
        py1 = py0 + pny

        icorners = {(px0, py0), (px1, py0), (px1, py1), (px1, py0)}

        # Find esafile (most of this copy-pasted from compare_wcs_fs)
        sltname = sci_header["SLTNAME"]
        _, raw_data_root_file = auxfunc.get_modeused_and_rawdatrt_PTT_cfg_file()
        specifics = [sltname]
        # check if ESA data is not in the regular directory tree
        NIDs = ["30055", "30055", "30205", "30133", "30133"]
        special_cutout_files = ["NRSSMOS-MOD-G1H-02-5344031756_1_491_SE_2015-12-10T03h25m56.fits",
                                "NRSSMOS-MOD-G1H-02-5344031756_1_492_SE_2015-12-10T03h25m56.fits",
                                "NRSSMOS-MOD-G2M-01-5344191938_1_491_SE_2015-12-10T19h29m26.fits",
                                "NRSSMOS-MOD-G3H-02-5344120942_1_491_SE_2015-12-10T12h18m25.fits",
                                "NRSSMOS-MOD-G3H-02-5344120942_1_492_SE_2015-12-10T12h18m25.fits"]
        if raw_data_root_file in special_cutout_files:
            nid = NIDs[special_cutout_files.index(raw_data_root_file)]
        else:
            nid = None
        esafile = auxfunc.get_esafile(esa_files_path, raw_data_root_file, "FS", specifics, nid=nid)

        if esafile == "ESA file not found":
            print(" * validate_wcs_extract2d is exiting because the corresponding ESA file was not found.")
            print("   -> The extract_2d test is now set to skip. ")
            result = "skip"
            return result

        if not isinstance(esafile, list):
            esafile_list = [esafile]
        else:
            esafile_list = esafile

        slit = sltname.replace("S", "")
        if "A" in slit:
            slit = "_"+slit.split("A")[0]+"_"
        if "B" in slit:
            slit = "_"+slit.split("B")[0]+"_"

        # choose corresponding esa file
        for esafile in esafile_list:
            if slit in esafile:
                print ("Using this ESA file: \n", esafile)
                esahdulist = fits.open(esafile)
                break

        # Find corners from ESA file
        esahdr1 = esahdulist["DATA1"].header
        enext = []
        for ext in esahdulist:
            enext.append(ext)
        if detector == "NRS1":
            eflux = fits.getdata(esafile, "DATA1")
        if detector == "NRS2":
            try:
                eflux = fits.getdata(esafile, "DATA2")
            except:
                IndexError
                print(" * Exiting extract_2d test because there are no extensions that match detector NRS2 in the ESA file.")
                print("   -> The extract_2d test is now set to skip. ")
                result = "skip"
                continue
        esahdulist.close()

        # get the origin of the subwindow
        ney, nex = eflux.shape
        if detector == "NRS1":
            ex0 = esahdr1["CRVAL1"] - esahdr1["CRPIX1"] + 1
            ey0 = esahdr1["CRVAL2"] - esahdr1["CRPIX2"] + 1
        else:
            ex0 = 2048.0 - esahdr1["CRPIX1"] + esahdr1["CRVAL1"]
            ey0 = 2048.0 - esahdr1["CRVAL2"] + esahdr1["CRPIX2"]

        ex1 = ex0 + nex
        ey1 = ey0 + ney

        ecorners = {(ex0, ey0), (ex1, ey0), (ex1, ey1), (ex0, ey1)}

        result[sltname] = ecorners == icorners

    return result


def find_MOSwindowcorners(infile_name, msa_conf_name, esa_files_path):
    """
	Find the slitlet corners of pipeline and ESA files and determine if they match.
	Args:
		infile_name - string, name of the input pipeline fits file
		msa_conf_name - string, path and name of the MSA configuration file
		esa_files_path - string, path to locate esa files
	Retrurs:
		result - string or dictionary
				if string, test is set to skip
				if dictionary, contains a boolean per slit test
	"""

    result = OrderedDict

    # Grab initial metadata
    primary_header = fits.getheader(infile_name, ext=0)
    msametfl = primary_header["MSAMETFL"]
    detector = primary_header["DETECTOR"]

    # check that shutter configuration file in header is the same as given in PTT_config file
    if msametfl != os.path.basename(msa_conf_name):
        print ("* WARNING! MSA config file name given in PTT_config file does not match the MSAMETFL keyword in main header.\n")

    # Get shutter info from metadata
    shutter_info = fits.getdata(msa_conf_name, extname="SHUTTER_INFO") # this is generally ext=2
    slit_ids = shutter_info.field("slitlet_id")
    quad = shutter_info.field("shutter_quadrant")
    row = shutter_info.field("shutter_row")
    col = shutter_info.field("shutter_column")

    # copy the MSA shutter configuration file into the pytest directory
    try:
        subprocess.run(["cp", msa_conf_name, "."])
    except FileNotFoundError:
        print(" * PTT is not able to locat the MSA shutter configuration file. Please make sure that the msa_conf_name variable in")
        print("   the PTT_config.cfg file is pointing exactly to where the fits file exists. ")
        print("   -> The extract_2d test is now set to skip. ")
        result = "skip"
        return result

    # Identify the science extensions
    sci_ext_list = auxfunc.get_sci_extensions(infile_name)

    # Iterate over sci extensions
    for i, s_ext in enumerate(sci_ext_list):
        sci_header = fits.getheader(infile_name, ext=s_ext)
        name = sci_header['SLTNAME']
        print ("\nWorking with slit: ", name)

        slitlet_idx = slit_ids.tolist().index(int(name))

        # grab corners of extracted subwindow
        px0 = sci_header["SLTSTRT1"] + primary_header["SUBSTRT1"] - 1
        py0 = sci_header["SLTSTRT2"] + primary_header["SUBSTRT2"] - 1
        pnx = sci_header["SLTSIZE1"]
        pny = sci_header["SLTSIZE2"]

        px1 = px0 + pnx
        py1 = py0 + pny

        icorners = {(px0, py0), (px1, py0), (px1, py1), (px1, py0)}

        # Identify the associated ESA file
        _, raw_data_root_file = auxfunc.get_modeused_and_rawdatrt_PTT_cfg_file()
        print("Using this raw data file to find the corresponding ESA file: ", raw_data_root_file)
        q, r, c = quad[slitlet_idx], row[slitlet_idx], col[slitlet_idx]
        print("Pipeline shutter info:   quadrant=", q, "   row=", r, "   col=", c)
        specifics = [q, r, c]
        esafile = auxfunc.get_esafile(esa_files_path, raw_data_root_file, "MOS", specifics)

        # skip the test if the esafile was not found
        if esafile == "ESA file not found":
            result[name] = "skip"

        # Open esafile and grab subarray coordinates
        with fits.open(esafile) as esahdulist:
            esahdr1 = esahdulist["DATA1"].header
            esa_shutter_i = esahdulist[0].header['SHUTTERI']
            esa_shutter_j = esahdulist[0].header['SHUTTERJ']
            esa_quadrant = esahdulist[0].header['QUADRANT']
            # first check if ESA shutter info is the same as pipeline
            print("For slitlet", name)
            if q == esa_quadrant:
                print("\n -> Same quadrant for pipeline and ESA data: ", q)
            else:
                print("\n -> Missmatch of quadrant for pipeline and ESA data: ", q, esa_quadrant)
            if r == esa_shutter_i:
                    print("\n -> Same row for pipeline and ESA data: ", r)
            else:
                print("\n -> Missmatch of row for pipeline and ESA data: ", r, esa_shutter_i)
            if c == esa_shutter_j:
                print("\n -> Same column for pipeline and ESA data: ", c, "\n")
            else:
                print("\n -> Missmatch of column for pipeline and ESA data: ", c, esa_shutter_j, "\n")

            if detector == "NRS1":
                esahdr = esahdulist['DATA1'].header
                ney, nex = esahdulist['DATA1'].data.shape
                ex0 = esahdr["CRVAL1"] - esahdr["CRPIX1"] + 1
                ey0 = esahdr["CRVAL2"] - esahdr["CRPIX2"] + 1
            else:
                esahdr = esahdulist['DATA2'].header
                ney, nex = esahdulist['DATA2'].data.shape
                ex0 = 2048.0 - esahdr1["CRPIX1"] + esahdr1["CRVAL1"]
                ey0 = 2048.0 - esahdr1["CRVAL2"] + esahdr1["CRPIX2"]

            ex1 = ex0 + nex
            ey1 = ey0 + ney

            ecorners = {(ex0, ey0), (ex1, ey0), (ex1, ey1), (ex0, ey1)}

            result[name] = ecorners == icorners

    return result

