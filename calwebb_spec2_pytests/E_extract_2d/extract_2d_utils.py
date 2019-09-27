import subprocess
import os
from collections import OrderedDict
from astropy.io import fits
from ..auxiliary_code import auxiliary_functions as auxfunc


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.2"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# Jan 2019 - Version 1.1: Maria modified and added Gray's code for validation tests
# Apr 2019 - Version 1.2: implemented logging capability


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
    log_msgs = []

    #iterate over slits
    sci_dict = auxfunc.get_sci_extensions(infile_name)

    primary_header = fits.getheader(infile_name, ext=0)
    detector = primary_header["DETECTOR"]

    for i, s_ext in enumerate(sci_dict):
        s_ext_number = sci_dict[s_ext]
        sci_header = fits.getheader(infile_name, ext=s_ext_number)

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
            msg1 = " * validate_wcs_extract2d is exiting because the corresponding ESA file was not found."
            msg2 = "   -> The extract_2d test is now set to skip. "
            print(msg1)
            print(msg2)
            log_msgs.append(msg1)
            log_msgs.append(msg2)
            result = "skip"
            return result, log_msgs

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
                with fits.open(esafile) as esahdulist:
                    # Find corners from ESA file
                    if "NRS1" in esafile or "491" in esafile:
                        dat = "DATA1"
                    else:
                        dat = "DATA2"
                    esahdr1 = esahdulist[dat].header
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
                            msg1 = " * Exiting extract_2d test because there are no extensions that match detector NRS2 in the ESA file."
                            msg2 = "   -> The extract_2d test is now set to skip. "
                            print(msg1)
                            print(msg2)
                            log_msgs.append(msg1)
                            log_msgs.append(msg2)
                            result = "skip"
                            continue

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

                    # Pytest pass/fail criterion: if esa corners match pipeline corners then True, else False
                    result[sltname] = ecorners == icorners

                    msg1 = "* Corners for slit "+sltname+":"
                    msg2 = "         ESA corners: "+repr(ecorners)
                    msg3 = "    Pipeline corners: "+repr(icorners)
                    print(msg1)
                    print(msg2)
                    print(msg3)
                    log_msgs.append(msg1)
                    log_msgs.append(msg2)
                    log_msgs.append(msg3)
                    if result[sltname]:
                        msg = "* Pytest PASSED "
                        print(msg)
                        log_msgs.append(msg)
                    else:
                        msg = "* Pytest FAILED "
                        print(msg)
                        log_msgs.append(msg)

    # If all tests passed then pytest will be marked as PASSED, else it will be FAILED
    for t, v in result.items():
        if not v:
            FINAL_TEST_RESULT = False
            break
        else:
            FINAL_TEST_RESULT = True

    if FINAL_TEST_RESULT:
        msg = "\n *** Final result for extract_2d test will be reported as PASSED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "All slits PASSED extract_2d test."
        print(result_msg)
        log_msgs.append(msg)
    else:
        msg = "\n *** Final result for extract_2d test will be reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "One or more slits FAILED extract_2d test."
        print(result_msg)
        log_msgs.append(msg)

    return FINAL_TEST_RESULT, log_msgs


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
		log_msgs - list of print statements to in the log file
	"""

    result = OrderedDict()
    log_msgs = []

    # Grab initial metadata
    primary_header = fits.getheader(infile_name, ext=0)
    msametfl = primary_header["MSAMETFL"]
    detector = primary_header["DETECTOR"]

    # check that shutter configuration file in header is the same as given in PTT_config file
    if msametfl != os.path.basename(msa_conf_name):
        msg = "* WARNING! MSA config file name given in PTT_config file does not match the MSAMETFL keyword in main header.\n"
        print(msg)
        log_msgs.append(msg)

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
        msg1 = " * PTT is not able to locat the MSA shutter configuration file. Please make sure that the msa_conf_name variable in"
        msg2 = "   the PTT_config.cfg file is pointing exactly to where the fits file exists. "
        msg3 = "   -> The extract_2d test is now set to skip. "
        print(msg1)
        print(msg2)
        print(msg3)
        log_msgs.append(msg1)
        log_msgs.append(msg2)
        log_msgs.append(msg3)
        result = "skip"
        return result, log_msgs

    # Identify the science extensions
    sci_ext_dict = auxfunc.get_sci_extensions(infile_name)

    # Iterate over sci extensions
    for i, s_ext in enumerate(sci_ext_dict):
        s_ext_number = sci_ext_dict[s_ext]
        sci_header = fits.getheader(infile_name, ext=s_ext_number)
        name = sci_header['SLTNAME']
        msg = "\nWorking with slit: "+name
        print(msg)
        log_msgs.append(msg)

        slitlet_idx = slit_ids.tolist().index(int(name))

        # grab corners of extracted subwindow
        px0 = sci_header["SLTSTRT1"] + primary_header["SUBSTRT1"] - 1
        py0 = sci_header["SLTSTRT2"] + primary_header["SUBSTRT2"] - 1

        pnx = sci_header["SLTSIZE1"]
        pny = sci_header["SLTSIZE2"]

        px1 = px0 + pnx
        py1 = py0 + pny

        icorners = {(px0, py0), (px1, py0), (px1, py1), (px0, py1)}

        # Identify the associated ESA file
        _, raw_data_root_file = auxfunc.get_modeused_and_rawdatrt_PTT_cfg_file()
        msg = "Using this raw data file to find the corresponding ESA file: "+raw_data_root_file
        print(msg)
        log_msgs.append(msg)
        q, r, c = quad[slitlet_idx], row[slitlet_idx], col[slitlet_idx]
        msg = "Pipeline shutter info:   quadrant="+str(q)+"   row="+str(r)+"   col="+str(c)
        print(msg)
        log_msgs.append(msg)
        specifics = [q, r, c]
        esafile = auxfunc.get_esafile(esa_files_path, raw_data_root_file, "MOS", specifics)
        if len(esafile) == 2:
            if len(esafile[-1]) == 0:
                esafile = esafile[0]
        msg = "Using this ESA file: \n"+esafile
        print(msg)
        log_msgs.append(msg)

        # skip the test if the esafile was not found
        if "ESA file not found" in esafile:
            result[name] = "skip"
            msg = "ESA file not found, skipping test for this shutter."
            print(msg)
            log_msgs.append(msg)
            continue

        # Open esafile and grab subarray coordinates
        with fits.open(esafile) as esahdulist:
            if "NRS1" in esafile  or  "491" in esafile:
                dat = "DATA1"
            else:
                dat = "DATA2"
            esahdr = esahdulist[dat].header
            esa_shutter_i = esahdulist[0].header['SHUTTERI']
            esa_shutter_j = esahdulist[0].header['SHUTTERJ']
            esa_quadrant = esahdulist[0].header['QUADRANT']
            # first check if ESA shutter info is the same as pipeline
            msg = "For slitlet"+name
            print(msg)
            log_msgs.append(msg)
            if q == esa_quadrant:
                msg = "\n -> Same quadrant for pipeline and ESA data: "+str(q)
                print(msg)
                log_msgs.append(msg)
            else:
                msg = "\n -> Missmatch of quadrant for pipeline and ESA data: "+str(q)+esa_quadrant
                print(msg)
                log_msgs.append(msg)
            if r == esa_shutter_i:
                msg = "\n -> Same row for pipeline and ESA data: "+str(r)
                print(msg)
                log_msgs.append(msg)
            else:
                msg = "\n -> Missmatch of row for pipeline and ESA data: "+str(r)+esa_shutter_i
                print(msg)
                log_msgs.append(msg)
            if c == esa_shutter_j:
                msg = "\n -> Same column for pipeline and ESA data: "+str(c)+"\n"
                print(msg)
                log_msgs.append(msg)
            else:
                msg = "\n -> Missmatch of column for pipeline and ESA data: "+str(c)+esa_shutter_j+"\n"
                print(msg)
                log_msgs.append(msg)

            if detector == "NRS1":
                ney, nex = esahdulist['DATA1'].data.shape
                ex0 = int(esahdr["CRVAL1"] - esahdr["CRPIX1"] + 1)
                ey0 = int(esahdr["CRVAL2"] - esahdr["CRPIX2"] + 1)
            else:
                ney, nex = esahdulist['DATA2'].data.shape
                ex0 = int(2048 - esahdr["CRPIX1"] + esahdr["CRVAL1"])
                ey0 = int(2048 - esahdr["CRVAL2"] + esahdr["CRPIX2"])

            ex1 = ex0 + nex
            ey1 = ey0 + ney

            ecorners = {(ex0, ey0), (ex1, ey0), (ex1, ey1), (ex0, ey1)}

            result[name] = ecorners == icorners

            print('ESA slit size      = ', nex, ney)
            print('Pipeline slit size = ', pnx, pny)

            msg1 = "\n* Corners for slitlet "+name+":"
            msg2 = "         ESA corners: "+repr(ecorners)
            msg3 = "    Pipeline corners: "+repr(icorners)
            print(msg1)
            print(msg2)
            print(msg3)
            log_msgs.append(msg1)
            log_msgs.append(msg2)
            log_msgs.append(msg3)
            if result[name]:
                msg = "* Pytest PASSED "
                print(msg)
                log_msgs.append(msg)
            else:
                msg = "* Pytest FAILED "
                print(msg)
                log_msgs.append(msg)

    # If all tests passed then pytest will be marked as PASSED, else it will be FAILED
    FINAL_TEST_RESULT = False
    for t, v in result.items():
        if not v:
            FINAL_TEST_RESULT = False
            break
        else:
            FINAL_TEST_RESULT = True

    if FINAL_TEST_RESULT:
        msg = "\n *** Final result for extract_2d test will be reported as PASSED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "All slitlets PASSED extract_2d test."
        print(result_msg)
        log_msgs.append(msg)
    else:
        msg = "\n *** Final result for extract_2d test will be reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "One or more slitlets FAILED extract_2d test."
        print(result_msg)
        log_msgs.append(msg)

    return FINAL_TEST_RESULT, log_msgs

