import subprocess
import os
import math
import argparse
import sys
from collections import OrderedDict
from astropy.io import fits

from gwcs import wcstools
from jwst import datamodels

from ..auxiliary_code import auxiliary_functions as auxfunc


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Jul 2020 - Version 1.0: initial version completed


"""
Check that the corners after the extract_2d step match the truth (benchmark) data.
"""


def find_FSwindowcorners(infile_name, truth_file=None, esa_files_path=None, extract_2d_threshold_diff=4):
    """
    Find the slits corners of pipeline and truth file and determine if they match.
    Args:
        infile_name - string, name of the input pipeline fits file
        truth_file - string or None, name of truth (or benchmark) file to compare to
        esa_files_path - string or None, path to locate esa files
        extract_2d_threshold_diff - integer, maximum allowed difference tolerance in pixels

    Returns:
        result - string or dictionary
                if string, test is set to skip
                if dictionary, contains a boolean per slit test
    """

    result = OrderedDict()
    large_diff_corners_dict = OrderedDict()
    log_msgs = []

    if esa_files_path is None:
        # get the model from the "truth" (or comparison) file
        truth_hdul = fits.open(truth_file)
        print("Information from the 'truth' (or comparison) file ")
        print(truth_hdul.info())
        truth_hdul.close()
        truth_img = datamodels.ImageModel(truth_file)

    # get grating and filter info from the pipeline file or datamodel
    if isinstance(infile_name, str):
        # get the datamodel from the assign_wcs output file
        img = datamodels.ImageModel(infile_name)
        msg = 'infile_name='+infile_name
        print(msg)
        log_msgs.append(msg)
    else:
        img = infile_name

    detector = img.meta.instrument.detector
    lamp = img.meta.instrument.lamp_state
    grat = img.meta.instrument.grating
    filt = img.meta.instrument.filter
    msg = "from assign_wcs file  -->     Detector: " + detector + "   Grating: " + grat + "   Filter: " + \
          filt + "   Lamp: " + lamp
    print(msg)
    log_msgs.append(msg)

    # To get the open and projected on the detector slits of the pipeline processed file
    for slit in img.slits:
        pipeslit = slit.name
        msg = "\nWorking with slit: "+pipeslit
        print(msg)
        log_msgs.append(msg)

        # grab corners of extracted sub-window
        px0 = slit.xstart + img.meta.subarray.xstart - 1
        py0 = slit.ystart + img.meta.subarray.ystart - 1
        pnx = slit.xsize
        pny = slit.ysize

        px1 = px0 + pnx
        py1 = py0 + pny

        pipeline_corners = [px0, py0, px1, py1]

        compare_to_esa_data = False
        if esa_files_path is not None:
            compare_to_esa_data = True

            # Find esafile (most of this copy-pasted from compare_wcs_fs)
            _, raw_data_root_file = auxfunc.get_modeused_and_rawdatrt_PTT_cfg_file(infile_name)
            specifics = [pipeslit]
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
            # print('got the ESA file:', esafile)

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
                if isinstance(esafile, tuple):
                    esafile_list = []
                    for ef in esafile:
                        if ef:
                            esafile_list.append(ef)
                else:
                    esafile_list = [esafile]
            else:
                esafile_list = esafile

            slit = pipeslit.replace("S", "")
            if "A" in slit:
                slit = "_"+slit.split("A")[0]+"_"
            if "B" in slit:
                slit = "_"+slit.split("B")[0]+"_"

            # choose corresponding esa file
            for esafile in esafile_list:
                if "not found" in esafile:
                    msg = " * No ESA file was found. Test will be skipped because there is no file to compare with."
                    print(msg)
                    log_msgs.append(msg)
                    return "skip", log_msgs

                if slit in esafile:
                    print("Using this ESA file: \n", esafile)
                    with fits.open(esafile) as esahdulist:
                        # Find corners from ESA file
                        # print(esahdulist.info())
                        if detector == "NRS1" or "NRS1" in infile_name  or  "491" in infile_name:
                            dat = "DATA1"
                        else:
                            dat = "DATA2"
                        try:
                            esahdr = esahdulist[dat].header
                            print("For detector=", detector, " reading data of extension name=", dat)
                        except KeyError:
                            print("This file contains no information for extension name=", dat)
                            continue

                        enext = []
                        for ext in esahdulist:
                            enext.append(ext)
                        if detector == "NRS1":
                            eflux = fits.getdata(esafile, "DATA1")
                        if detector == "NRS2":
                            try:
                                eflux = fits.getdata(esafile, "DATA2")
                            except IndexError:
                                msg1 = " * Exiting extract_2d test because there are no extensions that match " \
                                       "detector NRS2 in the ESA file."
                                msg2 = "   -> The extract_2d test is now set to skip. "
                                print(msg1)
                                print(msg2)
                                log_msgs.append(msg1)
                                log_msgs.append(msg2)
                                result = "skip"
                                continue

                        # get the origin of the subwindow
                        ney, nex = eflux.shape
                        truth_slit_size = [nex, ney]
                        if detector == "NRS1":
                            ex0 = int(esahdr["CRVAL1"]) - int(esahdr["CRPIX1"]) + 1
                            ey0 = int(esahdr["CRVAL2"]) - int(esahdr["CRPIX2"]) + 1
                        else:
                            ex0 = 2049 - int(esahdr["CRPIX1"]) - int(esahdr["CRVAL1"])
                            ey0 = 2049 - int(esahdr["CRVAL2"]) - int(esahdr["CRPIX2"])

                        ex1 = ex0 + nex
                        ey1 = ey0 + ney

                        truth_corners = [ex0, ey0, ex1, ey1]

        # In case we are NOT comparing to ESA data
        if not compare_to_esa_data:
            # determine if the current open slit is also open in the truth file
            continue_with_wcs_test = False
            for truth_slit in truth_img.slits:
                if truth_slit.name == pipeslit:
                    continue_with_wcs_test = True
                    break

            if not continue_with_wcs_test:
                msg1 = "\n * Script check_corners_extract_2d.py is exiting because open slit " + pipeslit + \
                       " is not open in truth file."
                msg2 = "   -> The extract_2d test is now set to FAILED. \n"
                print(msg1)
                print(msg2)
                log_msgs.append(msg1)
                log_msgs.append(msg2)
                FINAL_TEST_RESULT = "FAILED"
                return FINAL_TEST_RESULT, log_msgs

            else:
                # grab corners of extracted sub-window from the truth file
                tx0 = truth_slit.xstart + truth_img.meta.subarray.xstart - 1
                ty0 = truth_slit.ystart + truth_img.meta.subarray.ystart - 1
                tnx = truth_slit.xsize
                tny = truth_slit.ysize

                tx1 = tx0 + tnx
                ty1 = ty0 + tny

                truth_corners = [tx0, ty0, tx1, ty1]
                truth_slit_size = [tnx, tny]

        # pass/fail criterion: if truth corners match pipeline corners then True, else False
        result[pipeslit], msgs = truth_corners_in_pipeline_corners(truth_corners, pipeline_corners,
                                                                   extract_2d_threshold_diff=extract_2d_threshold_diff)
        for msg in msgs:
            log_msgs.append(msg)

        print('    Truth slit size = ', truth_slit_size)
        print(' Pipeline slit size = ', pnx, pny)

        msg1 = "Corners for slit " + pipeslit + ":  [x0, y0, x1, y1]"
        msg2 = "   Truth corners: " + repr(truth_corners)
        msg3 = "    Pipeline corners: " + repr(pipeline_corners)
        print(msg1)
        print(msg2)
        print(msg3)
        log_msgs.append(msg1)
        log_msgs.append(msg2)
        log_msgs.append(msg3)
        if result[pipeslit]:
            msg = "* Test PASSED: All corners match within the threshold."
            print(msg)
            log_msgs.append(msg)
        else:
            msg = "* Test FAILED: One or more corners have a difference larger than threshold."
            print(msg)
            log_msgs.append(msg)
            # Record all the corners that have points larger than threshold
            large_diff_corners_dict[pipeslit] = [truth_corners, pipeline_corners]

    # If all tests passed then test will be marked as PASSED, else it will be FAILED
    print('\nSummary of test results: \n', result)
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
        result_msg = "All slits PASSED extract_2d test."
        print(result_msg)
        log_msgs.append(msg)
    else:
        msg1 = '\n\n *** These slits have corners with differences larger than threshold of ' + \
               repr(extract_2d_threshold_diff)+' pixels: '
        msg2 = '   ESA corners                   Pipeline corners'
        log_msgs.append(msg1)
        log_msgs.append(msg2)
        print(msg1)
        print(msg2)
        for sltname, corners in large_diff_corners_dict.items():
            msg1 = 'slit '+sltname
            msg2 = ''+ repr(corners)
            log_msgs.append(msg1)
            log_msgs.append(msg2)
            print(msg1)
            print(msg2)
        msg = "\n *** Final result for extract_2d test will be reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)

    return FINAL_TEST_RESULT, log_msgs


def truth_corners_in_pipeline_corners(truth_corners, pipeline_corners, extract_2d_threshold_diff):
    """
    This function tests if the corners of the truth file are present in the corners of the pipeline.
    Args:
        truth_corners: list, integers for x0, y0, x1, y1
        pipeline_corners: list, integers for x0, y0, x1, y1
        extract_2d_threshold_diff: integer, maximum absolute tolerance in pixels

    Returns:
        results: list of booleans
    """
    result = False
    msgs = []
    for bi, pi in zip(truth_corners, pipeline_corners):
        # absolute difference <= threshold, then pass
        diff_within_threshold = math.isclose(bi, pi, abs_tol=extract_2d_threshold_diff)
        if diff_within_threshold:
            result = True
        else:
            # if difference of pipeline - truth is positive, then pass because pipeline extraction
            # includes truth extraction
            diff_within_threshold = pi - bi
            if diff_within_threshold > 0.0:
                msg = "* WARNING: Difference of pipeline - truth extraction window is larger than threshold BUT " \
                      "the truth extraction is fully within the pipeline extraction window."
                msgs.append(msg)
                result = True
            else:
                result = False
    return result, msgs


def find_MOSwindowcorners(infile_name, msa_conf_name, truth_file=None, esa_files_path=None,
                          extract_2d_threshold_diff=4):
    """
    Find the slitlet corners of pipeline and ESA files and determine if they match.
    Args:
        infile_name - string, name of the input pipeline fits file
        msa_conf_name - string, path and name of the MSA configuration file
        truth_file - string or None, name of truth (or benchmark) file to compare to
        esa_files_path - string or None, path to locate esa files
        extract_2d_threshold_diff - integer, maximum allowed difference tolerance in pixels
    Retrurs:
        result - string or dictionary
                if string, test is set to skip
                if dictionary, contains a boolean per slit test
        log_msgs - list of print statements to in the log file
    """

    result = OrderedDict()
    large_diff_corners_dict = OrderedDict()
    log_msgs = []

    if esa_files_path is None:
        # get the model from the "truth" (or comparison) file
        truth_hdul = fits.open(truth_file)
        print("Information from the 'truth' (or comparison) file ")
        print(truth_hdul.info())
        truth_hdul.close()
        truth_img = datamodels.ImageModel(truth_file)

    # get grating and filter info from the pipeline file or datamodel
    if isinstance(infile_name, str):
        # get the datamodel from the assign_wcs output file
        img = datamodels.ImageModel(infile_name)
        msg = 'infile_name='+infile_name
        print(msg)
        log_msgs.append(msg)
    else:
        img = infile_name

    detector = img.meta.instrument.detector
    lamp = img.meta.instrument.lamp_state
    grat = img.meta.instrument.grating
    filt = img.meta.instrument.filter
    msg = "from assign_wcs file  -->     Detector: " + detector + "   Grating: " + grat + "   Filter: " + \
          filt + "   Lamp: " + lamp
    print(msg)
    log_msgs.append(msg)

    # get name of shutter configuration file
    msametfl = img.meta.instrument.msa_metadata_file

    # check that shutter configuration file in header is the same as given in PTT_config file
    if msametfl != os.path.basename(msa_conf_name):
        msg = " *** WARNING! MSA config file name given in PTT_config file does not match the MSAMETFL " \
              "keyword in main header.\n"
        print(msg)
        log_msgs.append(msg)

    # Get shutter info from shutter configuration file
    shutter_info = fits.getdata(msa_conf_name, extname="SHUTTER_INFO")
    slit_ids = shutter_info.field("slitlet_id")
    quad = shutter_info.field("shutter_quadrant")
    row = shutter_info.field("shutter_row")
    col = shutter_info.field("shutter_column")

    # Iterate over the slits of the pipeline processed file
    for slit in img.slits:
        pipeslit = slit.name
        msg = "\nWorking with slitlet: "+pipeslit
        print(msg)
        log_msgs.append(msg)

        slitlet_idx = slit_ids.tolist().index(int(pipeslit))

        # grab corners of extracted sub-window
        px0 = slit.xstart + img.meta.subarray.xstart - 1
        py0 = slit.ystart + img.meta.subarray.ystart - 1
        pnx = slit.xsize
        pny = slit.ysize

        px1 = px0 + pnx
        py1 = py0 + pny

        pipeline_corners = [px0, py0, px1, py1]

        compare_to_esa_data = False
        if esa_files_path is not None:
            compare_to_esa_data = True

            # Identify the associated ESA file
            _, raw_data_root_file = auxfunc.get_modeused_and_rawdatrt_PTT_cfg_file(infile_name)
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
                result[pipeslit] = "skip"
                msg = "ESA file not found, skipping test for this shutter."
                print(msg)
                log_msgs.append(msg)
                continue

            # Open esafile and grab subarray coordinates
            with fits.open(esafile) as esahdulist:
                if "NRS1" in detector or "491" in detector:
                    dat = "DATA1"
                else:
                    dat = "DATA2"
                try:
                    esahdr = esahdulist[dat].header
                    print("For detector=", detector, " reading data of extension name=", dat)
                except KeyError:
                    print("This file contains no information for extension name=", dat)
                    continue
                esa_shutter_i = esahdulist[0].header['SHUTTERI']
                esa_shutter_j = esahdulist[0].header['SHUTTERJ']
                esa_quadrant = esahdulist[0].header['QUADRANT']
                # first check if ESA shutter info is the same as pipeline
                msg = "For slitlet" + pipeslit
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
                    ex0 = int(esahdr["CRVAL1"]) - int(esahdr["CRPIX1"]) + 1
                    ey0 = int(esahdr["CRVAL2"]) - int(esahdr["CRPIX2"]) + 1
                else:
                    ney, nex = esahdulist['DATA2'].data.shape
                    ex0 = 2049 - int(esahdr["CRPIX1"]) - int(esahdr["CRVAL1"])
                    ey0 = 2049 - int(esahdr["CRVAL2"]) - int(esahdr["CRPIX2"])

                ex1 = ex0 + nex
                ey1 = ey0 + ney

                truth_corners = [ex0, ey0, ex1, ey1]
                truth_slit_size = [nex, ney]

        # In case we are NOT comparing to ESA data
        if not compare_to_esa_data:
            # determine if the current open slit is also open in the truth file
            continue_with_wcs_test = False
            for truth_slit in truth_img.slits:
                if truth_slit.name == pipeslit:
                    continue_with_wcs_test = True
                    break

            if not continue_with_wcs_test:
                msg1 = "\n * Script check_corners_extract_2d.py is exiting because open slit " + pipeslit + \
                       " is not open in truth file."
                msg2 = "   -> The extract_2d test is now set to FAILED. \n"
                print(msg1)
                print(msg2)
                log_msgs.append(msg1)
                log_msgs.append(msg2)
                FINAL_TEST_RESULT = "FAILED"
                return FINAL_TEST_RESULT, log_msgs

            else:
                # grab corners of extracted sub-window from the truth file
                tx0 = truth_slit.xstart + truth_img.meta.subarray.xstart - 1
                ty0 = truth_slit.ystart + truth_img.meta.subarray.ystart - 1
                tnx = truth_slit.xsize
                tny = truth_slit.ysize

                tx1 = tx0 + tnx
                ty1 = ty0 + tny

                truth_corners = [tx0, ty0, tx1, ty1]
                truth_slit_size = [tnx, tny]

        result[pipeslit], msgs = truth_corners_in_pipeline_corners(truth_corners, pipeline_corners,
                                                                   extract_2d_threshold_diff=extract_2d_threshold_diff)
        for msg in msgs:
            log_msgs.append(msg)

        print('Truth slit size    = ', truth_slit_size)
        print('Pipeline slit size = ', pnx, pny)

        msg1 = "\n Corners for slitlet " + pipeslit + ":  [x0, y0, x1, y1]"
        msg2 = "       Truth corners: " + repr(truth_corners)
        msg3 = "    Pipeline corners: " + repr(pipeline_corners)
        print(msg1)
        print(msg2)
        print(msg3)
        log_msgs.append(msg1)
        log_msgs.append(msg2)
        log_msgs.append(msg3)
        if result[pipeslit]:
            msg = "* Pytest PASSED: All corners match within the threshold."
            print(msg)
            log_msgs.append(msg)
        else:
            msg = "* Pytest FAILED: One or more corners have a difference larger than threshold."
            print(msg)
            log_msgs.append(msg)
            # Record all the corners that have points larger than threshold
            large_diff_corners_dict[pipeslit] = [truth_corners, pipeline_corners]

    # If all tests passed then the test will be marked as PASSED, else it will be FAILED
    FINAL_TEST_RESULT = False
    for t, v in result.items():
        if not v:
            FINAL_TEST_RESULT = False
            break
        else:
            FINAL_TEST_RESULT = True

    if FINAL_TEST_RESULT:
        msg = "\n\n *** Final result for extract_2d test will be reported as PASSED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "All slitlets PASSED extract_2d test."
        print(result_msg)
        log_msgs.append(msg)
    else:
        msg1 = '\n\n *** These slitlets have corners with differences larger than threshold of ' + \
               repr(extract_2d_threshold_diff)+' pixels: '
        msg2 = '     Truth corners            Pipeline corners'
        log_msgs.append(msg1)
        log_msgs.append(msg2)
        print(msg1)
        print(msg2)
        for sltname, corners in large_diff_corners_dict.items():
            msg1 = 'slitlet '+sltname
            msg2 = ''+repr(corners)
            log_msgs.append(msg1)
            log_msgs.append(msg2)
            print(msg1)
            print(msg2)
        msg = "\n *** Final result for extract_2d test will be reported as FAILED *** \n"
        print(msg)
        log_msgs.append(msg)
        result_msg = "One or more slitlets FAILED extract_2d test."
        print(result_msg)
        log_msgs.append(msg)

    return FINAL_TEST_RESULT, log_msgs


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("infile_name",
                        action='store',
                        default=None,
                        help='Name of input fits file prior to assign_wcs step, i.e. blah_rate.fits')
    parser.add_argument("-b",
                        dest="truth_files_path",
                        action='store',
                        default=None,
                        help='Use the -b flag to provide the path were to locate the "truth" (or benchmark) '
                             'file to compare to')
    parser.add_argument("-e",
                        dest="esa_files_path",
                        action='store',
                        default=None,
                        help='Use flag -e to provide path were to locate the new benchmark files for comparison'
                             ' (to create new truth files)')
    parser.add_argument("-m",
                        dest="msa_conf_name",
                        action='store',
                        default=None,
                        help='Use the -m flag to provide name of the MSA shutter configuration file in '
                             'pipeline format, e.g. blah_msa.fits')
    parser.add_argument("-t",
                        dest="extract_2d_threshold_diff",
                        action='store',
                        default=None,
                        help='Use flag -t to provide a new threshold to use. Default is 4 pixels.')
    args = parser.parse_args()

    # Set the variables input from the command line
    infile_name = args.infile_name
    truth_file_path = args.truth_file_path
    esa_files_path = args.esa_files_path
    msa_conf_name = args.msa_conf_name
    extract_2d_threshold_diff = args.extract_2d_threshold_diff

    # choose to run FS or MOS check
    if msa_conf_name is not None:
        mos = True
    else:
        mos = False

    # convert threshold to integer if provided a new value
    if extract_2d_threshold_diff is not None:
        extract_2d_threshold_diff = int(extract_2d_threshold_diff)

    # Run the principal function of the script
    if mos:
        find_MOSwindowcorners(infile_name, msa_conf_name, truth_file=truth_file_path, esa_files_path=esa_files_path,
                              extract_2d_threshold_diff=extract_2d_threshold_diff)
    else:
        find_FSwindowcorners(infile_name, truth_file=truth_file_path, esa_files_path=esa_files_path,
                             extract_2d_threshold_diff=extract_2d_threshold_diff)


if __name__ == '__main__':
    sys.exit(main())


