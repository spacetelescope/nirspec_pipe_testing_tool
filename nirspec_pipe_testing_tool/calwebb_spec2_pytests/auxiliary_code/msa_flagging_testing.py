import numpy as np
import os
from astropy.io import fits
import json
import datetime
from astropy.visualization import (ImageNormalize, AsinhStretch)
import matplotlib.pyplot as plt
import time
import argparse
import sys

from jwst.assign_wcs import nirspec
from jwst import datamodels
from jwst.assign_wcs.assign_wcs_step import AssignWcsStep

"""
This script tests the pipeline msa flagging step output for MOS data. It is the scripted and generalized version of the 
jupyter notebook test_MSA_flagging.ipynb written by James Muzerolle in July of 2020.
"""

# HEADER
__author__ = "M. A. Pena-Guerrero & J. Muzerolle"
__version__ = "1.0"

# HISTORY
# Jul 2020 - Version 1.0: initial version completed


def create_metafile_fopens(outfile, allcols, allrows, quads, stellarity, failedopens,
                           save_fig=False, show_fig=False, debug=False):
    """
    Create the metafile for the failed open shutters according to msa failed open - operability - reference file.
    :param outfile: string, new fits msa meta file
    :param allcols: array
    :param allrows: array
    :param quads: array
    :param stellarity: array
    :param failedopens: array
    :param save_fig: boolean
    :param show_fig: boolean
    :param debug: boolean
    :return: nothing
    """
    now = datetime.datetime.now()
    # set up metafile table structure
    # SHUTTER_INFO table
    # slitlet IDs - use one per shutter
    slitlets = np.arange(0, len(allcols)) + 1
    # source IDs - arbitrary for ground data, use slitlet ID
    sources = slitlets
    # metadata IDs - arbitrary for ground data
    metaids = np.full_like(allcols, 1)
    # background shutter?  default to all "N" for ground data
    bkgd = np.full(len(allcols), 'N', dtype=str)
    # shutter state - "OPEN", by definition
    state = np.full(len(allcols), 'OPEN', dtype="<U8")
    # source position in shutter - N/A for ground data, assume centered
    srcx = np.full(len(allcols), 0.)
    srcy = srcx
    dithptind = np.full(len(allcols), 0)
    psrc = np.full(len(allcols), 'Y')
    tabcol1 = fits.Column(name='SLITLET_ID', format='I', array=slitlets)
    tabcol2 = fits.Column(name='MSA_METADATA_ID', format='I', array=metaids)
    tabcol3 = fits.Column(name='SHUTTER_QUADRANT', format='I', array=quads)
    tabcol4 = fits.Column(name='SHUTTER_ROW', format='I', array=allrows)
    tabcol5 = fits.Column(name='SHUTTER_COLUMN', format='I', array=allcols)
    tabcol6 = fits.Column(name='SOURCE_ID', format='I', array=sources)
    tabcol7 = fits.Column(name='BACKGROUND', format='A', array=bkgd)
    tabcol8 = fits.Column(name='SHUTTER_STATE', format='4A', array=state)
    tabcol9 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_X', format='E', array=srcx)
    tabcol10 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_Y', format='E', array=srcy)
    tabcol11 = fits.Column(name='DITHER_POINT_INDEX', format='I', array=dithptind)
    tabcol12 = fits.Column(name='PRIMARY_SOURCE', format='A', array=psrc)
    hdu2 = fits.BinTableHDU.from_columns(
        [tabcol1, tabcol2, tabcol3, tabcol4, tabcol5, tabcol6, tabcol7, tabcol8, tabcol9, tabcol10, tabcol11, tabcol12],
        name='SHUTTER_INFO')
    # SOURCE_INFO table
    # program ID - arbitrary
    program = np.full_like(allcols, 1)
    # source name - arbitrary
    name = np.full(len(allcols), 'lamp', dtype="<U8")
    # source alias - arbitrary
    alias = np.full(len(allcols), 'foo', dtype="<U8")
    # catalog ID - arbitrary
    catalog = np.full(len(allcols), 'foo', dtype="<U8")
    # RA, DEC - N/A for ground data
    ra = np.full(len(allcols), 0.)
    dec = ra
    # preimage file name - N/A
    preim = np.full(len(allcols), 'foo.fits', dtype="<U8")
    stellarity_arr = np.full(len(allcols), stellarity)
    tabcol1 = fits.Column(name='PROGRAM', format='I', array=program)
    tabcol2 = fits.Column(name='SOURCE_ID', format='I', array=sources)
    tabcol3 = fits.Column(name='SOURCE_NAME', format='4A', array=name)
    tabcol4 = fits.Column(name='ALIAS', format='3A', array=alias)
    tabcol5 = fits.Column(name='RA', format='D', array=ra)
    tabcol6 = fits.Column(name='DEC', format='D', array=dec)
    tabcol7 = fits.Column(name='PREIMAGE_ID', format='8A', array=preim)
    tabcol8 = fits.Column(name='STELLARITY', format='E', array=stellarity_arr)
    hdu3 = fits.BinTableHDU.from_columns([tabcol1, tabcol2, tabcol3, tabcol4, tabcol5, tabcol6, tabcol7, tabcol8],
                                         name='SOURCE_INFO')

    # create image of shutter status for first extension of the metafile
    # do each quadrant first, set value corresponding to each shutter depending on its status
    # (=0 if closed, 1 if open), then concatenate
    q1all = np.zeros(shape=(171, 365))
    q1all[[foid[2] for foid in failedopens if foid[0] == 1], [foid[1] for foid in failedopens if foid[0] == 1]] = 1
    q2all = np.zeros(shape=(171, 365))
    q2all[[foid[2] for foid in failedopens if foid[0] == 2], [foid[1] for foid in failedopens if foid[0] == 2]] = 1
    q3all = np.zeros(shape=(171, 365))
    q3all[[foid[2] for foid in failedopens if foid[0] == 3], [foid[1] for foid in failedopens if foid[0] == 3]] = 1
    q4all = np.zeros(shape=(171, 365))
    q4all[[foid[2] for foid in failedopens if foid[0] == 4], [foid[1] for foid in failedopens if foid[0] == 4]] = 1
    if debug:
        print("Quadrant 1 shape: ", q1all.shape)
        print("open shutters in Quadrant 1: ", np.where(q1all == 1))
        print("open shutters in Quadrant 2: ", np.where(q2all == 1))
        print("open shutters in Quadrant 3: ", np.where(q3all == 1))
        print("open shutters in Quadrant 4: ", np.where(q4all == 1))
    im1 = np.concatenate((q1all, q2all))
    im2 = np.concatenate((q3all, q4all))
    image = np.concatenate((im1, im2), axis=1)

    # plotting
    fig = plt.figure(figsize=(9, 9))
    norm = ImageNormalize(image, vmin=0., vmax=1., stretch=AsinhStretch())
    plt.imshow(image, norm=norm, aspect=1.0, origin='lower', cmap='viridis')
    if save_fig:
        datadir = os.path.dirname(outfile)
        plt_name = "FailedOpen_shutters" + ".png"
        plt_name = os.path.join(datadir, plt_name)
        plt.savefig(plt_name)
        print('Figure saved as: ', plt_name)
    if show_fig:
        plt.show()
    plt.close()

    # create fits file
    hdu0 = fits.PrimaryHDU()
    # add necessary keywords to primary header
    hdr = hdu0.header
    hdr.set('ORIGIN', 'STScI', 'institution responsible for creating FITS file')
    hdr.set('TELESCOP', 'JWST', 'telescope used to acquire data')
    hdr.set('INSTRUME', 'NIRSPEC', 'identifier for instrument used to acquire data')
    hdr.set('DATE', now.isoformat())
    hdr.set('FILENAME', outfile, 'name of file')
    hdr.set('PPSDBVER', 'PPSDB999',
            'version of PPS database used')  # N/A for non-OSS ground data, using arbitrary number
    hdr.set('PROGRAM', '1', 'program number')  # arbitrary
    hdr.set('VISIT', '1', 'visit number')  # arbitrary
    hdr.set('OBSERVTN', '1', 'observation number')  # arbitrary
    hdr.set('VISIT_ID', '1', 'visit identifier')  # arbitrary
    hdr.set('PNTG_SEQ', 1, 'pointing sequence number')  # arbitrary
    hdr.set('MSACFG10', 1, 'base 10 nirspec msa_at_pointing.msa_config_id')  # arbitrary
    hdr.set('MSACFG36', '01', 'base 36 version of MSACFG10')  # arbitrary
    hdu1 = fits.ImageHDU(image, name='SHUTTER_IMAGE')
    hdu_all = fits.HDUList([hdu0, hdu1, hdu2, hdu3])
    hdu_all.writeto(outfile)


def run_msa_flagging_testing(input_file, msa_flagging_threshold=99.5, stellarity=None, operability_ref=None, 
                             source_type=None, save_figs=False, show_figs=True, debug=False):
    """
    This is the validation function for the msa flagging step.
    :param input_file: string, fits file output from the msa_flagging step
    :param msa_flagging_threshold: float, percentage for all slits with more than 100 pixels
    :param stellarity: float, stellarity number fro 0.0 to 1.0
    :param operability_ref: string, msa failed open - operability - reference file
    :param source_type: string, options are point, extended, unknown
    :param save_figs: boolean
    :param show_figs: boolean
    :param debug: boolean
    :return:
        FINAL_TEST_RESULT: boolean, True if smaller than or equal to threshold
        result_msg: string, message with reason for passing, failing, or skipped
        log_msgs: list, diagnostic strings to be printed in log

    """
    # start the list of messages that will be added to the log file
    log_msgs = []

    # start the timer
    msa_flagging_test_start_time = time.time()

    # get the data model
    msaflag = datamodels.open(input_file)
    if debug:
        print('got MSA flagging datamodel!')

    # plot full image
    fig = plt.figure(figsize=(9, 9))
    norm = ImageNormalize(msaflag.data, vmin=0., vmax=50., stretch=AsinhStretch())
    plt.imshow(msaflag.data, norm=norm, aspect=1.0, origin='lower', cmap='viridis')
    # Show and/or save figures
    file_basename = os.path.basename(input_file.replace("_msa_flagging.fits", ""))
    datadir = os.path.dirname(input_file)
    detector = fits.getval(input_file, 'DETECTOR')
    if save_figs:
        t = (file_basename, "MSA_flagging_full_detector.png")
        plt_name = "_".join(t)
        plt_name = os.path.join(datadir, plt_name)
        plt.savefig(plt_name)
        print('Figure saved as: ', plt_name)
    if show_figs:
        plt.show()
    plt.close()

    # read in DQ flags from MSA_flagging product
    # find all pixels that have been flagged by this step -- DQ = 536870912
    dq_flag = 536870912
    msaflag_1d = msaflag.dq.flatten()
    index_opens = np.squeeze(np.asarray(np.where(msaflag_1d & dq_flag)))
    if debug:
        print("DQ array at 167, 1918: ", msaflag.dq[167, 1918])
        print("Index where Failed Open shutters exist: ", index_opens)

    # execute script that creates an MSA metafile for the failed open shutters
    # read operability reference file
    crds_path = os.environ.get('CRDS_PATH')
    references_json_file = "references/jwst/jwst_nirspec_msaoper_0001.json"

    if operability_ref is None:
        op_ref_file = os.path.join(crds_path, references_json_file)
    else:
        op_ref_file = operability_ref

    if not os.path.isfile(op_ref_file):
        result_msg = "Skipping msa_flagging test because the operability reference file does not exist: " + op_ref_file
        print(result_msg)
        log_msgs.append(result_msg)
        result = 'skip'
        return result, result_msg, log_msgs

    if debug:
        print("Using this operability reference file: ", op_ref_file)

    with open(op_ref_file) as f:
        msaoper_dict = json.load(f)
    msaoper = msaoper_dict["msaoper"]

    # find the failed open shutters
    failedopens = [(c["Q"], c["x"], c["y"]) for c in msaoper if c["state"] == 'open']
    if debug:
        print("Failed Open shutters: ", failedopens)

    # unpack the list of tuples into separate lists for MSA quadrant, row, column locations
    quads, allrows, allcols = zip(*failedopens)

    # stellarity -- internal lamps are uniform illumination, so set to 0
    # ** if considering a point source, need to change this to 1, or actual value if known
    if source_type is None:
        srctyapt = fits.getval(input_file, 'SRCTYAPT')
    else:
        srctyapt = source_type.upper()
    if stellarity is None:
        if "POINT" in srctyapt:
            stellarity = 1.0
        else:
            stellarity = 0.0
    else:
        stellarity = float(stellarity)

    # create MSA metafile with F/O shutters
    fometafile = os.path.join(datadir, 'fopens_metafile_msa.fits')
    if not os.path.isfile(fometafile):
        create_metafile_fopens(fometafile, allcols, allrows, quads, stellarity, failedopens,
                               save_fig=save_figs, show_fig=show_figs, debug=debug)

    # run assign_wcs on the science exposure using F/O metafile
    # change MSA metafile name in header to match the F/O metafile name
    rate_file = input_file.replace("msa_flagging", "rate")
    if not os.path.isfile(rate_file):
        # if a _rate.fits file does not exist try the usual name
        rate_file = os.path.join(datadir, 'final_output_caldet1_'+detector+'.fits')
        if not os.path.isfile(rate_file):
            result_msg = "Skipping msa_flagging test because no rate fits file was found in directory: " + datadir
            print(result_msg)
            log_msgs.append(result_msg)
            result = 'skip'
            return result, result_msg, log_msgs
    if debug:
        print("Will run assign_wcs with new Failed Open fits file on this file: ", rate_file)
    rate_mdl = datamodels.ImageModel(rate_file)
    if debug:
        print("MSA metadata file in initial rate file: ", rate_mdl.meta.instrument.msa_metadata_file)
    rate_mdl.meta.instrument.msa_metadata_file = fometafile
    if debug:
        print("New MSA metadata file in rate file: ", rate_mdl.meta.instrument.msa_metadata_file)
    # run assign_wcs; use +/-0.45 for the y-limits because the default is too big (0.6 including buffer)
    stp = AssignWcsStep()
    awcs_fo = stp.call(rate_mdl, slit_y_low=-0.45, slit_y_high=0.45)

    # get the slits from the F/O processing run
    slits_list = awcs_fo.meta.wcs.get_transform('gwa', 'slit_frame').slits

    # prepare arrays to hold info needed for validation test
    allsizes = np.zeros(len(slits_list))
    allchecks = np.zeros(len(slits_list))

    # loop over the slits and compare pixel bounds with the flagged pixels from the original product
    for i, slit in enumerate(slits_list):
        name = slit.name
        print("\nWorking with slit: ", name)
        print("Slit min and max in y direction: ", slit.ymin, slit.ymax)
        # get the WCS object for this particular slit
        wcs_slice = nirspec.nrs_wcs_set_input(awcs_fo, name)
        # get the bounding box for the 2D subwindow, round to nearest integer, and convert to integer
        bbox = np.rint(wcs_slice.bounding_box)
        bboxint = bbox.astype(int)
        print("bounding box rounded to next integer: ", bboxint)
        i1 = bboxint[0, 0]
        i2 = bboxint[0, 1]
        i3 = bboxint[1, 0]
        i4 = bboxint[1, 1]
        # make array of pixel locations within bounding box
        x, y = np.mgrid[i1:i2, i3:i4]
        index_1d = np.ravel_multi_index([[y], [x]], (2048, 2048))
        # get the slity WCS parameter to find which pixels are located in the actual spectrum
        det2slit = wcs_slice.get_transform('detector', 'slit_frame')
        slitx, slity, _ = det2slit(x, y)
        print("Max value in slity array (ignoring NANs): ", np.nanmax(slity))
        index_trace = np.squeeze(index_1d)[~np.isnan(slity)]
        n_overlap = np.sum(np.isin(index_opens, index_trace))
        if debug:
            print("Size of index_trace= ", index_trace.size)
            print("Size of index_opens=", index_opens.size)
            print("Sum of values found in index_opens and index_trace=", n_overlap)
        msg = 'percentage of F/O trace that was flagged: ' + repr(n_overlap/index_trace.size*100.)
        print(msg)
        log_msgs.append(msg)
        allchecks[i] = n_overlap/index_trace.size*100.
        allsizes[i] = index_trace.size

        # show 2D cutouts, with flagged pixels overlaid
        # calculate wavelength, slit_y values for the subwindow
        det2slit = wcs_slice.get_transform('detector', 'slit_frame')
        slitx, slity, _ = det2slit(x, y)

        # extract & display the F/O 2d subwindows from the msa_flagging sci image
        fig = plt.figure(figsize=(19, 19))
        subwin = msaflag.data[i3:i4, i1:i2].copy()
        # set all pixels outside of the nominal shutter length to 0, inside to 1
        subwin[np.isnan(slity.T)] = 0
        subwin[~np.isnan(slity.T)] = 1
        # find the pixels flagged by the msaflagopen step; set them to 1 and everything else to 0 for ease of display
        subwin_dq = msaflag.dq[i3:i4, i1:i2].copy()
        mask = np.zeros(subwin_dq.shape, dtype=bool)
        mask[np.where(subwin_dq & 536870912)] = True
        subwin_dq[mask] = 1
        subwin_dq[~mask] = 0
        # plot the F/O traces
        vmax = np.max(msaflag.data[i3:i4, i1:i2])
        norm = ImageNormalize(msaflag.data[i3:i4, i1:i2], vmin=0., vmax=vmax, stretch=AsinhStretch())
        plt.imshow(msaflag.data[i3:i4, i1:i2], norm=norm, aspect=10.0, origin='lower', cmap='viridis')
        plt.imshow(subwin, aspect=20.0, origin='lower', cmap='Reds', alpha=0.3)
        # overplot the flagged pixels in translucent grayscale
        plt.imshow(subwin_dq, aspect=20.0, origin='lower', cmap='gray', alpha=0.3)
        if save_figs:
            t = (file_basename, "FailedOpen_detector", detector, "slit", repr(name) + ".png")
            plt_name = "_".join(t)
            plt_name = os.path.join(datadir, plt_name)
            plt.savefig(plt_name)
            print('Figure saved as: ', plt_name)
        if show_figs:
            plt.show()
        plt.close()

    # validation: overlap should be >= msa_flagging_threshold percent for all slits with more than 100 pixels
    FINAL_TEST_RESULT = False
    if (allchecks[allsizes >= 100] >= msa_flagging_threshold).all():
        FINAL_TEST_RESULT = True
    if FINAL_TEST_RESULT:
        result_msg = "\n *** Final result for msa_flagging test will be reported as PASSED *** \n"
        print(result_msg)
        log_msgs.append(result_msg)
    else:
        result_msg = "\n *** Final result for msa_flagging test will be reported as FAILED *** \n"
        print(result_msg)
        log_msgs.append(result_msg)

    # end the timer
    msa_flagging_test_end_time = time.time() - msa_flagging_test_start_time
    if msa_flagging_test_end_time >= 60.0:
        msa_flagging_test_end_time = msa_flagging_test_end_time/60.0  # in minutes
        msa_flagging_test_tot_time = "* MSA flagging validation test took ", repr(msa_flagging_test_end_time) + \
                                     " minutes to finish."
        if msa_flagging_test_end_time >= 60.0:
            msa_flagging_test_end_time = msa_flagging_test_end_time/60.  # in hours
            msa_flagging_test_tot_time = "* MSA flagging validation test took ", repr(msa_flagging_test_end_time) + \
                                         " hours to finish."
    else:
        msa_flagging_test_tot_time = "* MSA flagging validation test took ", repr(msa_flagging_test_end_time) + \
                                  " seconds to finish."
    print(msa_flagging_test_tot_time)
    log_msgs.append(msa_flagging_test_tot_time)

    return FINAL_TEST_RESULT, result_msg, log_msgs


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("input_file",
                        action='store',
                        default=None,
                        help='Name of fits output file for the msa_flagging step, i.e. blah_msa_flagging.fits')
    parser.add_argument("-t",
                        dest="msa_flagging_threshold",
                        action='store',
                        default=99.5,
                        type=float,
                        help='Use flag -t to change the default threshold (currently set to 99.5%), e.g. -t=92.5')
    parser.add_argument("-s",
                        dest="stellarity",
                        action='store',
                        default=None,
                        help='Use flag -s to provide a specific stellarity value, e.g. -s=0.7')
    parser.add_argument("-n",
                        dest="save_figs",
                        action='store_false',
                        default=True,
                        help='Use flag -n to NOT save the figures.')
    parser.add_argument("-f",
                        dest="show_figs",
                        action='store_true',
                        default=False,
                        help='Use flag -f to show the figures.')
    parser.add_argument("-r",
                        dest="operability_ref",
                        action='store',
                        default=None,
                        help='Use flag -r to give a new operability reference file.')
    parser.add_argument("-o",
                        dest="object_type",
                        action='store',
                        default=None,
                        help='Use flag -o to provide the object (source) type, e.g. -o=point.')
    parser.add_argument("-d",
                        dest="debug",
                        action='store_true',
                        default=False,
                        help='Use flag -d to turn on debug mode.')
    args = parser.parse_args()

    # Set the variables input from the command line
    input_file = args.input_file
    msa_flagging_threshold = args.msa_flagging_threshold
    stellarity = args.stellarity
    save_figs = args.save_figs
    show_figs = args.show_figs
    operability_ref = args.operability_ref
    object_type = args.object_type
    debug = args.debug

    # Run the principal function of the script
    run_msa_flagging_testing(input_file, msa_flagging_threshold=msa_flagging_threshold, stellarity=stellarity,
                             operability_ref=operability_ref, save_figs=save_figs, show_figs=show_figs,
                             source_type=object_type, debug=debug)


if __name__ == '__main__':
    sys.exit(main())

