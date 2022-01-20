# create MOS metafile from MSA configuration

from astropy.io import fits, ascii
import numpy as np
import datetime
import argparse
import sys
import os
import urllib
import json
import tempfile
import collections

"""
This script will create the shutter configuration file for MOS data, while the pipeline is not automatically doing so.

 Usage:
  - Terminal
    CREATE NEW MSA SUTTER CONFIGURATION FILES
    In a terminal, in the directory where the testing tool lives, and within the pipeline environment, type:
    $ nptt_create_metafile CB10-GD-B.msa.fits  
    Output will be:
        CB10-GD-B_metafile_msa.fits
    
    FIX AN OLD SHUTTER CONFIGURATION FILE
    In a terminal, in the directory where the testing tool lives, and within the pipeline environment, type the command
    with the fix flag, -f:
    $ nptt_create_metafile V9621500100101_metafile_msa.fits -f
    Output will be:
        V9621500100101_metafile_msa_fixed.fits
    
    CREATE A SHUTTER CONFIGURATION FILE FROM SIMULATIONS OR DITHERS  
    In a terminal, in the directory where the testing tool lives, and within the pipeline environment, type the command
    with the dithers flag, -d, equal to the list of pointings (use a comma and no spaces bewteen file names), AND
    provide the number of shutters per slitlet with the flag -s:
    $ nptt_create_metafile V9621500100101.msa.fits -d=obs1.csv,obs2.csv,obs3.csv -s=3
    Output will be:
        V9621500100101_metafile_msa.fits
 
 - As a module
    # imports
    import nirspec_pipe_testing_tool as nptt
    
    # set the required variables 
    config_binary_file = something_msa.fits
    fix_old_config_file = False  # Only to be used if an older version of the file exists
    targ_file_list = ['something1.csv', 'something2.csv', 'something3.csv'] 
    shutters_in_slitlet = 5
    
    # set optional variables
    operability_ref = '/path_to_new/operability_ref.json' 
    flat_metafile = False  # set to True if this is for a flat metafile
    verbose = False
    
    # Run the function
    outfile = nptt.calwebb_spec2_pytests.auxiliary_code.create_metafile.run_create_metafile(config_binary_file, 
                                    fix_old_config_file, targ_file_list, shutters_in_slitlet, 
                                    operability_ref=operability_ref, flat_metafile=flat_metafile, 
                                    verbose=verbose)
    
"""

# HEADER
__author__ = "M. A. Pena-Guerrero & James Muzerolle"
__version__ = "1.2"

# HISTORY
# Oct 2019 - Version 1.0: initial version completed
# Mar 2020 - Version 1.1: minor bug fixes
# Dec 2021 - Version 1.2: added a condition to to create metafiles for flats

now = datetime.datetime.now()


def create_metafile(config_binary_file):
    """
    This function creates a new msa shutter configuration file starting from the .msa.fits files from APT.
    Args:
        config_binary_file: string, MSA configuration fits binary table

    Returns:
        outfile = string, the msa shutter configuration file with the right format for the pipeline
    """
    # read in configuration (binary fits format)
    hdul = fits.open(config_binary_file)
    # name of output metafile
    config_binary_file_list = config_binary_file.split(".")
    outfile = config_binary_file_list[0] + '_metafile_msa.fits'

    # create image for first extension
    q1all = hdul['Q1'].data.field('STATUS').reshape((171, 365))
    q2all = hdul['Q2'].data.field('STATUS').reshape((171, 365))
    q3all = hdul['Q3'].data.field('STATUS').reshape((171, 365))
    q4all = hdul['Q4'].data.field('STATUS').reshape((171, 365))
    # put quads together
    im1 = np.concatenate((q2all, q1all))
    im2 = np.concatenate((q3all, q4all))
    image = np.concatenate((im2, im1), axis=1)

    # find the open shutters (where status=1)
    # quad 1
    q1 = hdul['Q1'].data
    q1open = np.array(np.nonzero(q1.field('STATUS')))
    # convert to rows/cols
    q1col = (q1open - 1) // 365 + 1
    q1row = q1open - (q1col - 1) * 365
    # quad 2
    q2 = hdul['Q2'].data
    q2open = np.array(np.nonzero(q2.field('STATUS')))
    # convert to rows/cols
    q2col = (q2open - 1) // 365 + 1
    q2row = q2open - (q2col - 1) * 365
    # quad 3
    q3 = hdul['Q3'].data
    q3open = np.array(np.nonzero(q3.field('STATUS')))
    # convert to rows/cols
    q3col = (q3open - 1) // 365 + 1
    q3row = q3open - (q3col - 1) * 365
    # quad 4
    q4 = hdul['Q4'].data
    q4open = np.array(np.nonzero(q4.field('STATUS')))
    # convert to rows/cols
    q4col = (q4open - 1) // 365 + 1
    q4row = q4open - (q4col - 1) * 365

    hdul.close()

    # set up metafile table structure
    # SHUTTER_INFO table
    allcol = np.squeeze(np.concatenate((q1col, q2col, q3col, q4col), axis=1))
    allrow = np.squeeze(np.concatenate((q1row, q2row, q3row, q4row), axis=1))
    quads = np.squeeze(np.concatenate((np.full_like(q1col, 1), np.full_like(q2col, 2), np.full_like(q3col, 3),
                                       np.full_like(q4col, 4)), axis=1))
    # slitlet IDs - use one per shutter
    slitlets = np.arange(0, len(allcol)) + 1
    # source IDs - arbitrary for ground data, use slitlet ID
    sources = slitlets
    # metadata IDs - arbitrary for ground data
    metaids = np.full_like(allcol, 1)
    # background shutter?  default to all "N" for ground data
    bkgd = np.full(len(allcol), 'N', dtype=str)
    # shutter state - "OPEN", by definition
    state = np.full(len(allcol), 'OPEN', dtype="<U8")
    # source position in shutter - N/A for ground data, assume centered
    srcx = np.full(len(allcol), 0.)
    srcy = srcx
    # add dither point index and primary source columns
    nrows = len(bkgd)
    dithind = np.full(nrows, 1, dtype=int)
    psrc = np.full(nrows, 'Y', dtype=str)
    dither_point_index = fits.Column(name='DITHER_POINT_INDEX', format='I', array=dithind)
    primary_source = fits.Column(name='PRIMARY_SOURCE', format='1A', array=psrc)
    tabcol1 = fits.Column(name='SLITLET_ID', format='I', array=slitlets)
    tabcol2 = fits.Column(name='MSA_METADATA_ID', format='I', array=metaids)
    tabcol3 = fits.Column(name='SHUTTER_QUADRANT', format='I', array=quads)
    tabcol4 = fits.Column(name='SHUTTER_ROW', format='I', array=allrow)
    tabcol5 = fits.Column(name='SHUTTER_COLUMN', format='I', array=allcol)
    tabcol6 = fits.Column(name='SOURCE_ID', format='I', array=sources)
    tabcol7 = fits.Column(name='BACKGROUND', format='A', array=bkgd)
    tabcol8 = fits.Column(name='SHUTTER_STATE', format='4A', array=state)
    tabcol9 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_X', format='E', array=srcx)
    tabcol10 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_Y', format='E', array=srcy)
    hdu2 = fits.BinTableHDU.from_columns([tabcol1, tabcol2, tabcol3, tabcol4, tabcol5, tabcol6, tabcol7, tabcol8,
                                          tabcol9, tabcol10, dither_point_index, primary_source], name='SHUTTER_INFO')

    # SOURCE_INFO table
    # program ID - arbitrary
    program = np.full_like(allcol, 1)
    # source name - arbitrary
    name = np.full(len(allcol), 'lamp', dtype="<U8")
    # source alias - arbitrary
    alias = np.full(len(allcol), 'foo', dtype="<U8")
    # catalog ID - arbitrary
    catalog = np.full(len(allcol), 'foo', dtype="<U8")
    # RA, DEC - N/A for ground data
    ra = np.full(len(allcol), 0.)
    dec = ra
    # preimage file name - N/A
    preim = np.full(len(allcol), 'foo.fits', dtype="<U8")
    # stellarity -- internal lamps are uniform illumination, so set to 0
    # ** if considering a point source, need to change this to 0, or actual value if known
    stellarity = np.full(len(allcol), 0.)
    tabcol1 = fits.Column(name='PROGRAM', format='I', array=program)
    tabcol2 = fits.Column(name='SOURCE_ID', format='I', array=sources)
    tabcol3 = fits.Column(name='SOURCE_NAME', format='4A', array=name)
    tabcol4 = fits.Column(name='ALIAS', format='3A', array=alias)
    # tabcol5 = fits.Column(name='CATALOG_ID', format='3A', array=catalog)
    tabcol6 = fits.Column(name='RA', format='D', array=ra)
    tabcol7 = fits.Column(name='DEC', format='D', array=dec)
    tabcol8 = fits.Column(name='PREIMAGE_ID', format='8A', array=preim)
    tabcol9 = fits.Column(name='STELLARITY', format='E', array=stellarity)
    hdu3 = fits.BinTableHDU.from_columns([tabcol1, tabcol2, tabcol3, tabcol4, tabcol5, tabcol6, tabcol7,
                                          tabcol8, tabcol9], name='SOURCE_INFO')

    # create fits file
    hdu0 = fits.PrimaryHDU()
    # add necessary keywords to primary header
    hdr = hdu0.header
    hdr.set('ORIGIN', 'STScI', 'institution responsible for creating FITS file')
    hdr.set('TELESCOP', 'JWST', 'telescope used to acquire data')
    hdr.set('INSTRUME', 'NIRSPEC', 'identifier for instrument used to acquire data')
    hdr.set('DATE', now.isoformat())
    hdr.set('FILENAME', outfile, 'name of file')
    hdr.set('PPSDBVER', 'PPSDB999', 'version of PPS database used')  # N/A for non-OSS ground data,
    # using arbitrary number
    hdr.set('PROGRAM', '1', 'program number')  # arbitrary
    hdr.set('VISIT', '1', 'visit number')  # arbitrary
    hdr.set('OBSERVTN', '1', 'observation number')  # arbitrary
    hdr.set('VISIT_ID', '1', 'visit identifier')  # arbitrary
    # hdr.set('PNTG_SEQ', 1, 'pointing sequence number') # arbitrary
    hdr.set('MSACFG10', 1, 'base 10 nirspec msa_at_pointing.msa_config_id')  # arbitrary
    hdr.set('MSACFG36', '01', 'base 36 version of MSACFG10')  # arbitrary
    hdu1 = fits.ImageHDU(image, name='SHUTTER_IMAGE')
    hdu_all = fits.HDUList([hdu0, hdu1, hdu2, hdu3])
    hdu_all.writeto(outfile)
    return outfile


def fix_metafile(infile):
    """
    This function fixes the shutter configuration files for the new format (of build 7.3)
    Args:
        infile: string, old MSA configuration fits

    Returns:
        outfile = string, the msa shutter configuration file with the right format for the pipeline
    """
    # read in old metafile and set output file name
    hdul1 = fits.open(infile)
    outfile = infile.replace('.fits', '_fixed.fits')

    # SHUTTER_INFO table
    bkg = hdul1[2].data['BACKGROUND']
    nrows = len(bkg)
    # add dither point index and primary source columns
    dithind = np.full(nrows, 1, dtype=int)
    psrc = np.full(nrows, 'Y', dtype=str)
    dither_point_index = fits.Column(name='DITHER_POINT_INDEX', format='I', array=dithind)
    primary_source = fits.Column(name='PRIMARY_SOURCE', format='1A', array=psrc)
    new_columns = hdul1[2].columns + dither_point_index + primary_source
    hdu2 = fits.BinTableHDU.from_columns(new_columns, name="SHUTTER_INFO")

    # SOURCE_INFO table
    # separate original columns; need to leave out "catalog_id"
    nsrcs = len(hdul1[3].data['PROGRAM'])
    newcols = [hdul1[3].columns["PROGRAM"], hdul1[3].columns["SOURCE_ID"], hdul1[3].columns["SOURCE_NAME"],
               hdul1[3].columns["ALIAS"], hdul1[3].columns["RA"], hdul1[3].columns["DEC"],
               hdul1[3].columns["PREIMAGE_ID"], hdul1[3].columns["STELLARITY"]]
    # create new output hdu
    hdu3 = fits.BinTableHDU.from_columns(newcols, name="SOURCE_INFO", nrows=int(nsrcs))

    # create new fits file
    # get primary header from one of the inputs
    hdr = hdul1[0].header
    hdu0 = hdul1[0]
    del hdr["PNTG_SEQ"]  # delete this unnecessary keyword
    hdr["DATE"] = str(now)
    hdu0.header = hdr
    hdu1 = fits.ImageHDU(hdul1[1].data, name='SHUTTER_IMAGE')
    hdu_all = fits.HDUList([hdu0, hdu1, hdu2, hdu3])
    hdu_all.writeto(outfile)
    return outfile


def create_metafile_dither(config_binary_file, targ_file_list, shutters_in_slitlet, operability_ref=None,
                           output_dir=None, flat_metafile=False, verbose=False):
    """
    This function creates the shutter configuration file for dithers, using the csv files
    Args:
        config_binary_file: string, MSA configuration fits binary table
        targ_file_list: string, dither file or list of dither files
        shutters_in_slitlet: integer, number of shutters per slitlet
        operability_ref: string, MSA operability file
        output_dir: string, path to place the output file - if None output will be in same dir as input
        verbose: boolean

    Returns:
        outfile = string, the msa shutter configuration file with the right format for the pipeline
    """
    if shutters_in_slitlet is None:
        print('(create_metafile.create_metafile_dither:) Variable shutters_in_slitlet was not defined. Re-run the '
              'script with the flag -s set to the number of shutters per slitlet, e.g. -s=3.')
        print('                                          Exiting script.')
        exit()

    # get the MSA operability file
    crds_path = "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    op_ref_file = "jwst_nirspec_msaoper_0001.json"
    temp_dir = None
    if operability_ref is None:
        ref_file = os.path.join(crds_path, op_ref_file)
        temp_dir = tempfile.TemporaryDirectory()
        op_ref_file = os.path.join(temp_dir.name, op_ref_file)
        urllib.request.urlretrieve(ref_file, op_ref_file)
    else:
        op_ref_file = operability_ref
    if "https:" not in op_ref_file:
        if not os.path.isfile(op_ref_file):
            print("(create_metafile.create_metafile_dither:) WARNING! operability reference file does not exist: ",
                  op_ref_file)
            print('                                          Exiting script.')
            exit()
    print("(create_metafile.create_metafile_dither:) Using this operability reference file: ", op_ref_file)
    with open(op_ref_file) as f:
        msaoper_dict = json.load(f)
    msaoper = msaoper_dict["msaoper"]
    # find the failed open shutters
    failedopens = [(c["Q"], c["x"], c["y"]) for c in msaoper if c["state"] == 'open']
    if verbose:
        print("(create_metafile.create_metafile_dither:) Failed Open shutters: ", failedopens)
    # erase the local copy of the reference operability file
    if temp_dir is not None:
        temp_dir.cleanup()

    # read in configuration (binary fits format)
    hdul = fits.open(config_binary_file)
    # name of output metafile
    if '.msa.fits' in config_binary_file:
        outfile = config_binary_file.replace('.msa.fits', '') + '_metafile_msa.fits'
    else:
        outfile = config_binary_file[:-5] + '_metafile_msa.fits'
    if output_dir is not None:
        outfile = os.path.join(output_dir, os.path.basename(outfile))
    # check if outfile already exists, delete if it does
    if os.path.isfile(outfile):
        os.remove(outfile)

    # get number of nods from the number of input target files
    if isinstance(targ_file_list, str):
        targ_files = targ_file_list.split(',')
    else:
        targ_files = targ_file_list
    nnods = len(targ_files)
    print('(create_metafile.create_metafile_dither:) Number of nods =', nnods)

    # read in the csv target files
    iter = True
    for targfile in targ_files:
        print('(create_metafile.create_metafile_dither:) extracting target info from:', targfile)
        targinfo = ascii.read(targfile)
        if iter:
            source_id = targinfo['ID'].data
            source_ra = targinfo['Source RA (Degrees)'].data
            source_dec = targinfo['Source Dec (Degrees)'].data
            source_quad = targinfo['Quadrant'].data
            source_row = targinfo['Row'].data
            source_col = targinfo['Column'].data
            source_xpos = targinfo['Offset (x)'].data
            source_ypos = targinfo['Offset (y)'].data
            iter = False
        else:
            source_id = np.concatenate(([source_id, targinfo['ID'].data]))
            source_ra = np.concatenate(([source_ra, targinfo['Source RA (Degrees)'].data]))
            source_dec = np.concatenate(([source_dec, targinfo['Source Dec (Degrees)'].data]))
            source_quad = np.concatenate(([source_quad, targinfo['Quadrant'].data]))
            source_row = np.concatenate(([source_row, targinfo['Row'].data]))
            source_col = np.concatenate(([source_col, targinfo['Column'].data]))
            source_xpos = np.concatenate(([source_xpos, targinfo['Offset (x)'].data]))
            source_ypos = np.concatenate(([source_ypos, targinfo['Offset (y)'].data]))

    # create image for first extension
    q1all = hdul['Q1'].data.field('STATUS').reshape((171, 365))
    q2all = hdul['Q2'].data.field('STATUS').reshape((171, 365))
    q3all = hdul['Q3'].data.field('STATUS').reshape((171, 365))
    q4all = hdul['Q4'].data.field('STATUS').reshape((171, 365))
    # put quads together
    im1 = np.concatenate((q2all, q1all))
    im2 = np.concatenate((q3all, q4all))
    image = np.concatenate((im2, im1), axis=1)

    # find the open shutters (where status=1)
    # quad 1
    q1 = hdul['Q1'].data
    q1open = np.array(np.nonzero(q1.field('STATUS')))
    # convert to rows/cols
    q1col = (q1open - 1) // 365 + 1
    q1row = q1open - (q1col - 1) * 365 + 1
    # quad 2
    q2 = hdul['Q2'].data
    q2open = np.array(np.nonzero(q2.field('STATUS')))
    # convert to rows/cols
    q2col = (q2open - 1) // 365 + 1
    q2row = q2open - (q2col - 1) * 365 + 1
    # quad 3
    q3 = hdul['Q3'].data
    q3open = np.array(np.nonzero(q3.field('STATUS')))
    # convert to rows/cols
    q3col = (q3open - 1) // 365 + 1
    q3row = q3open - (q3col - 1) * 365 + 1
    # quad 4
    q4 = hdul['Q4'].data
    q4open = np.array(np.nonzero(q4.field('STATUS')))
    # convert to rows/cols
    q4col = (q4open - 1) // 365 + 1
    q4row = q4open - (q4col - 1) * 365 + 1
    hdul.close()

    # set up metafile table structure and fill

    # SHUTTER_INFO table
    # combine quadrants and copy rows of each shutter location array to account for nods
    allcols = np.repeat(np.squeeze(np.concatenate((q1col, q2col, q3col, q4col), axis=1)), nnods)
    allrows = np.repeat(np.squeeze(np.concatenate((q1row, q2row, q3row, q4row), axis=1)), nnods)
    quads = np.repeat(np.squeeze(
        np.concatenate((np.full_like(q1col, 1), np.full_like(q2col, 2), np.full_like(q3col, 3), np.full_like(q4col, 4)),
                       axis=1)), nnods)

    # initialize arrays for table columns
    # slitlet IDs - arbitrary numbering
    # tab_slitlet_id = np.full(allcols.size, 0)  # replaced by dictionary per source ID
    # metadata ID (arbitrary)
    tab_metadata_id = np.full(allcols.size, 1)
    # source IDs - from target file (based on original MPT catalog)
    tab_source_id = np.full(allcols.size, 0)
    # background shutter status
    tab_background = np.full(allcols.size, 'Y', dtype=str)
    # shutter state - "OPEN", by definition
    tab_shutter_state = np.full(allcols.size, 'OPEN', dtype="<U8")
    # source position in shutter (NaN for all except shutters containing the source
    tab_source_xpos = np.full(allcols.size, np.nan, dtype=float)
    tab_source_ypos = np.full(allcols.size, np.nan, dtype=float)
    tab_ditherpoint_index = np.full(allcols.size, 1)
    tab_primary_source = np.full(allcols.size, 'N', dtype="<U2")

    # loop over open shutters

    nshut = len(quads)
    # nod pattern - table rows that contain the source
    if nnods == 1:
        nod = [0]
    if nnods == 2:
        nod = [0, 3]
    if nnods == 3:
        nod = [1, 3, 8]
    if nnods > 3:
        print("(create_metafile.create_metafile_dither:) WARNING! nods > 3 not currently supported")
        print("                                          Exiting script.")
        sys.exit()

    obs_failedopen, all_source_tsid = 0, []
    slitlet_id_dict = collections.OrderedDict()
    for i in range(nshut):
        for fail_open_id in failedopens:
            # check if this is a failed open shutter
            if quads[i] == fail_open_id[0] and allrows[i] == fail_open_id[1] and allcols[i] == fail_open_id[2]:
                if verbose:
                    print('(create_metafile.create_metafile_dither:) Failed open shutter at ',
                          quads[i], allrows[i], allcols[i])
                obs_failedopen += 1
                break

        match = np.intersect1d(np.where(quads[i] == source_quad),
                               np.intersect1d(np.where(allrows[i] == source_row),
                                              np.where(allcols[i] == source_col)))
        match_length = len(match)
        if match_length != 0:
            if match_length > 1:
                match = match[0]

            # if there was a match AND this is not a failed open shutter, add source to the table
            tab_source_id[i] = source_id[match]

            if tab_source_id[i] not in slitlet_id_dict:
                tab_slitlet_id = int(i / nnods / shutters_in_slitlet) + 1
                while True:
                    if tab_slitlet_id in all_source_tsid:
                        tab_slitlet_id += 1
                    else:
                        break
                all_source_tsid.append(tab_slitlet_id)
                slitlet_id_dict[tab_source_id[i]] = {'nodsourceid': [],
                                                     'count_repetitions': 1,
                                                     'source_tsid': tab_slitlet_id}
            else:
                slitlet_id_dict[tab_source_id[i]]['count_repetitions'] += 1
                tab_slitlet_id = slitlet_id_dict[tab_source_id[i]]['source_tsid']

            # set the position of the source in the slitlet configuration
            count_repetitions = slitlet_id_dict[tab_source_id[i]]['count_repetitions']
            pos = count_repetitions - 1

            if verbose:
                print('(create_metafile.create_metafile_dither:) Found target, adding ID shutter ',
                      tab_slitlet_id)
                print('  Additional information:  shutter ID = ', tab_slitlet_id,
                      'quadrant = ', source_quad[match], '  row = ', source_row[match],
                      '  column = ', source_col[match], '  source ID = ', source_id[match],
                      '  source xpos = ', source_xpos[match], '  source ypos = ', source_ypos[match])

            if pos in nod:
                # some values only get filled if the source is in that shutter for that nod
                tab_background[i] = 'N'
                tab_source_xpos[i] = source_xpos[match]
                tab_source_ypos[i] = source_ypos[match]
                tab_primary_source[i] = 'Y'

            if flat_metafile:
                tab_source_xpos[i], tab_source_ypos[i] = 0.0, 0.0
            tab_ditherpoint_index[i] = np.mod(i, nnods) + 1

            # track the repetitions for source ID in the shutter configuration file
            nodsourceid = [tab_slitlet_id, tab_metadata_id[i], tab_source_id[i], tab_background[i],
                           tab_shutter_state[i], tab_source_xpos[i], tab_source_ypos[i],
                           tab_ditherpoint_index[i], tab_primary_source[i]]
            # store in the dictionary per source ID
            slitlet_id_dict[tab_source_id[i]]['nodsourceid'].append(nodsourceid)

        else:
            if verbose:
                print("(create_metafile.create_metafile_dither:) No source corresponding to shutter: ",
                      quads[i], allrows[i], allcols[i])

    # create and fill arrays to print into metafile
    tab_slitlet_id, tab_metadata_id, tab_source_id = np.array([]), np.array([]), np.array([])
    tab_shutter_state, tab_source_xpos, tab_source_ypos = np.array([]), np.array([]), np.array([])
    tab_background, tab_ditherpoint_index, tab_primary_source = np.array([]), np.array([]), np.array([])

    for sid, sid_dict in slitlet_id_dict.items():
        for nid in range(len(sid_dict['nodsourceid'])):
            tab_slitlet_id = np.append(tab_slitlet_id, sid_dict['nodsourceid'][nid][0])
            tab_metadata_id = np.append(tab_metadata_id, sid_dict['nodsourceid'][nid][1])
            tab_source_id = np.append(tab_source_id, sid_dict['nodsourceid'][nid][2])
            tab_background = np.append(tab_background, sid_dict['nodsourceid'][nid][3])
            tab_shutter_state = np.append(tab_shutter_state, sid_dict['nodsourceid'][nid][4])
            tab_source_xpos = np.append(tab_source_xpos, sid_dict['nodsourceid'][nid][5])
            tab_source_ypos = np.append(tab_source_ypos, sid_dict['nodsourceid'][nid][6])
            tab_ditherpoint_index = np.append(tab_ditherpoint_index, sid_dict['nodsourceid'][nid][7])
            tab_primary_source = np.append(tab_primary_source, sid_dict['nodsourceid'][nid][8])

    # collect all unique instances of the sources
    unique_source_id, match, match2 = np.intersect1d(source_id, tab_source_id, return_indices=True)
    # RA, DEC, rows, cols, and quadrants
    ra, dec = source_ra[match], source_dec[match]
    rows, cols, quadrants = source_row[match], source_col[match], source_quad[match]

    if verbose:
        for i in range(len(ra)):
            print('shutter ID =', tab_slitlet_id[i], '  quadrant =', quadrants[i], '  row =', rows[i],
                  '  column =', cols[i], '  source ID =', tab_source_id[i],
                  '  source xpos =', tab_source_xpos[i], '  source ypos =', tab_source_ypos[i],
                  '  sRA =', ra[i], '  sDEC =', dec[i])

    """
    # obtain the background slitlets
    remaining_open_shutters = len(tab_slitlet_id[mask]) / sourceid_repetitions
    obs_total_open = nshut - obs_failedopen
    #print(nshut, obs_failedopen, obs_total_open, remaining_open_shutters,
    # len(tab_slitlet_id[mask]), tab_slitlet_id[mask])
    # Potential code to add background slits
    background_id_starting_point = 900
    for bi in range(backgrounds):
        # assuming there should be 3 lines per background shutter
        for _ in range(3):
            tab_slitlet_id.append(background_id_starting_point+bi)
            tab_metadata_id.append(????)
            * determine location in quads, allrows, allcols
            tab_source_id.append(????)
            tab_background.append('Y')
            tab_shutter_state.append(???)
            tab_source_xpos.append('NULL')
            tab_source_ypos.append('NULL')
            tab_ditherpoint_index.append(1,2,3)   # representing the number of nods
            tab_primary_source.append('N')
    """

    # create and fill in fits table columns
    tabcol1 = fits.Column(name='SLITLET_ID', format='I', array=tab_slitlet_id)
    tabcol2 = fits.Column(name='MSA_METADATA_ID', format='I', array=tab_metadata_id)
    tabcol3 = fits.Column(name='SHUTTER_QUADRANT', format='I', array=quadrants)
    tabcol4 = fits.Column(name='SHUTTER_ROW', format='I', array=rows)
    tabcol5 = fits.Column(name='SHUTTER_COLUMN', format='I', array=cols)
    tabcol6 = fits.Column(name='SOURCE_ID', format='J', array=tab_source_id)
    tabcol7 = fits.Column(name='BACKGROUND', format='A', array=tab_background)
    tabcol8 = fits.Column(name='SHUTTER_STATE', format='4A', array=tab_shutter_state)
    tabcol9 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_X', format='E', array=tab_source_xpos)
    tabcol10 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_Y', format='E', array=tab_source_ypos)
    tabcol11 = fits.Column(name='DITHER_POINT_INDEX', format='I', array=tab_ditherpoint_index)
    tabcol12 = fits.Column(name='PRIMARY_SOURCE', format='1A', array=tab_primary_source)
    hdu2 = fits.BinTableHDU.from_columns(
        [tabcol1, tabcol2, tabcol3, tabcol4, tabcol5, tabcol6, tabcol7, tabcol8, tabcol9, tabcol10, tabcol11, tabcol12],
        name='SHUTTER_INFO')

    # SOURCE_INFO table
    nsources = len(unique_source_id)
    # program ID - arbitrary
    program = np.full(nsources, '111', dtype="<U8")
    # source name - arbitrary
    source_name = np.core.defchararray.add('111_', unique_source_id.astype(str))
    # source alias - arbitrary
    alias = np.full(nsources, 'foo', dtype="<U8")
    # preimage file name - arbitrary
    preim = np.full(nsources, 'foo_pre-image.fits', dtype="<U18")
    # stellarity -- assuming perfect point sources, so set to 1
    # ** if considering an extended source, need to change this to 0, or actual value if known
    stellarity = np.full(nsources, 1.)
    tabcol1 = fits.Column(name='PROGRAM', format='J', array=program)
    tabcol2 = fits.Column(name='SOURCE_ID', format='J', array=unique_source_id)
    tabcol3 = fits.Column(name='SOURCE_NAME', format='11A', array=source_name)
    tabcol4 = fits.Column(name='ALIAS', format='5A', array=alias)
    tabcol5 = fits.Column(name='RA', format='D', array=ra)
    tabcol6 = fits.Column(name='DEC', format='D', array=dec)
    tabcol7 = fits.Column(name='PREIMAGE_ID', format='18A', array=preim)
    tabcol8 = fits.Column(name='STELLARITY', format='E', array=stellarity)
    hdu3 = fits.BinTableHDU.from_columns([tabcol1, tabcol2, tabcol3, tabcol4, tabcol5, tabcol6, tabcol7, tabcol8],
                                         name='SOURCE_INFO')

    # create fits file
    hdu0 = fits.PrimaryHDU()
    # add necessary keywords to primary header
    hdr = hdu0.header
    hdr.set('ORIGIN', 'STScI', 'institution responsible for creating FITS file')
    hdr.set('TELESCOP', 'JWST', 'telescope used to acquire data')
    hdr.set('INSTRUME', 'NIRSPEC', 'identifier for instrument used to acquire data')
    hdr.set('DATE', now.isoformat())
    hdr.set('FILENAME', outfile, 'name of file')
    hdr.set('PPSDBVER', 'PPSDB999', 'version of PPS database used')  # using arbitrary number
    hdr.set('PROGRAM', '111', 'program number')  # arbitrary
    hdr.set('VISIT', '1', 'visit number')  # arbitrary
    hdr.set('OBSERVTN', '1', 'observation number')  # arbitrary
    hdr.set('VISIT_ID', '1', 'visit identifier')  # arbitrary
    hdr.set('MSACFG10', 1, 'base 10 nirspec msa_at_pointing.msa_config_id')  # arbitrary
    hdr.set('MSACFG36', '01', 'base 36 version of MSACFG10')  # arbitrary
    hdu1 = fits.ImageHDU(image.T, name='SHUTTER_IMAGE')
    hdu_all = fits.HDUList([hdu0, hdu1, hdu2, hdu3])
    hdu_all.writeto(outfile)
    return outfile


def run_create_metafile(config_binary_file, fix_old_config_file, targ_file_list, shutters_in_slitlet,
                        operability_ref=None, output_dir=None, flat_metafile=False, verbose=False):
    """
    This function is a wrapper for all cases.
    Args:
        config_binary_file: string, MSA configuration fits binary table
        fix_old_config_file: boolean
        targ_file_list: string, nod set of .csv files divided by commas (no spaces)
        shutters_in_slitlet: integer, number of shutters per slitlet
        operability_ref: string, MSA operability file
        output_dir: string, path to place the output file - if None output will be in same dir as input
        verbose: boolean

    Returns:
        outfile = string, the msa shutter configuration file with the right format for the pipeline
    """
    if targ_file_list is None:
        if not fix_old_config_file:
            outfile = create_metafile(config_binary_file)
        else:
            outfile = fix_metafile(config_binary_file)
    else:
        outfile = create_metafile_dither(config_binary_file, targ_file_list, shutters_in_slitlet,
                                         operability_ref=operability_ref, output_dir=output_dir,
                                         flat_metafile=flat_metafile, verbose=verbose)
    return outfile


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("config_binary_file",
                        action='store',
                        default=None,
                        help='Name MSA configuration fits binary table, i.e. CB10-GD-B.msa.fits')
    parser.add_argument("-d",
                        dest="targ_file_list",
                        action='store',
                        default=None,
                        help='For dithers only, use this flag to provide csv output from APT; list each one '
                             'in a nod set, e.g. -s=obs1.csv,obs1.csv')
    parser.add_argument("-f",
                        dest="fix_old_config_file",
                        action='store_true',
                        default=False,
                        help='If an old version of the shutter configuration file exists, use the -f flag '
                             'to fix it. In this case, the input to the script is the old config file.fits.')
    parser.add_argument("-s",
                        dest="shutters_in_slitlet",
                        action='store',
                        default=None,
                        help='Use the flag -s to provide the number of shutters per slitlet, e.g. -s=3.')
    parser.add_argument("-r",
                        dest="operability_ref",
                        action='store',
                        default=None,
                        help='Use the flag -r to provide the MSA operability file to use, '
                             'e.g. -o=jwst_nirspec_msaoper_0001.json. If None, the default in CRDS will be used.')
    parser.add_argument("-o",
                        dest="output_dir",
                        action='store',
                        default=None,
                        help='Use the flag -o to provide the output directory.')
    parser.add_argument("-mf",
                        dest="flat_metafile",
                        action='store_true',
                        default=False,
                        help='Use the flag -fm if building a metafile for a flat.')
    parser.add_argument("-v",
                        dest="verbose",
                        action='store_true',
                        default=False,
                        help='Use the flag -v to print a series of messages throughout the script.')
    args = parser.parse_args()

    # Set the variables
    config_binary_file = args.config_binary_file
    targ_file_list = args.targ_file_list
    fix_old_config_file = args.fix_old_config_file
    shutters_in_slitlet = args.shutters_in_slitlet
    operability_ref = args.operability_ref
    output_dir = args.output_dir
    flat_metafile = args.flat_metafile
    verbose = args.verbose

    if shutters_in_slitlet is not None:
        shutters_in_slitlet = int(shutters_in_slitlet)

    # Run the function
    run_create_metafile(config_binary_file, fix_old_config_file, targ_file_list, shutters_in_slitlet,
                        operability_ref=operability_ref, output_dir=output_dir,
                        flat_metafile=flat_metafile, verbose=verbose)

    print("\nScript create_metafile.py done.\n")


if __name__ == '__main__':
    sys.exit(main())
