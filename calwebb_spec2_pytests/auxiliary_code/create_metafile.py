# create MOS metafile from MSA configuration

from astropy.io import fits, ascii
import numpy as np
import datetime
import argparse
import sys
import os

"""
This script will create the shutter configuration file for MOS data, while the pipeline is not automatically doing so.

 Usage:
    CREATE NEW MSA SUTTER CONFIGURATION FILES
    In a terminal, in the directory where the testing tool lives, and within the pipeline environment, type:
    > python create_metafile.py CB10-GD-B.msa.fits  
    Output will be:
        CB10-GD-B_metafile_msa.fits
    
    FIX AN OLD SHUTTER CONFIGURATION FILE
    In a terminal, in the directory where the testing tool lives, and within the pipeline environment, type the command
    with the fix flag:
    > python create_metafile.py V9621500100101_metafile_msa.fits -f
    Output will be:
        V9621500100101_metafile_msa_fixed.fits
    
    CREATE A SHUTTER CONFIGURATION FILE FROM SIMULATIONS        
    In a terminal, in the directory where the testing tool lives, and within the pipeline environment, type the command
    with the simulations flag equal to the list of pointings (use a comma and no spaces bewteen file names):
    > python create_metafile.py V9621500100101.msa.fits -s=obs1.csv,obs2.csv,obs3.csv
    Output will be:
        V9621500100101_metafile_msa.fits
"""


# HEADER
__author__ = "M. A. Pena-Guerrero & James Muzerolle"
__version__ = "1.0"

# HISTORY
# Oct 2019 - Version 1.0: initial version completed

now = datetime.datetime.now()


def create_metafile(config_binary_file):
    """
    This function creates a new msa shutter configuration file starting from the .msa.fits files from APT.
    Args:
        config_binary_file: string, MSA configuration fits binary table

    Returns:
        The msa shutter configuration file with the right format for the pipeline
    """
    # read in configuration (binary fits format)
    hdul = fits.open(config_binary_file)
    # name of output metafile
    config_binary_file_list = config_binary_file.split(".")
    outfile = config_binary_file_list[0]+'_metafile_msa.fits'

    # create image for first extension
    q1all = hdul['Q1'].data.field('STATUS').reshape((171,365))
    q2all = hdul['Q2'].data.field('STATUS').reshape((171,365))
    q3all = hdul['Q3'].data.field('STATUS').reshape((171,365))
    q4all = hdul['Q4'].data.field('STATUS').reshape((171,365))
    # put quads together
    im1 = np.concatenate((q2all,q1all))
    im2 = np.concatenate((q3all,q4all))
    image = np.concatenate((im2,im1),axis=1)

    # find the open shutters (where status=1)
    # quad 1
    q1 = hdul['Q1'].data
    q1open = np.array(np.nonzero(q1.field('STATUS')))
    # convert to rows/cols
    q1col = (q1open-1)//365+1
    q1row = q1open-(q1col-1)*365
    # quad 2
    q2 = hdul['Q2'].data
    q2open = np.array(np.nonzero(q2.field('STATUS')))
    # convert to rows/cols
    q2col = (q2open-1)//365+1
    q2row = q2open-(q2col-1)*365
    # quad 3
    q3 = hdul['Q3'].data
    q3open = np.array(np.nonzero(q3.field('STATUS')))
    # convert to rows/cols
    q3col = (q3open-1)//365+1
    q3row = q3open-(q3col-1)*365
    # quad 4
    q4 = hdul['Q4'].data
    q4open = np.array(np.nonzero(q4.field('STATUS')))
    # convert to rows/cols
    q4col = (q4open-1)//365+1
    q4row = q4open-(q4col-1)*365

    hdul.close()

    # set up metafile table structure
    # SHUTTER_INFO table
    allcol = np.squeeze(np.concatenate((q1col,q2col,q3col,q4col),axis=1))
    allrow = np.squeeze(np.concatenate((q1row,q2row,q3row,q4row),axis=1))
    quads = np.squeeze(np.concatenate((np.full_like(q1col,1),np.full_like(q2col,2),np.full_like(q3col,3),np.full_like(q4col,4)),axis=1))
    # slitlet IDs - use one per shutter
    slitlets = np.arange(0,len(allcol))+1
    # source IDs - arbitrary for ground data, use slitlet ID
    sources = slitlets
    # metadata IDs - arbitrary for ground data
    metaids = np.full_like(allcol,1)
    # background shutter?  default to all "N" for ground data
    bkgd = np.full(len(allcol),'N',dtype=str)
    # shutter state - "OPEN", by definition
    state = np.full(len(allcol),'OPEN',dtype="<U8")
    # source position in shutter - N/A for ground data, assume centered
    srcx = np.full(len(allcol),0.)
    srcy = srcx
    # add dither point index and primary source columns
    nrows = len(bkgd)
    dithind = np.full(nrows, 1, dtype=int)
    psrc = np.full(nrows, 'Y', dtype=str)
    dither_point_index = fits.Column(name='DITHER_POINT_INDEX', format='I', array=dithind)
    primary_source = fits.Column(name='PRIMARY_SOURCE', format='1A', array=psrc)
    tabcol1 = fits.Column(name='SLITLET_ID',format='I',array=slitlets)
    tabcol2 = fits.Column(name='MSA_METADATA_ID',format='I',array=metaids)
    tabcol3 = fits.Column(name='SHUTTER_QUADRANT',format='I',array=quads)
    tabcol4 = fits.Column(name='SHUTTER_ROW',format='I',array=allrow)
    tabcol5 = fits.Column(name='SHUTTER_COLUMN',format='I',array=allcol)
    tabcol6 = fits.Column(name='SOURCE_ID',format='I',array=sources)
    tabcol7 = fits.Column(name='BACKGROUND',format='A',array=bkgd)
    tabcol8 = fits.Column(name='SHUTTER_STATE',format='4A',array=state)
    tabcol9 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_X',format='E',array=srcx)
    tabcol10 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_Y',format='E',array=srcy)
    hdu2 = fits.BinTableHDU.from_columns([tabcol1, tabcol2, tabcol3, tabcol4, tabcol5, tabcol6, tabcol7, tabcol8,
                                          tabcol9, tabcol10, dither_point_index, primary_source],name='SHUTTER_INFO')

    # SOURCE_INFO table
    # program ID - arbitrary
    program = np.full_like(allcol,1)
    # source name - arbitrary
    name = np.full(len(allcol),'lamp',dtype="<U8")
    # source alias - arbitrary
    alias = np.full(len(allcol),'foo',dtype="<U8")
    # catalog ID - arbitrary
    catalog = np.full(len(allcol),'foo',dtype="<U8")
    # RA, DEC - N/A for ground data
    ra = np.full(len(allcol),0.)
    dec = ra
    # preimage file name - N/A
    preim = np.full(len(allcol),'foo.fits',dtype="<U8")
    # stellarity -- internal lamps are uniform illumination, so set to 0
    # ** if considering a point source, need to change this to 0, or actual value if known
    stellarity = np.full(len(allcol),0.)
    tabcol1 = fits.Column(name='PROGRAM',format='I',array=program)
    tabcol2 = fits.Column(name='SOURCE_ID',format='I',array=sources)
    tabcol3 = fits.Column(name='SOURCE_NAME',format='4A',array=name)
    tabcol4 = fits.Column(name='ALIAS',format='3A',array=alias)
    #tabcol5 = fits.Column(name='CATALOG_ID',format='3A',array=catalog)
    tabcol6 = fits.Column(name='RA',format='D',array=ra)
    tabcol7 = fits.Column(name='DEC',format='D',array=dec)
    tabcol8 = fits.Column(name='PREIMAGE_ID',format='8A',array=preim)
    tabcol9 = fits.Column(name='STELLARITY',format='E',array=stellarity)
    hdu3 = fits.BinTableHDU.from_columns([tabcol1,tabcol2,tabcol3,tabcol4,tabcol5,tabcol6,tabcol7,tabcol8,tabcol9],name='SOURCE_INFO')

    # create fits file
    hdu0 = fits.PrimaryHDU()
    # add necessary keywords to primary header
    hdr = hdu0.header
    hdr.set('ORIGIN', 'STScI', 'institution responsible for creating FITS file')
    hdr.set('TELESCOP', 'JWST', 'telescope used to acquire data')
    hdr.set('INSTRUME', 'NIRSPEC', 'identifier for instrument used to acquire data')
    hdr.set('DATE', now.isoformat())
    hdr.set('FILENAME', outfile, 'name of file')
    hdr.set('PPSDBVER', 'PPSDB999', 'version of PPS database used') # N/A for non-OSS ground data, using arbitrary number
    hdr.set('PROGRAM', '1', 'program number') # arbitrary
    hdr.set('VISIT', '1', 'visit number') # arbitrary
    hdr.set('OBSERVTN', '1', 'observation number') # arbitrary
    hdr.set('VISIT_ID', '1', 'visit identifier') # arbitrary
    #hdr.set('PNTG_SEQ', 1, 'pointing sequence number') # arbitrary
    hdr.set('MSACFG10', 1, 'base 10 nirspec msa_at_pointing.msa_config_id') # arbitrary
    hdr.set('MSACFG36', '01', 'base 36 version of MSACFG10') #arbitrary
    hdu1 = fits.ImageHDU(image,name='SHUTTER_IMAGE')
    hdu_all = fits.HDUList([hdu0,hdu1,hdu2,hdu3])
    hdu_all.writeto(outfile)


def fix_metafile(infile):
    """
    This function fixes the shutter configuration files for the new format (of build 7.3)
    Args:
        infile: string, old MSA configuration fits

    Returns:
        The msa shutter configuration file with the right format for the pipeline
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


def create_metafile_sim(config_binary_file, targ_file_list):
    """
    This function creates the shutter configuration file for the simulations, using the nod set csv files
    Args:
        config_binary_file: string, MSA configuration fits binary table
        fix_old_config_file: boolean

    Returns:
        The msa shutter configuration file with the right format for the pipeline
    """
    # read in configuration (binary fits format)
    hdul = fits.open(config_binary_file)
    # name of output metafile
    if '.msa.fits' in config_binary_file:
        outfile = config_binary_file.replace('.msa.fits', '') + '_metafile_msa.fits'
    else:
        outfile = config_binary_file[:-5] + '_metafile_msa.fits'
    # check if outfile already exists, delete if it does
    if os.path.isfile(outfile):
        os.remove(outfile)

    # get number of nods from the number of input target files
    targ_files = targ_file_list.split(',')
    nnods = len(targ_files)
    print('# nods =', nnods)

    # read in the csv target files
    targ_files = targ_file_list.split(',')
    iter = True
    for targfile in targ_files:
        print('extracting target info from:', targfile)
        targinfo = ascii.read(targfile)
        if iter:
            source_id = targinfo['ID'].data
            source_ra = targinfo['Source RA (Degrees)'].data
            source_dec = targinfo['Source Dec (Degrees)'].data
            #		source_ra = targinfo['SourceRA(Degrees)'].data
            #		source_dec = targinfo['SourceDec(Degrees)'].data
            source_quad = targinfo['Quadrant'].data
            source_row = targinfo['Row'].data
            source_col = targinfo['Column'].data
            source_xpos = targinfo['Offset (x)'].data
            source_ypos = targinfo['Offset (y)'].data
            #		source_xpos = targinfo['Offset(x)'].data
            #		source_ypos = targinfo['Offset(y)'].data
            iter = False
        else:
            source_id = np.concatenate(([source_id, targinfo['ID'].data]))
            source_ra = np.concatenate(([source_ra, targinfo['Source RA (Degrees)'].data]))
            source_dec = np.concatenate(([source_dec, targinfo['Source Dec (Degrees)'].data]))
            #                source_ra = np.concatenate(([source_ra, targinfo['SourceRA(Degrees)'].data]))
            #                source_dec = np.concatenate(([source_dec, targinfo['SourceDec(Degrees)'].data]))
            source_quad = np.concatenate(([source_quad, targinfo['Quadrant'].data]))
            source_row = np.concatenate(([source_row, targinfo['Row'].data]))
            source_col = np.concatenate(([source_col, targinfo['Column'].data]))
            source_xpos = np.concatenate(([source_xpos, targinfo['Offset (x)'].data]))
            source_ypos = np.concatenate(([source_ypos, targinfo['Offset (y)'].data]))
    #                source_xpos = np.concatenate(([source_xpos, targinfo['Offset(x)'].data]))
    #                source_ypos = np.concatenate(([source_ypos, targinfo['Offset(y)'].data]))

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
    tab_slitlet_id = np.full(allcols.size, 0)
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

    nsources = len(source_id)
    nshut = len(quads)
    # nod pattern - table rows that contain the source
    if nnods == 1:
        nod = [0]
    if nnods == 2:
        nod = [0, 3]
    if nnods == 3:
        nod = [1, 3, 8]
    if nnods > 3:
        print("nods > 3 not currently supported")
        sys.exit()

    for i in range(nshut):
        #	print(quads[i],allrows[i],allcols[i])
        match = np.intersect1d(np.where(quads[i] == source_quad),
                               np.intersect1d(np.where(allrows[i] == source_row), np.where(allcols[i] == source_col)))
        if match.size != 0:
            #		print(source_quad[match],source_row[match],source_col[match],source_id[match],source_xpos[match])
            tab_slitlet_id[i] = int(i / nnods / nnods) + 1
            #		print(tab_slitlet_id[i])
            tab_source_id[i] = source_id[match]
            pos = np.mod(i, nnods * nnods)
            if pos in nod:
                # some values only get filled if the source is in that shutter for that nod
                tab_background[i] = 'N'
                tab_source_xpos[i] = source_xpos[match]
                tab_source_ypos[i] = source_ypos[match]
                tab_primary_source[i] = 'Y'
            tab_ditherpoint_index[i] = np.mod(i, nnods) + 1
        else:
            print("no source corresponding to shutter:", quads[i], allrows[i], allcols[i])
    # make a mask to ignore slitlets without sources (shouldn't normally happen)
    mask = tab_slitlet_id != 0
    # create and fill in fits table columns
    tabcol1 = fits.Column(name='SLITLET_ID', format='I', array=tab_slitlet_id[mask])
    tabcol2 = fits.Column(name='MSA_METADATA_ID', format='I', array=tab_metadata_id[mask])
    tabcol3 = fits.Column(name='SHUTTER_QUADRANT', format='I', array=quads[mask])
    tabcol4 = fits.Column(name='SHUTTER_ROW', format='I', array=allrows[mask])
    tabcol5 = fits.Column(name='SHUTTER_COLUMN', format='I', array=allcols[mask])
    tabcol6 = fits.Column(name='SOURCE_ID', format='J', array=tab_source_id[mask])
    tabcol7 = fits.Column(name='BACKGROUND', format='A', array=tab_background[mask])
    tabcol8 = fits.Column(name='SHUTTER_STATE', format='4A', array=tab_shutter_state[mask])
    tabcol9 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_X', format='E', array=tab_source_xpos[mask])
    tabcol10 = fits.Column(name='ESTIMATED_SOURCE_IN_SHUTTER_Y', format='E', array=tab_source_ypos[mask])
    tabcol11 = fits.Column(name='DITHER_POINT_INDEX', format='I', array=tab_ditherpoint_index[mask])
    tabcol12 = fits.Column(name='PRIMARY_SOURCE', format='1A', array=tab_primary_source[mask])
    hdu2 = fits.BinTableHDU.from_columns(
        [tabcol1, tabcol2, tabcol3, tabcol4, tabcol5, tabcol6, tabcol7, tabcol8, tabcol9, tabcol10, tabcol11, tabcol12],
        name='SHUTTER_INFO')

    # SOURCE_INFO table
    # collect all unique instances of the sources
    unique_source_id, match, match2 = np.intersect1d(source_id, tab_source_id[mask], return_indices=True)
    nsources = len(unique_source_id)
    # program ID - arbitrary
    program = np.full(nsources, '111', dtype="<U8")
    # source name - arbitrary
    source_name = np.core.defchararray.add('111_', unique_source_id.astype(str))
    # source alias - arbitrary
    alias = np.full(nsources, 'foo', dtype="<U8")
    # RA, DEC
    ra = source_ra[match]
    dec = source_dec[match]
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


def run_create_metafile(config_binary_file, fix_old_config_file, targ_file_list):
    """
    This function is a wrapper for all cases.
    Args:
        config_binary_file: string, MSA configuration fits binary table
        fix_old_config_file: boolean
        targ_file_list: string, nod set of .csv files divided by commas (no spaces)

    Returns:
        The msa shutter configuration file with the right format for the pipeline
    """
    if targ_file_list is None:
        if not fix_old_config_file:
            create_metafile(config_binary_file)
        else:
            fix_metafile(config_binary_file)
    else:
        create_metafile_sim(config_binary_file, targ_file_list)



if __name__ == '__main__':

    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("config_binary_file",
                        action='store',
                        default=None,
                        help='Name MSA configuration fits binary table, i.e. CB10-GD-B.msa.fits')
    parser.add_argument("-s",
                        dest="targ_file_list",
                        action='store',
                        default=None,
                        help='For simulations only, use this flag to provide csv output from APT; list each one in a nod set, e.g. -s=obs1.csv,obs1.csv')
    parser.add_argument("-f",
                        dest="fix_old_config_file",
                        action='store_true',
                        default=False,
                        help='If an old version of the shutter configuration file exists, use the -f flag to fix it. In this case, the input to the script is the old config file.fits.')
    args = parser.parse_args()

    # Set the variables
    config_binary_file = args.config_binary_file
    targ_file_list = args.targ_file_list
    fix_old_config_file = args.fix_old_config_file

    # Run the function
    run_create_metafile(config_binary_file, fix_old_config_file, targ_file_list)

    print("\nScript create_metafile.py done.\n")


