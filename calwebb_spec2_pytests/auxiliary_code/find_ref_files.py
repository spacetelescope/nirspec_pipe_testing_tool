
"""
This script uses the input raw fits file name to look for the right reference file that the pipeline should have
selected.
"""

def get_rf_list(mode, raw_fits_file):
    '''
    This function reads the text file ref_files.txt and collects the reference files for the given raw file.
    Args:
        mode: str, MOS, FS, or IFU
        raw_fits_file: str, base name of the raw data file

    Returns:
        rawff_reffiles: list, reference files that the pipeline should have used for that file
    '''
    reffilestxt = 'ref_files.txt'
    rawff_reffiles = []
    get_rf, collect_reffiles = False, False
    with open(reffilestxt) as rftxt:
        for line in rftxt.readlines():
            if mode in line:
                get_rf = True
            if get_rf:
                if raw_fits_file in line:
                    print ('Found the raw data file. Collecting the reference files...')
                    collect_reffiles = True
                if collect_reffiles:
                    if 'jwst' in line:
                        rf = line.replace('\n', '')
                        rawff_reffiles.append(rf)
                        if 'sflat' in line:
                            get_rf, collect_reffiles = False, False
    return rawff_reffiles




if __name__ == '__main__':

    # Set variables
    raw_fits_file = 'NRSV84600010001P0000000002101_4_492_SE_2016-01-17T17h34m08.fits'
    mode = 'MOS'

    # get reference files
    rawff_reffiles = get_rf_list(mode, raw_fits_file)