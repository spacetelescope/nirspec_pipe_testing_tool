from astropy.io import fits
import argparse
import subprocess
import sys
import pysiaf
import pysiaf.utils.rotations as rotations

"""
This script is meant to be called after the header has been fixed by either of the keywd_check scripts.

Usage:
    You will need the values for the target's RA, DEC, V2, and V3, as well as the aperture
    position angle: ra_targ, dec_targ, v2_targ, v3_targ, aper_angle

    Sample values are: ra_targ = 53.16199112, dec_targ = -27.79127312, v2_targ = 393.86285,
    v3_targ = -424.00329, aper_angle = 45.0

    - From a terminal type:
     $ nptt_fix_pointing blah.fits 53.16199112 -27.79127312 393.86285 -424.00329 45.0

     If the data is IFU add the flag -ifu at the end of the command. The output will be the updated file.
     To create a new updated file add flag -nf to the above command.
    
    - Within a script:
    import nirspec_pipe_testing_tool as nptt
    stsci_pipe_ready_file = '/somewhere/blah.fits'
    ra_targ = 53.16199112
    dec_targ = -27.79127312
    v2_targ = 393.86285,
    v3_targ = -424.00329
    aper_angle = 45.0
    ifu_used = False   # only True for IFU data
    nptt.utils.fix_pointing.fix_header_pointing(stsci_pipe_ready_file,
                                                RA_target, DEC_target,
                                                v2_target, v3_target,
                                                aper_angle, ifu=ifu_used)
    
"""

# HEADER
__author__ = "M. A. Pena-Guerrero & James Muzerolle"
__version__ = "1.0"


# HISTORY
# Mar 2020 - Version 1.0: initial version completed


def fix_header_pointing(infile, ra_targ, dec_targ, v2_targ, v3_targ, apa=45.0, ifu=False):
    """
        This function takes input FS, MOS, or IFU count rate file and adds/corrects the spacecraft
        pointing header information

    :param infile:
    :param ra_targ: float, RA of the target (the MSA target pointing; from an APT pointing file)
    :param dec_targ: float, DEC of the target (the MSA target pointing; from an APT pointing file)
    :param v2_targ: float, V2 of the target (the MSA target pointing; from an APT pointing file)
    :param v3_targ: float, V3 of the target (the MSA target pointing; from an APT pointing file)
    :param apa: float, aperture position angle [deg], e.g. -999.0
    :param ifu: boolean, True if IFU was used
    :return: input file with fixed pointing
    """
    # use pysiaf to convert MSA fiducial pointing to spacecraft pointing
    NIRSpec = pysiaf.Siaf('NIRSpec')
    if not ifu:
        msa_aper = NIRSpec['NRS_FULL_MSA']
    else:
        msa_aper = NIRSpec['NRS_FULL_IFU']
    v3pa = apa - msa_aper.V3IdlYAngle
    msa_attitude_matrix = rotations.attitude(v2_targ, v3_targ, ra_targ, dec_targ, v3pa)
    ra_tel, dec_tel = rotations.pointing(msa_attitude_matrix, 0., 0.)  # telescope pointing
    ra_msa, dec_msa = rotations.pointing(msa_attitude_matrix, msa_aper.V2Ref, msa_aper.V3Ref)  # MSA pointing

    # load input file
    primary_hdr = fits.getheader(infile, 'PRIMARY')
    sci_hdr = fits.getheader(infile, 'SCI')

    # replace the relevant keywords in the primary header
    for key, _ in primary_hdr.items():
        new_value = None
        # target coords in primary header (probably unneeded)
        if key == 'TARG_RA':
            new_value = ra_targ
            after_key = 'TARGTYPE'
        if key == 'TARG_DEC':
            new_value = dec_targ
            after_key = 'TARGTYPE'
        if key == 'PROP_RA':
            new_value = ra_targ
            after_key = 'MU_EPOCH'
        if key == 'PROP_DEC':
            new_value = dec_targ
            after_key = 'PROP_RA'

        # modify keyword value in input file
        if new_value is not None:
            fits.setval(infile, key, 0, value=new_value, after=after_key)

    # replace the relevant keywords in the primary header
    for key, _ in sci_hdr.items():
        new_value = None
        # SCI header
        if key == 'PA_APER':
            new_value = apa
            after_key = 'DEC_V1'
        if key == 'RA_V1':
            new_value = ra_tel
            after_key = 'PA_V3'
        if key == 'DEC_V1':
            new_value = dec_tel
            after_key = 'RA_V1'
        if key == 'PA_V3':
            new_value = v3pa
            after_key = 'COORDSYS'
        if key == 'RA_REF':
            new_value = ra_msa
            after_key = 'V3I_YANG'
        if key == 'DEC_REF':
            new_value = dec_msa
            after_key = 'RA_REF'
        if key == 'ROLL_REF':
            new_value = v3pa
            after_key = 'DEC_REF'

        # the following 3 are here just in case the default is not populated already
        if key == 'V2_REF':
            new_value = msa_aper.V2Ref
            after_key = 'SPORDER'
        if key == 'V3_REF':
            new_value = msa_aper.V3Ref
            after_key = 'V2_REF'
        if key == 'V3I_YANG':
            new_value = msa_aper.V3IdlYAngle
            after_key = 'VPARITY'

        # modify keyword value in input file
        if new_value is not None:
            fits.setval(infile, key, 1, value=new_value, after=after_key)


def main():
    # Get arguments to run script
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("input_fits_file",
                        action='store',
                        default=None,
                        help='Name of input fits file, i.e. blah.fits')
    parser.add_argument("-ifu",
                        dest="ifu_used",
                        action='store_true',
                        default=False,
                        help='Use flag -ifu if observation mode used was IFU.')
    parser.add_argument("ra_targ",
                        action='store',
                        default=None,
                        help='RA MSA pointing from APT file')
    parser.add_argument("dec_targ",
                        action='store',
                        default=None,
                        help='DEC MSA pointing from APT file')
    parser.add_argument("v2_targ",
                        action='store',
                        default=None,
                        help='V2 MSA pointing from APT file')
    parser.add_argument("v3_targ",
                        action='store',
                        default=None,
                        help='V3 MSA pointing from APT file')
    parser.add_argument("apa",
                        action='store',
                        default=None,
                        help='Aperture position angle [deg], e.g. 45.0')
    parser.add_argument("-nf",
                        dest="new_file",
                        action='store_true',
                        default=False,
                        help='Use -nf to create a new file with updated header.')
    args = parser.parse_args()

    # Set the variables
    input_fits_file = args.input_fits_file
    ifu_used = args.ifu_used
    ra_targ = float(args.ra_targ)
    dec_targ = float(args.dec_targ)
    v2_targ = float(args.v2_targ)
    v3_targ = float(args.v3_targ)
    apa = float(args.apa)
    new_file = args.new_file

    if new_file:
        input_file = input_fits_file.replace('.fits', '_updated_pointing.fits')
        subprocess.run(["cp", input_fits_file, input_file])
    else:
        input_file = input_fits_file

    # Fix the pointing in the file
    fix_header_pointing(input_file, ra_targ, dec_targ, v2_targ, v3_targ, apa, ifu=ifu_used)

    print('\n * Script  fix_pointing.py  finished * \n')


if __name__ == '__main__':
    sys.exit(main())
