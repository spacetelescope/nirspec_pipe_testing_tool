#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 12:50:23 2018

@author: gkanarek
"""

import sys
from jwst.datamodels import open as dmodel_open
from astropy.io import fits
from astropy.time import Time

#List of metadata to compare:
comparison_meta = {
        'exp_type': 'EXPOSURE.TYPE',
        'detector': 'INSTRUMENT.DETECTOR',
        'obs_filter': 'INSTRUMENT.FILTER',
        'grating': 'INSTRUMENT.GRATING',
        'instrument': 'INSTRUMENT.NAME',
        'subarray': 'SUBARRAY.NAME',
        'sub_xsize': 'SUBARRAY.XSIZE',
        'sub_xstart': 'SUBARRAY.XSTART',
        'sub_ysize': 'SUBARRAY.YSIZE',
        'sub_ystart': 'SUBARRAY.YSTART',
        'reftype': 'REFTYPE'}

def check_meta(input_meta, ref_meta, meta):
    try:
        in_val = input_meta[meta]
        ref_val = ref_meta[meta]
    except (KeyError, AttributeError):
        return None
    if ref_val == "N/A":
        return None
    #direct comparison for numeric types; "in" comparison for string types
    if isinstance(in_val, str):
        return in_val in ref_val
    return in_val == ref_val
    
def reffile_test(path_to_input_file, pipeline_step,
                        logfile=None,
                        useafter=True, 
                        readpattern=True,
                        exp_type=True,
                        detector=True,
                        obs_filter=True,
                        grating=True,
                        instrument=True,
                        subarray=True,
                        sub_xsize=True,
                        sub_ysize=True,
                        sub_xstart=True,
                        sub_ystart=True,
                        reftype=True,
                        **other_meta
                    ):
    """
    This is a function which compares the metadata values from the comparison 
    list (default being the "BESTREFS header" from the pipeline) between an
    input file and the reference file which was used for a given pipeline step.
    
    To omit any of the metadata tests from the default, set the corresponding
    keyword value to False when calling the function.
    
    To test any metadata not already included in the default list, include
    test names and metadata strings for each as keyword arguments.
    """
    
    #identify output streams
    if logfile is None:
        errstream = sys.stderr
        logstream = sys.stdout
    else:
        errstream = logstream = open(logfile, 'a')
    
    #Convert pipeline step to a header keyword if necessary
    if pipeline_step.upper().startswith("R_"):
        step_key = pipeline_step.upper()
    else:
        if len(pipeline_step) >= 6:
            step_key = "R_" + pipeline_step.upper()[:6]
        else:
            step_key = "R_" + pipeline_step.upper()
    
    #Identify the reference file
    try:
        reffile_name = fits.getval(path_to_input_file, step_key)
    except KeyError:
        print("Invalid pipeline step", file=errstream)
        return None
    reffile_path = reffile_name.replace('crds://', 
                                        '/grp/crds/jwst/references/jwst/')
    
    #Is there a reference file for this step? If not, PASS
    if reffile_path == "N/A":
        print("No reference file for this pipeline step.", file=errstream)
        return ""
    
    #Grab metadata from the input and reference files
    print("Loading input file...", file=logstream)
    try:
        input_file_meta = dmodel_open(path_to_input_file).meta
    except ValueError:
        import pdb; pdb.set_trace()
    print("Loading reference file...", file=logstream)
    try:
        reffile_meta = dmodel_open(reffile_path).meta
    except ValueError:
        import pdb; pdb.set_trace()
    
    #Go through the various check keywords
    keywords = {'exp_type':exp_type,
                'detector':detector,
                'obs_filter':obs_filter,
                'grating':grating,
                'instrument':instrument,
                'subarray':subarray,
                'sub_xsize':sub_xsize,
                'sub_ysize':sub_ysize,
                'sub_xstart':sub_xstart,
                'sub_ystart':sub_ystart,
                'reftype':reftype}
    
    tests = {}
    
    if useafter: #dates require special handling
        if reffile_meta.useafter is None:
            tests['useafter'] = True
        else:
            input_date = input_file_meta.observation.date
            input_time = input_file_meta.observation.time
            input_obstime = Time(input_date + "T" + input_time)
            ref_useafter = Time(reffile_meta.useafter)
            
            tests['useafter'] = input_obstime >= ref_useafter
        #Note that this does NOT check whether there is a more recent
        #(but still valid) reference file that could have been selected
    
    if readpattern: #read pattern requires special handling
        try: #first test P_READPA, then try READPATT if that doesn't work
            input_pattern = input_file_meta.exposure.readpatt
            ref_pattern = reffile_meta.exposure.preadpatt
            if ref_pattern is None:
                ref_pattern = reffile_meta.exposure.readpatt
        except (KeyError, AttributeError):
            tests['readpattern'] = None
        else:
            if ref_pattern == "N/A":
                tests['readpattern'] = None
            else:
                if "ALL" in ref_pattern:
                    tests['readpattern'] = ("IRS2" in ref_pattern) == ("IRS2" in input_pattern)
                else:
                    tests['readpattern'] = input_pattern in ref_pattern
    
    for meta in keywords:
        if keywords[meta]:
            try:
                meta_key = comparison_meta[meta].lower()
            except KeyError:
                tests[meta] = None
            else:
                tests[meta] = check_meta(input_file_meta, reffile_meta, meta_key)
        
    
    #Next, check any other pieces of metadata which were passed
    for meta, meta_key in other_meta:
        tests[meta] = check_meta(input_file_meta, reffile_meta, meta_key)
    
    final = all([x or x is None for x in tests.values()])
    
    failures = []
    failmsg = "{}: reffile value {}, input value {}"
    
    #Finally, print out the results of the tests
    print("REFERENCE FILE SELECTION TEST", file=logstream)
    print("  Input file: {}".format(path_to_input_file), file=logstream)
    print("  Pipeline step: {}".format(pipeline_step), file=logstream)
    print("  Header keyword: {}".format(step_key), file=logstream)
    print("  Reference file selected: {}".format(reffile_name), file=logstream)
    print("  **Metadata tests performed:**", file=logstream)
    rescode = {None: "N/A", True: "PASS", False: "FAIL"}
    for meta in sorted(tests):
        result = tests[meta]
        print("    {}: {}".format(meta, rescode[result]), file=logstream)
        if rescode[result] == "FAIL":
            if meta == "readpattern":
                ival = input_file_meta.exposure.readpatt
                rval = reffile_meta.exposure.preadpatt or reffile_meta.exposure.readpatt
            elif meta == "useafter":
                ival = input_obstime
                rval = ref_useafter
            else:
                ival = input_file_meta[comparison_meta[meta].lower()]
                rval = reffile_meta[comparison_meta[meta].lower()]
            failures.append(failmsg.format(meta, rval, ival))
            print("      Input file value: {}".format(ival), file=logstream)
            print("      Reference file value: {}".format(rval), file=logstream)
    print("  Final result: {}".format(rescode[final]), file=logstream)
    
    #Close the output stream if necessary
    if logfile is not None:
        logstream.close()
    
    return "\n".join(failures)

def create_rfile_test(step, doc_insert):
    """
    A factory to create wrappers for testing correct reference files.
    """
    
    def rfile_test_step(output_hdul):
        output_file = output_hdul[1]
        return reffile_test(output_file, step)
    
    rfile_test_step.__doc__ = """
    This function determines if the reference file for the {} matches the expected one.
    Args:
        output_hdul: output from the output_hdul function

    Returns:
        result: boolean, true if the reference file matches expected value
    """.format(doc_insert)
    
    return rfile_test_step

def check_all_reffiles(path_to_input_file, logfile=None):
    """
    A wrapper around reffile_test to test every reference file in the input
    file's header. A file path may be included to redirect output to a log.
    """
    
    all_steps = list(fits.getval(path_to_input_file, "R_*"))
    
    if logfile is not None: #erase existing log, since we'll be appending later
        with open(logfile, 'w'):
            pass
    
    for step in all_steps:
        reffile_test(path_to_input_file, step, logfile=logfile)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Test to see if the correct '
                                     'reference file was selected for the given'
                                     ' input file(s).')
    parser.add_argument('input_file', nargs='+', help="Paths to one or more input files")
    parser.add_argument('-l', '--log', help="Path to a desired output log file")
    
    args = parser.parse_args()
    
    for input_file in args.input_file:
        check_all_reffiles(input_file, logfile=args.log)
    