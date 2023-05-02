#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 12:50:23 2018

@author: gkanarek
"""

import sys, os
import argparse
from jwst.datamodels import open as dmodel_open
from astropy.io import fits
from astropy.time import Time
from crds.matches import find_match_paths_as_dict as ref_matches
from crds import getrecommendations


def check_meta(input_file, match_key, match_val):
    input_val = input_file[match_key.lower()]
    # direct comparison for numeric types; "in" comparison for string types
    if match_val in ["GENERIC", "ALL", "N/A", "ANY"]:
        return True
    if isinstance(input_val, str):
        return input_val in match_val
    return input_val == type(input_val)(match_val)  # have to make sure because match_val is always a string


def get_streams(logfile=None):
    """
    Utility function to identify the output streams.
    """
    if logfile is None:
        errstream = sys.stderr
        logstream = sys.stdout
    else:
        errstream = logstream = open(logfile, 'a')
    return logstream, errstream


def load_input_file(path_to_input_file, logstream=None):
    """
    Utility function as a wrapper around dmodel_open.
    """
    if not isinstance(path_to_input_file, str):
        return path_to_input_file
    if logstream is not None:
        print("Loading input file...", file=logstream)
    try:
        input_file = dmodel_open(path_to_input_file)
    except ValueError:
        import pdb; pdb.set_trace()
    
    return input_file


def reffile_test(path_to_input_file, pipeline_step, logfile=None,
                 input_file=None):
    """
    This is a new version of reffile_test which uses crds.matches instead of
    working with the reference file metadata directly. That way, if the rmap
    was updated manually on CRDS (to avoid redelivering files for a minor
    keyword change), this will test the actual match criteria.
    """
    log_msgs = []

    logstream, errstream = get_streams(logfile=logfile)

    # Convert pipeline step to a header keyword if necessary
    if pipeline_step.upper().startswith("R_"):
        step_key = pipeline_step.upper()
    else:
        if len(pipeline_step) >= 6:
            step_key = "R_" + pipeline_step.upper()[:6]
        else:
            step_key = "R_" + pipeline_step.upper()

    # Identify the context
    context = fits.getval(path_to_input_file, "CRDS_CTX")

    # Identify the reference file
    try:
        reffile_name = fits.getval(path_to_input_file, step_key)
    except KeyError:
        print("Invalid pipeline step", file=errstream)
        log_msgs.append("Invalid pipeline step")
        return None

    reffile_name = reffile_name.replace('crds://', '')

    # Is there a reference file for this step? If not, PASS
    if reffile_name == "N/A":
        print("No reference file for step {}.".format(pipeline_step), file=errstream)
        log_msgs.append("No reference file for step {}.".format(pipeline_step))
        return ""

    # Grab metadata from the input and reference files
    if input_file is None:
        input_file = load_input_file(path_to_input_file, logstream=logstream)
    print("Grabbing CRDS match criteria...", file=logstream)
    log_msgs.append("Grabbing CRDS match criteria...")
    try:
        match_criteria = ref_matches(context, reffile_name)
        if match_criteria:
            match_criteria = match_criteria[0]
        else:
            msg = 'Unable to find match criteria... skipping test'
            log_msgs.append(msg)
            return None
    except ValueError:
        import pdb
        pdb.set_trace()

    tests = {}  # store all the tests in a single dictionary

    # add instrument name in the expected keyword
    match_criteria['META.INSTRUMENT.NAME'] = 'NIRSPEC'

    # make sure that the subarray keyword is correct for the size of the data
    # subarray = get_subarray(path_to_input_file)
    subarray = fits.getval(path_to_input_file, 'SUBARRAY', 0)
    match_criteria['META.SUBARRAY.NAME'] = subarray

    # Test whether the recommended reference file was actually selected
    recommended_reffile = getrecommendations(match_criteria,
                                             reftypes=[pipeline_step],
                                             context=context,
                                             fast=True)

    if isinstance(recommended_reffile, str):
        recommended_reffile = os.path.basename(recommended_reffile)  # remove path, only want to test filename
        tests['RECOMMENDATION'] = recommended_reffile == reffile_name
    else:
        msg1 = '* WARNING: Unable to find recommendation for the reference file:'
        msg2 = '        Match criteria determined by pipeline to find reference file: '+repr(match_criteria)
        msg3 = '        Recommendation dictionary = '+repr(recommended_reffile)
        log_msgs.append(msg1)
        log_msgs.append(msg2)
        log_msgs.append(msg3)
        print(msg1)
        print(msg2)
        print(msg3)

    # Remove irrelevant match criteria
    del match_criteria['observatory']
    del match_criteria['instrument']
    del match_criteria['filekind']

    # Useafter dates require special handling
    if "META.OBSERVATION.DATE" not in match_criteria:
        tests['USEAFTER'] = True
    else:
        input_date = input_file.meta.observation.date
        input_time = input_file.meta.observation.time
        input_obstime = Time(input_date + "T" + input_time)
        ref_date = match_criteria.pop("META.OBSERVATION.DATE")
        ref_time = match_criteria.pop("META.OBSERVATION.TIME")
        ref_useafter = Time(ref_date + "T" + ref_time)
        tests["USEAFTER"] = input_obstime >= ref_useafter
        # Note that this does NOT check whether there is a more recent
        # (but still valid) reference file that could have been selected

    # Loop over the rest of the matching criteria
    for criterion, value in match_criteria.items():
        tests[criterion] = check_meta(input_file, criterion, value)

    final = all([x or x is None for x in tests.values()])

    failures = []
    failmsg = "{}: reffile value {}, input value {}"

    # Finally, print out the results of the tests
    print("REFERENCE FILE SELECTION TEST", file=logstream)
    print("  Input file: {}".format(path_to_input_file), file=logstream)
    print("  Pipeline step: {}".format(pipeline_step), file=logstream)
    print("  Header keyword: {}".format(step_key), file=logstream)
    print("  Reference file selected: {}".format(reffile_name), file=logstream)
    print("  **Metadata tests performed:**", file=logstream)
    log_msgs.append("REFERENCE FILE SELECTION TEST")
    log_msgs.append("  Input file: {}".format(path_to_input_file))
    log_msgs.append("  Pipeline step: {}".format(pipeline_step))
    log_msgs.append("  Header keyword: {}".format(step_key))
    log_msgs.append("  Reference file selected: {}".format(reffile_name))
    log_msgs.append("  **Metadata tests performed:**")
    rescode = {None: "N/A", True: "PASS", False: "FAIL"}
    for meta in sorted(tests):
        result = tests[meta]
        print("    {}: {}".format(meta, rescode[result]), file=logstream)
        if rescode[result] == "FAIL":
            if meta == "USEAFTER":
                ival = input_obstime
                rval = ref_useafter
            else:
                ival = input_file[meta.lower()]
                rval = match_criteria[meta]
            failures.append(failmsg.format(meta, rval, ival))
            print("      Input file value: {}".format(ival), file=logstream)
            print("      Reference file value: {}".format(rval), file=logstream)
            log_msgs.append("      Input file value: {}".format(ival))
            log_msgs.append("      Reference file value: {}".format(rval))

    print("  Final result: {}".format(rescode[final]), file=logstream)
    log_msgs.append("  Final result: {}".format(rescode[final]))

    # Close the output stream if necessary
    if logfile is not None:
        logstream.close()

    return "\n".join(failures), log_msgs


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
    
    if logfile is not None:   # erase existing log, since we'll be appending later
        with open(logfile, 'w'):
            pass

    # Only want to load the input file once
    input_file = load_input_file(path_to_input_file)

    failures = {}

    for step in all_steps:
        res = reffile_test(path_to_input_file, step, logfile=logfile, 
                           input_file=input_file)
        if res:
            failures[step] = res

    return failures


def main():

    parser = argparse.ArgumentParser(description='Test to see if the correct '
                                     'reference file was selected for the given'
                                     ' input file(s).')
    parser.add_argument('input_file', nargs='+', help="Paths to one or more input files")
    parser.add_argument('-l', '--log', help="Path to a desired output log file")

    args = parser.parse_args()

    for input_file in args.input_file:
        check_all_reffiles(input_file, logfile=args.log)


if __name__ == '__main__':
    sys.exit(main())
