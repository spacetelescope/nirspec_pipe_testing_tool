"""
This is a *FAILED* attempt to make the dictionary. BUT it contains all the info.
"""

CV3_testdata_dict["FS"]["PRISM"] = {
    "filter" : "CLEAR",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FFB",
    "NID" : 9852,
    "CV3filename" : ["NRSSLIT-COMBO-077_1_491_SE_2013-01-20T01h09m52.fits",
                     "NRSSLIT-COMBO-077_1_492_SE_2013-01-20T01h10m04.fits"],
    "level1Bfilenames" : ["jwtest1019001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1020001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["PRISM"] = {
    "filter" : "F100LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FFB",
    "NID" : 9845,
    "CV3filename" : ["NRSSLIT-COMBO-070_1_491_SE_2013-01-20T00h44m57.fits",
                     "NRSSLIT-COMBO-070_1_492_SE_2013-01-20T00h45m16.fits"],
    "level1Bfilenames" : ["jwtest1017001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1018001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G140M"] = {
    "filter" : "F070LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : ["CLS/FFV", "CLS/ARGON"],
    "NID" : [9814, 9812],
    "CV3filename" : [["NRSSLIT-COMBO-039_1_491_SE_2013-01-19T22h20m46.fits",
                     "NRSSLIT-COMBO-039_1_492_SE_2013-01-19T22h21m04.fits"],
                     ["NRSSLIT-COMBO-037_1_491_SE_2013-01-19T22h14m57.fits",
                      "NRSSLIT-COMBO-037_1_492_SE_2013-01-19T22h15m20.fits"]],
    "level1Bfilenames" : [["jwtest1011001_01101_00001_NRS1_uncal_mod.fits",
                           "jwtest1012001_01101_00001_NRS2_uncal_mod.fits"],
                          ["jwtest1009001_01101_00001_NRS1_uncal_mod.fits",
                           "jwtest1010001_01101_00001_NRS2_uncal_mod.fits"]]
}
CV3_testdata_dict["FS"]["G140H"] = {
    "filter" : "F070LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/ARGON",
    "NID" : 9776,
    "CV3filename" : ["NRSSLIT-COMBO-001_1_491_SE_2013-01-19T18h27m13.fits",
                     "NRSSLIT-COMBO-001_1_492_SE_2013-01-19T18h27m30.fits"],
    "level1Bfilenames" : ["jwtest1001001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1002001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G140H"] = {
    "filter" : "F100LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FF1",
    "NID" : 9791,
    "CV3filename" : ["NRSSLIT-COMBO-016_1_491_SE_2013-01-19T20h12m14.fits",
                     "NRSSLIT-COMBO-016_1_492_SE_2013-01-19T20h14m05.fits"],
    "level1Bfilenames" : ["jwtest1003001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1004001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G235M"] = {
    "filter" : "F170LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FF2",
    "NID" : 9830,
    "CV3filename" : ["NRSSLIT-COMBO-055_1_491_SE_2013-01-19T23h33m01.fits",
                     "NRSSLIT-COMBO-055_1_492_SE_2013-01-19T23h33m37.fits"],
    "level1Bfilenames" : ["jwtest1013001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1014001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G235H"] = {
    "filter" : "F290LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FF3",
    "NID" : 9799,
    "CV3filename" : ["NRSSLIT-COMBO-024_1_491_SE_2013-01-19T20h54m11.fits",
                     "NRSSLIT-COMBO-024_1_492_SE_2013-01-19T20h54m51.fits"],
    "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1006001_01101_00001_NRS2_uncal_mod.fits"]
}
CV3_testdata_dict["FS"]["G395H"] = {
    "filter" : "F290LP",
    "subarray" : "ALLSLITS",
    "CAA_lamp" : "CLS/FF3",
    "NID" : 9801,
    "CV3filename" : ["NRSSLIT-COMBO-026_1_491_SE_2013-01-19T21h12m50.fits",
                     "NRSSLIT-COMBO-026_1_492_SE_2013-01-19T21h14m03.fits"],
    "level1Bfilenames" : ["jwtest1007001_01101_00001_NRS1_uncal_mod.fits",
                          "jwtest1008001_01101_00001_NRS2_uncal_mod.fits"]
}


# MOS test data
CV3_testdata_dict["MOS"]["G140H"] = {
    "filter" : "OPAQUE",
    "CAA_lamp" : "LINE1",
    "NID" : 41598,
    "CV3filename" : ["NRSV96214001001P0000000002105_1_491_SE_2016-01-24T01h59m01.fits ",
                     "NRSV96214001001P0000000002105_1_492_SE_2016-01-24T01h59m01.fits"],
    "level1Bfilenames" : ["jwtest1015001_01101_00001_NRS1_uncal.fits", "jwtest1016001_01101_00001_NRS1_uncal.fits"],
    "MSA_config" : "V9621400100101"
}

