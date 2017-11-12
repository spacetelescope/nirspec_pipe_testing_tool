"""
This script contains a dictionary of the CV3 data used for testing build 7, as well
as the corresponding ESA intermediary products.
"""

CV3_testdata_dict = {}

# Fixed Slit test data
# ALLSLITS test data, this is in subarray mode 2048x256
CV3_testdata_dict["FS"] = {
    "NID" : {
        "9852":{
            "grism" : "PRISM",
            "filter" : "CLEAR",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/FFB",
            "CV3filename" : ["NRSSLIT-COMBO-077_1_491_SE_2013-01-20T01h09m52.fits",
                             "NRSSLIT-COMBO-077_1_492_SE_2013-01-20T01h10m04.fits"],
            "level1Bfilenames" : ["jwtest1019001_01101_00001_NRS1_uncal_mod.fits",
                                  "jwtest1020001_01101_00001_NRS2_uncal_mod.fits"]
        },
        "9845" : {
            "grism" : "PRISM",
            "filter" : "F100LP",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/FFB",
            "CV3filename" : ["NRSSLIT-COMBO-070_1_491_SE_2013-01-20T00h44m57.fits",
                             "NRSSLIT-COMBO-070_1_492_SE_2013-01-20T00h45m16.fits"],
            "level1Bfilenames" : ["jwtest1017001_01101_00001_NRS1_uncal_mod.fits",
                                  "jwtest1018001_01101_00001_NRS2_uncal_mod.fits"]
        },
        "9814": {
            "grism" : "G140M",
            "filter" : "F070LP",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/FFV",
            "CV3filename" : ["NRSSLIT-COMBO-039_1_491_SE_2013-01-19T22h20m46.fits",
                             "NRSSLIT-COMBO-039_1_492_SE_2013-01-19T22h21m04.fits"],
            "level1Bfilenames" : ["jwtest1011001_01101_00001_NRS1_uncal_mod.fits",
                                   "jwtest1012001_01101_00001_NRS2_uncal_mod.fits"],
        },
        "9812" :{
            "grism" : "G140M",
            "filter" : "F070LP",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/ARGON",
            "CV3filename" : ["NRSSLIT-COMBO-037_1_491_SE_2013-01-19T22h14m57.fits",
                              "NRSSLIT-COMBO-037_1_492_SE_2013-01-19T22h15m20.fits"],
            "level1Bfilenames" : ["jwtest1009001_01101_00001_NRS1_uncal_mod.fits",
                                   "jwtest1010001_01101_00001_NRS2_uncal_mod.fits"]
        },
        "9776" : {
            "grism" : "G140H",
            "filter" : "F070LP",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/ARGON",
            "CV3filename" : ["NRSSLIT-COMBO-001_1_491_SE_2013-01-19T18h27m13.fits",
                             "NRSSLIT-COMBO-001_1_492_SE_2013-01-19T18h27m30.fits"],
            "level1Bfilenames" : ["jwtest1001001_01101_00001_NRS1_uncal_mod.fits",
                                  "jwtest1002001_01101_00001_NRS2_uncal_mod.fits"]
        },
        "9791" : {
            "grism" : "G140H",
            "filter" : "F100LP",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/FF1",
            "CV3filename" : ["NRSSLIT-COMBO-016_1_491_SE_2013-01-19T20h12m14.fits",
                             "NRSSLIT-COMBO-016_1_492_SE_2013-01-19T20h14m05.fits"],
            "level1Bfilenames" : ["jwtest1003001_01101_00001_NRS1_uncal_mod.fits",
                                  "jwtest1004001_01101_00001_NRS2_uncal_mod.fits"]
        },
        "9830" : {
            "grism" : "G235M",
            "filter" : "F170LP",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/FF2",
            "CV3filename" : ["NRSSLIT-COMBO-055_1_491_SE_2013-01-19T23h33m01.fits",
                             "NRSSLIT-COMBO-055_1_492_SE_2013-01-19T23h33m37.fits"],
            "level1Bfilenames" : ["jwtest1013001_01101_00001_NRS1_uncal_mod.fits",
                                  "jwtest1014001_01101_00001_NRS2_uncal_mod.fits"]

        },
        "9799" : {
            "grism" : "G235H",
            "filter" : "F290LP",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/FF3",
            "CV3filename" : ["NRSSLIT-COMBO-024_1_491_SE_2013-01-19T20h54m11.fits",
                             "NRSSLIT-COMBO-024_1_492_SE_2013-01-19T20h54m51.fits"],
            "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal_mod.fits",
                                  "jwtest1006001_01101_00001_NRS2_uncal_mod.fits"]
       },
        "9801" : {
            "grism" : "G395H",
            "filter" : "F290LP",
            "subarray" : "ALLSLITS",
            "CAA_lamp" : "CLS/FF3",
            "CV3filename" : ["NRSSLIT-COMBO-026_1_491_SE_2013-01-19T21h12m50.fits",
                             "NRSSLIT-COMBO-026_1_492_SE_2013-01-19T21h14m03.fits"],
            "level1Bfilenames" : ["jwtest1007001_01101_00001_NRS1_uncal_mod.fits",
                                  "jwtest1008001_01101_00001_NRS2_uncal_mod.fits"]
        }
    }
}

# MOS test data
CV3_testdata_dict["MOS"] = {
    "NID":{
        "41598": {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "CAA_lamp" : "LINE1",
            "CV3filename" : ["NRSV96214001001P0000000002105_1_491_SE_2016-01-24T01h59m01.fits ",
                             "NRSV96214001001P0000000002105_1_492_SE_2016-01-24T01h59m01.fits"],
            "level1Bfilenames" : ["jwtest1015001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1016001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "V9621400100101"
        },
        #"37286" : {
        #    "grism" : "G140H",
        #    "filter" : "OPAQUE",
        #    "CAA_lamp" : "LINE1",
        #    "CV3filename" : ["NRSV00300060001P0000000002103_1_491_SE_2016-01-06T04h04m49.fits",
        #                     "NRSV00300060001P0000000002103_1_492_SE_2016-01-06T04h04m49.fits"],
        #    "level1Bfilenames" : ["", ""],
        #    "MSA_config" : "V0030006000101"
        #},
        "39547" : {
            "grism" : "G140M",
            "filter" : "OPAQUE",
            "CAA_lamp" : "LINE1",
            "CV3filename" : ["NRSV84600010001P0000000002101_4_491_SE_2016-01-17T17h34m08.fits",
                             "NRSV84600010001P0000000002101_4_492_SE_2016-01-17T17h34m08.fits"],
            "level1Bfilenames" : ["jwtest1001001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1002001_01101_00001_NRS2_uncal.fits"],
            "MSA_config" : "V8460001000101"
        },
        "39553" : {
            "grism" : "G235M",
            "filter" : "OPAQUE",
            "CAA_lamp" : "LINE2",
            "CV3filename" : ["NRSV84600011001P0000000002101_2_491_SE_2016-01-17T18h18m48.fits",
                             "NRSV84600011001P0000000002101_2_492_SE_2016-01-17T18h18m48.fits"],
            "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1006001_01101_00001_NRS2_uncal.fits"],
            "MSA_config" : "V8460001100101"
        },
        "41543" : {
            "grism" : "G395M",
            "filter" : "OPAQUE",
            "CAA_lamp" : "LINE3",
            "CV3filename" : ["NRSV96215001001P0000000002103_1_491_SE_2016-01-24T01h25m07.fits",
                             "NRSV96215001001P0000000002103_1_492_SE_2016-01-24T01h25m07.fits"],
            "level1Bfilenames" : ["jwtest1010001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1009001_01101_00001_NRS2_uncal.fits"],
            "MSA_config" : "V9621500100101"
        },
        "37328" : {
            "grism" : "PRISM",
            "filter" : "OPAQUE",
            "CAA_lamp" : "LINE4",
            "CV3filename" : ["NRSV00300060001P000000000210T_1_491_SE_2016-01-06T06h27m34.fits",
                             "NRSV00300060001P000000000210T_1_492_SE_2016-01-06T06h27m34.fits"],
            "level1Bfilenames" : ["jwtest1013001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1014001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "V0030006000104"
        }
    }
}

# IFU test data
CV3_testdata_dict["IFU"] = {
    "NID":{
        "37668" : {
            "grism" : "PRISM",
            "filter" : "CLEAR",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSSIMA-QUAL-04-B-6007022859_1_491_SE_2016-01-07T02h37m13.fits",
                              "NRSSIMA-QUAL-04-B-6007022859_1_492_SE_2016-01-07T02h37m13.fits"],
            "level1Bfilenames" : ["jwtest1001001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1002001_01101_00001_NRS2_uncal.fits"],
            "notes" : ["no external source, use for subtraction for NID37669,msa_config=ARDCLOSED",
                       "very bright external source, use with NID37668,msa_config=ARDCLOSED"],
        },
        "37669" : {
            "grism" : "PRISM",
            "filter" : "CLEAR",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSSIMA-QUAL-04-2-6007023323_1_491_SE_2016-01-07T02h41m22.fits",
                             "NRSSIMA-QUAL-04-2-6007023323_1_492_SE_2016-01-07T02h41m22.fits"],
            "level1Bfilenames" : ["jwtest1003001_01101_00001_NRS2_uncal.fits",
                                  "jwtest1004001_01101_00001_NRS2_uncal.fits"],
            "notes" : ["no external source, use for subtraction for NID37669,msa_config=ARDCLOSED",
                       "very bright external source, use with NID37668,msa_config=ARDCLOSED"],
        },
        "30192" : {
            "grism" : "G140M",
            "filter" : "OPAQUE",
            "CAA_lamp" : "LINE1",
            "CV3filename" : ["NRSSMOS-MOD-G1M-17-5344175105_1_491_SE_2015-12-10T18h00m06.fits",
                             "NRSSMOS-MOD-G1M-17-5344175105_1_492_SE_2015-12-10T18h00m05.fits"],
            "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1006001_01101_00001_NRS1_uncal.fits"],
            "notes" : ""
        },
        "30227" : {
            "grism" : "G235M",
            "filter" : "OPAQUE",
            "CAA_lamp" : "LINE2",
            "CV3filename" : ["NRSSMOS-MOD-G2M-17-5344211451_1_491_SE_2015-12-10T21h23m36.fits",
                             "NRSSMOS-MOD-G2M-17-5344211451_1_492_SE_2015-12-10T21h23m37.fits"],
            "level1Bfilenames" : ["jwtest1007001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1008001_01101_00001_NRS1_uncal.fits"],
            "notes" : ""
        },
        "30273" : {
            "grism" : "G395M",
            "filter" : "OPAQUE",
            "CAA_lamp" : "LINE3",
            "CV3filename" : ["NRSSMOS-MOD-G3M-17-5345014854_1_491_SE_2015-12-11T01h57m10.fits",
                             "NRSSMOS-MOD-G3M-17-5345014854_1_492_SE_2015-12-11T01h57m10.fits"],
            "level1Bfilenames" : ["jwtest1009001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1010001_01101_00001_NRS1_uncal.fits"],
            "notes" : ""
        }
    }
}

# IRS2 test data
# IRS2 data has a raw physical size of 2048x3200
# In some places IRS2 considered subarray data only because it is not 2048x2048.
CV3_testdata_dict["IRS2"] = {
    "NID" : {
        #"30055": {
        #    "grism" : "G140H",
        #    "filter" : "OPAQUE",
        #    "mode" : "MOS",
        #    "CAA_lamp" : "LINE1",
        #    "CV3filename" : ["NRSSMOS-MOD-G1H-02-5344031756_1_491_SE_2015-12-10T03h25m56.fits",
        #                     "NRSSMOS-MOD-G1H-02-5344031756_1_492_SE_2015-12-10T03h25m56.fits"],
        #    "level1Bfilenames" : ["", ""],
        #    "MSA_config" : ""
        #},
        "35373" : {
            "grism" : "G235H",
            "filter" : "OPAQUE",
            "mode" : "MOS",
            "CAA_lamp" : "LINE2",
            "CV3filename" : ["NRSV00300010001P0000000002109_1_491_SE_2016-01-02T19h18m49.fits",
                             "NRSV00300010001P0000000002109_1_492_SE_2016-01-02T19h18m49.fits"],
            "level1Bfilenames" : ["jwtest1011001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1012001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "V0030001000101"
        }
    }
}

# DARKS
# the IRS2 data was taken with readout pattern NRSIRS2RAPID
CV3_testdata_dict["DARK"] = {
    "NID" : {
        "30487" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345184144_1_491_SE_2015-12-11T19h03m38.fits",
                             "NRSDET-DARK-IRS2-5345184144_1_492_SE_2015-12-11T19h03m36.fits"],
            "level1Bfilenames" : ["jwtest1001001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1002001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30489" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345184144_2_491_SE_2015-12-11T19h20m47.fits",
                             "NRSDET-DARK-IRS2-5345184144_2_492_SE_2015-12-11T19h20m46.fits"],
            "level1Bfilenames" : ["jwtest1003001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1004001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30491" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345184144_3_491_SE_2015-12-11T19h38m37.fits",
                             "NRSDET-DARK-IRS2-5345184144_3_492_SE_2015-12-11T19h38m37.fits"],
            "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1006001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30494" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345184144_4_491_SE_2015-12-11T19h55m00.fits",
                             "NRSDET-DARK-IRS2-5345184144_4_492_SE_2015-12-11T19h54m59.fits"],
            "level1Bfilenames" : ["jwtest1007001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1008001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30496" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345184144_5_491_SE_2015-12-11T20h12m29.fits",
                             "NRSDET-DARK-IRS2-5345184144_5_492_SE_2015-12-11T20h12m29.fits"],
            "level1Bfilenames" : ["jwtest1009001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1010001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30500" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345184144_6_491_SE_2015-12-11T20h32m39.fits",
                             "NRSDET-DARK-IRS2-5345184144_6_492_SE_2015-12-11T20h32m40.fits"],
            "level1Bfilenames" : ["jwtest1011001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1012001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30506" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345214649_1_491_SE_2015-12-11T22h09m20.fits",
                             "NRSDET-DARK-IRS2-5345214649_1_492_SE_2015-12-11T22h09m19.fits"],
            "level1Bfilenames" : ["jwtest1013001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1014001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30507" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345214649_2_491_SE_2015-12-11T22h26m19.fits",
                             "NRSDET-DARK-IRS2-5345214649_2_492_SE_2015-12-11T22h26m20.fits"],
            "level1Bfilenames" : ["jwtest1015001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1016001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30508" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" :["NRSDET-DARK-IRS2-5345214649_3_491_SE_2015-12-11T22h43m19.fits",
                            "NRSDET-DARK-IRS2-5345214649_3_492_SE_2015-12-11T22h43m20.fits"],
            "level1Bfilenames" : ["jwtest1017001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1018001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30513" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345214649_4_491_SE_2015-12-11T23h00m49.fits",
                             "NRSDET-DARK-IRS2-5345214649_4_492_SE_2015-12-11T23h00m50.fits"],
            "level1Bfilenames" : ["jwtest1019001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1020001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30514" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5345214649_5_491_SE_2015-12-11T23h19m32.fits",
                             "NRSDET-DARK-IRS2-5345214649_5_492_SE_2015-12-11T23h19m31.fits"],
            "level1Bfilenames" : ["jwtest1021001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1022001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30577" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_10_491_SE_2015-12-12T10h55m20.fits",
                             "NRSDET-DARK-IRS2-5346073733_10_492_SE_2015-12-12T10h55m20.fits"],
            "level1Bfilenames" : ["jwtest1024001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1023001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30555" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_1_491_SE_2015-12-12T07h59m04.fits",
                             "NRSDET-DARK-IRS2-5346073733_1_492_SE_2015-12-12T07h59m04.fits"],
            "level1Bfilenames" : ["jwtest1025001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1026001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30560" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_2_491_SE_2015-12-12T08h19m25.fits",
                             "NRSDET-DARK-IRS2-5346073733_2_492_SE_2015-12-12T08h19m24.fits"],
            "level1Bfilenames" : ["jwtest1027001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1028001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30562" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_3_491_SE_2015-12-12T08h36m35.fits",
                             "NRSDET-DARK-IRS2-5346073733_3_492_SE_2015-12-12T08h36m34.fits"],
            "level1Bfilenames" : ["jwtest1029001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1030001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30564" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_4_491_SE_2015-12-12T08h54m35.fits",
                             "NRSDET-DARK-IRS2-5346073733_4_492_SE_2015-12-12T08h54m35.fit"],
            "level1Bfilenames" : ["jwtest1031001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1032001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30565" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_5_491_SE_2015-12-12T09h11m25.fits",
                             "NRSDET-DARK-IRS2-5346073733_5_492_SE_2015-12-12T09h11m24.fits"],
            "level1Bfilenames" : ["jwtest1033001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1034001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30567" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_6_491_SE_2015-12-12T09h30m45.fits",
                             "NRSDET-DARK-IRS2-5346073733_6_492_SE_2015-12-12T09h30m44.fits"],
            "level1Bfilenames" : ["jwtest1035001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1036001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30568" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_7_491_SE_2015-12-12T09h47m05.fits",
                             "NRSDET-DARK-IRS2-5346073733_7_492_SE_2015-12-12T09h47m05.fits"],
            "level1Bfilenames" : ["jwtest1037001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1038001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30573" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_8_491_SE_2015-12-12T10h04m02.fits",
                             "NRSDET-DARK-IRS2-5346073733_8_492_SE_2015-12-12T10h04m03.fits"],
            "level1Bfilenames" : ["jwtest1039001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1040001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30576" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346073733_9_491_SE_2015-12-12T10h21m33.fits",
                             "NRSDET-DARK-IRS2-5346073733_9_492_SE_2015-12-12T10h21m33.fits"],
            "level1Bfilenames" : ["jwtest1041001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1042001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30624" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346161822_1_491_SE_2015-12-12T16h48m34.fits",
                             "NRSDET-DARK-IRS2-5346161822_1_492_SE_2015-12-12T16h48m33.fits"],
            "level1Bfilenames" : ["jwtest1043001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1044001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30625" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346161822_2_491_SE_2015-12-12T17h05m54.fits",
                             "NRSDET-DARK-IRS2-5346161822_2_492_SE_2015-12-12T17h05m54.fits"],
            "level1Bfilenames" : ["jwtest1045001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1046001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30626" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346161822_3_491_SE_2015-12-12T17h21m54.fits",
                             "NRSDET-DARK-IRS2-5346161822_3_492_SE_2015-12-12T17h21m54.fits"],
            "level1Bfilenames" : ["jwtest1047001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1048001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30628" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346161822_4_491_SE_2015-12-12T17h39m43.fits",
                             "NRSDET-DARK-IRS2-5346161822_4_492_SE_2015-12-12T17h39m44.fits"],
            "level1Bfilenames" : ["jwtest1049001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1050001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30631" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSIRS2RAPID",
            "mode" : "IRS2",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-IRS2-5346161822_5_491_SE_2015-12-12T18h58m53.fits",
                             "NRSDET-DARK-IRS2-5346161822_5_492_SE_2015-12-12T18h58m53.fits"],
            "level1Bfilenames" : ["jwtest1051001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1052001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },

        # Change in readout mode
        "30419" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_10_491_SE_2015-12-11T13h10m09.fits",
                             "NRSDET-DARK-TRAD-5345123434_10_492_SE_2015-12-11T13h10m09.fits"],
            "level1Bfilenames" : ["jwtest1001001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1002001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30421" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_11_491_SE_2015-12-11T13h15m19.fits",
                             "NRSDET-DARK-TRAD-5345123434_11_492_SE_2015-12-11T13h15m19.fits"],
            "level1Bfilenames" : ["jwtest1003001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1004001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30422" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_12_491_SE_2015-12-11T13h17m39.fits",
                             "NRSDET-DARK-TRAD-5345123434_12_492_SE_2015-12-11T13h17m39.fits"],
            "level1Bfilenames" : ["jwtest1005001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1006001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30423" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_13_491_SE_2015-12-11T13h19m49.fits",
                             "NRSDET-DARK-TRAD-5345123434_13_492_SE_2015-12-11T13h19m49.fits"],
            "level1Bfilenames" : ["jwtest1007001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1008001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30424" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_14_491_SE_2015-12-11T13h23m29.fits",
                             "NRSDET-DARK-TRAD-5345123434_14_492_SE_2015-12-11T13h23m29.fits"],
            "level1Bfilenames" : ["jwtest1009001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1010001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30425" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_15_491_SE_2015-12-11T13h26m19.fits",
                             "NRSDET-DARK-TRAD-5345123434_15_492_SE_2015-12-11T13h26m19.fits"],
            "level1Bfilenames" : ["jwtest1011001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1012001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30426" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_16_491_SE_2015-12-11T13h29m19.fits",
                             "NRSDET-DARK-TRAD-5345123434_16_492_SE_2015-12-11T13h29m19.fits"],
            "level1Bfilenames" : ["jwtest1013001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1014001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30427" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_17_491_SE_2015-12-11T13h31m39.fits",
                             "NRSDET-DARK-TRAD-5345123434_17_492_SE_2015-12-11T13h31m39.fits"],
            "level1Bfilenames" : ["jwtest1015001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1016001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30429" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_18_491_SE_2015-12-11T13h34m09.fits",
                             "NRSDET-DARK-TRAD-5345123434_18_492_SE_2015-12-11T13h34m09.fits"],
            "level1Bfilenames" : ["jwtest1017001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1018001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30430" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_19_491_SE_2015-12-11T13h36m29.fits",
                             "NRSDET-DARK-TRAD-5345123434_19_492_SE_2015-12-11T13h36m29.fits"],
            "level1Bfilenames" : ["jwtest1019001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1020001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30409" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_1_491_SE_2015-12-11T12h49m00.fits",
                             "NRSDET-DARK-TRAD-5345123434_1_492_SE_2015-12-11T12h49m00.fits"],
            "level1Bfilenames" : ["jwtest1021001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1022001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30431" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_20_491_SE_2015-12-11T13h39m02.fits",
                             "NRSDET-DARK-TRAD-5345123434_20_492_SE_2015-12-11T13h39m02.fits"],
            "level1Bfilenames" : ["jwtest1023001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1024001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30432" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_21_491_SE_2015-12-11T13h41m52.fits",
                             "NRSDET-DARK-TRAD-5345123434_21_492_SE_2015-12-11T13h41m51.fits"],
            "level1Bfilenames" : ["jwtest1025001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1026001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30433" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_22_491_SE_2015-12-11T13h43m51.fits",
                             "NRSDET-DARK-TRAD-5345123434_22_492_SE_2015-12-11T13h43m51.fits"],
            "level1Bfilenames" : ["jwtest1027001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1028001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30434" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_23_491_SE_2015-12-11T13h46m11.fits",
                             "NRSDET-DARK-TRAD-5345123434_23_492_SE_2015-12-11T13h46m11.fits"],
            "level1Bfilenames" : ["jwtest1029001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1030001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30435" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_24_491_SE_2015-12-11T13h48m42.fits",
                             "NRSDET-DARK-TRAD-5345123434_24_492_SE_2015-12-11T13h48m42.fits"],
            "level1Bfilenames" : ["jwtest1031001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1032001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30441" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_25_491_SE_2015-12-11T15h07m32.fits",
                             "NRSDET-DARK-TRAD-5345123434_25_492_SE_2015-12-11T15h07m32.fits"],
            "level1Bfilenames" : ["jwtest1033001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1034001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30410" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_2_491_SE_2015-12-11T12h51m29.fits",
                             "NRSDET-DARK-TRAD-5345123434_2_492_SE_2015-12-11T12h51m29.fits"],
            "level1Bfilenames" : ["jwtest1035001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1036001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30411" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_3_491_SE_2015-12-11T12h53m39.fits",
                             "NRSDET-DARK-TRAD-5345123434_3_492_SE_2015-12-11T12h53m39.fits"],
            "level1Bfilenames" : ["jwtest1037001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1038001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30413" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_4_491_SE_2015-12-11T12h55m59.fits",
                             "NRSDET-DARK-TRAD-5345123434_4_492_SE_2015-12-11T12h55m59.fits"],
            "level1Bfilenames" : ["jwtest1039001_01101_00001_NRS1_uncal.fits", "jwtest1040001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30414" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_5_491_SE_2015-12-11T12h58m29.fits",
                             "NRSDET-DARK-TRAD-5345123434_5_492_SE_2015-12-11T12h58m29.fits"],
            "level1Bfilenames" : ["jwtest1041001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1042001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30415" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_6_491_SE_2015-12-11T13h00m39.fits",
                             "NRSDET-DARK-TRAD-5345123434_6_492_SE_2015-12-11T13h00m39.fits"],
            "level1Bfilenames" : ["jwtest1043001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1044001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30416" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_7_491_SE_2015-12-11T13h02m59.fits",
                             "NRSDET-DARK-TRAD-5345123434_7_492_SE_2015-12-11T13h02m59.fit"],
            "level1Bfilenames" : ["jwtest1045001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1046001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30417" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_8_491_SE_2015-12-11T13h05m19.fits",
                             "NRSDET-DARK-TRAD-5345123434_8_492_SE_2015-12-11T13h05m19.fits"],
            "level1Bfilenames" : ["jwtest1047001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1048001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        },
        "30418" : {
            "grism" : "G140H",
            "filter" : "OPAQUE",
            "readpatt": "NRSRAPID",
            "mode" : "MOS",
            "CAA_lamp" : "NO_LAMP",
            "CV3filename" : ["NRSDET-DARK-TRAD-5345123434_9_491_SE_2015-12-11T13h07m39.fits",
                             "NRSDET-DARK-TRAD-5345123434_9_492_SE_2015-12-11T13h07m39.fit"],
            "level1Bfilenames" : ["jwtest1049001_01101_00001_NRS1_uncal.fits",
                                  "jwtest1050001_01101_00001_NRS1_uncal.fits"],
            "MSA_config" : "ADRCLOSED"
        }
    }
}


if __name__ == '__main__':
    # example of how to extract info from the dictionary
    detector = "491"
    d = CV3_testdata_dict
    for nid, nid_dict_key in d["MOS"]["NID"].items():
        if nid_dict_key["grism"] == "G140H":
            if nid_dict_key["filter"] == "OPAQUE":
                print "NID =", nid
                print "CV3filename = ", nid_dict_key["CV3filename"][0] # the 0 is because of detector 491
