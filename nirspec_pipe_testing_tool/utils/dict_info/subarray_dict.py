
# HEADER
__author__ = "M. A. Pena-Guerrero and E. Puga"
__version__ = "2.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed
# May 2019 - Version 2.0: corrected version


# Subarray dictionary  ->  values taken from ESA-JWST-SCI-NRS-TN-2018-002 by S. Birkmann, N. Luetzgendorf

subarray_dict = {}
subarray_dict["FULL-FRAME"] = {"substrt1" : 1,
                               "substrt2" : {# grat  NRS1, NRS2
                                            "G140H": (1,   1),
                                            "G140M": (1,   1),
                                            "G235H": (1,   1),
                                            "G235M": (1,   1),
                                            "G395H": (1,   1),
                                            "G395M": (1,   1),
                                            "PRISM": (1,   1),
                                            "MIRROR":(1,   1)},
                               "subsize1" : 2048,
                               "subsize2" : 2048}

subarray_dict["ALLSLITS"] = {"substrt1" : 1,
                             "substrt2" : {# grat  NRS1, NRS2
                                           "G140H":(897, 895),
                                           "G140M":(897, 895),
                                           "G235H":(897, 895),
                                           "G235M":(897, 895),
                                           "G395H":(897, 895),
                                           "G395M":(897, 895),
                                           "PRISM":(897, 895),
                                           "MIRROR":(897,895)},
                             "subsize1" : 2048,
                             "subsize2" : 256}

subarray_dict["S200A1"]   = {"substrt1" : 1,
                             "substrt2" :{# grat    NRS1, NRS2
                                           "G140H":(1054, 930),
                                           "G140M":(1054, 930),
                                           "G235H":(1054, 930),
                                           "G235M":(1054, 930),
                                           "G395H":(1054, 930),
                                           "G395M":(1041, 945),
                                           "PRISM":(1054, 930),
                                           "MIRROR":(1054,930)},
                             "subsize1" : 2048,
                             "subsize2" : 64}

subarray_dict["S200A2"]   = {"substrt1" : 1,
                             "substrt2" : {# grat  NRS1,  NRS2
                                           "G140H":(1018, 966),
                                           "G140M":(1018, 966),
                                           "G235H":(1018, 966),
                                           "G235M":(1018, 966),
                                           "G395H":(1018, 966),
                                           "G395M":(1005, 981),
                                           "PRISM":(1018, 966),
                                           "MIRROR":(1018,966)},
                             "subsize1" : 2048,
                             "subsize2" : 64}

subarray_dict["S200B1"]   = {"substrt1" : 1,
                             "substrt2": {# grat   NRS1, NRS2
                                          "G140H": (918, 1068),
                                          "G140M": (918, 1068),
                                          "G235H": (918, 1068),
                                          "G235M": (918, 1068),
                                          "G395H": (918, 1068),
                                          "G395M": (905, 1083 ),
                                          "PRISM": (918, 1068),
                                          "MIRROR":(918, 1068)},
                             "subsize1" : 2048,
                             "subsize2" : 64}

subarray_dict["S400A1"]   = {"substrt1" : 1,
                             "substrt2" : {# grat  NRS1, NRS2
                                           "G140H":(978, 1006),
                                           "G140M":(978, 1006),
                                           "G235H":(978, 1006),
                                           "G235M":(978, 1006),
                                           "G395H":(978, 1006),
                                           "G395M":(965, 1021),
                                           "PRISM":(978, 1006),
                                           "MIRROR":(978,1006)},
                             "subsize1" : 2048,
                             "subsize2" : 64}

subarray_dict["SUB1024A"] = {"substrt1" : 1,
                             "substrt2" : {# grat  NRS1, NRS2
                                           "G140H":(963, 1055),
                                           "G140M":(963, 1055),
                                           "G235H":(967, 1051),
                                           "G235M":(963, 1055),
                                           "G395H":(963, 1055),
                                           "G395M":(946, 1072),
                                           "PRISM":(957, 1059),
                                           "MIRROR":(975, 975)},
                             "subsize1" : 1024,
                             "subsize2" : 32}

subarray_dict["SUB1024B"] = {"substrt1" : 1025,
                             "substrt2" : {# grat  NRS1, NRS2
                                           "G140H":(963, 1055),
                                           "G140M":(963, 1055),
                                           "G235H":(967, 1051),
                                           "G235M":(963, 1055),
                                           "G395H":(963, 1055),
                                           "G395M":(946, 1072),
                                           "PRISM":(957, 1059),
                                           "MIRROR":(975, 975)},
                             "subsize1" : 1024,
                             "subsize2" : 32}

subarray_dict["SUB2048"]  = {"substrt1" : 1,
                             "substrt2" : {# grat  NRS1, NRS2
                                           "G140H":(963, 1055),
                                           "G140M":(963, 1055),
                                           "G235H":(967, 1051),
                                           "G235M":(963, 1055),
                                           "G395H":(963, 1055),
                                           "G395M":(946, 1072),
                                           "PRISM":(957, 1059),
                                           "MIRROR":(975, 975)},
                             "subsize1" : 2048,
                             "subsize2" : 32}

subarray_dict["SUB32"]    = {"substrt1" : 1399,
                             "substrt2" : {# grat  NRS1,  NRS2
                                           "G140H":(None, None),
                                           "G140M":(None, None),
                                           "G235H":(None, None),
                                           "G235M":(None, None),
                                           "G395H":(None, None),
                                           "G395M":(None, None),
                                           "PRISM":(None, None),
                                           "MIRROR":(975, 975)},
                             "subsize1" : 32,
                             "subsize2" : 32}

subarray_dict["SUB512"]   = {"substrt1" : 1025,
                             "substrt2" : {# grat  NRS1, NRS2
                                           "G140H":(963, 1055),
                                           "G140M":(963, 1055),
                                           "G235H":(967, 1051),
                                           "G235M":(963, 1055),
                                           "G395H":(963, 1055),
                                           "G395M":(946, 1072),
                                           "PRISM":(957, 1059),
                                           "MIRROR":(975, 975)},
                             "subsize1" : 512,
                             "subsize2" : 32}

subarray_dict["SUB512S"]  = {"substrt1" : 1025,
                             "substrt2" : {# grat  NRS1, NRS2
                                           "G140H":(971, 1063),
                                           "G140M":(971, 1063),
                                           "G235H":(975, 1059),
                                           "G235M":(971, 1063),
                                           "G395H":(971, 1063),
                                           "G395M":(954, 1080),
                                           "PRISM":(965, 1067),
                                           "MIRROR":(983, 983)},
                             "subsize1" : 512,
                             "subsize2" : 16}
