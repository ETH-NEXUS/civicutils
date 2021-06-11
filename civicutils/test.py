#!/usr/bin/env python

import sys
import os
import re
import json

# from read_and_write import get_dict_aminoacids
# dict_codes = get_dict_aminoacids()
# print(json.dumps(dict_codes,indent=1))

# from read_and_write import readInSnvs
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v1.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v2.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v3.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v4.txt"
# (rawData,snvData) = readInSnvs(inFile)
# print("rawData:")
# print(json.dumps(rawData,indent=1))
# print("snvData:")
# print(json.dumps(snvData,indent=1))

# from read_and_write import readInCnvs
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_cnv.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_cnv.v1.txt"
# (rawData,cnvData) = readInCnvs(inFile)
# print("rawData:")
# print(json.dumps(rawData,indent=1))
# print("snvData:")
# print(json.dumps(cnvData,indent=1))

from read_and_write import readInExpr
inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_expr.txt"
(rawData,exprData) = readInExpr(inFile)
print("rawData:")
print(json.dumps(rawData,indent=1))
print("snvData:")
print(json.dumps(exprData,indent=1))
