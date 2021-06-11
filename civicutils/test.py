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

# from read_and_write import readInExpr
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_expr.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_expr.v1.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_expr.v2.txt"
# (rawData,exprData) = readInExpr(inFile)
# print("rawData:")
# print(json.dumps(rawData,indent=1))
# print("snvData:")
# print(json.dumps(exprData,indent=1))


# from read_and_write import readInSnvs
# from query import query_civic

# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v1.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v2.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v3.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v4.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v5.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v6.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v7.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v8.txt"
# (rawData,snvData) = readInSnvs(inFile)
# print("rawData:")
# print(json.dumps(rawData,indent=1))
# print("snvData:")
# print(json.dumps(snvData,indent=1))

# allGenes = list(snvData.keys())
# allGenes = "a"
# allGenes = []
# print("allGenes:")
# print(len(allGenes))

# print("Querying CIVIC...")
# varMap = query_civic(allGenes, identifier_type="entrez_symbol")
# varMap = query_civic(allGenes, identifier_type="ENTREZ_SYMBOL")
# varMap = query_civic(allGenes, identifier_type="dummy")
# varMap = query_civic(allGenes, identifier_type="entrez_id")
# varMap = query_civic(allGenes, identifier_type="civic_id")
# print("varMap:")
# print(json.dumps(varMap,indent=1))


# from utils import check_is_list,check_is_dict,check_is_str,check_tier_selection

# thisList = []
# thisList = [2,3,"q"]
# thisDict = {}
# thisStr = ""

# check_is_list(thisList,"thisList")
# check_is_list(thisDict,"thisDict")
# check_is_list(thisStr,"thisStr")

# check_is_dict(thisList,"thisList")
# check_is_dict(thisDict,"thisDict")
# check_is_dict(thisStr,"thisStr")

# check_is_str(thisList,"thisList")
# check_is_str(thisDict,"thisDict")
# check_is_str(thisStr,"thisStr")

# sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]

# select_tier = []
# select_tier = None
# select_tier = {"2":None}
# select_tier = "x"
# select_tier = "all"
# select_tier = "highest"
# select_tier = ["tier_3","highest"]
# select_tier = ["tier_3","tier_3","tier_2","tier_1b","tier_4"]
# select_tier = ["tier_3","tier_3","tier_2","tier_1b","tier_4","tier_1"]

# print("select_tier (before):")
# print(select_tier)
# select_tier = check_tier_selection(select_tier,sorted_tiers)
# print("select_tier (after):")
# print(select_tier)


from read_and_write import readInSnvs
from query import query_civic
from match import match_in_civic
from filtering import filter_civic

inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.txt"
(rawData,snvData) = readInSnvs(inFile)
print("rawData:")
print(json.dumps(rawData,indent=1))
print("snvData:")
print(json.dumps(snvData,indent=1))

print("\nQuerying CIVIC...")
allGenes = list(snvData.keys())
varMap = query_civic(allGenes, identifier_type="entrez_symbol")
print("varMap:")
print(json.dumps(varMap,indent=1))

print("\nFiltering CIVIC...")
varMap = filter_civic(varMap, source_status_in = ['ACCEPTED'], var_origin_not_in = ['GERMLINE'], output_empty=False)
print("varMap:")
print(json.dumps(varMap,indent=1))

data_type = "SNV"
tierSelect = "all"
print("\nMatching CIVIC...")
# (matchMap,matchedIds,varMap) = match_in_civic(snvData, data_type, identifier_type="entrez_symbol", select_tier=tierSelect, varMap=None)
(matchMap,matchedIds,varMap) = match_in_civic(snvData, data_type, identifier_type="entrez_symbol", select_tier="all", varMap=varMap)
print("matchMap:")
print(json.dumps(matchMap,indent=1))
print("matchedIds:")
print(len(matchedIds))
print(matchedIds)
print("varMap:")
print(json.dumps(varMap,indent=1))


