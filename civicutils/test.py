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


from read_and_write import readInSnvs,get_dict_support,write_to_json,write_to_yaml,write_match
from query import query_civic
from filtering import filter_civic
from match import match_in_civic,filter_matches,annotate_ct,process_drug_support,filter_ct,add_ct,add_match

# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v1.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v2.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v3.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v4.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v9.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v10.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v11.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_snv.v12.txt"
# inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/complete_snv.txt"
inFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/complete_snv.v1.txt"
(rawData,snvData,extraHeader) = readInSnvs(inFile)
print("rawData:")
print(json.dumps(rawData,indent=1))
print("snvData:")
print(json.dumps(snvData,indent=1))
print("extraHeader:")
print(extraHeader)

# print("\nQuerying CIVIC...")
# allGenes = list(snvData.keys())
# varMap = query_civic(allGenes, identifier_type="entrez_symbol")
# origMap = query_civic(allGenes, identifier_type="entrez_symbol")
# print("varMap:")
# print(json.dumps(varMap,indent=1))
# print("origMap:")
# print(json.dumps(origMap,indent=1))

# print("\nFiltering CIVIC...")
# varMap = filter_civic(varMap, evidence_status_in = ['ACCEPTED'], var_origin_not_in = ['GERMLINE'], output_empty=False)
# print("varMap:")
# print(json.dumps(varMap,indent=1))

# print("\nAnnotating CT...")
# blackList = []
# whiteList = ["leukemia"]
# altList = []
# annotMap = annotate_ct(varMap,blackList,whiteList,altList)
# annotMap = annotate_ct(origMap,blackList,whiteList,altList)
# print("annotMap:")
# print(json.dumps(annotMap,indent=1))

data_type = "SNV"
tierSelect = "highest"
# tierSelect = "all"
# tierSelect = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]
print("\nMatching CIVIC...")
(matchMap,matchedIds,varMap) = match_in_civic(snvData, data_type, identifier_type="entrez_symbol", select_tier=tierSelect, varMap=None)
# (matchMap,matchedIds,varMap) = match_in_civic(snvData, data_type, identifier_type="entrez_symbol", select_tier=tierSelect, varMap=varMap)
# (matchMap,matchedIds,varMap) = match_in_civic(snvData, data_type, identifier_type="entrez_symbol", select_tier=tierSelect, varMap=annotMap)
# (matchMap,matchedIds,varMap) = match_in_civic(snvData, data_type, identifier_type="entrez_symbol", select_tier=tierSelect, varMap=origMap)
print("matchMap:")
print(json.dumps(matchMap,indent=1))
print("matchedIds:")
print(len(matchedIds))
print(matchedIds)
print("varMap:")
print(json.dumps(varMap,indent=1))

# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output.json"
# write_to_json(matchedIds,outFile,indent=1)
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output.yaml"
# write_to_yaml(data_type,outFile)

# print("\nFiltering match...")
# tierSelect="highest"
# (matchMap,matchedIds) = filter_matches(matchMap,tierSelect)
# (annotMatch,matchedIds) = filter_matches(annotMatch,tierSelect)
# print("matchMap:")
# print(json.dumps(matchMap,indent=1))
# print("matchedIds:")
# print(len(matchedIds))
# print(matchedIds)

# print("\nFiltering CIVIC...")
# fMatch = filter_civic(varMap, var_id_in=["3109"], output_empty=False)
# fMatch = filter_civic(varMap, var_id_in=["3109"], output_empty=True)
# print("fMatch:")
# print(json.dumps(fMatch,indent=1))

# (test,testIds,testMap) = match_in_civic(snvData, data_type, identifier_type="entrez_symbol", select_tier=tierSelect, varMap=fMatch)
# print("test:")
# print(json.dumps(test,indent=1))


# print("\nAnnotating CT...")
# blackList = []
# blackList = ["melanoma"]
# whiteList = ["breast","ovarian"]
# altList = ["cancer"]
# annotMap = annotate_ct(varMap,blackList,whiteList,altList)
# print("annotMap:")
# print(json.dumps(annotMap,indent=1))
# annotOrig = annotate_ct(origMap,blackList,whiteList,altList)
# print("annotOrig:")
# print(json.dumps(annotOrig,indent=1))


# print("\nFiltering CT...")
# ctSelect = "highest"
# annotMap = filter_ct(annotMap,ctSelect)
# print("annotMap:")
# print(json.dumps(annotMap,indent=1))


# print("\nProcessing drug support...")
# Get dict of drug support
supportDict = get_dict_support()
# print("supportDict:")
# print(json.dumps(supportDict,indent=1))

# annotMatch = process_drug_support(matchMap,annotMap,supportDict)
# print("annotMatch:")
# print(json.dumps(annotMatch,indent=1))

# annotMatch2 = process_drug_support(matchMap,annotOrig,supportDict)
# print("annotMatch2:")
# print(json.dumps(annotMatch2,indent=1))


# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v1.txt"
# write_match(matchMap,varMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False)
# write_match(matchMap,varMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False)
# write_match(matchMap,varMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=True, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(matchMap,varMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=False, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(matchMap,varMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=False, writeCt=False, writeSupport=False, writeComplete=True)

# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v2.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v3.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v4.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v5.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v6.txt"
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=True, hasCt=True, writeCt=True, writeSupport=True, writeComplete=False)
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=True, writeCt=True, writeSupport=True, writeComplete=False)
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=True, writeCt=True, writeSupport=False, writeComplete=False)
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=True, writeCt=True, writeSupport=False, writeComplete=False)
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=True, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(matchMap,varMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=False, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=True, writeCt=False, writeSupport=False, writeComplete=True)
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile,hasSupport=False, hasCt=True, writeCt=True, writeSupport=False, writeComplete=True)

# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v7.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v8.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v9.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v10.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v11.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v12.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v13.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.v14.txt"
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=False, hasCt=False, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=False, hasCt=True, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(matchMap,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(annotMatch,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(annotMatch,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=True, writeSupport=False, writeComplete=False)

# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.r1.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.r2.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.r3.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.r4.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.r5.txt"
# write_match(annotMatch,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=True, writeSupport=True, writeComplete=False)
# write_match(matchMap,fMatch,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=False, hasCt=False, writeCt=False, writeSupport=False, writeComplete=False)

# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.r3b.txt"
# write_match(annotMatch,annotOrig,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False)

# write_match(annotMatch,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False)
# write_match(annotMatch,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=False, writeSupport=False, writeComplete=False)
# write_match(annotMatch,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=False, writeSupport=False, writeComplete=True)
# write_match(annotMatch,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=True)
# write_match(annotMatch,annotMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=True, hasCt=True, writeCt=True, writeSupport=True, writeComplete=True)

# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.nv1.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.nv2.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.nv3.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.nv4.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.nv5.txt"
# outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.nv6.txt"
outFile = "/cluster/work/nexus/lourdes/civic_query/tests/new/test_output_snv.nv7.txt"
write_match(matchMap,varMap,rawData,extraHeader,dataType="SNV",outfile=outFile, hasSupport=False, hasCt=False, writeCt=False, writeSupport=False, writeComplete=False)


sys.exit(1)



# gene = "ERBB4"
# variant = "310"
# variant = "1673"
# evidence_type = "PREDICTIVE"
# newMap = {}
# newMap[gene] = {}
# newMap[gene][variant] = {}
# newMap[gene][variant]['evidence_items'] = {}
# newMap[gene][variant]['evidence_items'][evidence_type] = {}
# print("newMap:")
# print(json.dumps(newMap,indent=1))

# diseases = ["melanoma"]
# diseases = ["dummy","melanoma"]
# diseases = ["MELANOMA"]
# diseases = ["HEAD AND NECK SQUAMOUS CELL CARCINOMA"]
# ct = "ct"
# newMap = add_ct(diseases,ct,gene,variant,evidence_type,newMap,varMap,isAnnot=False)
# newMap = add_ct(diseases,ct,gene,variant,evidence_type,newMap,varMap,isAnnot=True)
# newMap = add_ct(diseases,ct,gene,variant,evidence_type,newMap,annotMap,isAnnot=False)
# newMap = add_ct(diseases,ct,gene,variant,evidence_type,newMap,annotMap,isAnnot=True)
# newMap = add_ct(diseases,ct,gene,variant,evidence_type,newMap,origMap,isAnnot=False)
# newMap = add_ct(diseases,ct,gene,variant,evidence_type,newMap,origMap,isAnnot=True)
# print("newMap:")
# print(json.dumps(newMap,indent=1))

# matchMap = process_drug_support(matchMap,varMap,supportDict)
# matchMap = process_drug_support(matchMap,annotMap,supportDict)
# annotMatch = process_drug_support(matchMap,annotMap,supportDict)
# annotMatch = process_drug_support(annotMatch,annotMap,supportDict)
# print("annotMatch:")
# print(json.dumps(annotMatch,indent=1))
# print("matchMap:")
# print(json.dumps(matchMap,indent=1))
# print("matchedIds:")
# print(len(matchedIds))
# print(matchedIds)

# gene = "ERBB4"
# variant = "c.83-126671T>A|.|sequence_feature|1/27|1"
# match = matchMap[gene][variant]
# match = {}
# match[gene] = {}
# match[gene][variant] = matchMap[gene][variant]
# match[gene][variant] = {}
# print("match:")
# print(json.dumps(match,indent=1))

# annot_match = process_drug_support(match,varMap,supportDict)
# annot_match = process_drug_support(match,annotMap,supportDict)
# print("annot_match:")
# print(json.dumps(annot_match,indent=1))

# print("match:")
# print(json.dumps(annot_match[gene][variant],indent=1))

# print("matchMap:")
# print(json.dumps(matchMap,indent=1))

# print("\nAdding match...")
# (matchMap2,matchedIds) = add_match(annotMatch,gene,variant,match)
# gene = "DUMMY"
# variant = "310"
# (matchMap2,matchedIds) = add_match(matchMap,gene,variant,match)
# (matchMap2,matchedIds) = add_match(matchMap,gene,variant,annot_match[gene][variant])
# print("matchMap2:")
# print(json.dumps(matchMap2,indent=1))
# print("matchedIds:")
# print(len(matchedIds))

# print("\nFiltering match...")
# tierSelect=["tier_1","tier_1b","tier_4"]
# tierSelect=["tier_4"]
# tierSelect="highest"
# tierSelect="all"
# (matchMap,matchedIds) = filter_matches(matchMap,tierSelect)
# print("matchMap:")
# print(json.dumps(matchMap,indent=1))
# print("matchedIds:")
# print(len(matchedIds))
# print(matchedIds)

# testMap = process_drug_support(matchMap,varMap,supportDict)
# testMap = process_drug_support(matchMap,annotMap,supportDict)
# print("testMap:")
# print(json.dumps(matchMap,indent=1))
# print("matchedIds:")
# print(len(matchedIds))
# print(matchedIds)

# print("\nAnnotating CT...")
# blackList = []
# whiteList = ["leukemia"]
# altList = []
# annotMap = annotate_ct(varMap,blackList,whiteList,altList)
# print("annotMap:")
# print(json.dumps(annotMap,indent=1))

# blackList = ["leukemia"]
# whiteList = ["melanoma"]
# altList = []
# annotMap = annotate_ct(varMap,blackList,whiteList,altList)
# annotMap = annotate_ct(annotMap,blackList,whiteList,altList)
# print("annotMap:")
# print(json.dumps(annotMap,indent=1))

# print("\nFiltering CT...")
# ctSelect = ["nct"]
# annotMap = filter_ct(annotMap,ctSelect)
# print("annotMap:")
# print(json.dumps(annotMap,indent=1))
# testMap = filter_ct(varMap,ctSelect)
# print("testMap:")
# print(json.dumps(testMap,indent=1))

# print("\nFiltering CIVIC...")
# testMap = filter_civic(varMap, evidence_status_in = ['ACCEPTED'], var_origin_not_in = ['GERMLINE'], output_empty=False)
# print("testMap:")
# print(json.dumps(testMap,indent=1))
# testMap = filter_civic(annotMap, evidence_status_in = ['ACCEPTED'], var_origin_not_in = ['GERMLINE'], output_empty=False)
# print("testMap:")
# print(json.dumps(testMap,indent=1))

