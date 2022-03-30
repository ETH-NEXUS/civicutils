#!/usr/bin/env python

'''
Process (and optionally) combine CIVIC variant info
Lourdes Rosano, Feb 2022
'''

import sys
import os
import argparse
import re
import copy

## Dictionary that allows mapping of special cases where drugs are referred to in CIVIC using synonym names or terms
## Structure of the dictionary follows the form: {CIVIC SYNONYM (key): DRUG (value)}
drugSyn_dict = {'DOVITINIB DILACTIC ACID (TKI258 DILACTIC ACID)':'DOVITINIB', '5-FLUOROURACIL':'FLUOROURACIL', '5-FU':'FLUOROURACIL', 'ADO-TRASTUZUMAB EMTANSINE':'TRASTUZUMAB EMTANSINE', 'PD0325901':'PD-0325901', 'PD173074':'PD-173074', 'BGJ-398':'INFIGRATINIB', 'BGJ398':'INFIGRATINIB'}

## Dictionary that allows prioritization of CIVIC support categories (in oncoprint) when >1 gene has CIVIC info for a given drug+sample
supportPrior = {'civic_support':1, 'civic_resistance':1, 'civic_conflict':2, 'civic_unknown':3, 'civic_unspecific_vars':4, 'dgidb_only':5}

'''
Functions
'''

# Given a header already split by tabs, return the position of a given column name
# Throw an error if provided column name is not found
def get_column_position(column_name, header_split):
    if column_name in header_split:
        pos = header_split.index(column_name)
    else:
        print("Error! Column '%s' could not be found in header!" %(column_name))
        sys.exit(1)
    return pos


## Do no prioritize evidence levels but account for all supporting evidence items across all levels
def parseInputFile(sampleFile, sampleName, clinInfoDict):
# def parseInputFile(sampleFile,sampleName,clinInfoDict,drugsToVariants):
    print("Sample %s. File: %s" %(sampleName, sampleFile))

    infile = open(sampleFile,'r')
    header = infile.readline()
    header_split = header.strip().split("\t")
    # Retrieve column positions for the required columns
    tier_pos = get_column_position("CIViC_Tier", header_split)
    civic_score_pos = get_column_position("CIViC_Score", header_split)
    drug_supp_pos = get_column_position("CIViC_Drug_Support", header_split)
    drug_pos = get_column_position("CIViC_PREDICTIVE", header_split)
    diag_pos = get_column_position("CIViC_DIAGNOSTIC", header_split)
    prog_pos = get_column_position("CIViC_PROGNOSTIC", header_split)
    pred_pos = get_column_position("CIViC_PREDISPOSING", header_split)

    ## Keep track of number of lines with available CIVIC information,
    ## number of parsed drugs and number of parsed genes per sample
    civicLines = 0
    parsedDrugs = []
    parsedGenes = []

    n_tier_1 = 0
    n_tier_1b = 0
    n_tier_1_agg = 0
    n_tier_2 = 0
    n_tier_3 = 0
    n_tier_4 = 0
    matched_variants_mapping = {}
    matched_variants = 0
    n_drug_info = 0

    for line in infile:
        lineSplit = line.strip().split("\t")
        tier = str(lineSplit[tier_pos].strip())
        ## Tier=4 should be skipped as no information was found on CIVIC for the current gene
        ## Special case Tier='SNV_only' should also be skipped. This only occurs for CNV CIVIC data and corresponds to a special tier=3 case
        ## Also skip line if Tier=1,2 or 3, but still no drug information was available on CIVIC for the current gene (eg. prognostic info only)
        if tier=="4":
            n_tier_4 += 1
            continue

        ## Only process further lines which have CIVIC info available
        civicLines += 1
        if tier=="3":
            n_tier_3 += 1
        if tier=="2":
            n_tier_2 += 1
        if tier=="1b":
            n_tier_1b += 1
        if tier=="1":
            n_tier_1 += 1
        if tier=="1" or tier=="1b":
            n_tier_1_agg += 1

        ## 1. TODO: check number of variants matched per line (and keep track of associated tier)
        # Keep track of number of matched variants, keep counts per tier
# n_matched["tier_1"] = # matches, sum up for each variant encountered which is tagged as tier_1, and so on. Then, to comput mean simply divide the final count of # matches, by the total number of variants tagged as tier_1 (n_tier_1)
# we can compute the mean number of matched variants across the whole variant set (patient) by simply aggregating the match count across all tiers, and then dividing by the total number of variants with civic info

        # NOTE: take into account CIVIC variant names containing ":", e.g. 'NTRK1:LMNA::NTRK1 E11-E10:18.5;'
        civic_infos = lineSplit[civic_score_pos].strip()
        if civic_infos == ".":
            if not (tier=="3" or tier=="4"):
                raise ValueError("Encountered unexpected case of variant with tier!=3 and tier!=4 but no associated variant matches in CIVIC in line %s" %(line.strip()))
            n_variants = 0
        else:
            civic_infos_list = civic_infos.split(';')
            n_variants = len(civic_infos_list)
            if n_variants == 0:
                raise ValueError("Encountered unexpected case of no associated variant matches in CIVIC for line %s" %(line.strip()))

        # tier -> # CIVIC variants matched
        if tier not in matched_variants_mapping.keys():
            matched_variants_mapping[tier] = 0
        matched_variants_mapping[tier] += n_variants
        matched_variants += n_variants


        # TODO: parse column containing civic scores (there, all matched CIVIC variants are listed), and simply count the different variants (double check no duplicates here?). NOTE: for SNVs, all variant matches are ensure to originate for a single variant annotation + gene, while for CNVs, several variant matches arising from different genes are taken into account within the same line


        ## 2. TODO: check across all 4 evidence columns, and retrieve all unique instances of disease_name + ct classification
        # TODO: parse all diseases + ct available for each variant line (and associated tier) -> report # diseases for the whole variant set (patient) as well as counts per tier, and also report # diseases per ct classification for the whole variant set (patient) as well as counts per tier
        # TODO: do not count unique diseases + ct per variant but rather, look at the whole set and retrieve the unique diseases+ct, and look at the tier subsets and retrieve the unique diseases+ct
        for tmp_pos in [drug_pos, diag_pos, prog_pos, pred_pos]:
            evidence_infos = lineSplit[tmp_pos].strip()
            if evidence_infos == ".":
                continue
            evidence_infos_list = evidence_infos.split(";")
            # NOTE: take into account CIVIC variant names containing ":", e.g. 'NTRK1:LMNA::NTRK1 E11-E10:18.5;'
            for evidence_info in evidence_infos_list:
                tmp_evidence_split = evidence_info.strip().split("(")
                evidence_split = tmp_evidence_split[0].strip().split("|")
                if tmp_pos == drug_pos:
                    # drug = evidence_split[
                else:
                    ct_type = evidence_split[len(evidence_split)-1]


        ## 3. TODO: check in the consensus drug support column (double check in the corresponding Predictive column?)
        # TODO: for each variant line, select only drugs associated with the highest available ct classification (if only nct, then all will be taken into account) -> 
        # TODO: count # of unique drugs (already filtered for only highest ct) predicted for the whole variant set (patient), as well as counts per tier
        # TODO: also, per ct classification, report # unique drugs, and also report # unique drugs per tier
        # TODO: taking into account all available drugs + associated ct, for each unique drug, report # of ct classifications associated with it across the file, also, report counts per tier (then do mean across all drugs for the patient, as well as mean across all variants assigned to each tier)

        drug_infos = lineSplit[drug_supp_pos].strip()
        if drug_infos != ".":
            n_drug_info += 1
            # E.g.: 'ENTRECTINIB:NCT:CIVIC_SUPPORT;LAROTRECTINIB:NCT:CIVIC_SUPPORT;..'
            drug_infos_list = drug_infos.split(';')
            if not drug_infos_list:
                raise ValueError("Encountered unexpected case of no associated variant matches in CIVIC for line %s" %(line.strip()))

            ct_mapping = {}

            # Format: 'DRUG_NAME:CT:CIVIC_SUPPORT'
            for drug_info in drug_infos_list:
                tmp_drug_split = drug_info.strip().split(":")
                if len(tmp_drug_split) != 3:
                    # TODO
                drug = tmp_drug_split[0].strip()
                ct_type = tmp_drug_split[1].strip()
                civic_support = tmp_drug_split[2].strip()
                if ct_type not in ct_mapping.keys():
                    ct_mapping[ct_type] = []

        ## 4. TODO: check in the consensus drug support column: report, out of all unique drugs, how many (%) have associated: a. civic_support, civic_resistance, etc. Also, take into account situations where the same drug is associated to different classifications depending on the ct, variant, etc. -> how to handle?




    infile.close()
    print("Parsed %s lines with available CIVIC information associated to %s genes and %s drugs" %(civicLines,len(parsedGenes),len(parsedDrugs)))

    return (clinInfoDict)


'''
Script
'''

parser = argparse.ArgumentParser(description='Combine SNV, CNV, DRS and CIVIC drug predictions.')
parser.add_argument('--inputTable', dest='inputTable', required=True, help='Input file with combined SNV+DRS+CNV gene-drug predictions.')
parser.add_argument('--inputDir_civic_snv', dest='inputDir_civic_snv', required=True, help='Input directory with SNV CIVIC files of format [sample].civic_snv.txt.')
parser.add_argument('--inputDir_civic_cnv', dest='inputDir_civic_cnv', required=True, help='Input directory with CNV CIVIC files of format [sample].civic_cnv.txt.')
parser.add_argument('--fileEnding_snv', dest='fileEnding_snv', required=True, help='To retrieve correct input files, specify the desired file ending ("civic_snv.txt" for SNV data).')
parser.add_argument('--fileEnding_cnv', dest='fileEnding_cnv', required=True, help='To retrieve correct input files, specify the desired file ending ("civic_cnv.txt" for CNV data).')
parser.add_argument('--outFileTag', dest='outFileTag', required=True, help='Name prefix of the output files.')

args = parser.parse_args()


### Process input files for CIVIC (pool SNV+CNV data together)

# Global dictionary (pools SNV+CNV CIVIC info) to match drugs to the underlying target genes and CIVIC support in each sample
## TODO: Beware of inconsistency between tier3 in SNV and tier3 in CNV (all vs. CNV-only CIVIC records)
##       Because of this, when pooling and same drug+sample+gene has tier3 for both, CNV entries will be duplicated because they are also in SNV results

clinInfoDict = {}

## 4a) Process input files for SNV CIVIC
seenSamples_snv = []
for file in os.listdir(args.inputDir_civic_snv):
    sampleFile_snv = "%s%s" %(args.inputDir_civic_snv,os.path.basename(file))
    if os.path.isfile(sampleFile_snv) and sampleFile_snv.endswith(args.fileEnding_snv):
        name_snv = os.path.basename(file).split("_")[0].split("-")[1]
        if name_snv not in seenSamples_snv:
            seenSamples_snv.append(name_snv)
        else:
            print("Error! Multiple SNV CIVIC files for %s!" %(name_snv))
            sys.exit(1)
        (clinInfoDict) = parseInputFile(sampleFile_snv,name_snv,clinInfoDict)

## 4b) Process input files for CNV CIVIC
seenSamples_cnv = []
for file in os.listdir(args.inputDir_civic_cnv):
    sampleFile_cnv = "%s%s" %(args.inputDir_civic_cnv,os.path.basename(file))
    if os.path.isfile(sampleFile_cnv) and sampleFile_cnv.endswith(args.fileEnding_cnv):
        name_cnv = os.path.basename(file).split("_")[0].split("-")[1]
        if name_cnv not in seenSamples_cnv:
            seenSamples_cnv.append(name_cnv)
        else:
            print("Error! Multiple CNV CIVIC files for %s!" %(name_cnv))
            sys.exit(1)
        (clinInfoDict) = parseInputFile(sampleFile_cnv,name_cnv,clinInfoDict)

## SNV and CNV samples must be the same (n=412)
if (set(seenSamples_snv) != set(seenSamples_cnv)):
    print("Error! Sample names in input directories %s and %s do not match." %(args.inputDir_civic_snv,args.inputDir_civic_cnv))
    sys.exit(1)

## TODO: For now, ignore DRS information, ie. only DGIDB+CIVIC
infile = open(args.inputTable,'r')
outfile_det = open(args.outFileTag + ".details.tsv",'w') # Details table
outfile_onco = open(args.outFileTag + ".condensed.tsv",'w') # Oncoprint table
outfile_infos = open(args.outFileTag + ".civicInfos.tsv",'w') # CIVIC details table

onco_header = "Drug\tNumberSamples_genomic\tNumberSamples_drs\tNumberSamples_civic\tCancer_related\tContained_in_snv\tContained_in_drs\tContained_in_cnv\tContained_in_civic"
infos_header = "GENE\tDRUG\tCONTAINED_IN_CIVIC_RESULTS\tMUTATED_SAMPLES\tCIVIC_SAMPLES\tPOSITIVE_SAMPLES\tNEGATIVE_SAMPLES\tUNKNOWN_CONFLICT_SAMPLES\tUNKNOWN_DNS_SAMPLES\tUNKNOWN_BLANK_SAMPLES\tUNKNOWN_NOVARMATCH_SAMPLES"

# Dictionary to retrieve the corresponding sample name given a position (ie. column) in the line
dict_predPos = {}
headerPred = infile.readline().strip().split('\t')
startSampleCols = 7
for posHead in range(startSampleCols,len(headerPred)):
    ## Add current sample to header of output files
    onco_header += "\t" + headerPred[posHead]
    dict_predPos[posHead] = headerPred[posHead]

## Write complete header to output files
outfile_onco.write(onco_header + "\n")
outfile_det.write(onco_header + "\n")  # header is the same for both the oncoprint and details tables

## Classification of evidences will be counted for each gene and drug (based on structure of supportDict)
## Generate all combinations of direction + clinical significance
countDict = {}
## Use the following dictionary to classify each combination of direction+clinSignf as POSITIVE or NEGATIVE and use this information in the output header column names
colNameDict = {}
for direction in supportDict:
    for clinSignf in supportDict[direction].keys():
        supportType = direction + "_" + clinSignf
        countDict[supportType] = 0
        colNameDict[supportType] = supportDict[direction][clinSignf]
## Add 'UNKNOWN_BLANK' support as well (not present in supportDict)
## Also, this type is not classified as POSITIVE or NEGATIVE for the column names (see below)
countDict['UNKNOWN_BLANK'] = 0

# Dictionary to retrieve the corresponding support count given a position (ie. column) in the line
dict_infosPos ={}
## Generate and write header to output file
for posInfos,supportType in enumerate(countDict.keys()):
    if supportType in colNameDict.keys():
        ## Change UNKNOWN_DNS (due to DOES_NOT_SUPPORT) to simply UNKNOWN
        if colNameDict[supportType] == 'UNKNOWN_DNS':
            infos_header += "\t" + 'UNKNOWN_' + supportType + "_TOTAL"
        else:
            infos_header += "\t" + colNameDict[supportType] + '_' + supportType + "_TOTAL"
    ## For UNKNOWN (due to blanks), there is no classification as POSITIVE OR NEGATIVE, so do not add in column name
    else:
        infos_header += "\t" + supportType + "_TOTAL"
    dict_infosPos[posInfos] = supportType

## Write complete header to output file
outfile_infos.write(infos_header + "\n")


seenDrugs = []
## Iterate drug table (ie. each line corresponds to unique drug)
for line in infile:
    ## Initialize dictionary for counting classification of evidences for current drug
    drugCountMap = {}
    ## Initialize dictionary for counting classification of samples for each drug-gene interaction relating to current drug
    ## array[0] = mutSamples, array[1] = posSamples, array[2] = negSamples, array[3] = unkSamples, array[4] = unkSamples_noMatchVar,
    ## array[5] = civicSamples, array[6] = dnsSamples, array[7] = conflictSamples
    drugCountSampleMap = {}

    lineSplit = line.strip().split('\t')
    ## Retrieve drug name for current line
    predDrug = lineSplit[0].upper()
    if predDrug not in seenDrugs:
        seenDrugs.append(predDrug)
    else:
        print("Warning! Drug %s contained twice in input file." %(predDrug))
        continue

    sampleNum_genomic = lineSplit[1]
    sampleNum_drs = lineSplit[2]
    cancer_rel = lineSplit[3]
    matchedSNV = lineSplit[4]
    matchedDRS = lineSplit[5]
    matchedCNV = lineSplit[6]

    matchedCIVIC = 'n'
    ## Pooled CIVIC data for SNV and CNV
    if predDrug in clinInfoDict.keys():
        matchedCIVIC = 'y'

    ## Drug and NumberSamples columns will be added at the end
    outString_det = cancer_rel + '\t' + matchedSNV + '\t' + matchedDRS + '\t' + matchedCNV + '\t' + matchedCIVIC
    outString_onco = cancer_rel + '\t' + matchedSNV + '\t' + matchedDRS + '\t' + matchedCNV + '\t' + matchedCIVIC

    ## Keep track of total no. samples that had support for CIVIC (any type)
    samples_civic = 0

    ## Iterate samples and combine predictions for current drug
    for pos in range(startSampleCols,len(lineSplit)):
        ## Retrieve corresponding sample name for current position
        sample = dict_predPos[pos]
        predDGIDB = lineSplit[pos]
        predSupport = "no_mutation"
        ctSupport = ''
        tierSupport = '4'

        ## If at least one gene was found, this means there is a DGIDB prediction
        if 'score' in predDGIDB:
            predSupport = "dgidb_only"

        ## Parse gene(s) predicted by DGIDB for currently evaluated sample and drug (if any), and combine CIVIC support for each (if available)
        ## The input table in this case corresponds to a combination of SNV+DRS+CNV predictions, in the specified order and separated with '|'
        ## Format eg.: 'PIK3CA (multiple variants, score: 9) | ERBB2;MIR4728 (score: 5) | -0.0183 | TP53 (score: 30) | ERBB2 (score: 7)'
        predSplit = predDGIDB.split(' | ')
        ## There can be several DGIDB gene predictions within one data type (ie. snv or cnv)
        predGenes = []
        genesWithSupport = []
        predStrings = []
        skip_conflict = False
        ## Only parse DGIDB predicionts and skip DRS; anyway only possible to retrieve necessary gene info from DGIDB
        for singlePred in predSplit:
            tempPred = singlePred
            ## Skip accounting for current gene because it was already seen (due to snv+cnv mutation info)
            skip_count = False  
            if 'score' in singlePred:
                gene = singlePred.split('(')[0].strip()
                ## Check for format containing variant annotations
                ## Format: eg. CDKN2A:c.355G>A|p.Glu119Lys;c.365G>A|p.Arg122Gln;c.*245G>A| (multiple variants, score: 4)
                gene_and_var = ''
                if ":" in gene:
                    gene_and_var = singlePred.split('(')[0].strip()
                    gene = gene_and_var.split(":")[0].strip()
                if gene not in predGenes:
                    predGenes.append(gene)
                else:
                    skip_count = True

                ## Combine DGIDB prediction with CIVIC information (if available)
                ## Since many genes can be predicted for a single drug+sample, iterate all and report best CIVIC support
                tempSupport = "dgidb_only"
                ctGene = ''         # whether CIVIC info derives from ct/gt/nct
                tierGene = ''       # whether CIVIC info derives from tier 1/2/3
                ## Initialize dictionaries for counting classification of evidences and samples, for current gene (and drug)
                if gene not in drugCountMap.keys():
                    ## One entry per support type (ie. combination of direction+clinical significance), including 'UNKNOWN'
                    ## If gene has no support for CIVIC, ie. dgidb_only, then all entries will be 0
                    drugCountMap[gene] = copy.deepcopy(countDict)
                    ## For each gene-drug interaction of current drug, keep track of support across samples
                    ##  array[0] = mutSamples, array[1] = posSamples, array[2] = negSamples, array[3] = unkSamples, array[4] = unkSamples_noMatchVar,
                    ##  array[5] = civicSamples, array[6] = dnsSamples, array[7] = conflictSamples
                    ## If gene has no support for CIVIC, dgidb_only, then all entries will be 0 except first one (dgidb)
                    drugCountSampleMap[gene] = [0,0,0,0,0,0,0,0]

                ## Account this gene as predicted by DGIDB for the currently evaluated drug
                if not skip_count:
                    drugCountSampleMap[gene][0] += 1
                geneSupport = []
                ## Check whether there is any CIVIC info available for current DGIDB prediciton
                noCIVIC = False
                if predDrug not in clinInfoDict.keys():
                    noCIVIC = True
                elif sample not in clinInfoDict[predDrug].keys():
                    noCIVIC = True
                elif gene not in clinInfoDict[predDrug][sample].keys():
                    noCIVIC = True

                ## If not, then continue iteration of predictions and append (unchanged) prediction details string for reporting in final output
                if noCIVIC:
                    predStrings.append(tempPred)
                    continue

                ## Add sanity check for gene having already been evaluated for CIVIC support (ie. when same gene is mutated in SNV and CNV)
                ## Required as we will only add 'MULTIPLE' tag when several **different** genes have been evaluated for CIVIC support
                if gene not in genesWithSupport:
                    genesWithSupport.append(gene)
                ## Account this gene as having CIVIC info available for the currently evaluated drug
                if not skip_count:
                    drugCountSampleMap[gene][5] += 1

                ## In dict clinInfoDict, only available tiers are 1,2,3 (4 is excluded since no clinical info in this case)
                allTiers = list(clinInfoDict[predDrug][sample][gene].keys())
                ## Just report clinical info for a single tier (if many are available). Prioritize available tiers for this (1>2>3)
                ## Only combination of drug+sample+gene that has >1 tier -> AFATINIB + A6B1 + ERBB3 (tier 2, tier 3)
                if (len(allTiers)>1):
                    tier = min(allTiers)
                    tierGene = tier
                else:
                    tier = allTiers[0]
                    tierGene = tier

                ## Classify sample as unknown when only tier 3 is available for drug+sample+gene and no other
                if tier == '3':
                    geneSupport = ['UNKNOWN_TIER3']
                else:
                    ## Do not ignore cancer specificity but generate several 'supports' depending on the ct
                    geneSupport_ct = []
                    geneSupport_gt = []
                    geneSupport_nct = []
                    for ct in clinInfoDict[predDrug][sample][gene][tier].keys():
                        for direction in clinInfoDict[predDrug][sample][gene][tier][ct].keys():
                            for clinSig in clinInfoDict[predDrug][sample][gene][tier][ct][direction]:
                                ## Unknown is assigned when there is at least one 'N/A' or 'NULL for direction or clinical significance
                                if ('NULL' in direction) or ('N/A' in direction) or ('NULL' in clinSig) or ('N/A' in clinSig):
                                    thisSupport = 'UNKNOWN_BLANK'
                                    ## Increase by 1 the counter for unknown significance
                                    if not skip_count:
                                        drugCountMap[gene][thisSupport] += 1
                                else:
                                    if direction not in supportDict.keys():
                                        print("Error! Could not find direction %s in support dictionary." %(direction))
                                        sys.exit(1)
                                    if clinSig not in supportDict[direction].keys():
                                        print("Error! Could not find clinical significance %s in support dictionary." %(clinSig))
                                        sys.exit(1)
                                    thisSupport = supportDict[direction][clinSig]
                                    ## Increase by 1 the counter for current combination of direction+clinical significance
                                    if not skip_count:
                                        drugCountMap[gene][direction + "_" + clinSig] += 1
                                ## Split gene support by cancer-type specificity
                                ## Keep track of all available support evidence
                                if ct == 'ct':
                                    geneSupport_ct.append(thisSupport)
                                elif ct == 'gt':
                                    geneSupport_gt.append(thisSupport)
                                elif ct == 'nct':
                                    geneSupport_nct.append(thisSupport)

                    ## Once all clinical info has been iterated, choose what supporting evidence to report for current gene based on available cancer-specificity
                    ## Apply hierarchy: ct > gt > nct
                    if geneSupport_ct:
                        geneSupport = geneSupport_ct
                        ctGene = 'ct'
                    elif geneSupport_gt:
                        geneSupport = geneSupport_gt
                        ctGene = 'gt'
                    elif geneSupport_nct:
                        geneSupport = geneSupport_nct
                    else:
                        print("Error! No clinical information was found for gene %s." %(gene))
                        sys.exit(1)

                ## Counters for available CIVIC evidence items
                count_pos = 0
                count_neg = 0
                count_unk = 0
                count_dns = 0

                ## Classify sample (for given gene) depending on supporting evidence
                ## Tier 3 case: collected evidence is too disperse to be reported
                if ('UNKNOWN_TIER3' in geneSupport):
                    tempSupport = "civic_unspecific_vars"
                    if not skip_count:
                        drugCountSampleMap[gene][4] += 1

                ## Remaining classification includes POSITIVE, NEGATIVE, UNKNOWN_BLANK (ie. due to blank info) and UNKNOWN_DNS (ie. due to DOES_NOT_SUPPORT evidence)
                else:
                    count_pos = geneSupport.count('POSITIVE')
                    count_neg = geneSupport.count('NEGATIVE')
                    count_unk = geneSupport.count('UNKNOWN_BLANK')
                    count_dns = geneSupport.count('UNKNOWN_DNS')
                    ## Pool UNKNOWN_BLANK and UNKNOWN_DNS together (as both result in unknown CIVIC support)
                    count_total_unk = count_unk + count_dns
                    ## Sanity check that there is at least some support
                    if (count_pos == 0) and (count_neg == 0) and (count_total_unk == 0):
                        print("Error! Unexpected support case for gene %s." %(gene))
                        sys.exit(1)
                    ## Resolve contradicting evidence (if any) by majority vote
                    ## For this, pool UNKNOWN_BLANK and UNKNOWN_DNS together
                    if (count_total_unk > count_pos) and (count_total_unk > count_neg):
                        ## In oncoprint, simply classify samples as 'civic_unknown'
                        tempSupport = "civic_unknown"
                        ## However, for detailed table, distinguish between UNKNOWN_BLANK and UNKNOWN_DNS
                        ## If count is equal for both types of unknowns, classify as UNKNOWN_DNS
                        if not skip_count:
                            if count_dns >= count_unk:
                                drugCountSampleMap[gene][6] += 1
                            else:
                                drugCountSampleMap[gene][3] += 1
                    elif count_pos == count_neg:
                        tempSupport = "civic_conflict"
                        if not skip_count:
                            drugCountSampleMap[gene][7] += 1
                    elif (count_pos > count_neg) and (count_pos >= count_total_unk):
                        tempSupport = "civic_support"
                        if not skip_count:
                            drugCountSampleMap[gene][1] += 1
                    elif (count_neg > count_pos) and (count_neg >= count_total_unk):
                        tempSupport = "civic_resistance"
                        if not skip_count:
                            drugCountSampleMap[gene][2] += 1
                    else:
                        print("Error! Unexpected support case for drug %s, sample %s and gene %s." %(predDrug,sample,gene))
                        sys.exit(1)

                ## Decide final support to report based on the type of ct and achieved civic support
                ##  - First, decide on cancer specificity: ct > gt > nct (always take the best possible ct no matter the associated CIVIC support)
                ##  - Then, decide on civic support (use dictionary 'supportPrior' for this):
                ##      support = resistance > civic_conflict > civic_unknown > civic_unspecific_vars (tier3) > dgidb_only
                if tempSupport != "dgidb_only":
                    ## If available evidence was found to be conflicting across the mutated genes (conflicting signals), skip further classification and report this information
                    ## However, if available evidence was found to be conflicting WITHIN a gene (conflicting evidences), then consider it for prioritization of CIVIC support
                    if not skip_conflict:
                        ## Change final civic support if current gene has same or better tier
                        if (tierGene <= tierSupport):
                            ## Case of better tier, always keep it
                            if (tierGene < tierSupport):
                                ctSupport = ctGene
                                tierSupport = tierGene
                                predSupport = tempSupport
                            ## Case of same tier, change final civic support if current gene has same or better ct
                            elif (ctGene == ctSupport) or (ctGene == 'ct' and ctSupport == 'gt') or (ctGene and not ctSupport):
                                ## Case of better ct, always keep it
                                if (ctGene != ctSupport):
                                    ctSupport = ctGene
                                    tierSupport = tierGene
                                    predSupport = tempSupport
                                ## Case of same ct, decide on civic support based on priority of categories
                                else:
                                    predSupport_score = supportPrior[predSupport]
                                    tempSupport_score = supportPrior[tempSupport]
                                    ## Case of both support and resistance: conflicting signals
                                    if (predSupport_score==1 and tempSupport_score==1) and (predSupport != tempSupport):
                                        ctSupport = ctGene
                                        tierSupport = tierGene
                                        predSupport = "civic_conflict"
                                        # Do not prioritize support any further from now on
                                        skip_conflict = True
                                    ## Case of better support, always keep it
                                    elif (predSupport_score > tempSupport_score):
                                        ctSupport = ctGene
                                        tierSupport = tierGene
                                        predSupport = tempSupport
                                    ## Same support will stay the same anyway

                ## Finally, add CIVIC support info to current gene prediction details
                ## First, retrieve DGIDB score (really, whole content inside brackets) for currently evaluated DGIDB prediction
                predContent = singlePred.split('(')[1].strip()
                predContent = predContent.split(')')[0].strip()
                ctString = ctGene
                if not ctString:
                    ctString = 'nct'
                ## Format of CIVIC details string: tier | ct | no.positives | no.negatives | no.unknown_dns | no.unknown_blank
                civicString = tierGene + '|' + ctString + '|' + str(count_pos) + '|' + str(count_neg) + '|' + str(count_dns) + '|' + str(count_unk)
                ## Reconstruct initial DGIDB prediction details string, now also adding CIVIC info
                predContent += ', civic: ' + civicString
                ## Here, tempPred value is changed (if no CIVIC details were available, then temPred would correspond to original value in table)
                ## Check for format containing variant annotations
                if gene_and_var:
                    tempPred = gene_and_var + ' (' + predContent + ')'
                else:
                    tempPred = gene + ' (' + predContent + ')'

                #### **End of genes iteration (within single DGIDB prediction)

            ## Append prediction details string for reporting in final output
            ## Note that 'tempPred' may have changed due to addition to CIVIC details (if any)
            predStrings.append(tempPred)

        #### **End of DGIDB predictions iteration**

        ## Add support for current sample to output file
        ## Include whether current support was derived from ct/gt information or not
        if ctSupport:
            ctSupport = '_' + ctSupport
        outString_onco += '\t' + str(predSupport + ctSupport)
        ## Add details string (containing CIVIC support details as well when available) for current sample to output file
        outString_det += '\t' + ' | '.join(predStrings)

        ## Account for sample having some kind of CIVIC support (any) for current drug
        if len(genesWithSupport) > 0:
            samples_civic += 1
        ## Sanity check of whether >1 **different** gene per drug+patient has CIVIC support
        ## Cases of same gene being mutated in SNV and CNV do not get MULTIPLE tag
        if len(genesWithSupport) > 1:
            ## Include tag for the oncoprint to indicate that >1 gene had CIVIC support for the given drug+sample
            outString_onco += ';MULTIPLE'
            print("Warning! Found >1 gene with CIVIC support for drug %s and sample %s." %(predDrug,sample))

    #### **End of sample iteration**

    ## Add drug name and sample count for current drug (line)
    ## Write line for current drug (all samples) to output files
    outfile_onco.write("%s\t%s\t%s\t%s\t%s\n" %(predDrug,sampleNum_genomic,sampleNum_drs,samples_civic,outString_onco))
    outfile_det.write("%s\t%s\t%s\t%s\t%s\n" %(predDrug,sampleNum_genomic,sampleNum_drs,samples_civic,outString_det))

    ## Write to output aggregated info across all samples and associated genes, for current drug
    for gene in drugCountMap.keys():
        geneCounts = drugCountMap[gene]
        ## array[0] = mutSamples, array[1] = posSamples, array[2] = negSamples, array[3] = unkSamples, array[4] = unkSamples_noMatchVar,
        ## array[5] = civicSamples, array[6] = dnsSamples, array[7] = conflictSamples
        geneSampleCounts = drugCountSampleMap[gene]
        ## If there is at least one sample with CIVIC information,
        ## that means current gene+drug combination is present in CIVIC
        contained_CIVIC = 'n'
        if geneSampleCounts[5]>0:
            contained_CIVIC = 'y'
        outString_infos = gene + '\t' + predDrug + '\t' + contained_CIVIC + '\t' + str(geneSampleCounts[0]) + '\t' + str(geneSampleCounts[5]) + '\t' + str(geneSampleCounts[1]) + '\t' + str(geneSampleCounts[2]) + '\t' + str(geneSampleCounts[7]) + '\t' + str(geneSampleCounts[6]) + '\t' + str(geneSampleCounts[3]) + '\t' +  str(geneSampleCounts[4])
        for i in range(0,len(dict_infosPos)):
            supportType = dict_infosPos[i]
            supportCount = geneCounts[supportType]
            outString_infos += '\t' + str(supportCount)

        ## Add line for support of current gene+drug to output file
        outfile_infos.write(outString_infos + '\n')

outfile_onco.close()
outfile_det.close()
outfile_infos.close()
print("\nParsed %s predicted drugs for %s samples." %(len(seenDrugs),len(dict_predPos.keys())))
