#!/usr/bin/env python

'''
Use CIViCutils to query and process CIViC variant information
Requires newer format of variant annotations, i.e.: GENE:variant|variant;GENE:variant;...
Version queries offline cache provided in civicpy module (updating cache is managed by Vipin)

Lourdes Rosano, Feb 2022
'''

import sys
import argparse

### Load relevant functions from CIViCutils package

sys.path.append('/cluster/work/nexus/lourdes/civic_query/git/civicutils/')

from read_and_write import get_dict_support,write_evidences,write_output_line
from query import query_civic
from filtering import filter_civic
from match import match_in_civic,annotate_ct,filter_ct,process_drug_support,reprocess_drug_support_across_selected_variants
from utils import check_match_before_writing,check_keys,check_keys_not,check_data_type,check_dict_entry


'''
Functions
'''

## Given the header, find the input columns containing genes and variants/cnvs values (use provided arguments directly)
## For SNVs, also retrieve the input columns containing variant impacts and exons
def getColumnIndices(header, dataType):
    header_split = header.split('\t')
    index_geneCol = -1
    index_varCol = -1
    index_impCol = -1
    index_exonCol = -1
    for pos in range(0,len(header_split)):
        # Avoid mismatches due to case by always using uppercase
        if args.colName_gene.upper() == header_split[pos].upper():
            index_geneCol = pos
        if args.colName_data.upper() == header_split[pos].upper():
            index_varCol = pos
        if dataType == 'SNV':
            if args.colName_varImpact.upper() == header_split[pos].upper():
                index_impCol = pos
            if args.colName_exon.upper() == header_split[pos].upper():
                index_exonCol = pos
    if (index_geneCol == -1) or (index_varCol == -1):
        print("\nError! Could not match all input columns in header %s." %(header))
        sys.exit(1)
    if dataType == 'SNV':
        if (index_impCol == -1) or (index_exonCol == -1):
            print("\nError! Could not match all input columns in header %s." %(header))
            sys.exit(1)
    return (index_geneCol, index_varCol, index_impCol, index_exonCol)


def readInSnvsMultipleGenes(input_snv_file):
    # lineNumber -> [[var, gene,dna,prot,impact,exon],[var, gene,dna,prot,impact,exon],..]
    rawData = {}
    # gene -> variant (dna|prot|impact|exon|lineNumber) -> null
    snvData = {}

    infile = open(input_snv_file,'r')
    firstInputLine = infile.readline().strip()
    ## Retrieve indices for all relevant columns (ie. gene and variant annotation; for SNV also variant impact and exon info)
    ## If not all are found, script will exit with an error
    ## When dataType is CNV, index_impCol=-1 and index_exonCol=-1 since they will be disregared
    (index_geneCol,index_varCol,index_impCol,index_exonCol) = getColumnIndices(firstInputLine, dataType)

    # Each parsed line corresponds to a different genomic variant (which can have several annotations available)
    for nLine,line in enumerate(infile):
        lineSplit = line.strip().split("\t")
        # FIXME: do uppercase? -> ".strip().upper().split(";")"?
        all_variants = lineSplit[index_varCol].strip().split(";")
        all_impacts = lineSplit[index_impCol].strip().split(";")
        all_exons = lineSplit[index_exonCol].strip().split(";")

        # Sanity check the the 3 columns contain the same number of variants annotations, sorted in the same order (1-1 correspondance of positions)
        # Adapted from: https://stackoverflow.com/questions/35791051/better-way-to-check-if-all-lists-in-a-list-are-the-same-length
        check_lists = [all_variants, all_impacts, all_exons]
        it = iter(check_lists)
        n_variants = len(next(it))
        if not all(len(l) == n_variants for l in it):
            raise ValueError("Encountered different number of available annotations in columns '%s', '%s' and '%s'!" %(args.colName_data, args.colName_varImpact, args.colName_exon))

        # lineNumber -> [var1,..,varN]}
        rawData[str(nLine)] = []

        # Loop all available items (same number), and retrieve for each variant corresponding info from each list
        # For each element, retrieve single gene (separated by ':'), retrieve single c./g./n. annotation, and retrieve single p. annotation when available
        # For each element, retrieve the corresponding variant impact (separated from gene by ':') and exon info (separated from gene and variant by ':')
        for n_variant in range(0, n_variants):
            variant = all_variants[n_variant]
            variant_split = variant.strip().split(":")
            # FIXME: do uppercase? -> ".strip().upper()"?
            gene = variant_split[0].strip()
            hgvs_split = variant_split[1].strip().split("|")
            cVar = hgvs_split[0].strip()
            pVar = hgvs_split[1].strip()

            # Sanity check required information is never empty
            # NOTE: pVar is considered optional
            if not gene:
                raise ValueError("Encountered empty gene annotation in column '%s'!" %(args.colName_data))
            if not cVar:
                raise ValueError("Encountered empty variant annotation in column '%s'!" %(args.colName_data))

            impact_tag = all_impacts[n_variant]
            impact_split = impact_tag.strip().split(":")
            if impact_split[0].strip() != gene:
                raise ValueError("Encountered different genes in positional annotations of columns '%s' and '%s'!" %(args.colName_data, args.colName_varImpact))
            impact = impact_split[1].strip()

            exon_tag = all_exons[n_variant]
            exon_split = exon_tag.strip().split(":")
            if exon_split[0].strip() != gene:
                raise ValueError("Encountered different genes in positional annotations of columns '%s' and '%s'!" %(args.colName_data, args.colName_exon))
            if exon_split[1].strip() != variant_split[1].strip():
                raise ValueError("Encountered different variants in positional annotations of columns '%s' and '%s'!" %(args.colName_data, args.colName_exon))
            exon = exon_split[2].strip()

            # Sanity check required information is never empty
            # NOTE: exon is considered optional
            if not impact:
                raise ValueError("Encountered empty variant impact annotation in column '%s'!" %(args.colName_varImpact))

            # NOTE: for now, skip variants which do not have a valid c.HGVS expression available (e.g. n.HGVS)
            # We can do this because we will use alternative custom function 'write_match_multiple_annotations' to write results to output, and these variant annotations will not be removed from the affected columns
            if not cVar.startswith("c."):
                continue

            # Process rawData to have gene-centered dict
            # Returns dict of gene -> [var1,var2,..,varN], where a given var="dna|prot|impact|exon|lineNumber"
            if gene not in snvData.keys():
                snvData[gene] = {}
            # Collapse variant info separated with "|"
            # Keep track of what line each variant comes from
            variant = cVar + "|" + pVar + "|" + impact + "|" + exon + "|" + str(nLine)
        # NOTE: Variants can never be duplicated because of the different row numbers assigned
        # if variant in snvData[gene].keys():
        #     print("Found duplicated variant '%s|%s' for gene '%s' in line '%s'!" %(cVar,pVar,gene,str(nLine)))
        #     sys.exit(1)
            snvData[gene][variant] = None
            rawData[str(nLine)].append([variant, gene, cVar, pVar, impact, exon])
    infile.close()
    return (rawData, snvData, firstInputLine)


def readInCnvsMultipleGenes(input_cnv_file):
    # lineNumber -> {details -> cnv, genes -> [gene1,..,geneN]}
    rawData = {}
    # gene -> variant (cnv|lineNumber) -> null
    cnvData = {}

    infile = open(input_cnv_file,'r')
    firstInputLine = infile.readline().strip()
    ## Retrieve indices for all relevant columns (ie. gene and variant annotation; for SNV also variant impact and exon info)
    ## If not all are found, script will exit with an error
    ## When dataType is CNV, index_impCol=-1 and index_exonCol=-1 since they will be disregared
    (index_geneCol,index_varCol,index_impCol,index_exonCol) = getColumnIndices(firstInputLine, dataType)

    # Each parsed line corresponds to a different CNV (which can have several genes annotated)
    for nLine,line in enumerate(infile):
        lineSplit = line.strip().split("\t")
        cnv = lineSplit[index_varCol].strip()
        # lineNumber -> {details -> cnv, genes -> [gene1,..,geneN]}
        rawData[str(nLine)] = {}
        rawData[str(nLine)]["details"] = cnv
        rawData[str(nLine)]["genes"] = []
        # NOTE: in this case, the same CNV can contain several annotated genes separated with ';' in the Gene column
        genes = lineSplit[index_geneCol].strip().split(";")
        # Collapse variant info separated with "|"
        # Keep track of what line each variant comes from
        variant = cnv + "|" + str(nLine)
        for gene in genes:
            rawData[str(nLine)]["genes"].append(gene)
            # Process rawData to have gene-centered dict
            # Returns dict of gene -> [var1,var2,..,varN], where a given var="cnv|lineNumber"
            if gene not in cnvData.keys():
                cnvData[gene] = {}
        # NOTE: Variants can never be duplicated because of the different row numbers assigned
        # if variant in cnvData[gene].keys():
        #     print("Found duplicated variant '%s' for gene '%s' in line '%s'!" %(cnv,gene,str(nLine)))
        #     sys.exit(1)
            cnvData[gene][variant] = None
    infile.close()
    return (rawData, cnvData, firstInputLine)


def check_tier_and_matches(gene, combId, matchMap, hasSupport=True):
    sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]
    # Check if matchMap contains the provided input variants
    if gene not in matchMap.keys():
        raise ValueError("Provided gene '%s' is not contained in provided 'matchMap'." %(gene))
    if combId not in matchMap[gene].keys():
        raise ValueError("Provided variant '%s' is not contained in supplied 'matchMap' for gene '%s'." %(combId, gene))

    all_tiers = list(matchMap[gene][combId].keys())
    selected_tier = "tier_4"
    # Pick highest tier available
    for tmp_tier in all_tiers:
        # Only look into tiers that have at least 1 match in CIVIC
        skip_tier = False
        if hasSupport:
            if not matchMap[gene][combId][tmp_tier]["matched"]:
                skip_tier = True
        else:
            if not matchMap[gene][combId][tmp_tier]:
                skip_tier = True
        if skip_tier:
            continue

        # Compare the currently evaluated tier with the highest one encountered (so far)
        # NOTE: use the order of priority already defined in 'sorted_tiers'
        tmp_tier_index = sorted_tiers.index(tmp_tier)
        tier_index = sorted_tiers.index(selected_tier)
        if tmp_tier_index < tier_index:
            selected_tier = tmp_tier

    # Return all variants matched for the selected tier and the provided gene+variant
    allVariants = []
    if selected_tier != "tier_4":
        if hasSupport:
            for tmpVar in matchMap[gene][combId][selected_tier]["matched"]:
                allVariants.append(tmpVar)
        else:
            for tmpVar in matchMap[gene][combId][selected_tier]:
                allVariants.append(tmpVar)
    return (selected_tier, allVariants)


def write_match_multiple_annotations(matchMap, varMap, rawMap, input_file, outfile, dataType="SNV", hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False):
    # NOTE: uppercase is critical for matching!
    sorted_evidence_types = ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']
    evidenceType = "PREDICTIVE"
    special_cases = ["NON_SNV_MATCH_ONLY","NON_CNV_MATCH_ONLY","NON_EXPR_MATCH_ONLY"]
    sorted_cts = ["ct","gt","nct"]
    sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]

    # Use CIViCutils functionality to sanity check that all expected data and formats are correct
    check_match_before_writing(matchMap, varMap, rawMap, hasSupport, hasCt, writeCt, writeSupport, writeComplete)

    if writeSupport:
        if not hasSupport:
            raise ValueError("Option 'writeSupport' cannot be selected when 'hasSupport'=False!")
    if writeCt:
        if not hasCt:
            raise ValueError("Option 'writeCt' cannot be selected when 'hasCt'=False!")

    # Keep track of all matches and non-matches
    exact_matches = 0           # Tier 1
    syn_matches = 0             # Tier 1b
    pos_matches = 0             # Tier 2
    no_matches = 0              # Tier 3 (includes special cases when there is no match of the same variant type returned)
    gene_not_found = 0          # Tier 4

    # Iterate through the input table once more and simultaneously write results to output table
    # Each parsed line corresponds to a different variant
    infile = open(input_file,'r')
    ignore_header = infile.readline().strip()
    for nLine,line in enumerate(infile):
        nLine = str(nLine)
        if nLine not in rawMap.keys():
            raise ValueError("Line %s could not be found in provided 'rawMap'!" %(nLine))

        tier_to_write = -1
        # gene -> variant
        genes_and_variants_to_write = {}

        # Process and prioritize CIVIC results matched across all variants and genes annotated for the current SNV line
        # NOTE: when reporting results for a variant line, make sure to only report one set of results (i.e. corresponding to one single variant annotation c. + p., if avail). This is to ensure that the same CIVIC variant record is not reported multiple times for many variant annotations (as they all refer to the same variant, this would bias the available evidence)
        #   1. if none have CIVIC, then variant line is tier4
        #   2. if only one variant annotation has CIVIC, then report the corresponding tier (only 1 since we already applied highest priority filter)
        #   3. if >1 variant annotations have CIVIC:
        #    a. choose the one with the highest tier
        #    b. if >1 with the same tier, pick the one with the highest number of CIVIC variant matches
        #   c. if same number of CIVIC variants, choose the first variant annotation encountered (random choice)
        if dataType == "SNV":
            nested_variant_list = rawMap[nLine]
            selected_tier = "tier_4"
            n_variants = -1
            selected_variant = []
            for variant_list in nested_variant_list:
                if (len(variant_list) != 6):
                    raise ValueError("Must provide 6 elements to describe a SNV variant (even if some can be empty): varId,gene,dna,[prot],impact,[exon],..'")
                combId = variant_list[0]
                gene = variant_list[1]
                cVar = variant_list[2]
                pVar = variant_list[3]
                impact = variant_list[4]
                exon = variant_list[5]

                # NOTE: for now, skip variants which do not have a valid c.HGVS expression available (e.g. n.HGVS)
                # We can do this because we will use alternative custom function 'write_match_multiple_annotations' to write results to output, and these variant annotations will not be removed from the affected columns
                if not cVar.startswith("c."):
                    continue

                (tmp_tier, all_variants) = check_tier_and_matches(gene, combId, matchMap, hasSupport)
                # Compare the currently evaluated tier with the highest one encountered (so far)
                # NOTE: use the order of priority already defined in 'sorted_tiers'
                prioritize_tier = False
                tmp_tier_index = sorted_tiers.index(tmp_tier)
                tier_index = sorted_tiers.index(selected_tier)
                # Keep the highest tier available
                if tmp_tier_index < tier_index:
                    prioritize_tier = True
                # If both tiers are the same, keep the variant annotation with the highest number of CIVIC variant ids matched
                elif tmp_tier_index == tier_index:
                    if len(all_variants) > n_variants:
                        prioritize_tier = True
                # If both the tier and the number of matched variants are the same, then pick the first variant annotation encountered (random choice)
                if prioritize_tier:
                    selected_tier = tmp_tier
                    selected_variant = [gene, combId]
                    n_variants = len(all_variants)

# FIXME: sanity check potential issues due to skipping non c.HGVS variant annotations -> what would happen if all annotations for the line are n.HGVS? This would likely trigger errors, so we need to ensure that in these cases a tier4 is reported? what other alternative for the output table?

            if not selected_variant:
                raise ValueError("Could not select a single valid variant annotation for line %s" %(line.strip()))

            # In the end, we only select one single variant annotation per SNV line (even if >1 genes are annotated)
            # FIXME: improve check to pick one single variant annotation per available gene -> here, we are sure that the corresponding variants matched in CIVIC will never be duplicated (because of the different gene)
            tier_to_write = selected_tier
            final_gene = selected_variant[0]
            final_variant = selected_variant[1]
            # Keep track of the gene and variant annotation selected for the current line, as well as all the corresponding variant ids matched in CIVIC
            genes_and_variants_to_write[final_gene] = {}
            genes_and_variants_to_write[final_gene][final_variant] = []
            (tmp_tier, all_variants) = check_tier_and_matches(final_gene, final_variant, matchMap, hasSupport)
            for tmp_variant in all_variants:
                genes_and_variants_to_write[final_gene][final_variant].append(tmp_variant)


        # Process and prioritize CIVIC results matched across all genes annotated for the current CNV line
        # NOTE: when reporting results for a CNV line, sanity check there are no duplicated genes, and report all matched CIVIC information available for the highest tier (as they will always correspond to different genes, hence the retrieved CIVIC records will never be duplicated for a given CNV type and tier)
        if dataType == "CNV":
            if ("details" not in rawMap[nLine].keys()) or ("genes" not in rawMap[nLine].keys()):
                raise ValueError("Must provide 2 elements to describe a CNV variant: gene,cnv,..'")
            cnv = rawMap[nLine]["details"]
            combId = cnv + "|" + nLine
            all_genes = rawMap[nLine]["genes"]
            # Classify all annotated genes for the current CNV according to the corresponding tier from CIVICutils
            tmp_cnv_results = {}
            for gene in all_genes:
                (tmp_tier, all_variants) = check_tier_and_matches(gene, combId, matchMap, hasSupport)
                if tmp_tier not in tmp_cnv_results.keys():
                    tmp_cnv_results[tmp_tier] = []
                if gene in tmp_cnv_results[tmp_tier]:
                    raise ValueError("Encountered duplicated gene annotation in line %s" %(line.strip()))
                tmp_cnv_results[tmp_tier].append(gene)

            # Select the highest available tier for the current line
            all_available_tiers = list(tmp_cnv_results.keys())
            selected_tier = "tier_4"
            for tmp_avail_tier in all_available_tiers:
                # Compare the currently evaluated tier with the highest one encountered (so far)
                # NOTE: use the order of priority already defined in 'sorted_tiers'
                tmp_tier_index = sorted_tiers.index(tmp_avail_tier)
                tier_index = sorted_tiers.index(selected_tier)
                if tmp_tier_index < tier_index:
                    selected_tier = tmp_avail_tier

            # Once the highest (single) tier classification has been selected, keep track of all the annotated genes which have said tier, and keep track of their corresponding variant ids matched in CIVIC
            tier_to_write = selected_tier
            for final_gene in tmp_cnv_results[selected_tier]:
                if final_gene not in genes_and_variants_to_write.keys():
                    genes_and_variants_to_write[final_gene] = {}
                    genes_and_variants_to_write[final_gene][combId] = []
                (tier, all_variants) = check_tier_and_matches(final_gene, combId, matchMap, hasSupport)
                for tmp_variant in all_variants:
                    genes_and_variants_to_write[final_gene][combId].append(tmp_variant)


        ## At this point, the genes and variant annotations which will be reported for the current variant line have already been parsed and selected
        geneScores = []
        geneVarTypes = []
        drugSupport = []
        resultMap = {}
        for gene in genes_and_variants_to_write.keys():
            # Expectation is that there is only 1 variant available per gene
            for variant_annotation in genes_and_variants_to_write[gene].keys():
                all_civic_variants = genes_and_variants_to_write[gene][variant_annotation]
                if tier_to_write == "tier_4":
                    if all_civic_variants:
                        raise ValueError("Unexpectedly found matched variants for a line classified as 'tier_4': %s" %(all_civic_variants))
                elif hasSupport:
                    drug_support = matchMap[gene][variant_annotation][tier_to_write]["drug_support"]
                    if writeSupport:
                        for i in drug_support:
                            drugSupport.append(i.upper())

                for varId in all_civic_variants:
                    # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                    # This is a dummy tag and not an actual variant record from CIVICdb, so skip checking in varMap
                    if varId.upper() in special_cases:
                        # In this case, current line will be associated with tier3, but all columns will be empty with "."
                        for evidence_type in sorted_evidence_types:
                            if evidence_type not in resultMap.keys():
                                resultMap[evidence_type] = []
                        continue

                    civic_variant = varMap[gene][varId]['name']
                    geneScores.append(gene + ':' + civic_variant + ':' + str(varMap[gene][varId]['civic_score']))
                    geneVarTypes.append(gene + ':' + civic_variant + ':' + ','.join(varMap[gene][varId]['types']))
                    for evidence_type in sorted_evidence_types:
                        if evidence_type in varMap[gene][varId]["evidence_items"].keys():
                            if evidence_type not in resultMap.keys():
                                resultMap[evidence_type] = []
                            writeDrug = False
                            if evidence_type == evidenceType:
                                writeDrug=True
                            if hasCt:
                                check_keys(list(varMap[gene][varId]["evidence_items"][evidence_type].keys()), "varMap", sorted_cts, matches_all=True)
                                for ct in varMap[gene][varId]["evidence_items"][evidence_type].keys():
                                    if writeCt:
                                        results = write_evidences(varMap[gene][varId]["evidence_items"][evidence_type][ct], writeDrug=writeDrug, writeCt=ct, writeComplete=writeComplete)
                                    else:
                                        results = write_evidences(varMap[gene][varId]["evidence_items"][evidence_type][ct], writeDrug=writeDrug, writeCt=None, writeComplete=writeComplete)
                                    for x in results:
                                        resultMap[evidence_type].append(gene + ":" + civic_variant + ":" + x)
                            else:
                                check_keys_not(list(varMap[gene][varId]["evidence_items"][evidence_type].keys()),"varMap", sorted_cts)
                                if writeCt:
                                    raise ValueError("Option 'writeCt' cannot be selected when 'hasCt'=False!")
                                results = write_evidences(varMap[gene][varId]["evidence_items"][evidence_type], writeDrug=writeDrug, writeCt=None, writeComplete=writeComplete)
                                for x in results:
                                    resultMap[evidence_type].append(gene + ":" + civic_variant + ":" + x)

        # NOTE: in the above case, if CIVIC info is reported across multiple different variant annotations within the same line (e.g. case of the CNVs with >1 genes having CIVIC info with the same tier), then the drug support will need to be recomputed to aggregate the information across all reported variants

        # Use custom function to recompute final consensus drug support for the current line
        # Now, support is aggregated across all genes and variants for which CIVIC info is to be reported
        # NOTE: currently, for SNV data, only 1 single gene and variant annotation is selected per line
        # NOTE: currently, for CNV data, multiple genes annotated for the same CNV can be reported with CIVIC info (only those having the same, highest tier)
        drugSupport = []
        drugSupport = reprocess_drug_support_across_selected_variants(genes_and_variants_to_write, matchMap, varMap, supportDict, hasSupport=True)

        ### Write result to output
        outLine = write_output_line(tier_to_write, line.strip(), geneScores, geneVarTypes, drugSupport, resultMap, writeSupport)

        ### Once here, all genes and variants within the row have been classified into tiers
        if tier_to_write == "tier_4":
            gene_not_found += 1
        elif (tier_to_write == "tier_3") or (tier_to_write in special_cases):
            no_matches += 1
        elif tier_to_write == "tier_2":
            pos_matches += 1
        elif tier_to_write == "tier_1b":
            syn_matches += 1
        elif tier_to_write == "tier_1":
            exact_matches += 1

        outfile.write(outLine + "\n")

    infile.close()
    return (exact_matches, syn_matches, pos_matches, no_matches, gene_not_found)

# NOTE: for now, always write column for consensus drug support
def write_header(outfile, writeSupport=True):
    sorted_evidence_types = ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']
    if writeSupport:
        out_header = "%s\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\tCIViC_Drug_Support\t%s" %(header, "\t".join(["CIViC_" + x for x in sorted_evidence_types]))
    else:
        out_header = "%s\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\t%s" %(header, "\t".join(["CIViC_" + x for x in sorted_evidence_types]))
    outfile.write(out_header + "\n")
    return None


'''
Script
'''

parser = argparse.ArgumentParser(description='Query CIViC to retrieve drug information for snvs or cnvs.')
parser.add_argument('--inputTable', dest='inputFile', required=True, help='Input table with genes and variants, needs to be tab separated.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of the output file.')
parser.add_argument('--cancerTypeList', dest='cancerTypeList', required=True, help='Comma-separated list of accepted cancer types. Partial matches will be sought.')
parser.add_argument('--blackList', dest='blackList', required=True, help='Comma-separated list of not accepted cancer types. Partial matches will be sought, and those matched will be excluded from the list of considered evidence.')
parser.add_argument('--highLevelList', dest='highLevelList', required=True, help='Comma-separated list of high level cancer types (e.g. Cancer). Only exact matches will be sought. The user should be aware that results for high level cancer types will only be retrieved when no match for cancer specific types (ie. --cancerTypeList) is found.')
parser.add_argument('--colName_gene', dest='colName_gene', required=True, help='Name of column containing gene symbols.')
parser.add_argument('--colName_data', dest='colName_data', required=True, help='Name of column containing variant annotations (SNV) or cnv categories (CNV).')
parser.add_argument('--dataType', dest='dataType', required=True, choices=['snv','cnv'], default='snv', help='Type of data in inputTable. Possible options are: snv, cnv.')
parser.add_argument('--colName_varImpact', dest='colName_varImpact', required=False, help='Name of column containing variant impacts (for SNV only).')
parser.add_argument('--colName_exon', dest='colName_exon', required=False, help='Name of column containing variant exon information (for SNV only).')

args = parser.parse_args()

## Check input parameters
dataType = args.dataType.upper() # SNV or CNV

## Check expected data type using CIViCutils functionality
check_data_type(dataType)

## For SNV, two additional arguments are required (ie. column names for variant impact and exon information)
if (dataType=='SNV') and (args.colName_varImpact is None or args.colName_exon is None):
    parser.error("--dataType 'snv' requires --colName_varImpact and --colName_exon.")

print("\nParameters:\n inputTable: %s\n outFile: %s\n dataType: %s\n colName_gene: %s\n colName_data: %s\n colName_varImpact: %s\n colName_exon: %s\n" %(args.inputFile,args.outFile,args.dataType,args.colName_gene,args.colName_data,args.colName_varImpact,args.colName_exon))


cancerTypeList = args.cancerTypeList
blackList = args.blackList
highLevelList = args.highLevelList
cancerTypeList = [s.strip().upper() for s in cancerTypeList.split(',')]
blackList = [s.strip().upper() for s in blackList.split(',')]
highLevelList = [s.strip().upper() for s in highLevelList.split(',')]

print('\nWhite listed cancer types: {}'.format(','.join(cancerTypeList)))
print('Black listed cancer types: {}'.format(','.join(blackList)))
print('High level cancer types: {}'.format(','.join(highLevelList)))

## Sanity check for "empty" input lists (ie. in the form of [''])
## Turn them into "real" empty lists for implementation purposes
if cancerTypeList == ['']:
    cancerTypeList = []
if blackList == ['']:
    blackList = []
if highLevelList == ['']:
    highLevelList = []


# Read in and process file of input SNV variants
if dataType=='SNV':
    (rawData, varData, header) = readInSnvsMultipleGenes(args.inputFile)
# Read in and process file of input CNV variants
if dataType=='CNV':
    (rawData, varData, header) = readInCnvsMultipleGenes(args.inputFile)


# Already write new header into output file
outfile = open(args.outFile,'w')
write_header(outfile, writeSupport=True)

# Extract all genes parsed in the provided input file
genes = list(varData.keys())

# When no genes were found in the input file, write empty files (necessary for snakemake pipeline) and exit without error
# eg. empty input file containing only header because patient had no variants at all
if not genes:
    print("\nDid not find any genes in column \'{}\' of file {}".format(args.colName_gene, args.inputFile))
    outfile.close()
    sys.exit(0)


### Query CIViC for the genes of interest

print('\nTotal # genes to query: {}'.format(len(genes)))
print('\nRetrieving data from CIViC...')
# Query input genes in CIVIC
varMap = query_civic(genes, identifier_type="entrez_symbol")


# (varMap,retrieved_genes,no_variants,all_variants) = query_civic_genes(genes, identifier_type="entrez_symbol", dataType=dataType)
# retrieved = set(retrieved_genes)
# unmatched = list(set(genes) - retrieved)
# print("Found %s/%s genes in CIVIC associated to %s variants. Found %s CIVIC genes that had no variants available." %(len(retrieved_genes),len(genes),len(all_variants),len(no_variants)))
# print('\nGenes with no CIVIC data: {}'.format(','.join(unmatched)))

# Filter undesired evidences to avoid matching later on
varMap = filter_civic(varMap, evidence_status_in = ['ACCEPTED'], var_origin_not_in = ['GERMLINE'], output_empty=False)


# Match input SNV variants in CIVIC, pick highest tier available per input gene+variant
# Tier hierarchy: 1 > 1b > 2 > 3 > 4
(matchMap, matchedIds, varMap) = match_in_civic(varData, dataType=dataType, identifier_type="entrez_symbol", select_tier="highest", varMap=varMap)


# Annotate matched CIVIC evidences with cancer specificity of the associated diseases
disease_name_not_in = []
disease_name_in = ['bladder']
alt_disease_names = ['solid tumor']
annotMap = annotate_ct(varMap, disease_name_not_in, disease_name_in, alt_disease_names)

# Filter CIVIC evidences to pick only those for the highest cancer specificity available
# ct hierarchy: ct > gt > nct
annotMap = filter_ct(annotMap, select_ct="highest")

# Get custom dictionary of support from data.yml (provided within the package)
# This defines how each combination of evidence direction + clinical significance in CIVIC is classified in terms of drug support (eg. sensitivity, resistance, unknown, etc.)
supportDict = get_dict_support()

# Process drug support of the matched variants using the annotated CIVIC evidences
annotMatch = process_drug_support(matchMap, annotMap, supportDict)

# Write to output
# Parse input file again, now checking (and prioritizing, if necessary) the available CIVIC info per variant line
# Report the CT classification of each disease, and write column with the overall drug support of the match for each available CT class
(exact_matches, syn_matches, pos_matches, no_matches, gene_not_found) = write_match_multiple_annotations(annotMatch, annotMap, rawData, args.inputFile, outfile, dataType=dataType, hasSupport=True, hasCt=True, writeCt=True, writeSupport=True, writeComplete=True)

outfile.close()

print('\nTotal # exact matches: {}'.format(exact_matches))
print('Total # synonym matches: {}'.format(syn_matches))
print('Total # positional matches: {}'.format(pos_matches))
print('Total # gene-only matches: {}'.format(no_matches))
stringNotFound = 'Total # rows without variant data: {}'
if dataType == 'CNV':
    stringNotFound = 'Total # genes without variant data: {}'
print(stringNotFound.format(gene_not_found))
print('---------------------')
