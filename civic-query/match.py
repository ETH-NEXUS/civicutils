import sys
import os
import re

# TODO: import all functions
import .utils
from .query import query_civic#,reformat_civic?
# TODO
from .filter import


# TODO: expand framework to provide variant type (SNV or CNV) along with each individual variant, always generate matching strings for both types, and match to one depending on retrieved type -> would allow mixing SNV and CNV variants in the same file/input object
# Note that rowMap will have a different structure depending on dataType:
#  - For SNV: gene -> annotationType -> value (either SNV annotation, variant impact or variant exon)
#  - For CNV: gene -> value (CNV category)


# Given a CIVIC variant name and its corresponding hgvs expressions (if available),
# generate a list of potential strings that will be used to match the variant to our input
# For CNVs, variant matching is not based on HGVS, since input data is different
def generate_civic_matchStrings(varName, hgvsExpressions, dataType):
    matchStrings = []
    ## For CNVs, only the last step 5) is executed, since variant matching is not done at the HGVS level
    if dataType == 'SNV':
        # 1) First, remove reference from annotation (ie. 'transcriptID:')
        # This step will be skipped for CNVs
        for x in hgvsExpressions:
            # CIViC HGVS expressions are of the form reference + single annotation
            # eg. [NM_007313.2:c.1001C>T, NP_005148.2:p.Thr315Ile]
            new = x.split(':')[-1].upper()
            if new not in matchStrings:
                matchStrings.append(new)
                # 2) Second, generate modified HGVS strings that comply with input table format
                newAnnot = civic_hgvs_to_input(new)
                # Resulting annotation will be None unless something was modified (special cases only)
                if (newAnnot is not None) and (newAnnot not in matchStrings):
                    matchStrings.append(newAnnot)
        # 3) Third, generate a list of potential HGVS strings based on the actual CIVIC variant name
        potentialHGVS = civicName_to_hgvs(varName)
        for x in potentialHGVS:
            # Sometimes, duplicated HGVS annotations will get generated (eg. V600E already has p.Val600Glu)
            if x not in matchStrings:
                matchStrings.append(x)
        # 4) Last, add potential positional matches for already existing strings (will only affect p. annotations)
        for x in matchStrings:
            start = extract_p_start(x)
            # Only p. annotations will return a positional string (eg. p.Val600Glu -> p.Val600)
            if (start is not None) and (start not in matchStrings): 
                matchStrings.append(start)
    # 5) Also, add CIVIC variant name to allow for matching using input 'descriptional' strings (eg. EXON 15 MUTATION)
    ## Only this step is executed for CNV
    matchStrings.append(varName)

    return matchStrings



# Given a set of one or more annotations (SNV or CNV), generate a list of potential strings
# that will be used to match the variant in CIVIC eg. EXON 15 MUTATION, AMPLIFICATION
# Only data types allowed: "SNV","CNV" (strings generated are different in each case)
# Arguments impactAnnots and exonAnnots are optional for "SNV" argument. When provided, variant impacts will be used to generate additional potential record names. If exon annotations are provided, then impacts must also be provided, and there must be a 1-1 correspondance with this list (to determine if variant is intronic or exonic, from the impact)
def generate_input_matchStrings(varAnnotations, dataType, impactAnnots=[], exonAnnots=[]):
    matchStrings = []
    # Used for SNV: Keep track of whether a string corresponds to an exact (True) or positional match (False)
    isExact = []
    # Used for SNV: isTrueExact can be used together with isExact to distinguish between a true exact match
    # (isExact=True, isTrueExact=True) corresponding to input HGVS strings, or a more general match
    # (isExact=True, isTrueExact=False) corresponding to descriptive terms eg. EXON 1 MUTATION with lower preference
    isTrueExact = []

    ## a) Given a set of HGVS expressions, impacts and exon information for a single variant, generate a list of
    ## potential strings that will be used to match the variant in CIVIC eg. EXON 15 MUTATION, TRUNCATING FRAMESHIFT
    if dataType == 'SNV':
        # 1) First, add all (unique) HGVS annotations gathered from input table
        # If matched, they will correspond to exact matches (True in isExact)
        for varAnnot in varAnnotations:
            if varAnnot not in matchStrings:
                matchStrings.append(varAnnot)
                isExact.append(True)
                isTrueExact.append(True)
                # Special case for protein extensions (eg. p.Ter130Tyrext*?) in input table
                # Subset string to match CIVIC's format (p.Ter130Tyr)
                if re.match(r'(P\.TER[0-9]+[A-Z]+)EXT', varAnnot):
                    newVarAnnot = re.match(r'(P\.TER[0-9]+[A-Z]+)EXT', varAnnot).groups()[0]
                    if newVarAnnot not in matchStrings:
                        matchStrings.append(newVarAnnot)
                        isExact.append(True)
                        isTrueExact.append(True)
        # 2) Second, add potential positional matches for already existing strings (will only affect p. annotations)
        # If matched, they will correspond to positional matches (False in isExact)
        for x in matchStrings:
            start = extract_p_start(x)
            # Only p. annotations will return a positional string (eg. p.Val600Glu -> p.Val600)
            if (start is not None) and (start not in matchStrings):
                matchStrings.append(start)
                isExact.append(False)
                isTrueExact.append(False)

        # 3) Last, add potential variant synonyms for matching to CIVIC's record names (eg. EXON 15 MUTATION)
        # If matched, they will correspond to exact matches (isExact=True) but not to true exact matches (isTrueExact=False)
        ## Sanity check that both lists (impact and exon annotations) are the same length
        newTags = []
        ## Always include 'MUTATION' as a potential variant tag
        newTags.append('MUTATION')

        ## Iterate impact information to generate and include additional variant tags
        for impact in impactAnnots:
            ## Introduce sanity check on variant impact being available
            if not impact:
                continue
            ## Generate potential additional variant tags based on the variant impact
            ## (can be contained in impact eg. 'frameshift_variant&stop_gained')
            if re.search('3_PRIME_UTR_VARIANT',impact):
                newTags.append("3' UTR MUTATION")

            if re.search('5_PRIME_UTR_VARIANT',impact):
                newTags.append("5' UTR MUTATION")

            if re.search('STOP_GAINED',impact):
                newTags.append('TRUNCATING MUTATION')

            if re.search('FRAMESHIFT_VARIANT',impact):
                newTags.append('FRAMESHIFT MUTATION')

        ## When provided, exon info must be matched with impact information, to indicate if variant is exonic or intronic
        if exonAnnots:
            ## Sanity check that both lists (impact and exon annotations) are the same length
            if len(impactAnnots) != len(exonAnnots):
                print("\nError! A variant impact must be provided for each of the provided exon annotations (to guess whether exonic or intronic variant)")
                sys.exit(1)


        ## Iterate exon information to generate and include additional variant tags
        for i,exon in enumerate(exonAnnots):
            ## Introduce sanity check on variant exon being available
            ## eg. with impacts downstream_gene_variant, upstream_gene_variant, intergenic_region
            if not exon:
                continue
            ## Based on rank (exon/intron) information, generate an additional tag
            rank = exon.split('/')[0]
            ## After checking some examples, it seems that intron information is associated to variant impacts 'intron_variant'
            ## and 'sequence_feature' (can be contained in impact eg. 'splice_donor_variant&intron_variant')
            ## NOTE: any variant impact other than the above two will generate a tag 'EXON XX MUTATION',
            ## even 'synonymous' mutations
            if re.search('INTRON_VARIANT',impactAnnots[i]) or re.search('SEQUENCE_FEATURE',impactAnnots[i]):
                newTags.append('INTRON ' + rank + ' MUTATION')
            else:
                newTags.append('EXON ' + rank + ' MUTATION')
                ## As of 13/02/2019, only one case of variant name 'EXON XX FRAMESHIFT' (gene CALR)
                if re.search('FRAMESHIFT_VARIANT',impactAnnots[i]):
                    newTags.append('EXON ' + rank + ' FRAMESHIFT')

        ## Add unique tags and assign corresponding information about type of match
        for tag in newTags:
            if tag not in matchStrings:
                matchStrings.append(tag)
                isExact.append(True)
                isTrueExact.append(False)

    ## Opposite to the SNV case, for CNVs, varAnnotations should correspond to a single element (ie. CNV category)
    if dataType == 'CNV':
        newTags = []
        ## CIVIC seems to consider that GAIN and AMP are the same CNV
        if (varAnnotations == 'AMPLIFICATION') or (varAnnotations == 'AMP') or (varAnnotations == 'GAIN') or (varAnnotations == 'DUPLICATION') or (varAnnotations == 'DUP'):
            newTags.append('AMPLIFICATION')
        ## CIVIC seems to consider that DELETION and LOSS are the same CNV
        ## Both CNV categories 'DELETION' and 'LOSS' occur in CIVIC
        elif (varAnnotations == 'DELETION') or (varAnnotations == 'DEL') or (varAnnotations == 'LOSS'):
            newTags.append('DELETION')
            newTags.append('LOSS')
        newTags.append('COPY NUMBER VARIATION')

        for tag in newTags:
            if tag not in matchStrings:
                matchStrings.append(tag)
                isExact.append(True)
                isTrueExact.append(True)

    return (matchStrings,isExact,isTrueExact)



# FIXME
# Attempt match to a CIViC record using all available variant annotations for the given gene
# A tier and match will always be retrieved (either exact, positional or all gene variants when no match)
# FIXME
# Given a list of variant annotations for a gene in the input table (HGVS but also synonym descriptive terms eg. EXON 15 MUTATION and positional
# strings eg. p.Val600), and all CIVIC variant records for said gene, attempt mapping of any input annotations to one or more CIVIC records.
# Return tier of the match (exact, positional, no match) and corresponding CIVIC variant records that were matched

# Do exhaustive search as there could be >1 perfect match in cases of redundancy (eg. E55FS and E55RFSTER11
# both translate to p.Glu55fs) and there can be >1 positional matches (for SNV case, 'general' variants
# eg. V600 have preference over the rest). Also, for SNVs in case of having two types of exact matches (truly exact
# eg. V600E and descriptive term eg. EXON 15 MUTATION), give preference to the former as descriptive terms should only
# be used when no true exact match was found.

# For CNV, inputStrings will be limited (usually only 1 or 2), and all will always correspond to exact matches (no positional)
def match_variants_in_civic(gene, variants, varMap, dataType, impacts=None, exons=None):
    match = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}

# TODO: sanity check that arguments are of the correct type
    check_is_dict(varMap)

    # Function for generating additional synonym strings for input annotations (eg. EXON 15 MUTATION)
    # isExact is a list of equal length to inputStrings, indicating whether a given string corresponds
    # to an exact (True) or positional (False) match
    # allVariants must refer to the same variant, eg. "c." and "p." annotations of the same variant
    # FIXME: allImpacts and allExons can be empty
    (inputStrings,isExact,isTrueExact) = generate_input_matchStrings(variants, dataType, impacts, exons)

    ## If gene is in CIVIC, possible tier levels are 1,2,3
    if gene in varMap.keys():
        for var_id in varMap[gene].keys():
            check_variant_level_entry(varMap[gene][var_id],"name")
            variant_name = varMap[gene][var_id]["name"]
            check_variant_level_entry(varMap[gene][var_id],"hgvs")
            hgvs_expressions = varMap[gene][var_id]["hgvs"]

            # Generate list of strings that will be used to match input variants to CIVIC database
            # Returned list always has at least length=1 (in this case, containing only variant name)
            # For dataType=CNV, variant matching is not based on HGVS, so matchStrings will only contain the variant name
            civicStrings = generate_civic_matchStrings(variant_name, hgvs_expressions, dataType)

## TODO: in the future, change this if block to situation of 1) having impact and exon 2) or not

            # Iterate input annotations strings and attempt match to any CIVIC variant string
            # The position of each input string corresponds to the type of match (true exact, synonym exact or positional)
            for indx,inputAnnot in enumerate(inputStrings):
                if inputAnnot in civicStrings:
                    # Determine type of match using its position and store accordingly
                    # For CNV, all matches will be exact and true exact matches (either variant is in CIVIC or not)
                    if isExact[indx]:
                        # For SNV: give preference to truly exact matches over descriptive synonym terms
                        # eg. report V600E over EXON 15 MUTATION
                        if isTrueExact[indx]:
                            # For SNV: True exact match (eg. V600E)
                            # For CNV: all exact matches will also be true exact matches (either CNV is in CIVIC or not)
                            # There could be 0,1,>1 exact matches (eg. DELETION + LOSS + COPY NUMBER for cnvs)
                            if var_id not in match["tier_1"]:
                                match["tier_1"].append(var_id)
                        else:
                            # For SNV: descriptive term match (eg. EXON 15 MUTATION, FRAMESHIFT MUTATION)
                            # There could be 0,1,>1 synonym matches
                            if var_id not in match["tier_1b"]:
                                match["tier_1b"].append(var_id)
                    else:
                        # Positional match
                        # There could be 0,1,>1 positional matches
                        if var_id not in match["tier_2"]:
                            match["tier_2"].append(var_id)
                ## For CNV: When there is no exact match for a 'DELETION', also consider special CIVIC CNV records related to exons (these will be positional matches)
                else:
                    ## TODO CNV: Generalize for DELETION and AMPLIFICATION as well
                    if (dataType == 'CNV') and (inputAnnot == 'DELETION'):
                        ## In CNV case, civicStrings corresponds to the CIVIC variant name (single string)
                        for tempString in civicStrings:
                            ## Look for special cases like eg. 'EXON 5 DELETION', 'EXON 1-2 DELETION' or '3' EXON DELETION'
                            isExon = cnv_is_exon_string(tempString)
                            if isExon:
                                ## These special cases are accounted for as positional matches (tier 2)
                                # There could be 0,1,>1 positional matches
                                if var_id not in match["tier_2"]:
                                    match["tier_2"].append(var_id)

        ## Once all CIVIC variants have been iterated, determine final tier and corresponding matched variants
        ## For CNV: either exact or positional matches will occurr due to the implementation design

        # For CNV: positional matches will correspond to EXON records (eg. EXON 1-2 DELETION, 3' EXON DELETION..)
        # For SNV: if there are positional matches, check for preferential positional matches, ie. general variants (like V600)
        if match["tier_2"] and (dataType == 'SNV'):
            for var in match["tier_2"]:
                isGeneral = check_general_variant(var)
                if isGeneral:
                    # Stop as soon as a general variant is found and report only this
                    match["tier_2"] = [var]
                    break

        # If no match was found, then tier case is 3, and return all relevant variants
        if not (match["tier_1"] or match["tier_1b"] or match["tier_2"]):
            ## For SNV: when input variant could not be matched in CIVIC (tier3), return all CIVIC variants associated
            ## to the given gene but that do not correspond to a CNV or EXPRESSION related variant
            if dataType == 'SNV':
                match["tier_3"] = civic_return_all_snvs(varMap[gene])
            ## For CNV: when input cnv could not be matched in CIVIC (tier3), return all CIVIC cnvs associated the given gene (if any)
            ## ie. not all CIVIC records are returned but only those that are matched to a 'CNV' tags
            elif dataType == 'CNV':
                match["tier_3"] = civic_return_all_cnvs(varMap[gene])

    ## If gene is not in CIVIC, tier level is 4 and match will be empty
    else:
        # Sanity check that no other match should have been found
        if (match["tier_1"] or match["tier_1b"] or match["tier_2"] or match["tier_3"]):
            raise ValueError(
                f"Unexpected variant match condition was met.\n"
            )
        # No variants will be reported in this case (as the information is not available)
        match["tier_4"] = True

    return match


def add_civic_match(matchMap,gene,variant,match):
    check_is_dict(matchMap)
    check_is_dict(match)

    # Sanity check that gene and variant are keys of the dict
    if gene not in matchMap.keys():
        # TODO
    if variant not in matchMap[gene].keys():
        # TODO

    # Keep track of all the variant ids matched across tiers for the current gene + variant
    matchedIds = []
    # Add match results to the current gene + variant in matchMap
    for x in match.keys():
        if x == "tier_4":
            matchMap[gene][variant][x] = match[x]
        else:
            for varId in match[x]:
                matchedData[gene][variant][x].append(varId)
                # It is possible that the same CIVIC variant id has been matched for several tiers
                if varId not in matchedIds:
                    matchedIds.append(varId)

    return (matchedData,matchedIds)


# FIXME: currently, only for SNV
def match_in_civic(varData, dataType, identifier_type, select_tier="all", varMap=None):

# TODO: sanity check that arguments are of the correct type
    check_is_dict(varData)

    # Process and sanity check provided select_tier (expected format, valid values, etc.)
    select_tier = check_tier_selection(select_tier)

    matchMap = {}
    matchedIds = []     # keep track of all matched variant ids across genes and tiers

# TODO: add option to provide user-specified varMap (eg. when filters need to be applied) -> add NOTE or warning about anything not being provided will be interpreted as not available in CIVIC

    if varMap is None:
        # gene -> variant -> null
        all_genes = list(varData.keys())
        varMap = query_civic(all_genes, identifier_type)
    # else:
    # TODO: Check correct structure of varMap?

    # gene -> variant -> null
    # where variant -> var="lineNumber|dna|prot|impact|exon"
    for gene in varData.keys():
        if gene not in matchMap.keys():
            matchMap[gene] = {}
        for variant in varData[gene].keys():
            matchMap[gene][variant] = {}
            varArr = variant.split("|")

            variants = []
            impactArr = []
            exonArr = []
            # Format: var="lineNumber|dna|prot|[impact|exon]"
            if dataType == "SNV":
                cVars = varArr[1]
                pVars = varArr[2]
                impacts = varArr[3]
                exons = varArr[4]

# TODO: use a function for this parsing
                cVarArr = cVars.split(",")
                pVarArr = pVars.split(",")
                for cVar in cVarArr:
                    if cVar not in variants:
                        variants.append(cVar)
                for pVar in pVarArr:
                    if pVar not in variants:
                        variants.append(pVar)
                if impacts:
                    impactSplit = impacts.split(",")
                    for impact in impactSplit:
                        impactArr.append(impact)
                if exons:
                    exonSplit = exons.split(",")
                    for exon in exonSplit:
                        exonArr.append(exon)

            # Format: var="lineNumber|cnv"
            if dataType == "CNV":
                cnvVars = varArr[1]
                cnvArr = cnvVars.split(",")
                for cnv in cnvArr:
                    if cnv not in variants:
                        variants.append(cnv)

            # tier -> [matched_vars]
            match = match_variants_in_civic(gene, variants, varMap, dataType, impactAnnots=impactArr, exonAnnots=exonArr)
            # Avoid executing unnecessary filter when select_tier='all' (as resulting dict will be identical)
            if select_tier != "all":
                match = filter_match(match, select_tier)
            # Add the match to the current entry for gene + variant
            # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
            (matchMap,allIds) = add_civic_match(matchMap,gene,variant,match)

            # Keep track of all the variant ids matched across genes and tiers
            for varId in allIds:
                if varId not in matchedIds:
                    matchedIds.append(varId)

    # At this point, all input gene + variants in varData have been matched to CIVIC
    # Filter varMap used to match CIVIC info based on the matched variant ids (to avoid returning unnecessary records)
    varMap = filter_civic(varMap, var_id_in=matchedIds, output_empty=False)

    # Return matchMap, list of all matched variant ids, and associated variant records retrieved from CIVIC (or provided by user), already filtered for the matched variants
    return (matchMap,matchedIds,varMap)


def filter_match(match, select_tier):
    sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]
    # Sanity check that dict contains all expected keys (throws an error if not)
    check_keys(match,sorted_tiers)
    # Process and sanity check provided select_tier (expected format, valid values, etc.)
    select_tier = check_tier_selection(select_tier,sorted_tiers)

    new = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}
    keepTiers = []

    # When select_tier="all", filter is off (keep data for all tiers)
    if isinstance(select_tier, str) and (select_tier == "all"):
        for tmpTier in sorted_tiers:
            keepTiers.append(tmpTier)

    # When select_tier="highest", then keep data only for the highest tier 1>1b>2>3>4
    elif isinstance(select_tier, str) and (select_tier == "highest"):
        for tmpTier in sorted_tiers:
            if tmpTier != "tier_4":
                check_is_list(match[tmpTier])
                if match[tmpTier]:
                    keepTiers.append(tmpTier)
                    break
            else:
                # If we reached to this point, that means all other tiers were empty (tier_4=True)
                keepTiers.append(tmpTier)
                break

    # When select_tier is a list of tiers, then keep data only for those tiers
    elif isinstance(select_tier, list):
        for tmpTier in sorted_tiers:
            if tmpTier in select_tier:
                keepTiers.append(tmpTier)

    # Sanity check for any duplicated tiers
    keepTiers = list(set(keepTiers))

    # Keep variant matches only for the selected tiers (all other data will be excluded from output)
    for tier in sorted_tiers:
        if tier in keepTiers:
            if tier != "tier_4":
                for varId in match[tier]:
                    new[tier].append(varId)

    # Always check if current match corresponds to tier_4 after filtering
    # eg. if only tier_1 is selected and current match did not have any
    if not (new["tier_1"] or new["tier_1b"] or new["tier_2"] or new["tier_3"]):
        new["tier_4"] = True

    return new


def filter_matches(matchMap, select_tier):
    cleanMap = {}
    matchedIds = []

    # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
    # where variant -> var="lineNumber|dna|prot|impact|exon"
    for gene in matchMap.keys():
        if gene not in cleanMap.keys():
            cleanMap[gene] = {}
        for variant in matchMap[gene].keys():
            cleanMap[gene][variant] = {}
            new = filter_match(matchMap[gene][variant], select_tier)
            # Add the filtered match to the current entry for gene + variant
            (cleanMap,allIds) = add_civic_match(cleanMap,gene,variant,new)
            # Keep track of all the variant ids matched across genes and tiers
            for varId in allIds:
                if varId not in matchedIds:
                    matchedIds.append(varId)

    return (cleanMap,matchedIds)


def process_drug_support(matchMap, varMap):
    newMap = {}

    # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
    # where variant -> var="lineNumber|dna|prot|impact|exon"
    for gene in matchMap.keys():
        if gene not in varMap.keys():
            # TODO
        if gene not in newMap.keys():
            newMap[gene] = {}
        for variant in matchMap[gene].keys():
            newMap[gene][variant] = {}
            for tier in matchMap[gene][variant].keys():
                newMap[gene][variant][tier] = {}
                for varId in matchMap[gene][variant][tier]:
                    if varId not in varMap[gene].keys():
                        # TODO
                    # 





    return newMap


# Given a list of disease names, first filter out those matching the black-listed terms (if any), and from the remaining set, return either:
#   - Subset of diseases matching white-listed terms - partial match (eg. 'Melanoma')
#   - If previous is not available, subset of matching high-level diseases - exact match (eg. 'Cancer', 'Solid tumor')
#   - If previous is not available, return all available diseases (ie. so, original set except those matched to the black list, if any)
def classify_diseases(diseaseArr, disease_name_not_in, disease_name_in, alt_disease_names):

    ## Keep track of diseases that did not pass the black-list filter
    blackMatched = []
    ## Keep track of diseases that passed the black-list filter
    cleanSet = []
    ## Final list of matched diseases (associated to one of 3 categories depending on match type: 'ct', 'gt' or 'nct')
    matched = []
    ## Keep track of which type of disease specificy was matched in the end ('ct', 'gt' or 'nct')
    ctArr = []
    gtArr = []
    nctArr = []

    ## 1) First, remove diseases that partially match black-listed terms (if any are provided)
    # NOTE: PARTIAL matches to the black list are allowed! eg:
    #   - including 'small' will remove 'non-small cell lung cancer' and 'lung small cell carcinoma'
    #   - including 'non-small' will remove 'non-small cell lung cancer' but not 'lung small cell carcinoma'
    if disease_name_not_in:
        # Iterate available diseases and keep track of those that partially match to black list (there can be several)
        for disease in diseaseArr:
            # To find partial matches, it is necessary to iterate through both lists (input and civic)
            for blackListed in disease_name_not_in:
                # Search for partial match of INPUT disease in CIVIC disease e.g. 'Melanoma' (input list) in 'Skin Melanoma' (CIVIC) and not opposite
                if blackListed in disease:
                    if disease not in blackMatched:
                        blackMatched.append(disease)
        # Iterate available diseases once again to retrieve those that passed the black list filter
        for disease in diseaseArr:
            # Retrieve valid diseases only
            if disease in blackMatched:
                continue
            if disease not in cleanSet:
                cleanSet.append(disease)

    ## If no black list was provided, then all available diseases constitute the clean set
    else:
        cleanSet = diseaseArr

    ## 2) Now, iterate the list of "allowed" diseases (ie. passing the black list filter) and attempt to match to white-listed terms
    # NOTE: again, PARTIAL matches to the white list are allowed! eg:
    #   - including 'melanoma' will match 'melanoma', 'skin melanoma' and 'uveal melanoma', but not 'skin cancer' (hypothetical)
    #   - including 'uveal melanoma' will only match 'uveal melanoma'
    for cleanType in cleanSet:
        # To find partial matches, it is necessary to iterate through both lists (input and civic)
        for whiteListed in disease_name_in:
            # Search for partial match of INPUT disease in CIVIC disease e.g. 'Melanoma' (input list) in 'Skin Melanoma' (CIVIC) and not opposite
            if whiteListed in cleanType:
                # Keep track of diseases that passed the white list filter
                if cleanType not in matched:
                    matched.append(cleanType)
                    if cleanType not in ctArr:
                        ctArr.append(cleanType)

    ## 3) For those not already matched as ct, attempt a second-best match strategy to higher-level diseases
    ## In CIVIC, some 'general' diseases are included as disease, eg. 'cancer' or 'solid tumor'. Hence, must be exact matches since they are DB-specific
    # NOTE: here, only PERFECT matches are allowed!
    #   - including 'cancer' will only match 'cancer' and not 'lung cancer'
    for cleanType in cleanSet:
        if cleanType in alt_disease_names:
            if cleanType not in matched:
                matched.append(cleanType)
                if cleanType not in gtArr:
                    gtArr.append(cleanType)

    ## 4) For those not already matched as either ct or gt, return all 'allowed' (ie. not black-listed) diseases available in CIVIC
    ## These will be considered as non-specific diseases (ie. off label)
    for cleanType in cleanSet:
        if cleanType not in matched:
            matched.append(cleanType)
            if cleanType not in nctArr:
                nctArr.append(cleanType)

    ## Return diseases that did not match any blacklisted terms, classified and split into: ct, gt, nct
    return (ctArr,gtArr,nctArr)


def add_ct(diseases,ct,gene,variant,evidence_type,newMap,varMap):
    check_is_dict(newMap)
    check_is_dict(varMap)

    # Sanity check that gene and variant are keys of the dict
    if gene not in newMap.keys():
        # TODO
    if variant not in newMap[gene].keys():
        # TODO
    if evidence_type not in newMap[gene][variant]['evidence_items'].keys():
        # TODO
    if gene not in varMap.keys():
        # TODO
    if variant not in varMap[gene].keys():
        # TODO
    if evidence_type not in varMap[gene][variant]['evidence_items'].keys():
        # TODO

    if ct not in newMap[gene][variant]['evidence_items'][evidence_type].keys():
        newMap[gene][variant]['evidence_items'][evidence_type][ct] = {}
    for disease in diseases:
        if disease not in varMap[gene][variant]['evidence_items'][evidence_type].keys():
            # TODO
        newMap[gene][variant]['evidence_items'][evidence_type][ct][disease] = {}
        for drug in varMap[gene][variant]['evidence_items'][evidence_type][disease].keys():
            newMap[gene][variant]['evidence_items'][evidence_type][ct][disease][drug] = {}
            for evidence in varMap[gene][variant]['evidence_items'][evidence_type][disease][drug].keys():
                newMap[gene][variant]['evidence_items'][evidence_type][ct][disease][drug][evidence] = {}
                for evidence_level in varMap[gene][variant]['evidence_items'][evidence_type][disease][drug][evidence].keys():
                    newMap[gene][variant]['evidence_items'][evidence_type][ct][disease][drug][evidence][evidence_level] = []
                    for thisString in varMap[gene][variant]['evidence_items'][evidence_type][disease][drug][evidence][evidence_level]:
                        newMap[gene][variant]['evidence_items'][evidence_type][ct][disease][drug][evidence][evidence_level].append(thisString)

    return (newMap)


def annotate_ct(varMap, disease_name_not_in, disease_name_in, alt_disease_names):
    sorted_cts = ["ct","gt","nct"]
    newMap = {}

    # Iterate the complete varMap dict and reorganize it to classify diseases
    for gene in varMap.keys():
        newMap[gene] = {}
        for variant in varMap[gene].keys():
            newMap[gene][variant] = {}
# FIXME: This would currently fail: adapt to have strict option
            check_keys(varMap[gene][variant], ['name','civic_score','hgvs','types','n_evidence_items','evidence_items'])
            newMap[gene][variant]['name'] = varMap[gene][variant]['name']
            newMap[gene][variant]['civic_score'] = varMap[gene][variant]['civic_score']
            newMap[gene][variant]['hgvs'] = [a for a in varMap[gene][variant]['hgvs']]
            newMap[gene][variant]['types'] = [b for b in varMap[gene][variant]['types']]
            newMap[gene][variant]['n_evidence_items'] = varMap[gene][variant]['n_evidence_items']
            newMap[gene][variant]['evidence_items'] = {}
            for evidence_type in varMap[gene][variant]['evidence_items'].keys():
                newMap[gene][variant]['evidence_items'][evidence_type] = {}
                # Retrieve all disease names associated witht he current evidence type, and classify them into ct, gt, nct
                allDiseases = list(varMap[gene][variant]['evidence_items'][evidence_type].keys())
                (ctDis,gtDis,nctDis) = classify_diseases(allDiseases, disease_name_not_in, disease_name_in, alt_disease_names)
                # Add new layer of classification within the evidence types: specificity of the disease (ct, gt, nct)
                for ct in sorted_cts:
                    if ct == "ct":
                        newMap = add_ct(ctDis,ct,gene,variant,evidence_type,newMap,varMap)
                    if ct == "gt":
                        newMap = add_ct(gtDis,ct,gene,variant,evidence_type,newMap,varMap)
                    if ct == "nct":
                        newMap = add_ct(nctDis,ct,gene,variant,evidence_type,newMap,varMap)

    return newMap


def filter_ct(varMap, select_ct):
    sorted_cts = ["ct","gt","nct"]
    # Process and sanity check provided select_ct (expected format, valid values, etc.)
    select_ct = check_tier_selection(select_ct,sorted_cts)

    newMap = {}

    # When select_ct="all", filter is off (keep data for all 3 cts)
    # Do not copy dict again, simply return input dict untouched
    if isinstance(select_ct, str) and (select_ct == "all"):
        return varMap
    # When select_ct="highest", then keep data only for the highest specificity ct>gt>nct
    elif isinstance(select_ct, str) and (select_ct == "highest"):
        # Iterate the complete varMap dict and reorganize it to classify diseases
        for gene in varMap.keys():
            newMap[gene] = {}
            for variant in varMap[gene].keys():
                newMap[gene][variant] = {}
# FIXME: This would currently fail: adapt to have strict option
                check_keys(varMap[gene][variant], ['name','civic_score','hgvs','types','n_evidence_items','evidence_items'])
                newMap[gene][variant]['name'] = varMap[gene][variant]['name']
                newMap[gene][variant]['civic_score'] = varMap[gene][variant]['civic_score']
                newMap[gene][variant]['hgvs'] = [a for a in varMap[gene][variant]['hgvs']]
                newMap[gene][variant]['types'] = [b for b in varMap[gene][variant]['types']]
                newMap[gene][variant]['n_evidence_items'] = varMap[gene][variant]['n_evidence_items']
                newMap[gene][variant]['evidence_items'] = {}
                for evidence_type in varMap[gene][variant]['evidence_items'].keys():
                    newMap[gene][variant]['evidence_items'][evidence_type] = {}
                    check_keys(varMap[gene][variant]['evidence_items'][evidence_type], ['ct','gt','nct']
                    for ct in sorted_cts:
                        newMap[gene][variant]['evidence_items'][evidence_type][ct] = {}
                        ctArr = list(varMap[gene][variant]['evidence_items'][evidence_type][ct].keys())
                        if ctArr:
                            newMap = add_ct(ctArr,ct,gene,variant,evidence_type,newMap,varMap)
                            break

    # When select_ct is a list of cts, then keep data only for those selected
    elif isinstance(select_ct, list):
        # Iterate the complete varMap dict and reorganize it to classify diseases
        for gene in varMap.keys():
            newMap[gene] = {}
            for variant in varMap[gene].keys():
                newMap[gene][variant] = {}
# FIXME: This would currently fail: adapt to have strict option
                check_keys(varMap[gene][variant], ['name','civic_score','hgvs','types','n_evidence_items','evidence_items'])
                newMap[gene][variant]['name'] = varMap[gene][variant]['name']
                newMap[gene][variant]['civic_score'] = varMap[gene][variant]['civic_score']
                newMap[gene][variant]['hgvs'] = [a for a in varMap[gene][variant]['hgvs']]
                newMap[gene][variant]['types'] = [b for b in varMap[gene][variant]['types']]
                newMap[gene][variant]['n_evidence_items'] = varMap[gene][variant]['n_evidence_items']
                newMap[gene][variant]['evidence_items'] = {}
                for evidence_type in varMap[gene][variant]['evidence_items'].keys():
                    newMap[gene][variant]['evidence_items'][evidence_type] = {}
                    check_keys(varMap[gene][variant]['evidence_items'][evidence_type], ['ct','gt','nct']
                    for ct in sorted_cts:
                        newMap[gene][variant]['evidence_items'][evidence_type][ct] = {}
                        if (ct in select_ct):
                            ctArr = list(varMap[gene][variant]['evidence_items'][evidence_type][ct].keys())
                            newMap = add_ct(ctArr,ct,gene,variant,evidence_type,newMap,varMap)

    return newMap

