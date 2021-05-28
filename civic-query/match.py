import sys
import os
import re

from .utils import check_variant_level_entry,civic_hgvs_to_input,civicName_to_hgvs,extract_p_start,cnv_is_exon_string,civic_return_all_snvs,civic_return_all_cnvs,check_general_variant


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


# FIXME: make impactAnnots and exonAnnots optional!

# Given a set of one or more annotations (SNV or CNV), generate a list of potential strings
# that will be used to match the variant in CIVIC eg. EXON 15 MUTATION, AMPLIFICATION
def generate_input_matchStrings(varAnnotations, dataType, impactAnnots=None, exonAnnots=None):
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
        if len(impactAnnots) != len(exonAnnots):
            print("\nError! Number of variant impacts and exon annotations are not equal")
            sys.exit(1)
        newTags = []
        ## Always include 'MUTATION' as a potential variant tag
        newTags.append('MUTATION')



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
    matchMap = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}
    exactMatches = []
    synMatches = []
    posMatches = []

# TODO: sanity check that arguments are of the correct type
    if not (isinstance(varMap, dict):
        raise TypeError(
            f"'{varMap}' must be of type 'dict'.\n"
        )

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
                            # For CNV: all matches will be exact and true exact matches (either variant is in CIVIC or not)
                            if var_id not in exactMatches:
                                exactMatches.append(var_id)
                        else:
                            # For SNV: descriptive term match (eg. EXON 15 MUTATION, FRAMESHIFT MUTATION)
                            if var_id not in synMatches:
                                synMatches.append(var_id)
                    else:
                        # Positional match
                        if var_id not in posMatches:
                            posMatches.append(var_id)
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
                                if var_id not in posMatches:
                                    posMatches.append(var_id)

        ## Once all CIVIC variants have been iterated, determine final tier and corresponding matched variants
        ## For CNV: either exactMatches or posMatches will occurr due to the implementation design

        # There could be 0,1,>1 perfect matches
        # For CNV: all exact matches will also be true exact matches (either CNV is in CIVIC or not). Multiple matches could occur for a single CNV (eg. DELETION + LOSS + COPY NUMBER)
        matchMap["tier_1"] = exactMatches

        # There could be 0,1,>1 synonym matches
        matchMap["tier_1b"] = synMatches

        # There could be 0,1,>1 positional matches
        # For CNV: positional matches will correspond to EXON records (eg. EXON 1-2 DELETION, 3' EXON DELETION..)
        matchMap["tier_2"] = posMatches
        # For SNV: if there are positional matches, check for preferential positional matches, ie. general variants (like V600)
        if posMatches and (dataType == 'SNV'):
            for var in posMatches:
                isGeneral = check_general_variant(var)
                if isGeneral:
                    # Stop as soon as a general variant is found and report only this
                    matchMap["tier_2"] = [var]
                    break

        # If no match was found, then tier case is 3, and return all relevant variants
        if not (exactMatches or synMatches or posMatches):
            ## For SNV: when input variant could not be matched in CIVIC (tier3), return all CIVIC variants associated
            ## to the given gene but that do not correspond to a CNV or EXPRESSION related variant
            if dataType == 'SNV':
                matchMap["tier_3"] = civic_return_all_snvs(varMap[gene])
            ## For CNV: when input cnv could not be matched in CIVIC (tier3), return all CIVIC cnvs associated the given gene (if any)
            ## ie. not all CIVIC records are returned but only those that are matched to a 'CNV' tags
            elif dataType == 'CNV':
                matchMap["tier_3"] = civic_return_all_cnvs(varMap[gene])

    ## If gene is not in CIVIC, tier level is 4 and matchMap will be empty
    else:
        # No variants will be reported in this case (as the information is not available)
        matchMap["tier_4"] = True

    return matchMap


# Given a list of disease names, first filter out those matching the black-listed terms (if any), and from the remaining set, return either:
#   - Subset of diseases matching white-listed terms - partial match (eg. 'Melanoma')
#   - If previous is not available, subset of matching high-level diseases - exact match (eg. 'Cancer', 'Solid tumor')
#   - If previous is not available, return all available diseases (ie. so, original set except those matched to the black list, if any)
def match_and_classify_diseases(diseaseArr, disease_name_not_in, disease_name_in, alt_disease_names):

    ## Keep track of diseases that did not pass the black-list filter
    blackMatched = []
    ## Keep track of diseases that passed the black-list filter
    cleanSet = []
    ## Final list of matched diseases (associated to one of 3 categories depending on match type: 'ct', 'gt' or 'nct')
    matched = []
    ## Keep track of which type of disease specificy was matched in the end ('ct', 'gt' or 'nct')
    ct = ''

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
                    ct = 'ct'
                    matched.append(cleanType)

    ## 3) When nothing could be matched to the white-list terms, attempt a second-best match strategy to higher-level diseases
    ## In CIVIC, some 'general' diseases are included as disease, eg. 'cancer' or 'solid tumor'. Hence, must be exact matches since they are DB-specific
    # NOTE: here, only PERFECT matches are allowed!
    #   - including 'cancer' will only match 'cancer' and not 'lung cancer'
    if not matched:
        for cleanType in cleanSet:
            if cleanType in alt_disease_names:
                if cleanType not in matched:
                    ct = 'gt'
                    matched.append(cleanType)

    ## 4) If nothing could be matched to either the white-list or higher-level diseases, return all 'allowed' (ie. not black-listed) diseases available in CIVIC
    ## These will be considered as non-specific diseases (ie. off label)
    if not matched:
        for cleanType in cleanSet:
            if cleanType not in matched:
                ct = 'nct'
                matched.append(cleanType)

    ## Now, list 'matched' contains all cancer types for which their evidences items will be returned
    ## They can correspond to either 'ct' (from white-list), 'gt' (from high-level list) or 'nct' (if nothing could be matched, return all that is available)
    return (matched,ct)
