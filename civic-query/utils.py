import sys
import os
import re

# import read_and_write
# dict_coodes = get_dict_aminoacids()

def check_is_list(myList):
    # Check object is a list (even if empty)
    if not isinstance(myList, list):
        raise TypeError(
            f"'{myList}' is not of type 'list'.\n"
        )


def check_string_filter_arguments(field, myList):
    # Check expected class was provided for each argument
    # myList should be a list (even if empty)
    check_is_list(myList)
    # field should be a non-empty string
    if not isinstance(field, str):
        raise TypeError(
            f"'{field}' is not of type 'str'.\n"
        )
    if not field:
        raise TypeError(
            f"'{field}' cannot be an empty string.\n"
        )
    # Use uppercase and remove leading/trailing spaces, for consistency of strings
    newField = field.strip().upper()
    newList = []
    # Convert individual elements into strings, in case provided list contains numbers (eg. gene or variant ids)
    for tmpItem in myList:
        newItem = str(tmpItem)
        newItem = newItem.strip().upper()
        newList.append(newItem)
    return (newField,newList)


def apply_in_filter(field, inList, matchType="exact"):
    (field, inList) = check_string_filter_arguments(field, inList)
    match = False
    # Only perform check when inList is not empty
    if inList:
        # If a list was provided and field is empty, then filter is not passed
        if field != "NULL":
            if (matchType == "exact"):
                if field in inList:
                    match = True
            if (matchType == "partial"):
                for tmp in inList:
                    if tmp in field:
                        match = True
    else:
        match = True
    return match


def apply_not_in_filter(field, outList, matchType="exact"):
    (field, outList) = check_string_filter_arguments(field, outList)
    match = False
    # If no list was provided, even if field is empty, filter is always "passed" (ie. False to not remove variant)
    # This allows to filter record if a particular field is "NULL"
    if outList:
        if (matchType == "exact"):
            if field in outList:
                match = True
        if (matchType == "partial"):
            for tmp in outList:
                if tmp in field:
                    match = True
    return match


# Only keep if cutoff or more instances
# Ignore filter when cutoff=0
def apply_cutoff_filter(field, cutoff):

## TODO check for correct classes provided
# field should be numeric or float
# cutoff should be numeric or float

    match = True
    cutoff_f = float(cutoff)
    field_f = float(field)
    if (cutoff_f != float(0)):
        if (field_f < cutoff_f):
            match = False
    return match


def filter_civic_results(varMap, gene_id_in=[], gene_id_not_in=[], min_variants=0, var_id_in=[], var_id_not_in=[], var_name_in=[], var_name_not_in=[], min_civic_score=0, var_type_in=[], var_type_not_in=[], min_evidence_items=0, evidence_type_in=[], evidence_type_not_in=[], disease_in=[], disease_not_in=[], drug_name_in=[], drug_name_not_in=[], evidence_dir_in=[], evidence_dir_not_in=[], evidence_clinsig_in=[], evidence_clinsig_not_in=[], evidence_level_in=[], evidence_level_not_in=[], evidence_status_in=[], evidence_status_not_in=[], source_status_in=[], source_status_not_in=[], var_origin_in=[], var_origin_not_in=[], source_type_in=[], source_type_not_in=[], min_evidence_rating=0, output_empty=False):

## filters are applied in the same order as their corresponding arguments (as they are defined in the function, not the order specified in the function call)
## so, logic is always AND for all selected filters?

## if desired filter logic is not possible, then the function would need to be run several times, applyin the filters subsequently

## all array or cutoff filters work identically, except for disease-related arguments (ie. disease_in=[], disease_not_in=[], alt_diseases=[])

## Assume varMap has the following structure:
# TODO

## TODO: talk about arg output_empty -> if True, then empty entries can be returned (it will be empty at the level of the filter that failed)

## TODO: add framework for flexibly choosing what is filtered based of constructed argument?
## eg. feature_name (have a catalogue of available and throw error if not matched) + in/not_in + partial/exact (feature type would determine if expectation is list or cutoff value)

## TODO: add framework to apply filters based on order of arguments?
## eg. for each feature that is provided, filtering fuction is called, and the dict is subsequently filtered in this manner until all provided filters have been applied

    # sanity check that provided output_empty is logical
    if not (isinstance(output_empty, int) or isinstance(output_empty, float)):
        raise TypeError(
            f"'{output_empty}' must be of type 'int' or 'float'.\n"
        )

    # Iterate complete dict of results and apply selected filters to generate a new (filtered) dict
    cleanMap = {}
    for gene_id in varMap.keys():
        # gene id should always be available  (type of id is chosen during the query)
        gene_id_str = str(gene_id)
        keepGene = apply_in_filter(gene_id_str, gene_id_in, matchType="exact")
        if not keepGene:
            continue
        removeGene = apply_not_in_filter(gene_id_str, gene_id_not_in, matchType="exact")
        if removeGene:
            continue

        # allow number of min_variants to be 0
        n_variants = len(varMap[gene_id].keys())
        keepGene = apply_cutoff_filter(n_variants, min_variants)
        if not keepGene:
            continue

        # current gene has passed all gene-level filters
        # write to output (filtered) dict only when output_empty=True
        if output_empty:
            if gene_id not in cleanMap.keys():
                cleanMap[gene_id] = {}

        # variant id should always be available
        for var_id in varMap[gene_id].keys():
            var_id_str = str(var_id)
            keepVar = apply_in_filter(var_id_str, var_id_in, matchType="exact")
            if not keepVar:
                continue
            removeVar = apply_not_in_filter(var_id_str, var_id_not_in, matchType="exact")
            if removeVar:
                continue

            # variant name should always be available (never "NULL")
            variant = varMap[gene_id][var_id]["name"]
            keepVar = apply_in_filter(variant, var_name_in, matchType="partial")
            if not keepVar:
                continue
            removeVar = apply_not_in_filter(variant, var_name_not_in, matchType="partial")
            if removeVar:
                continue

            # civic score is always a number (can be 0)
            var_score = varMap[gene_id][var_id]["civic_score"]
            keepVar = apply_cutoff_filter(var_score, min_civic_score)
            if not keepVar:
                continue

            # variant types will never be empty list
            # use ["NULL"] when not available
            variant_types = varMap[gene_id][var_id]["types"]
            nKeep = 0
            nRemove = 0
            # when "NULL":
            #  - if any var_type_in was provided, then variant will fail this filter
            #  - if any var_type_not_in was provided, then variant will not fail this filter (given that the previous filter did not either)
            for var_type in variant_types:
                keepType = apply_in_filter(var_type, var_type_in, matchType="partial")
                if keepType:
                    nKeep += 1
                removeType = apply_not_in_filter(var_type, var_type_not_in, matchType="partial")
                if removeType:
                    nRemove += 1
            if (nKeep == 0):
                continue
            if (nRemove > 0):
                continue

            # allow number of evidence items to be 0
            n_evidence_items = varMap[gene_id][var_id]["n_evidence_items"]
            keepVar = apply_cutoff_filter(n_evidence_items, min_evidence_items)
            if not keepVar:
                continue

            # current variant has passed all variant-level filters
            # write to output (filtered) dict only when output_empty=True
            if output_empty:
                if var_id not in cleanMap[gene_id].keys():
                    cleanMap[gene_id][var_id] = {}
                    cleanMap[gene_id][var_id]["name"] = variant
                    cleanMap[gene_id][var_id]["civic_score"] = var_score
                    cleanMap[gene_id][var_id]["hgvs"] = varMap[gene_id][var_id]['hgvs']
                    cleanMap[gene_id][var_id]["types"] = variant_types
                    cleanMap[gene_id][var_id]["n_evidence_items"] = 0
                    cleanMap[gene_id][var_id]["evidence_items"] = {}

            n_evidence_items_after = 0
            evidence_types = list(varMap[gene_id][var_id]["evidence_items"].keys())
            allowed_evidence_types = []
            # evidence type should always be available
            for evidence_type in evidence_types:
                keepType = apply_in_filter(evidence_type, evidence_type_in, matchType="exact")
                if not keepType:
                    continue
                removeType = apply_not_in_filter(evidence_type, evidence_type_not_in, matchType="exact")
                if removeType:
                    continue
                allowed_evidence_types.append(evidence_type)

            # only work on subset of evidence types that passed filters (if any)
            for evidence_type in allowed_evidence_types:
                # current evidence type has passed filters
                # write to output (filtered) dict only when output_empty=True
                if output_empty:
                    if evidence_type not in cleanMap[gene_id][var_id]["evidence_items"].keys():
                        cleanMap[gene_id][var_id]["evidence_items"][evidence_type] = {}

                type_diseases = list(varMap[gene_id][var_id]["evidence_items"][evidence_type].keys())
                # disease names should always be available
                allowed_diseases = []
                for type_disease in type_diseases:
                    keepDisease = apply_in_filter(type_disease, disease_in, matchType="partial")
                    if not keepDisease:
                        continue
                    removeDisease = apply_not_in_filter(type_disease, disease_not_in, matchType="partial")
                    if removeDisease:
                        continue
                    allowed_diseases.append(type_disease)

                # only work on subset of diseases that passed filters (if any)
                for disease in allowed_diseases:
                    # current disease name has passed filters
                    # write to output (filtered) dict only when output_empty=True
                    if output_empty:
                        if disease not in cleanMap[gene_id][var_id]["evidence_items"][evidence_type].keys():
                            cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease] = {}

                    # all evidence types except 'PREDICTIVE' will have a single entry 'NULL' (because drugs are not associated to these evidence types)
                    # when "NULL", the drug-related filters do not apply (to avoid filtering evidences from types other than 'PREDICTIVE')

                    disease_drugs = list(varMap[gene_id][var_id]["evidence_items"][evidence_type][disease].keys())

                    allowed_drugs = []
                    if (evidence_type != "PREDICTIVE"):
                        allowed_drugs = disease_drugs   # ['NULL']
                    else:
                        # iterate existing drugs and apply filters
                        for disease_drug in disease_drugs:
                            keepDrug = apply_in_filter(disease_drug, drug_name_in, matchType="partial")
                            if not keepDrug:
                                continue
                            removeDrug = apply_not_in_filter(disease_drug, drug_name_not_in, matchType="partial")
                            if removeDrug:
                                continue
                            allowed_drugs.append(disease_drug)

                    # only work on subset of drugs that passed filters (if any)
                    for drug in allowed_drugs:
                        # current drug has passed filters ("NULL" whenever evidence type is not 'PREDICTIVE')
                        # write to output (filtered) dict only when output_empty=True
                        if output_empty:
                            if drug not in cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease].keys():
                                cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug] = {}

                        evidences = list(varMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug].keys())
                        allowed_evidences = []
                        for evidence in evidences:
                            # split the evidence string into direction and clinical significance
                            # format: 'DIRECTION:CLINSIGF'
                            evidenceArr = evidence.strip().split(':')
                            # sanity check for correct format of evidence string
                            if (len(evidenceArr) != 2):
                                raise ValueError(
# TODO
                                    f"'{disease_drugs}' .\n"
                                )

                            # check evidence direction filters
                            direction = evidenceArr[0]
                            keepDir = apply_in_filter(direction, evidence_dir_in, matchType="exact")
                            if not keepDir:
                                continue
                            removeDir = apply_not_in_filter(direction, evidence_dir_not_in, matchType="exact")
                            if removeDir:
                                continue
                            # check evidence clinical significance filters
                            clin_signf = evidenceArr[1]
                            keepClin = apply_in_filter(clin_signf, evidence_clinsig_in, matchType="exact")
                            if not keepClin:
                                continue
                            removeClin = apply_not_in_filter(clin_signf, evidence_clinsig_not_in, matchType="exact")
                            if removeClin:
                                continue
                            allowed_evidences.append(evidence)

                        # only work on subset of evidences that passed filters (if any)
                        for this_evidence in allowed_evidences:
                            # current evidence has passed filters
                            # write to output (filtered) dict only when output_empty=True
                            if output_empty:
                                if this_evidence not in cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug].keys():
                                    cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence] = {}

                            evidence_levels = list(varMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence].keys())
                            allowed_levels = []
                            for evidence_level in evidence_levels:
                                keepLevel = apply_in_filter(evidence_level, evidence_level_in, matchType="exact")
                                if not keepLevel:
                                    continue
                                removeLevel = apply_not_in_filter(evidence_level, evidence_level_not_in, matchType="exact")
                                if removeLevel:
                                    continue
                                allowed_levels.append(evidence_level)

                            # only work on subset of evidence levels that passed filters (if any)
                            for level in allowed_levels:
                                # current evidence level has passed filters
                                # write to output (filtered) dict only when output_empty=True
                                if output_empty:
                                    if level not in cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence].keys():
                                        cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence][level] = []

                                # Format of list: ['TYPE_ID1:EVIDENCESTATUS1:SOURCESTATUS1:VARORIGIN1:RATING1',..,'TYPE_IDN:EVIDENCESTATUSN:SOURCESTATUSN:VARORIGINN:RATINGN']
                                all_evidence_items = varMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence][level]
                                for evidence_item in all_evidence_items:
                                    # split the evidence item string into the 5 separate fields
                                    itemArr = evidence_item.strip().split(":")
                                    # sanity check for correct format of evidence item string
                                    if (len(itemArr) != 5):
                                        raise ValueError(
# TODO
                                            f"'{itemArr}' .\n"
                                        )
                                    tmpId = itemArr[0]
                                    evidence_status = itemArr[1]
                                    source_status = itemArr[2]
                                    var_origin = itemArr[3]
                                    rating = itemArr[4]

                                    # split the id string into the 2 separate fields
                                    # format: 'TYPE_ID'
                                    tmpIdArr = tmpId.split("_")
                                    # sanity check for correct format of evidence id string
                                    if (len(tmpIdArr) != 2):
                                        raise ValueError(
# TODO
                                            f"'{tmpIdArr}' .\n"
                                        )
                                    idType = tmpIdArr[0]
                                    this_id = tmpIdArr[1]

                                    # apply filters on evidence status
                                    keepStatus = apply_in_filter(evidence_status, evidence_status_in, matchType="exact")
                                    if not keepStatus:
                                        continue
                                    removeStatus = apply_not_in_filter(evidence_status, evidence_status_not_in, matchType="exact")
                                    if removeStatus:
                                        continue
                                    # apply filters on source status
                                    keepSource = apply_in_filter(source_status, source_status_in, matchType="partial")
                                    if not keepSource:
                                        continue
                                    removeSource = apply_not_in_filter(source_status, source_status_not_in, matchType="partial")
                                    if removeSource:
                                        continue
                                    # apply filters on variant origin
                                    keepOrigin = apply_in_filter(var_origin, var_origin_in, matchType="partial")
                                    if not keepOrigin:
                                        continue
                                    removeOrigin = apply_not_in_filter(var_origin, var_origin_not_in, matchType="partial")
                                    if removeOrigin:
                                        continue
                                    # apply filters on source type
                                    keepType = apply_in_filter(idType, source_type_in, matchType="exact")
                                    if not keepType:
                                        continue
                                    removeType = apply_not_in_filter(idType, source_type_not_in, matchType="exact")
                                    if removeType:
                                        continue
                                    # when rating in not available (ie. 'NULL'), the evidence will directly fail if corresponding filter is set
                                    if (rating == "NULL"):
                                        if (float(min_evidence_rating) != float(0)):
                                            continue
                                    # when rating is available, apply filters on evidence rating
                                    else:
                                        keepRating = apply_cutoff_filter(rating, min_evidence_rating)
                                        if not keepRating:
                                            continue

                                    # current evidence item has passed all filters
                                    # write complete entry (from gene to evidence item) to output (filtered) dict only when output_empty>0
                                    if not output_empty:
                                        if gene_id not in cleanMap.keys():
                                            cleanMap[gene_id] = {}
                                        if var_id not in cleanMap[gene_id].keys():
                                            cleanMap[gene_id][var_id] = {}
                                            cleanMap[gene_id][var_id]["name"] = variant
                                            cleanMap[gene_id][var_id]["civic_score"] = var_score
                                            cleanMap[gene_id][var_id]["hgvs"] = varMap[gene_id][var_id]['hgvs']
                                            cleanMap[gene_id][var_id]["types"] = variant_types
                                            cleanMap[gene_id][var_id]["n_evidence_items"] = 0
                                            cleanMap[gene_id][var_id]["evidence_items"] = {}
                                        if evidence_type not in cleanMap[gene_id][var_id]["evidence_items"].keys():
                                            cleanMap[gene_id][var_id]["evidence_items"][evidence_type] = {}
                                        if disease not in cleanMap[gene_id][var_id]["evidence_items"][evidence_type].keys():
                                            cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease] = {}
                                        if drug not in cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease].keys():
                                            cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug] = {}
                                        if this_evidence not in cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug].keys():
                                            cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence] = {}
                                        if level not in cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence].keys():
                                            cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence][level] = []

                                    cleanMap[gene_id][var_id]["evidence_items"][evidence_type][disease][drug][this_evidence][level].append(evidence_item)
                                    n_evidence_items_after += 1

            # At this point, all evidence items available for this variant have been parsed and filtered
            # Add number of evidence items remaining after filtering for this variant
            # check if output_empty=True or not
            if not output_empty:
                # In this case, corresponding gene and variant entries will only be available when n_items > 0 (avoid empty entries)
                if (n_evidence_items_after > 0):
                    cleanMap[gene_id][var_id]["n_evidence_items"] = n_evidence_items_after
            else:
                # In this case, corresponding gene and variant entries are always available
                cleanMap[gene_id][var_id]["n_evidence_items"] = n_evidence_items_after

    return cleanMap


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


# Translate a 1-letter aminoacid code (if it exists) into a 3-letter code
def translate_aa(oneLetter):
    aaNew = None
    if oneLetter.upper() in dict_codes,dict_codes:
        aaNew = dict_codes[oneLetter.upper()]
    return aaNew

# Given a single CIVIC variant name, extract potential HGVS annotations by
# parsing and modifying this string using knowledge on CIVIC naming conventions
# Function only applicable to SNV data
def civicName_to_hgvs(varName):
    # TODO: variant names like T193C (refer to DNA) would be translated to proteins
    hgvsStrings = []
    # 1) HGVS 1-letter protein code, including stop codons and general variants (eg. V600E, *600E, V600* or V600)
    # TODO: multiple cases, what to do? pos = re.findall
    if re.match(r'([A-Z*])(\d+)([A-Z*]?)($|\s\()', varName):
        pos = re.match(r'([A-Z*])(\d+)([A-Z*]?)($|\s\()', varName).groups()
        aa1 = pos[0]
        npos = pos[1]
        aa2 = pos[2]

        ## Generate the HGVS string by translating 1-letter aa codes into 3-letter

        # Special case of protein extensions where '*' is represented by 'Ter' in input table
        if aa1 == '*':
            aa1New = 'TER'
        else:
            # Translate 1-letter aa code to 3-letter code
            aa1New = translate_aa(aa1)
        # Check for general variants (second aa will be '')
        if aa2:
            # Special cases like p.Ter600Ter happen in input table, so in theory,
            # *600* should also be possible (however does not seem to happen in CIVIC)
            if (aa1New == 'TER') and (aa2=='*'):
                aa2New = 'TER'
            else:
                # Translate 1-letter aa code to 3-letter code
                # In cases where aa2='*', it will remain as '*' (using Ter is a special case when aa1=Ter)
                aa2New = translate_aa(aa2)
        else:
            aa2New = ''
        # Sanity check that translated strings correspond to valid aa codes
        if (aa1New is not None) and (aa2New is not None):
            # Construct new protein annotation
            newAnnotation = 'P.' + aa1New + npos + aa2New
            hgvsStrings.append(newAnnotation)

    # 2) Extract embedded c. annotation if available ie. '(c.XXX)'
    # Use re.search to allow for matches not in the string start
    if re.search(r'\((C\..+?)\)', varName):
        # TODO: at least two multiple cases, what to do?
        # 'A50A (c.150C>G); Splicing alteration (c.463-1G>T)'
        # r = re.findall(r'\((c\..+?)\)', varName)
        annot = re.search(r'\((C\..+?)\)', varName).groups()[0]
        hgvsStrings.append(annot)

    # 3) Frameshifts (eg. T157FS or T157MFS)
    if re.match(r'([A-Z])(\d+)([A-Z]?)FS', varName):
        pos = re.match(r'([A-Z])(\d+)([A-Z]?)FS', varName).groups()
        aa = pos[0]
        npos = pos[1]
        # Translate aa in 1st position and if ok, write frameshift in short form
        aaNew = translate_aa(aa)
        if aaNew is not None:
            annot = 'P.' + aaNew + npos + 'FS'
            hgvsStrings.append(annot)

        # TODO: extend cases?

    return hgvsStrings


# Given a single CIVIC hgvs expression, parse and modify it to ensure
# it complies with the HGVS format followed by the input table
# Function only applicable to SNV data
# TODO: Currently, only modification of p. annotations are supported
def civic_hgvs_to_input(civic_hgvs):
    ## The following "special hgvs" cases should be mutually exclusive
    ## ie. a single HGVS expression should match only 1 of the cases
    newAnnot = None

    ## Frameshift
    # Input table seems to use short form for frameshift annotations (ie. p.Glu55fs) while CIVIC does not
    if re.match(r'(P\.[A-Z]+[0-9]+)[A-Z]+FS.*', civic_hgvs):
        newAnnot = re.sub(r'(P\.[A-Z]+[0-9]+)[A-Z]+FS.*', r'\1FS', civic_hgvs)
        # Only return new HGVS expression if something changed
        if newAnnot != civic_hgvs:
            return newAnnot

    ## Nonsense mutations (gain of stop codon)
    # In general, CIViC seems to use 'Ter' for stop codons, only 1 exception so far (p.F76Lfs*56)
    # Input table uses '*' for stop codons; exception in protein extensions (eg. p.Ter370Tyrext*?)
    # '*' is also used for nucleotide numbering in c.,n., but not relevant here (only p.)
    if re.match(r'(P\.[A-Z]+[0-9]+)TER', civic_hgvs):
        newAnnot = re.sub(r'(P\.[A-Z]+[0-9]+)TER', r'\1*', civic_hgvs)
        # Only return new HGVS expression if something changed
        if newAnnot != civic_hgvs:
            return newAnnot

    ## Loss of stop codon
    # In CIVIC, a stop codon loss is expressed as p.Ter214Cys
    # In HGVS, this is enconded as an extension (eg. p.Ter370Tyrext*?)
    # In this case, we need to change annotation from input table (generate_civic_matchStrings)

    ## Silent mutations (no aa change)
    # CIVIC uses a '=' in the 2nd position to represent no change (eg. p.Pro61=)
    # In the input table, the 2nd aa is specified (eg. p.Pro61Pro)
    if re.match(r'P\.([A-Z]+)([0-9]+)=', civic_hgvs):
        newAnnot = re.sub(r'P\.([A-Z]+)([0-9]+)=', r'P.\1\2\1', civic_hgvs)
        # Only return new HGVS expression if something changed
        if newAnnot != civic_hgvs:
            return newAnnot

    return newAnnot


# Given a single HGVS p. annotation (of the form p.Pro61Cys...), extract the start of the string (ie. p.Pro61)
# Function only applicable to SNV data
def extract_p_start(pAnnotation):
    startAnnot = None
    if re.match(r'(P\.[A-Z]+[0-9]+)', pAnnotation):
        startAnnot = re.match(r'(P\.[A-Z]+[0-9]+)', pAnnotation).groups()[0]

    return startAnnot


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


# Given a single CIVIC variant name, return whether it corresponds to a CNV variant record related to exons
# For this, attemp to match the variant name to special CNV exon cases present in CIVIC (eg. EXON 1-2 DELETION, EXON 5 DELETION, 3' EXON DELETION...)
def cnv_is_exon_string(varName):
    exonStrings = ['^EXON [0-9-]+ DELETION$', '^[35\']+ EXON DELETION$', '^EXON [0-9-]+ SKIPPING MUTATION$']
    isExon = False
    for exonString in exonStrings:
        if re.search(exonString, varName):
            isExon = True

    return isExon

# Given a single CIVIC record name, return whether it corresponds to a EXPRESSION record related to exons
# For this, attemp to match the variant name to special EXPRESSION exon cases present in CIVIC (eg. EXON 1-2 EXPRESSION, EXON 5 OVEREXPRESSION...)
def expr_is_exon_string(varName):
    exonStrings = ['^EXON [0-9-]+ EXPRESSION$', '^EXON [0-9-]+ OVEREXPRESSION$', '^EXON [0-9-]+ UNDEREXPRESSION$']
    isExon = False
    for exonString in exonStrings:
        if re.search(exonString, varName):
            isExon = True

    return isExon


# Given a set of one or more annotations (SNV or CNV), generate a list of potential strings
# that will be used to match the variant in CIVIC eg. EXON 15 MUTATION, AMPLIFICATION
def generate_input_matchStrings(varAnnotations, dataType, impactAnnots, exonAnnots):
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

        ## Iterate exon and impact information to generate and include additional variant tags
        for i,impact in enumerate(impactAnnots):
            ## Introduce sanity check on variant impact being available
            if not impact:
                continue
            ## First, generate potential additional variant tags based on the variant impact
            ## (can be contained in impact eg. 'frameshift_variant&stop_gained')
            if re.search('3_PRIME_UTR_VARIANT',impact):
                newTags.append("3' UTR MUTATION")

            if re.search('5_PRIME_UTR_VARIANT',impact):
                newTags.append("5' UTR MUTATION")

            if re.search('STOP_GAINED',impact):
                newTags.append('TRUNCATING MUTATION')

            if re.search('FRAMESHIFT_VARIANT',impact):
                newTags.append('FRAMESHIFT MUTATION')

            ## Second, get rank (exon/intron) information to generate an additional tag
            rank = exonAnnots[i]
            ## Skip cases where exon information is not available
            ## eg. with impacts downstream_gene_variant, upstream_gene_variant, intergenic_region
            if not rank:
                continue
            ## Get the actual rank (eg. 2/30)
            rank = rank.split('/')[0]
            ## After checking some examples, it seems that intron information is associated to variant impacts 'intron_variant'
            ## and 'sequence_feature' (can be contained in impact eg. 'splice_donor_variant&intron_variant')
            ## NOTE: any variant impact other than the above two will generate a tag 'EXON XX MUTATION',
            ## even 'synonymous' mutations
            if re.search('INTRON_VARIANT',impact) or re.search('SEQUENCE_FEATURE',impact):
                newTags.append('INTRON ' + rank + ' MUTATION')
            else:
                newTags.append('EXON ' + rank + ' MUTATION')
                ## As of 13/02/2019, only one case of variant name 'EXON XX FRAMESHIFT' (gene CALR)
                if re.search('FRAMESHIFT_VARIANT',impact):
                    newTags.append('EXON ' + rank + ' FRAMESHIFT')

        ## Add unique tags and assign corresponding information about type of match
        for tag in newTags:
            if tag not in matchStrings:
                matchStrings.append(tag)
                isExact.append(True)
                isTrueExact.append(False)

    ## For data type CNV, follow a different procedure
    else:
        newTags = []

        ## b) Given a single CNV annotation (AMP, GAIN, DEL, LOSS), generate a list of potential strings that will be used to match CNV in CIVIC
        if dataType == 'CNV':
            ## Opposite to the SNV case, for CNVs, varAnnotations should correspond to a single element (ie. CNV category)
            ## CIVIC seems to consider that GAIN and AMP are the same CNV
            if (varAnnotations == 'AMP') or (varAnnotations == 'GAIN'):
                newTags.append('AMPLIFICATION')
            ## CIVIC seems to consider that DELETION and LOSS are the same CNV
            ## Both CNV categories 'DELETION' and 'LOSS' occur in CIVIC
            elif (varAnnotations == 'DEL') or (varAnnotations == 'LOSS'):
                newTags.append('DELETION')
                newTags.append('LOSS')
            newTags.append('COPY NUMBER VARIATION')

        for tag in newTags:
            if tag not in matchStrings:
                matchStrings.append(tag)
                isExact.append(True)
                isTrueExact.append(True)

    return matchStrings,isExact,isTrueExact


# Given all CIVIC variant records associated to a given gene, return all those variant names that should correspond to SNV records
# For this, remove the most common CNV and EXPRESSION names existing in CIVIC from the complete set of the gene's variant names
# Additionally, consider other special cases present in CIVIC (eg. EXON 1-2 DELETION/EXPRESSION, EXON 5 DELETION/OVEREXPRESSION, 3' EXON DELETION/UNDEREXPRESSION...)
def civic_return_all_snvs(geneData):
    # Common CNV record names in CIVIC
    cnvNames = ['AMPLIFICATION','DELETION','LOSS','COPY NUMBER VARIATION']
    # Common EXPRESSION record names in CIVIC
    exprNames = ['OVEREXPRESSION','UNDEREXPRESSION','EXPRESSION']

    # All variant names not matching a CNV or EXPRESSION will be returned
    matches = []
    for varName in list(geneData.keys()):
        ## Also attempt matching of variant name to common CNV and EXPRESSION names related to exons
        isExon_cnv = cnv_is_exon_string(varName)
        isExon_expr = expr_is_exon_string(varName)
        ## Skip variant names matching a CNV or EXPRESSION
        if (varName in cnvNames) or isExon_cnv:
            continue
        if (varName in exprNames) or isExon_expr:
            continue
        if varName not in matches:
            matches.append(varName)

    return matches


# Given all CIVIC variant records associated to a given gene, return all those variant names that correspond to CNV records
# For this, attemp to match the gene's variant names to the most common CNV names existing in CIVIC
# Additionally, consider other special CNV cases present in CIVIC (eg. EXON 1-2 DELETION, EXON 5 DELETION, 3' EXON DELETION...)
def civic_return_all_cnvs(geneData):
    # Common CNV record names in CIVIC
    cnvNames = ['AMPLIFICATION','DELETION','LOSS','COPY NUMBER VARIATION']

    # All matched variant names will be returned
    matches = []
    for varName in list(geneData.keys()):
        ## Also attempt matching of variant name to common CNV names related to exons
        isExon = cnv_is_exon_string(varName)
        if (varName in cnvNames) or isExon:
#             ## Special case for 'LOSS': only considered synonym of 'DELETION' when a certain variant type is present
#             if varName == 'LOSS' and ('TRANSCRIPT_ABLATION' not in geneData[varName]['types']):
#                 continue
            if varName not in matches:
                matches.append(varName)

    return matches


# Given a list of variant annotations for a gene in the input table (HGVS but also synonym descriptive terms eg. EXON 15 MUTATION and positional
# strings eg. p.Val600), and all CIVIC variant records for said gene, attempt mapping of any input annotations to one or more CIVIC records.
# Return tier of the match (exact, positional, no match) and corresponding CIVIC variant records that were matched

# Do exhaustive search as there could be >1 perfect match in cases of redundancy (eg. E55FS and E55RFSTER11
# both translate to p.Glu55fs) and there can be >1 positional matches (for SNV case, 'general' variants
# eg. V600 have preference over the rest). Also, for SNVs in case of having two types of exact matches (truly exact
# eg. V600E and descriptive term eg. EXON 15 MUTATION), give preference to the former as descriptive terms should only
# be used when no true exact match was found.

# For CNV, inputStrings will be limited (usually only 1 or 2), and all will always correspond to exact matches (no positional)
def matchVariants(inputStrings, isExact, isTrueExact, geneData, dataType):
    matchMap = {}
    exactMatches = []
    synMatches = []
    posMatches = []
    # Iterate CIVIC variants for given gene and attempt match of any input string to any CIVIC string
    for varName in geneData.keys():
        # Retrieve all available CIVIC strings (HGVS, manually generated, name, etc.) for the given CIVIC variant
        # For CNV, civicStrings will contain one single element corresponding to the CIVIC variant name
        civicStrings = geneData[varName]['match_hgvs']
        # Iterate input annotations strings and attempt match to any
        # The position of each input string corresponds to the type of match (true exact, synonym exact or positional)
        for indx,inputAnnot in enumerate(inputStrings):
            if inputAnnot in civicStrings:
#                 ## For CNV: especial case for cnv 'LOSS', depending on variant type we should disregard the match or not
#                 if (dataType == 'CNV') and (inputAnnot == 'LOSS'):
#                     ## CIVIC CNVs of type 'LOSS' will only be considered when they are associated to variant type 'transcript_ablation'
#                     if 'TRANSCRIPT_ABLATION' not in geneData[varName]['types']:
#                         ## Disregard incorrect CNV match
#                         continue

                # Determine type of match using its position and store accordingly
                # For CNV, all matches will be exact and true exact matches (either variant is in CIVIC or not)
                if isExact[indx]:
                    # For SNV: give preference to truly exact matches over descriptive synonym terms
                    # eg. report V600E over EXON 15 MUTATION
                    if isTrueExact[indx]:
                        # For SNV: True exact match (eg. V600E)
                        # For CNV: all matches will be exact and true exact matches (either variant is in CIVIC or not)
                        if varName not in exactMatches:
                            exactMatches.append(varName)
                    else:
                        # For SNV: descriptive term match (eg. EXON 15 MUTATION, FRAMESHIFT MUTATION)
                        if varName not in synMatches:
                            synMatches.append(varName)
                else:
                    # Positional match
                    if varName not in posMatches:
                        posMatches.append(varName)
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
                            if varName not in posMatches:
                                posMatches.append(varName)

    ## Once all CIVIC variants have been iterated, determine final tier and corresponding matched variants
    ## For CNV: either exactMatches or posMatches will occurr due to the implementation design

    # If there is at least one truly exact match, return this (highest preference)
    # For CNV: all exact matches will also be true exact matches (either CNV is in CIVIC or not). Multiple matches could occur for a single CNV (eg. DELETION + LOSS + COPY NUMBER)
    if exactMatches:
        tier = 1
        # There could be >1 perfect matches
        matchedVars = exactMatches
    # If there are no true exact matches but at least one synonym exact match, return this (still exact match)
    elif synMatches:
        tier = 1
        # There could be >1 perfect matches
        matchedVars = synMatches
    # If there are no exact matches of any kind but at least one positional, return this
    elif posMatches:
        tier = 2
        # For CNV: positional matches will correspond to EXON records (eg. EXON 1-2 DELETION, 3' EXON DELETION..)
        # There can be >1 positional matches
        matchedVars = posMatches
        # For SNV: check for presence of preferential positional matches, ie. general variants (like V600)
        if dataType == 'SNV':
            for var in posMatches:
                isGeneral = check_general_variant(var)
                if isGeneral:
                    # Stop as soon as a general variant is found and report only this
                    matchedVars = [var]
                    break

    # If no match was found, then tier case is 3, and return all relevant variants
    else:
        tier = 3
        ## For SNV: when input variant could not be matched in CIVIC (tier3), return all CIVIC variants associated
        ## to the given gene but that do not correspond to a CNV or EXPRESSION related variant
        if dataType == 'SNV':
#             matchedVars = list(geneData.keys())
            matchedVars = civic_return_all_snvs(geneData)
        ## For CNV: when input cnv could not be matched in CIVIC (tier3), return all CIVIC cnvs associated the given gene (if any)
        ## ie. not all CIVIC records are returned but only those that are matched to a 'CNV' tags
        elif dataType == 'CNV':
            matchedVars = civic_return_all_cnvs(geneData)

    return tier,matchedVars


# Check whether a variant name in CIVIC refers to a group of variants (eg. V600)
# Function only applicable to SNV data
def check_general_variant(varName):
    isGeneral = False
    if re.match(r'[A-Z]\d+($|\s\()', varName):
        isGeneral = True
    return isGeneral

