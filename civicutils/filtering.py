import sys
import os
import re

sys.path.insert(0, '/cluster/work/nexus/antoine/Projects/2022_12_integrate_Civicutils/git_files/civicutils/civicutils/')
from utils import check_string_filter_arguments,check_cutoff_filter_arguments,check_is_bool,check_keys,check_keys_not,check_is_none
sys.path.remove('/cluster/work/nexus/antoine/Projects/2022_12_integrate_Civicutils/git_files/civicutils/civicutils/')


def filter_in(field, fieldName, inList, listName, matchType="exact"):
    """
    Check if a given field is contained in a list of strings to be selected for during filtering.
    :param field:       A string to check in the provided list. If 'NULL', then match will not be attempted (ie. if list is non-empty, then 'match=False' will be returned by the function).
    :param fieldName:   The name of the field being checked (for error handling).
    :param inList:      A list of strings to be compared against field. If empty list is provided, then 'match=True' will be returned by the function.
    :param listName:    The name of the list being checked (for error handling).
    :param matchType:   ['exact', 'partial']
                        exact:   the exact provided field needs to be found in the list to declare 'match=True'.
                        partial: partial matches of field in the list is enough to declare 'match=True'.
                        Type of match to be sought.
                        This parameter defaults to 'exact'.
    :return:            Returns a boolean indicating if field should be kept (True) or not (None).
    """
    (field, inList) = check_string_filter_arguments(field, fieldName, inList, listName)
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


def filter_not_in(field, fieldName, outList, listName, matchType="exact"):
    """
    Check if a given field is contained in a list of strings to be removed during filtering.
    :param field:       A string to check in the provided list. Can be 'NULL' (match will still be attempted).
    :param fieldName:   The name of the field being checked (for error handling).
    :param outList:     A list of strings to be compared against field. If empty list is provided, then 'match=False' will be returned by the function.
    :param listName:    The name of the list being checked (for error handling).
    :param matchType:   ['exact', 'partial']
                        exact:   the exact provided field needs to be found in the list to declare 'match=True'.
                        partial: partial matches of field in the list is enough to declare 'match=True'.
                        Type of match to be sought.
                        This parameter defaults to 'exact'.
    :return:            Returns a boolean indicating if field should be excluded (True) or not (None).
    """
    (field, outList) = check_string_filter_arguments(field, fieldName, outList, listName)
    match = False
    # If no list was provided, even if field is empty, filter is always "passed" (ie. False to avoid removing variant)
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


def filter_cutoff(field, fieldName, cutoff, cutoffName):
    """
    Check if a given numeric field is equal or above a certain threshold.
    :param field:       A number to compare against the provided cutoff. When cutoff!=0, return 'match=True' when field => cutoff.
    :param fieldName:   The name of the numeric field being checked (for error handling).
    :param cutoff:      A numeric threshold to be compared against the given field. If cutoff=0, then 'match=True' will be returned by the function.
    :param cutoffName:  The name of the cutoff value being checked (for error handling).
    :return:            Returns a boolean indicating if field is equal/above (True) or below (None) the threshold.
    """
    # Sanity check for provided arguments (both should be numeric or floats)
    check_cutoff_filter_arguments(field, fieldName)
    check_cutoff_filter_arguments(cutoff, cutoffName)
    match = True
    cutoff_f = float(cutoff)
    field_f = float(field)
    if (cutoff_f != float(0)):
        if (field_f < cutoff_f):
            match = False
    return match


def filter_civic(varMap, gene_id_in=[], gene_id_not_in=[], min_variants=0, var_id_in=[], var_id_not_in=[], var_name_in=[], var_name_not_in=[], min_civic_score=0, var_type_in=[], var_type_not_in=[], min_evidence_items=0, evidence_type_in=[], evidence_type_not_in=[], disease_in=[], disease_not_in=[], drug_name_in=[], drug_name_not_in=[], evidence_dir_in=[], evidence_dir_not_in=[], evidence_clinsig_in=[], evidence_clinsig_not_in=[], evidence_level_in=[], evidence_level_not_in=[], evidence_status_in=[], evidence_status_not_in=[], source_status_in=[], source_status_not_in=[], var_origin_in=[], var_origin_not_in=[], source_type_in=[], source_type_not_in=[], min_evidence_rating=0, output_empty=False):
    """
    Filter a nested dictionary of results from querying genes in CIViCdb. Dictionary is assumed to have been generated via function 'reformat_civic()'. Provided filtering parameters are evaluated in the order in which they are listed in the function definition (ie. not in the order specified during the function call). The logic of any filters applied is always 'AND'. If the desired filter logic in not possible in one single function call, then the function needs to be applied several times subsequently.
    :param varMap:                  Nested dictionary of results from querying genes in CIViCdb. Provided dictionary cannot be annotated with disease specificity, and it is assumed to have been generated via 'reformat_civic()'.
    :param gene_id_in:              List of gene identifiers to filter for (must match the type of identifier used in the provided 'varMap'). Default is empty list. If non-empty, then only gene records matching the provided identifiers will be returned (exact string match is sought).
    :param gene_id_not_in:          List of gene identifiers to filter out (must match the type of identifier used in the provided 'varMap'). Default is empty list. If non-empty, then any gene records matching the provided identifiers will be excluded from the returned dictionary (exact string match is sought).
    :param min_variants:            Minimum number of variants the need to be available for the gene. Default is zero. If non-zero, then only genes having an equal or larger number of variants will be returned.
    :param var_id_in:               List of variant CIViC identifiers to filter for. Default is empty list. If non-empty, then only variant records matching the provided identifiers will be returned (exact string match is sought).
    :param var_id_not_in:           List of variant CIViC identifiers to filter out. Default is empty list. If non-empty, then variant records matching the provided identifiers will be excluded from the returned dictionary (exact string match is sought).
    :param var_name_in:             List of variant names to filter for. Default is empty list. If non-empty, then only variant records matching the provided identifiers will be returned (partial string match is sought).
    :param var_name_not_in:         List of variant names to filter out. Default is empty list. If non-empty, then variant records matching the provided identifiers will be excluded from the returned dictionary (partial string match is sought).
    :param min_civic_score:         Minimum CIViC score required for the variant record. Default is zero. If non-zero, then only variant records having an equal or larger score will be returned.
    :param var_type_in:             List of variant types to filter for. Default is empty list. If non-empty, then only variant records associated to at least one type matching the provided list will be returned (partial string match is sought). When there are no variant types available for a given variant record, then 'NULL' will used in this case. Filtering can still be applied if desired, however current filter will always fail for these records.
    :param var_type_not_in:         List of variant types to filter out. Default is empty list. If non-empty, then variant records associated to at least one type matching the provided list will be excluded from the returned dictionary (partial string match is sought). When there are no variant types available for a given variant record, then 'NULL' will used in this case. Filtering can still be applied if desired, however current filter will never fail for these records.
    :param min_evidence_items:      Minimum number of evidence items that need to be available for the variant. Default is zero. If non-zero, then only variant records having an equal or larger number of evidences will be returned. This number is computed summing evidences available across all downstream layers (ie. evidence types, directions, clinical significances, diseases, drugs, etc.). Filtering is applied directly on the original value from 'varMap' (ie. before the number of evidences is updated to reflect the records that passed the filter selection).
    :param evidence_type_in:        List of evidence types to filter for. Default is empty list. If non-empty, then only evidence records matching the provided types will be returned (exact string match is sought).
    :param evidence_type_not_in:    List of evidence types to filter out. Default is empty list. If non-empty, then evidence records matching the provided types will be excluded from the returned dictionary (exact string match is sought).
    :param disease_in:              List of disease names to filter for. Default is empty list. If non-empty, then only evidence records matching the provided diseases will be returned (partial string match is sought).
    :param disease_not_in:          List of disease names to filter out. Default is empty list. If non-empty, then evidence records matching the provided diseases will be excluded from the returned dictionary (partial string match is sought).
    :param drug_name_in:            List of drug names to filter for. Default is empty list. If non-empty, then only evidence records matching the provided drugs will be returned (partial string match is sought). Filtering at this level will only be applied for evidence records of type 'PREDICTIVE', otherwise 'NULL' is used and filtering is ignored. Drug name can still be 'NULL' for 'PREDICTIVE' evidences in some cases (eg. submitted evidence).
    :param drug_name_not_in:        List of drug names to filter out. Default is empty list. If non-empty, then evidence records matching the provided drugs will be excluded from the returned dictionary (partial string match is sought). Filtering at this level will only be applied for evidence records of type 'PREDICTIVE', otherwise 'NULL' is used and filtering is ignored. Drug name can still be 'NULL' for 'PREDICTIVE' evidences in some cases (eg. submitted evidence).
    :param evidence_dir_in:         List of evidence directions to filter for. Default is empty list. If non-empty, then only evidence records matching the provided directions will be returned (exact string match is sought).
    :param evidence_dir_not_in:     List of evidence directions to filter out. Default is empty list. If non-empty, then evidence records matching the provided directions will be excluded from the returned dictionary (exact string match is sought).
    :param evidence_clinsig_in:     List of clinical significances to filter for. Default is empty list. If non-empty, then only evidence records matching the provided significances will be returned (exact string match is sought).
    :param evidence_clinsig_not_in: List of clinical significances to filter out. Default is empty list. If non-empty, then evidence records matching the provided significances will be excluded from the returned dictionary (exact string match is sought).
    :param evidence_level_in:       List of evidence levels to filter for. Default is empty list. If non-empty, then only evidence records matching the provided levels will be returned (exact string match is sought).
    :param evidence_level_not_in:   List of evidence levels to filter out. Default is empty list. If non-empty, then evidence records matching the provided levels will be excluded from the returned dictionary (exact string match is sought).
    :param evidence_status_in:      List of evidence statuses to filter for. Default is empty list. If non-empty, then only evidence records matching the provided statuses will be returned (exact string match is sought).
    :param evidence_status_not_in:  List of evidence statuses to filter out. Default is empty list. If non-empty, then evidence records matching the provided statuses will be excluded from the returned dictionary (exact string match is sought).
    :param source_status_in:        List of source statuses to filter for. Default is empty list. If non-empty, then only evidence records matching the provided source statuses will be returned (partial string match is sought).
    :param source_status_not_in:    List of source statuses to filter out. Default is empty list. If non-empty, then evidence records matching the provided source statuses will be excluded from the returned dictionary (partial string match is sought).
    :param var_origin_in:           List of variant origins to filter for. Default is empty list. If non-empty, then only evidence records matching the provided variant origins will be returned (partial string match is sought).
    :param var_origin_not_in:       List of variant origins to filter out. Default is empty list. If non-empty, then evidence records matching the provided variant origins will be excluded from the returned dictionary (partial string match is sought).
    :param source_type_in:          List of types of sources to filter for. Default is empty list. If non-empty, then only evidence records matching the provided source types will be returned (exact string match is sought).
    :param source_type_not_in:      List of types of sources to filter out. Default is empty list. If non-empty, then evidence records matching the provided source types will be excluded from the returned dictionary (exact string match is sought).
    :param min_evidence_rating:     Minimum rating required for the evidence record. Default is zero. If non-zero, then only evidence records having an equal or larger rating will be returned. Some evidence records can have a 'NULL' rating, and in this case it will fail this filter if a cutoff!=0 is set.
    :param output_empty:            Boolean indicating if empty entries can be contained in the filtered dictionary returned by the function. It is recommended to only use this argument for checking at which level the different records contained in 'varMap' failed the applied filters. It is not recommended to use this argument when intention is to apply other functions from this packge on the returned dictionary (behavior can be unexpected).
    :return:                        Returns a new nested dictionary, containing the filtered results.
    """
# TODO: add framework for flexibly choosing what is filtered based of constructed argument?
# eg. feature_name (have a catalogue of available and throw error if not matched) + in/not_in + partial/exact (feature type would determine if expectation is list or cutoff value)
# TODO: add framework to apply filters based on order of arguments?
# eg. for each feature that is provided, filtering fuction is called, and the dict is subsequently filtered in this manner until all provided filters have been applied

    # Assume variant entries in varMap have the following structure of the first layer:
    varmap_entries = ['name','civic_score','hgvs','types','n_evidence_items','evidence_items']
    # For checking that provided nested dictionary is not annotated for disease specificity and has been generated by reformat_civic()
    sorted_cts = ["ct","gt","nct"]

    # Sanity check that provided output_empty is logical
    check_is_none(output_empty,"output_empty")
    check_is_bool(output_empty,"output_empty")

    # Iterate complete dict of results and apply selected filters to generate a new (filtered) dict
    cleanMap = {}
    for gene_id in varMap.keys():
        # gene id should always be available  (type of id is chosen during the query)
        gene_id_str = str(gene_id)
        keepGene = filter_in(gene_id_str, "gene_id", gene_id_in, "gene_id_in", matchType="exact")
        if not keepGene:
            continue
        removeGene = filter_not_in(gene_id_str, "gene_id", gene_id_not_in, "gene_id_not_in", matchType="exact")
        if removeGene:
            continue

        # allow number of min_variants to be 0
        n_variants = len(varMap[gene_id].keys())
        keepGene = filter_cutoff(n_variants, "n_variants", min_variants, "min_variants")
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
            keepVar = filter_in(var_id_str, "var_id", var_id_in, "var_id_in", matchType="exact")
            if not keepVar:
                continue
            removeVar = filter_not_in(var_id_str, "var_id", var_id_not_in, "var_id_not_in", matchType="exact")
            if removeVar:
                continue

            # Check that the expected entries are found in the dictionary
            check_keys(list(varMap[gene_id][var_id].keys()),"varMap",varmap_entries,matches_all=True)

            # variant name should always be available (never "NULL")
            variant = varMap[gene_id][var_id]["name"]
            keepVar = filter_in(variant, "variant", var_name_in, "var_name_in", matchType="partial")
            if not keepVar:
                continue
            removeVar = filter_not_in(variant, "variant", var_name_not_in, "var_name_not_in", matchType="partial")
            if removeVar:
                continue

            # civic score is always a number (can be 0)
            var_score = varMap[gene_id][var_id]["civic_score"]
            keepVar = filter_cutoff(var_score, "var_score", min_civic_score, "min_civic_score")
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
                keepType = filter_in(var_type, "var_type", var_type_in, "var_type_in", matchType="partial")
                if keepType:
                    nKeep += 1
                removeType = filter_not_in(var_type, "var_type", var_type_not_in, "var_type_not_in", matchType="partial")
                if removeType:
                    nRemove += 1
            if (nKeep == 0):
                continue
            if (nRemove > 0):
                continue

            # allow number of evidence items to be 0
            n_evidence_items = varMap[gene_id][var_id]["n_evidence_items"]
            keepVar = filter_cutoff(n_evidence_items, "n_evidence_items", min_evidence_items, "min_evidence_items")
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
                keepType = filter_in(evidence_type, "evidence type", evidence_type_in, "evidence_type_in", matchType="exact")
                if not keepType:
                    continue
                removeType = filter_not_in(evidence_type, "evidence type", evidence_type_not_in, "evidence_type_not_in", matchType="exact")
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
                # Check that provided varMap is not annotated with disease specificity info (ct/gt/nct)
                check_keys_not(type_diseases,"varMap",sorted_cts)

                # disease names should always be available
                allowed_diseases = []
                for type_disease in type_diseases:
                    keepDisease = filter_in(type_disease, "disease type", disease_in, "disease_in", matchType="partial")
                    if not keepDisease:
                        continue
                    removeDisease = filter_not_in(type_disease, "disease type", disease_not_in, "disease_not_in", matchType="partial")
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
                        allowed_drugs = [x for x in disease_drugs]   # ['NULL']
                    else:
                        # iterate existing drugs and apply filters
                        for disease_drug in disease_drugs:
                            keepDrug = filter_in(disease_drug, "drug", drug_name_in, "drug_name_in", matchType="partial")
                            if not keepDrug:
                                continue
                            removeDrug = filter_not_in(disease_drug, "drug", drug_name_not_in, "drug_name_not_in", matchType="partial")
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
                                raise ValueError("Unexpected format of evidence '%s'! Please provide string as 'EVIDENCE_DIRECTION:CLINICAL_SIGNIFICANCE'." %(evidence))

                            # check evidence direction filters
                            direction = evidenceArr[0]
                            keepDir = filter_in(direction, "evidence direction", evidence_dir_in, "evidence_dir_in", matchType="exact")
                            if not keepDir:
                                continue
                            removeDir = filter_not_in(direction, "evidence direction", evidence_dir_not_in, "evidence_dir_not_in", matchType="exact")
                            if removeDir:
                                continue
                            # check evidence clinical significance filters
                            clin_signf = evidenceArr[1]
                            keepClin = filter_in(clin_signf, "clinical significance", evidence_clinsig_in, "evidence_clinsig_in", matchType="exact")
                            if not keepClin:
                                continue
                            removeClin = filter_not_in(clin_signf, "clinical significance", evidence_clinsig_not_in, "evidence_clinsig_not_in", matchType="exact")
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
                                keepLevel = filter_in(evidence_level, "evidence level", evidence_level_in, "evidence_level_in",  matchType="exact")
                                if not keepLevel:
                                    continue
                                removeLevel = filter_not_in(evidence_level, "evidence level", evidence_level_not_in, "evidence_level_not_in", matchType="exact")
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
                                        raise ValueError("Unexpected format of evidence item '%s'! Please provide string as 'EVIDENCE_ID:EVIDENCE_STATUS:SOURCE_STATUS:VARIANT_ORIGIN:RATING'." %(evidence_item))
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
                                        raise ValueError("Unexpected format of evidence id '%s'! Please provide string as 'PUBMED_[ID]' or 'ASCO_[ID]'." %(tmpId))
                                    idType = tmpIdArr[0]
                                    this_id = tmpIdArr[1]

                                    # apply filters on evidence status
                                    keepStatus = filter_in(evidence_status, "evidence status", evidence_status_in, "evidence_status_in", matchType="exact")
                                    if not keepStatus:
                                        continue
                                    removeStatus = filter_not_in(evidence_status, "evidence status", evidence_status_not_in, "evidence_status_not_in", matchType="exact")
                                    if removeStatus:
                                        continue
                                    # apply filters on source status
                                    keepSource = filter_in(source_status, "source status", source_status_in, "source_status_in", matchType="partial")
                                    if not keepSource:
                                        continue
                                    removeSource = filter_not_in(source_status, "source status", source_status_not_in, "source_status_not_in", matchType="partial")
                                    if removeSource:
                                        continue
                                    # apply filters on variant origin
                                    keepOrigin = filter_in(var_origin, "variant origin", var_origin_in, "var_origin_in", matchType="partial")
                                    if not keepOrigin:
                                        continue
                                    removeOrigin = filter_not_in(var_origin, "variant origin", var_origin_not_in, "var_origin_not_in", matchType="partial")
                                    if removeOrigin:
                                        continue
                                    # apply filters on source type
                                    keepType = filter_in(idType, "source id type", source_type_in, "source_type_in", matchType="exact")
                                    if not keepType:
                                        continue
                                    removeType = filter_not_in(idType, "source id type", source_type_not_in, "source_type_not_in", matchType="exact")
                                    if removeType:
                                        continue
                                    # when rating in not available (ie. 'NULL'), the evidence will directly fail if corresponding filter is set
                                    if (rating == "NULL"):
                                        if (float(min_evidence_rating) != float(0)):
                                            continue
                                    # when rating is available, apply filters on evidence rating
                                    else:
                                        keepRating = filter_cutoff(rating, "rating", min_evidence_rating, "min_evidence_rating")
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
