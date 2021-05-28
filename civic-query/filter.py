import sys
import os
import re

from .utils import check_string_filter_arguments

def filter_in(field, inList, matchType="exact"):
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


def filter_not_in(field, outList, matchType="exact"):
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
def filter_cutoff(field, cutoff):

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


def filter_civic(varMap, gene_id_in=[], gene_id_not_in=[], min_variants=0, var_id_in=[], var_id_not_in=[], var_name_in=[], var_name_not_in=[], min_civic_score=0, var_type_in=[], var_type_not_in=[], min_evidence_items=0, evidence_type_in=[], evidence_type_not_in=[], disease_in=[], disease_not_in=[], drug_name_in=[], drug_name_not_in=[], evidence_dir_in=[], evidence_dir_not_in=[], evidence_clinsig_in=[], evidence_clinsig_not_in=[], evidence_level_in=[], evidence_level_not_in=[], evidence_status_in=[], evidence_status_not_in=[], source_status_in=[], source_status_not_in=[], var_origin_in=[], var_origin_not_in=[], source_type_in=[], source_type_not_in=[], min_evidence_rating=0, output_empty=False):

## filters are applied in the same order as their corresponding arguments (as they are defined in the function, not the order specified in the function call)
## so, logic is always AND for all selected filters?

## if desired filter logic is not possible, then the function would need to be run several times, applyin the filters subsequently

# FIXME: this is not true
## all array or cutoff filters work identically, except for disease-related arguments (ie. disease_in=[], disease_not_in=[], alt_diseases=[])

# FIXME: warn about very hard filters -> if variant matches any filter, then it would be removed from final results


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
        keepGene = filter_in(gene_id_str, gene_id_in, matchType="exact")
        if not keepGene:
            continue
        removeGene = filter_not_in(gene_id_str, gene_id_not_in, matchType="exact")
        if removeGene:
            continue

        # allow number of min_variants to be 0
        n_variants = len(varMap[gene_id].keys())
        keepGene = filter_cutoff(n_variants, min_variants)
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
            keepVar = filter_in(var_id_str, var_id_in, matchType="exact")
            if not keepVar:
                continue
            removeVar = filter_not_in(var_id_str, var_id_not_in, matchType="exact")
            if removeVar:
                continue

            # variant name should always be available (never "NULL")
            variant = varMap[gene_id][var_id]["name"]
            keepVar = filter_in(variant, var_name_in, matchType="partial")
            if not keepVar:
                continue
            removeVar = filter_not_in(variant, var_name_not_in, matchType="partial")
            if removeVar:
                continue

            # civic score is always a number (can be 0)
            var_score = varMap[gene_id][var_id]["civic_score"]
            keepVar = filter_cutoff(var_score, min_civic_score)
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
                keepType = filter_in(var_type, var_type_in, matchType="partial")
                if keepType:
                    nKeep += 1
                removeType = filter_not_in(var_type, var_type_not_in, matchType="partial")
                if removeType:
                    nRemove += 1
            if (nKeep == 0):
                continue
            if (nRemove > 0):
                continue

            # allow number of evidence items to be 0
            n_evidence_items = varMap[gene_id][var_id]["n_evidence_items"]
            keepVar = filter_cutoff(n_evidence_items, min_evidence_items)
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
                keepType = filter_in(evidence_type, evidence_type_in, matchType="exact")
                if not keepType:
                    continue
                removeType = filter_not_in(evidence_type, evidence_type_not_in, matchType="exact")
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
                    keepDisease = filter_in(type_disease, disease_in, matchType="partial")
                    if not keepDisease:
                        continue
                    removeDisease = filter_not_in(type_disease, disease_not_in, matchType="partial")
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
                            keepDrug = filter_in(disease_drug, drug_name_in, matchType="partial")
                            if not keepDrug:
                                continue
                            removeDrug = filter_not_in(disease_drug, drug_name_not_in, matchType="partial")
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
                            keepDir = filter_in(direction, evidence_dir_in, matchType="exact")
                            if not keepDir:
                                continue
                            removeDir = filter_not_in(direction, evidence_dir_not_in, matchType="exact")
                            if removeDir:
                                continue
                            # check evidence clinical significance filters
                            clin_signf = evidenceArr[1]
                            keepClin = filter_in(clin_signf, evidence_clinsig_in, matchType="exact")
                            if not keepClin:
                                continue
                            removeClin = filter_not_in(clin_signf, evidence_clinsig_not_in, matchType="exact")
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
                                keepLevel = filter_in(evidence_level, evidence_level_in, matchType="exact")
                                if not keepLevel:
                                    continue
                                removeLevel = filter_not_in(evidence_level, evidence_level_not_in, matchType="exact")
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
                                    keepStatus = filter_in(evidence_status, evidence_status_in, matchType="exact")
                                    if not keepStatus:
                                        continue
                                    removeStatus = filter_not_in(evidence_status, evidence_status_not_in, matchType="exact")
                                    if removeStatus:
                                        continue
                                    # apply filters on source status
                                    keepSource = filter_in(source_status, source_status_in, matchType="partial")
                                    if not keepSource:
                                        continue
                                    removeSource = filter_not_in(source_status, source_status_not_in, matchType="partial")
                                    if removeSource:
                                        continue
                                    # apply filters on variant origin
                                    keepOrigin = filter_in(var_origin, var_origin_in, matchType="partial")
                                    if not keepOrigin:
                                        continue
                                    removeOrigin = filter_not_in(var_origin, var_origin_not_in, matchType="partial")
                                    if removeOrigin:
                                        continue
                                    # apply filters on source type
                                    keepType = filter_in(idType, source_type_in, matchType="exact")
                                    if not keepType:
                                        continue
                                    removeType = filter_not_in(idType, source_type_not_in, matchType="exact")
                                    if removeType:
                                        continue
                                    # when rating in not available (ie. 'NULL'), the evidence will directly fail if corresponding filter is set
                                    if (rating == "NULL"):
                                        if (float(min_evidence_rating) != float(0)):
                                            continue
                                    # when rating is available, apply filters on evidence rating
                                    else:
                                        keepRating = filter_cutoff(rating, min_evidence_rating)
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
