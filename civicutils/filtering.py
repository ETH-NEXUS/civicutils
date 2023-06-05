import sys
import os
import re

from civicutils.utils import check_string_filter_arguments, check_cutoff_filter_arguments, check_is_bool, check_keys, check_keys_not, check_is_none


def filter_in(field, field_name, in_list, list_name, match_type="exact"):
    """
    Check if a given field is contained in a list of strings to be selected for during filtering.
    :param field:       String to check in the provided list. If 'NULL', then match will not be attempted (i.e. if list is non-empty, then 'match=False' will be returned by the function).
    :param field_name:  Name of the field being checked (for error handling).
    :param in_list:     List of strings to be compared against field. If empty list is provided, then 'match=True' will be returned by the function.
    :param list_name:   Name of the list being checked (for error handling).
    :param match_type:  ['exact', 'partial']
                        exact:   	the exact field that was provided needs to be found in the list of strings in order to declare 'match=True'
                        partial: 	partial matches of the provided field in the list is enough to declare 'match=True'
                        Type of match to search for. This parameter defaults to 'exact'.
    :return:            Returns a boolean indicating if field should be kept (True) or not (None).
    """
    (field, in_list) = check_string_filter_arguments(field, field_name, in_list, list_name)
    match = False
    # Only perform check when in_list is not empty
    if in_list:
        # If a list was provided and field is empty, then filter is not passed
        if field != "NULL":
            if (match_type == "exact"):
                if field in in_list:
                    match = True
            if (match_type == "partial"):
                for tmp in in_list:
                    if tmp in field:
                        match = True
    else:
        match = True

    return match


def filter_not_in(field, field_name, out_list, list_name, match_type="exact"):
    """
    Check if a given field is contained in a list of strings to be removed during filtering.
    :param field:       String to check in the provided list. Can be 'NULL' (match will still be attempted).
    :param field_name:  Name of the field being checked (for error handling).
    :param out_list:    List of strings to be compared against field. If empty list is provided, then 'match=False' will be returned by the function.
    :param list_name:   Name of the list being checked (for error handling).
    :param match_type:  ['exact', 'partial']
                        exact:   	the exact field that was provided needs to be found in the list of strings in order to declare 'match=True'
                        partial: 	partial matches of the provided field in the list is enough to declare 'match=True'
                        Type of match to search for. This parameter defaults to 'exact'.
    :return:            Returns a boolean indicating if field should be excluded (True) or not (None).
    """
    (field, out_list) = check_string_filter_arguments(field, field_name, out_list, list_name)
    match = False
    # If no list was provided, even if field is empty, filter is always "passed" (i.e. False to avoid removing variant)
    # This allows to filter record if a particular field is "NULL"
    if out_list:
        if (match_type == "exact"):
            if field in out_list:
                match = True
        if (match_type == "partial"):
            for tmp in out_list:
                if tmp in field:
                    match = True

    return match


def filter_cutoff(field, field_name, cutoff, cutoff_name):
    """
    Check if a given numeric field is equal or above a certain threshold.
    :param field:	Number to compare against the provided cutoff. When cutoff!=0, return 'match=True' when field => cutoff.
    :param field_name:	Name of the numeric field being checked (for error handling).
    :param cutoff:	Numeric threshold to be compared against the given field. If cutoff=0, then 'match=True' will be returned by the function.
    :param cutoff_name:	Name of the cutoff value being checked (for error handling).
    :return:            Returns a boolean indicating if field is equal/above (True) or below (None) the threshold.
    """
    # Sanity check for provided arguments (both should be numeric or floats)
    check_cutoff_filter_arguments(field, field_name)
    check_cutoff_filter_arguments(cutoff, cutoff_name)
    match = True
    cutoff_f = float(cutoff)
    field_f = float(field)
    if (cutoff_f != float(0)):
        if (field_f < cutoff_f):
            match = False

    return match


def filter_civic(var_map, gene_id_in=[], gene_id_not_in=[], min_variants=0, var_id_in=[], var_id_not_in=[], var_name_in=[], var_name_not_in=[], min_civic_score=0, var_type_in=[], var_type_not_in=[], min_evidence_items=0, evidence_type_in=[], evidence_type_not_in=[], disease_in=[], disease_not_in=[], drug_name_in=[], drug_name_not_in=[], evidence_dir_in=[], evidence_dir_not_in=[], evidence_clinsig_in=[], evidence_clinsig_not_in=[], evidence_level_in=[], evidence_level_not_in=[], evidence_status_in=[], evidence_status_not_in=[], source_status_in=[], source_status_not_in=[], var_origin_in=[], var_origin_not_in=[], source_type_in=[], source_type_not_in=[], min_evidence_rating=0, output_empty=False):
    """
    Filter a nested dictionary of results from querying genes in CIViC. Dictionary is assumed to have been generated via function 'reformat_civic()'.
    Provided filtering parameters are evaluated in the order in which they are listed in the function definition (i.e. not in the order specified during the function call).
    The logic of any filters applied is always 'AND'. If the desired filter logic in not possible in one single function call, then the function needs to be applied several times subsequently.
    :param var_map:                 Nested dictionary of results from querying genes in CIViC. Provided dictionary cannot be annotated with disease specificity, and it is assumed to have been generated via 'reformat_civic()'. See README for more details about the specific structure of 'var_map'.
    :param gene_id_in:              List of gene identifiers to filter for (must match the type of identifier used in the provided 'var_map'). Default is empty list. If non-empty, then only gene records matching the provided identifiers will be returned (exact string match is sought).
    :param gene_id_not_in:          List of gene identifiers to filter out (must match the type of identifier used in the provided 'var_map'). Default is empty list. If non-empty, then any gene records matching the provided identifiers will be excluded from the returned dictionary (exact string match is sought).
    :param min_variants:            Minimum number of variants the need to be available for the gene. Default is zero. If non-zero, then only genes having an equal or larger number of variants will be returned.
    :param var_id_in:               List of variant CIViC identifiers to filter for. Default is empty list. If non-empty, then only variant records matching the provided identifiers will be returned (exact string match is sought).
    :param var_id_not_in:           List of variant CIViC identifiers to filter out. Default is empty list. If non-empty, then variant records matching the provided identifiers will be excluded from the returned dictionary (exact string match is sought).
    :param var_name_in:             List of variant names to filter for. Default is empty list. If non-empty, then only variant records matching the provided identifiers will be returned (partial string match is sought).
    :param var_name_not_in:         List of variant names to filter out. Default is empty list. If non-empty, then variant records matching the provided identifiers will be excluded from the returned dictionary (partial string match is sought).
    :param min_civic_score:         Minimum CIViC score required for the variant record. Default is zero. If non-zero, then only variant records having an equal or larger score will be returned.
    :param var_type_in:             List of variant types to filter for. Default is empty list. If non-empty, then only variant records associated to at least one type matching the provided list will be returned (partial string match is sought). When there are no variant types available for a given variant record, then 'NULL' will used in this case. Filtering can still be applied if desired, however current filter will always fail for these records.
    :param var_type_not_in:         List of variant types to filter out. Default is empty list. If non-empty, then variant records associated to at least one type matching the provided list will be excluded from the returned dictionary (partial string match is sought). When there are no variant types available for a given variant record, then 'NULL' will used in this case. Filtering can still be applied if desired, however current filter will never fail for these records.
    :param min_evidence_items:      Minimum number of evidence items that need to be available for the variant. Default is zero. If non-zero, then only variant records having an equal or larger number of evidences will be returned. This number is computed summing evidences available across all downstream layers (i.e. evidence types, directions, clinical significances, diseases, drugs, etc.). Filtering is applied directly on the original value from 'var_map' (i.e. before the number of evidences is updated to reflect the records that passed the filter selection).
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
    :param output_empty:            Boolean indicating if empty entries can be contained in the filtered dictionary returned by the function. It is recommended to only use this argument for checking at which level the different records contained in 'var_map' failed the applied filters. It is not recommended to use this argument when intention is to apply other functions from this packge on the returned dictionary (behavior can be unexpected).
    :return:                        Returns a new nested dictionary, containing the filtered results.
    """

    # Assume variant entries in var_map have the following structure in the first and second layers:
    var_map_entries_variant = ["name", "hgvs", "types"]
    var_map_entries_molecular_profile = ["civic_score", "n_evidence_items", "evidence_items"]
    # For checking that provided nested dictionary is not annotated for disease specificity and has been generated by reformat_civic()
    sorted_cts = ["ct", "gt", "nct"]

    # Sanity check that provided output_empty is logical
    check_is_none(output_empty, "output_empty")
    check_is_bool(output_empty, "output_empty")

    # Iterate complete dict of results and apply selected filters to generate a new (filtered) dict
    clean_map = {}
    for gene_id in var_map.keys():
        # Gene ID should always be available  (type of id is chosen during the query)
        gene_id_str = str(gene_id)
        keep_gene = filter_in(gene_id_str, "gene_id", gene_id_in, "gene_id_in", match_type="exact")
        if not keep_gene:
            continue
        remove_gene = filter_not_in(gene_id_str, "gene_id", gene_id_not_in, "gene_id_not_in", match_type="exact")
        if remove_gene:
            continue

        # Allow number of min_variants to be 0
        n_variants = len(var_map[gene_id].keys())
        keep_gene = filter_cutoff(n_variants, "n_variants", min_variants, "min_variants")
        if not keep_gene:
            continue

        # Current gene has passed all gene-level filters
        # Write to output (filtered) dict only when output_empty=True
        if output_empty:
            if gene_id not in clean_map.keys():
                clean_map[gene_id] = {}

        # Variant ID should always be available
        for var_id in var_map[gene_id].keys():
            var_id_str = str(var_id)
            keep_var = filter_in(var_id_str, "var_id", var_id_in, "var_id_in", match_type="exact")
            if not keep_var:
                continue
            remove_var = filter_not_in(var_id_str, "var_id", var_id_not_in, "var_id_not_in", match_type="exact")
            if remove_var:
                continue

            # Check that the expected entries are found in the dictionary
            check_keys(list(var_map[gene_id][var_id].keys()), "var_map", var_map_entries_variant, matches_all=False)

            # Variant name should always be available (never "NULL")
            variant = var_map[gene_id][var_id]["name"]
            keep_var = filter_in(variant, "variant", var_name_in, "var_name_in", match_type="partial")
            if not keep_var:
                continue
            remove_var = filter_not_in(variant, "variant", var_name_not_in, "var_name_not_in", match_type="partial")
            if remove_var:
                continue

            # Variant types will never be empty list
            # Use ["NULL"] when not available
            variant_types = var_map[gene_id][var_id]["types"]
            n_keep = 0
            n_remove = 0
            # When "NULL":
            #  - If any var_type_in was provided, then variant will fail this filter
            #  - If any var_type_not_in was provided, then variant will not fail this filter (given that the previous filter did not either)
            for var_type in variant_types:
                keep_type = filter_in(var_type, "var_type", var_type_in, "var_type_in", match_type="partial")
                if keep_type:
                    n_keep += 1
                remove_type = filter_not_in(var_type, "var_type", var_type_not_in, "var_type_not_in", match_type="partial")
                if remove_type:
                    n_remove += 1
            if (n_keep == 0):
                continue
            if (n_remove > 0):
                continue

            # Current variant has passed all variant-level filters
            # Write to output (filtered) dict only when output_empty=True
            if output_empty:
                if var_id not in clean_map[gene_id].keys():
                    clean_map[gene_id][var_id] = {}
                    clean_map[gene_id][var_id]["name"] = variant
                    clean_map[gene_id][var_id]["hgvs"] = var_map[gene_id][var_id]["hgvs"]
                    clean_map[gene_id][var_id]["types"] = variant_types


            molecular_profile_ids = set(list(var_map[gene_id][var_id].keys())) ^ set(var_map_entries_variant)

            for molecular_profil_id in molecular_profile_ids:
                # Check that the expected entries are found in the dictionary
                check_keys(list(var_map[gene_id][var_id][molecular_profil_id].keys()), "var_map", var_map_entries_molecular_profile, matches_all=False)
                
                # Allow number of evidence items to be 0
                n_evidence_items = var_map[gene_id][var_id][molecular_profil_id]["n_evidence_items"]
                keep_mp = filter_cutoff(n_evidence_items, "n_evidence_items", min_evidence_items, "min_evidence_items")
                if not keep_mp:
                    continue

                # CIViC Score is always a number (can be 0)
                var_score = var_map[gene_id][var_id][molecular_profil_id]["civic_score"]
                keep_mp = filter_cutoff(var_score, "var_score", min_civic_score, "min_civic_score")
                if not keep_mp:
                    continue

                n_evidence_items_after = 0
                evidence_types = list(var_map[gene_id][var_id][molecular_profil_id]["evidence_items"].keys())
                allowed_evidence_types = []
                
                if output_empty:
                    if molecular_profil_id not in clean_map[gene_id][var_id].keys():
                        clean_map[gene_id][var_id][molecular_profil_id] = {}
                        clean_map[gene_id][var_id][molecular_profil_id]["civic_score"] = var_score
                        clean_map[gene_id][var_id][molecular_profil_id]["n_evidence_items"] = 0
                        clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"] = {}
                
                # Evidence type should always be available
                for evidence_type in evidence_types:
                    keep_type = filter_in(evidence_type, "evidence type", evidence_type_in, "evidence_type_in", match_type="exact")
                    if not keep_type:
                        continue
                    remove_type = filter_not_in(evidence_type, "evidence type", evidence_type_not_in, "evidence_type_not_in", match_type="exact")
                    if remove_type:
                        continue
                    allowed_evidence_types.append(evidence_type)

                # Only work on subset of evidence types that passed filters (if any)
                for evidence_type in allowed_evidence_types:
                    # Current evidence type has passed filters
                    # Write to output (filtered) dict only when output_empty=True
                    if output_empty:
                        if evidence_type not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"].keys():
                            clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type] = {}

                    type_diseases = list(var_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys())
                    # Check that provided var_map is not annotated with disease specificity info (ct/gt/nct)
                    check_keys_not(type_diseases, "var_map", sorted_cts)

                    # Disease names should always be available
                    allowed_diseases = []
                    for type_disease in type_diseases:
                        keep_disease = filter_in(type_disease, "disease type", disease_in, "disease_in", match_type="partial")
                        if not keep_disease:
                            continue
                        remove_disease = filter_not_in(type_disease, "disease type", disease_not_in, "disease_not_in", match_type="partial")
                        if remove_disease:
                            continue
                        allowed_diseases.append(type_disease)

                    # Only work on subset of diseases that passed filters (if any)
                    for disease in allowed_diseases:
                        # Current disease name has passed filters
                        # Write to output (filtered) dict only when output_empty=True
                        if output_empty:
                            if disease not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys():
                                clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease] = {}
                        # All evidence types except "PREDICTIVE" will have a single entry "NULL" (because drugs are not associated to these evidence types)
                        # When "NULL", the drug-related filters do not apply (to avoid filtering evidences from types other than "PREDICTIVE"
                        disease_drugs = list(var_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease].keys())
                        if (len(disease_drugs) == 0):
                            continue                        
                        
                        allowed_drugs = []
                        if (evidence_type != "PREDICTIVE"):
                            allowed_drugs = [x for x in disease_drugs]   # ["NULL"]
                        else:
                            # Iterate existing drugs and apply filters
                            for disease_drug in disease_drugs:
                                keep_drug = filter_in(disease_drug, "drug", drug_name_in, "drug_name_in", match_type="partial")
                                if not keep_drug:
                                    continue
                                remove_drug = filter_not_in(disease_drug, "drug", drug_name_not_in, "drug_name_not_in", match_type="partial")
                                if remove_drug:
                                    continue
                                allowed_drugs.append(disease_drug)

                        # Only work on subset of drugs that passed filters (if any)
                        for drug in allowed_drugs:
                            # Current drug has passed filters ("NULL" whenever evidence type is not "PREDICTIVE")
                            # Write to output (filtered) dict only when output_empty=True
                            if output_empty:
                                if drug not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease].keys():
                                    clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug] = {}

                            evidences = list(var_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug].keys())
                            allowed_evidences = []
                            for evidence in evidences:
                                # Split the evidence string into direction and clinical significance
                                # Format: "DIRECTION:CLINSIGF"
                                evidence_list = evidence.strip().split(":")
                                # Sanity check for correct format of evidence string
                                if (len(evidence_list) != 2):
                                    raise ValueError("Unexpected format of evidence '%s'! Please provide string as 'EVIDENCE_DIRECTION:CLINICAL_SIGNIFICANCE'." %(evidence))

                                # Check evidence direction filters
                                direction = evidence_list[0]
                                keep_dir = filter_in(direction, "evidence direction", evidence_dir_in, "evidence_dir_in", match_type="exact")
                                if not keep_dir:
                                    continue
                                remove_dir = filter_not_in(direction, "evidence direction", evidence_dir_not_in, "evidence_dir_not_in", match_type="exact")
                                if remove_dir:
                                    continue
                                # Check evidence clinical significance filters
                                clin_signf = evidence_list[1]
                                keep_clin = filter_in(clin_signf, "clinical significance", evidence_clinsig_in, "evidence_clinsig_in", match_type="exact")
                                if not keep_clin:
                                    continue
                                remove_clin = filter_not_in(clin_signf, "clinical significance", evidence_clinsig_not_in, "evidence_clinsig_not_in", match_type="exact")
                                if remove_clin:
                                    continue
                                allowed_evidences.append(evidence)

                            # Only work on subset of evidences that passed filters (if any)
                            for this_evidence in allowed_evidences:
                                # Current evidence has passed filters
                                # Write to output (filtered) dict only when output_empty=True
                                if output_empty:
                                    if this_evidence not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug].keys():
                                        clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence] = {}

                                evidence_levels = list(var_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence].keys())
                                allowed_levels = []
                                for evidence_level in evidence_levels:
                                    keep_level = filter_in(evidence_level, "evidence level", evidence_level_in, "evidence_level_in",  match_type="exact")
                                    if not keep_level:
                                        continue
                                    remove_level = filter_not_in(evidence_level, "evidence level", evidence_level_not_in, "evidence_level_not_in", match_type="exact")
                                    if remove_level:
                                        continue
                                    allowed_levels.append(evidence_level)

                                # Only work on subset of evidence levels that passed filters (if any)
                                for level in allowed_levels:
                                    # Current evidence level has passed filters
                                    # Write to output (filtered) dict only when output_empty=True
                                    if output_empty:
                                        if level not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence].keys():
                                            clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence][level] = []

                                    # Format of list: ["TYPE_ID1:EVIDENCESTATUS1:SOURCESTATUS1:VARORIGIN1:RATING1",..,"TYPE_IDN:EVIDENCESTATUSN:SOURCESTATUSN:VARORIGINN:RATINGN"]
                                    all_evidence_items = var_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence][level]
                                    for evidence_item in all_evidence_items:
                                        # Split the evidence item string into the 5 separate fields
                                        item_list = evidence_item.strip().split(":")
                                        # Sanity check for correct format of evidence item string
                                        if (len(item_list) != 5):
                                            raise ValueError("Unexpected format of evidence item '%s'! Please provide string as 'EVIDENCE_ID:EVIDENCE_STATUS:SOURCE_STATUS:VARIANT_ORIGIN:RATING'." %(evidence_item))
                                        tmp_id = item_list[0]
                                        evidence_status = item_list[1]
                                        source_status = item_list[2]
                                        var_origin = item_list[3]
                                        rating = item_list[4]

                                        # Split the id string into the 2 separate fields
                                        # Format: "TYPE_ID"
                                        tmp_id_list = tmp_id.split("_")
                                        # Sanity check for correct format of evidence id string
                                        if (len(tmp_id_list) != 2):
                                            raise ValueError("Unexpected format of evidence id '%s'! Please provide string as 'PUBMED_[ID]' or 'ASCO_[ID]'." %(tmp_id))
                                        id_type = tmp_id_list[0]
                                        this_id = tmp_id_list[1]

                                        # Apply filters on evidence status
                                        keep_status = filter_in(evidence_status, "evidence status", evidence_status_in, "evidence_status_in", match_type="exact")
                                        if not keep_status:
                                            continue
                                        remove_status = filter_not_in(evidence_status, "evidence status", evidence_status_not_in, "evidence_status_not_in", match_type="exact")
                                        if remove_status:
                                            continue
                                        # Apply filters on source status
                                        keep_source = filter_in(source_status, "source status", source_status_in, "source_status_in", match_type="partial")
                                        if not keep_source:
                                            continue
                                        remove_source = filter_not_in(source_status, "source status", source_status_not_in, "source_status_not_in", match_type="partial")
                                        if remove_source:
                                            continue
                                        # Apply filters on variant origin
                                        keep_origin = filter_in(var_origin, "variant origin", var_origin_in, "var_origin_in", match_type="partial")
                                        if not keep_origin:
                                            continue
                                        remove_origin = filter_not_in(var_origin, "variant origin", var_origin_not_in, "var_origin_not_in", match_type="partial")
                                        if remove_origin:
                                            continue
                                        # Apply filters on source type
                                        keep_type = filter_in(id_type, "source id type", source_type_in, "source_type_in", match_type="exact")
                                        if not keep_type:
                                            continue
                                        remove_type = filter_not_in(id_type, "source id type", source_type_not_in, "source_type_not_in", match_type="exact")
                                        if remove_type:
                                            continue
                                       # When rating in not available (ie. "NULL"), the evidence will directly fail if corresponding filter is set
                                        if (rating == "NULL"):
                                            if (float(min_evidence_rating) != float(0)):
                                                continue
                                        # When rating is available, apply filters on evidence rating
                                        else:
                                            keep_rating = filter_cutoff(rating, "rating", min_evidence_rating, "min_evidence_rating")
                                            if not keep_rating:
                                                continue

                                        # Current evidence item has passed all filters
                                        # Write complete entry (from gene to evidence item) to output (filtered) dict only when output_empty>0
                                        if not output_empty:
                                            if gene_id not in clean_map.keys():
                                                clean_map[gene_id] = {}
                                            if var_id not in clean_map[gene_id].keys():
                                                clean_map[gene_id][var_id] = {}
                                                clean_map[gene_id][var_id]["name"] = variant
                                                clean_map[gene_id][var_id]["hgvs"] = var_map[gene_id][var_id]["hgvs"]
                                                clean_map[gene_id][var_id]["types"] = variant_types
                                            if molecular_profil_id not in clean_map[gene_id][var_id].keys():
                                                clean_map[gene_id][var_id][molecular_profil_id] = {}
                                                clean_map[gene_id][var_id][molecular_profil_id]["civic_score"] = var_score
                                                clean_map[gene_id][var_id][molecular_profil_id]["n_evidence_items"] = 0
                                                clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"] = {}                                                
                                            if evidence_type not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"].keys():
                                                clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type] = {}
                                            if disease not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys():
                                                clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease] = {}
                                            if drug not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease].keys():
                                                clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug] = {}
                                            if this_evidence not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug].keys():
                                                clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence] = {}
                                            if level not in clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence].keys():
                                                clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence][level] = []

                                        clean_map[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][drug][this_evidence][level].append(evidence_item)
                                        n_evidence_items_after += 1

                # At this point, all evidence items available for this variant have been parsed and filtered
                # Add number of evidence items remaining after filtering for this variant
                # Check if output_empty=True or not
                if not output_empty:
                    # In this case, corresponding gene and variant entries will only be available when n_items > 0 (avoid empty entries)
                    if (n_evidence_items_after > 0):
                        clean_map[gene_id][var_id][molecular_profil_id]["n_evidence_items"] = n_evidence_items_after
                else:
                    # In this case, corresponding gene and variant entries are always available
                    clean_map[gene_id][var_id][molecular_profil_id]["n_evidence_items"] = n_evidence_items_after

    return clean_map
