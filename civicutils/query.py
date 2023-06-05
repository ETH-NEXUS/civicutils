import sys
import os
import re

from civicutils.utils import check_identifier_type, check_empty_field, check_argument, check_is_list


def query_civic(genes, identifier_type="entrez_symbol"):
    """
    Given a list of gene identifiers, query CIViC using 'civic.get_all_genes()' (from Python package CIViCpy) and return a nested dictionary summarizing the evidences from the database. Only a selection of the complete record is processed and returned.
    :param genes:    		List containing gene identifiers to query.
    :param identifier_type:     ['entrez_symbol', 'entrez_id', 'civic_id']
                        	entrez_symbol:  Entrez gene symbol
                        	entrez_id: 	Entrez gene identifier
                        	civic_id: 	CIViC internal identifier
                        	Type of gene identifier used in your query. This parameter defaults to 'entrez_symbol'.
    :return:            	Returns a nested dictionary from 'reformat_civic()'; see README for more details about the specific structure (i.e. var_map). 
    """
    # Check arguments
    check_argument(genes, "genes")
    check_is_list(genes, "genes")
    # Check that id type corresponds to one of the allowed options
    check_identifier_type(identifier_type)

    # Load CIViC offline cache file
    from civicpy import civic
    success = civic.load_cache(on_stale="ignore")
    if not success:
        raise ValueError("Could not load cache file from 'civicpy'!")

    # Offline cache of CIViC should already be loaded
    # NOTE: querying functionality provided by module "civicpy" can only use internal CIViC ids
    # NOTE: when using "civic.get_genes_by_ids()" with CIViC ids, if any of them is not contained in the cache file, then it will try to directly query the CIViC, causing this script to crash

    # Workaround: retrieve all gene records available in the offline cache and parse them to match to input genes
    # Avoid at all costs directly querying the CIViC: even for CIViC ids, we will match to available records from the offline cache
    all_results = civic.get_all_genes()

    # Iterate individual records and retrieve only those matching the provided gene ids
    results = []
    for gene_record in all_results:
        to_keep = False

        if (identifier_type == "civic_id"):
            gene_id = str(gene_record.id)                    # expectation is a single number
            if gene_id in genes:
                to_keep = True

        if (identifier_type == "entrez_id"):
            gene_id = str(gene_record.entrez_id)             # expectation is a single number
            if gene_id in genes:
                to_keep = True

        if (identifier_type == "entrez_symbol"):
            gene_id = gene_record.name.strip()          # expectation is single string
            aliases = gene_record.aliases               # expectation is list of strings
            # Perform union of all gene symbols available for the current record
            tmp_symbols = list(set([gene_id]) | set(aliases))
            # Try to match current gene record (any alias) to the provided gene symbols
            for tmp_symbol in tmp_symbols:
                # Use uppercase to ensure consistency of gene symbols
                this_symbol = tmp_symbol.strip().upper()
                if this_symbol in genes:
                    to_keep = True
                    break

        if to_keep:
            results.append(gene_record)

    # At this point, all CIViC results for queried genes have been retrieved in a list
    # Process gene records into a dictionary with structured format
    # gene -> variants -> evidence_items
    dict_results = reformat_civic(results, identifier_type)

    return dict_results


def reformat_civic(results, identifier_type="entrez_symbol"):
    """
    Given a list of results from querying genes in CIViC, reformat the returned records into a nested dictionary with a specific structure. Only a selection of the complete record from the database is processed and returned.
    :param results:    		List of objects of class 'civicpy.civic.Gene', returned from querying genes in CIViC via 'civic.get_all_genes()'. Can be empty.
    :param identifier_type:     ['entrez_symbol', 'entrez_id', 'civic_id']
                        	entrez_symbol:	Entrez gene symbol
                        	entrez_id: 	Entrez gene identifier
                        	civic_id: 	CIViC internal identifier
                        	Type of gene identifier used in your query. This parameter defaults to 'entrez_symbol'.
    :return:            	Returns a nested dictionary with a specific structure. See README for more details (i.e. var_map).
    """
    # Check arguments
    check_is_list(results, "results")
    # Check that id type corresponds to one of the allowed options
    check_identifier_type(identifier_type)

    # Returned nested dictionary
    var_map = {}

    # Iterate individual gene records to retrieve associated variants and evidence information
    for gene_record in results:
        # Retrieve all ids associated to this gene (CIViC, entrez or symbol)
        gene_civic_id = str(gene_record.id)
        gene_id = str(gene_record.entrez_id)
        # Use uppercase for consistency of the gene symbols
        gene_symbol = gene_record.name.strip().upper()
        # Retrieve variants associated to this gene (can be empty)
        gene_variants = gene_record.variants

        # Use the provided gene id (CIViC, entrez or symbol) to uniquely identify genes
        if identifier_type == "civic_id":
            gene_key = gene_civic_id
        if identifier_type == "entrez_id":
            gene_key = gene_id
        if identifier_type == "entrez_symbol":
            gene_key = gene_symbol

        # Keep track of genes records returned by the query, even if they have no variants associated
        if gene_key not in var_map.keys():
            var_map[gene_key] = {}

        # NOTE: it seems that only genes having at least 1 variant record available are included in the offline cache (e.g. gene ADORA1 has no variants and is found via API but not in the cache)
        # Skip genes that do not have any variants available in CIViC
        if not gene_variants:
            continue

        # Iterate variant records associated to the current gene
        # Retrieve all relevant info listed for each variant
        for variant_record in gene_variants:
            # Internal variant id in CIViC
            variant_id = str(variant_record.id)
            # Variant name in CIViC; use uppercase for consistency
            variant_name = variant_record.name.strip().upper()
            hgvs_expressions = variant_record.hgvs_expressions
            molecular_profiles = variant_record.molecular_profiles
            
            # Variant id should be unique even across genes
            if variant_id not in var_map[gene_key].keys():
                var_map[gene_key][variant_id] = {}
                var_map[gene_key][variant_id]["name"] = variant_name
                # Keep original HGVS annotations (empty list when nothing is available)
                # Use uppercase to avoid mismatches due to case
                var_map[gene_key][variant_id]["hgvs"] = [h.strip().upper() for h in hgvs_expressions]
                
                # Include associated variant types (sequence ontology terms). There can be multiple terms
                var_map[gene_key][variant_id]["types"] = []
                for vartype_record in variant_record.variant_types:
                    var_map[gene_key][variant_id]["types"].append(vartype_record.name.strip().upper())
                    # Account for empty variant types (can happen)
                    # "NULL" is introduced to distinguish from "N/A" tag
                if not var_map[gene_key][variant_id]["types"]:
                    var_map[gene_key][variant_id]["types"] = ["NULL"]
                
                # Iterate molecular profile associated to the current variant
                # Retrieve all relevant info listed for each molecular profile
                for molecular_profile in molecular_profiles:
                    
                    molecular_profile_name = molecular_profile.name.strip().upper()
                    # Internal molecular profile id in CIVIC
                    molecular_profile_id = str(molecular_profile.id)
                    # Score to assess the accumulation of evidence for each molecular profile (quantity and quality)
                    civic_score = molecular_profile.molecular_profile_score
                    # Sanity check for empty scores
                    if (civic_score is None) or (not civic_score):
                        civic_score = 0.0
                    # List of evidence records available for the current molecular profile (can be empty)
                    evidence_items = molecular_profile.evidence_items
                    n_items = len(evidence_items)

                    var_map[gene_key][variant_id][molecular_profile_id] = {}
                    var_map[gene_key][variant_id][molecular_profile_id]["name"] = molecular_profile_name
                    var_map[gene_key][variant_id][molecular_profile_id]["civic_score"] = civic_score

                    # Keep track of number of evidence items associated with the current molecular profile
                    var_map[gene_key][variant_id][molecular_profile_id]["n_evidence_items"] = n_items
                    var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"] = {}

                    # Iterate through the listed evidence items and store relevant information for this variant
                    # Variants which do not have any clinical data associated to them will be directly skipped
                    for evidence_record in evidence_items:
                        
                        # Use uppercase for consistency of the tags (except id which should be unique)
                        # These fields are expected to never be empty or None
                        evidence_type = evidence_record.evidence_type.strip().upper()
                        # Check if disease name is available for the current evidence record
                        # NOTE: newly introduced evidence type "FUNCTIONAL" can be associated to records without any disease name assigned to it. In this case, disease has no attribute "name" and it is an empty dict
                        if hasattr(evidence_record.disease, "name"):
                            disease = evidence_record.disease.name.strip().upper()
                        else:
                            disease = "NULL"
                        evidence_status = evidence_record.status.strip().upper()
                        source_type = evidence_record.source.source_type.strip().upper()
                        source_status = evidence_record.status.strip().upper()
                        evidence_id = str(evidence_record.source.citation_id).strip()        # expected to be entirely numeric

                        if evidence_type not in var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"].keys():
                            var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type] = {}
                        if disease not in var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type].keys():
                            var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type][disease] = {}

                        # Sanity check for fields that can be empty
                        # "NULL" is introduced to distinguish from "N/A" tag
                        evidence_level = check_empty_field(evidence_record.evidence_level).strip().upper()                     # string (just in case, should never be None)
                        variant_origin = check_empty_field(evidence_record.variant_origin).strip().upper()                     # string, can be None
                        evidence_direction = check_empty_field(evidence_record.evidence_direction).strip().upper()             # string, can be None
                        clinical_significance = check_empty_field(evidence_record.significance).strip().upper()                # string, can be None
                        evidence_rating = check_empty_field(evidence_record.rating)                                            # numeric, can be None
                        # Sanity check for cases when rating=0 (check will return "NULL")
                        if (isinstance(evidence_record.rating, int) or isinstance(evidence_record.rating, float)) and (evidence_rating == "NULL"):
                            evidence_rating = 0.0

                        # Combine the direction and significance of the evidence in one term
                        evidence = evidence_direction + ":" + clinical_significance

                        drugs = []
                        evidence_drugs = evidence_record.therapies
                        for evidence_drug in evidence_drugs:
                            drug_name = evidence_drug.name.strip().upper()
                            if drug_name not in drugs:
                                drugs.append(drug_name)

                        # When more than 1 drug are listed for the same evidence item, "drug_interaction_type" is not null and defines the nature of this multiple drug entry
                        drug_interaction = evidence_record.therapy_interaction_type
                        if drug_interaction is not None:
                            drug_interaction = drug_interaction.strip().upper()
                            # "Substitutes" indicates that drugs can be considered individually
                            if drug_interaction != "SUBSTITUTES":
                                # Remaining terms ("Sequential" and "Combination") indicate that drugs should be considered together, so join their names into a single tag
                                # Sort drugs alphabetically to ensure that their order in the combination treatment is always the same
                                drugs.sort()
                                drugs = ["+".join(drugs)]

                        if not drugs:
                            # Only non-Predictive evidences and Predictive ones without drugs will have this dummy level
                            # Introduced for consistency purposes within the varMap structure
                            drugs = ["NULL"]

                        # Sanity checks that only "PREDICTIVE" evidences have drugs associated
                        # Submitted evidence items can fulfill having "PREDICTIVE" evidence type and no drugs ("NULL")
                        if (evidence_type != "PREDICTIVE") and (drugs != ["NULL"]):
                            raise ValueError("Only evidences of type 'PREDICTIVE' can have drugs associated!")
                            
                # Iterate through drugs to add evidences associated to them
                #   For non-Predictive evidences or Predictive with empty drugs, drugs=["NULL"]
                #   For Predictive and interaction=None, len(drugs) = 1
                #   For Predictive and interaction="Substitutes", len(drugs)>1
                #   For Predictive and interaction!="Substitutes", len(drugs)=1 (combiantion of several using "+")
                for drug in drugs:
                    if drug not in var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type][disease].keys():
                        var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type][disease][drug] = {}
                    if evidence not in var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type][disease][drug].keys():
                        var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type][disease][drug][evidence] = {}
                    if evidence_level not in var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type][disease][drug][evidence].keys():
                        var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type][disease][drug][evidence][evidence_level] = []

                    # Group all publications associated to the same level
                    # Keep track of associated info: source type, source id, evidence status, publication status, variant origin, evidence rating
                    # Format: "TYPE_ID:EVIDENCESTATUS:SOURCESTATUS:VARORIGIN:RATING"
                    var_map[gene_key][variant_id][molecular_profile_id]["evidence_items"][evidence_type][disease][drug][evidence][evidence_level].append(source_type + "_" + str(evidence_id) + ":" + evidence_status + ":" + source_status + ":" + variant_origin + ":" + str(evidence_rating))

    return var_map
