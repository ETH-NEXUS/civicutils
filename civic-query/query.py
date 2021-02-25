import sys
import os
import re


def check_identifier_type(identifier_type):
    """Check that a given identifier type is valid.

## TODO

    """

    if identifier_type not in ["entrez_id", "entrez_symbol", "civic_id"]:
        raise ValueError(
            f"'{identifier_type}' is not a valid identifier type.\n"
            f"Please, provide one of the following identifier types: 'entrez_id','entrez_symbol','civic_id'.\n"
        )


def reformat_results(results, identifier_type="entrez_symbol"):
    """Reformat records returned from CIVIC query into a dictionary.

    Arguments
    ---------
    results: list
        List of CIVIC records
    identifier_type: `str` (default: `entrez_symbol`)
        Type of gene identifier to be used

    Returns
    -------

## TODO

    """

# TODO: for now, keep the same filters as we had in the previous version
# TODO: to be done soon; do not filter anything but retrieve all the gene-variant records unchanged. we will apply filtering using functions at a later step

    varMap = {}
    retrieved_genes = []    # keep track of genes that could be retrieved from CIVIC
    no_variants = []        # keep track of genes retrieved from CIVIC but with no variants available
    all_variants = []       # keep track of all variants retrieved from CIVIC

    # Check that id type corresponds to one of the allowed options
    ignore = check_identifier_type(identifier_type)

    # Iterate individual gene records to retrieve associated variants and evidence information
    for gene_record in results:
        # Retrieve all ids associated to this gene (CIVIC, entrez, symbol)
        gene_civic_id = str(gene_record.id)
        gene_id = str(gene_record.entrez_id)
        # Use uppercase for consistency of the gene symbols
        gene_symbol = gene_record.name.strip().upper()
        # Retrieve variants associated to this gene (can be empty)
        gene_variants = gene_record.variants

        # Use the provided gene id (civic, entrez or symbol) to uniquely identify genes
        if identifier_type == "civic_id":
            gene_key = gene_civic_id
        if identifier_type == "entrez_id":
            gene_key = gene_id
        if identifier_type == "entrez_symbol":
            gene_key = gene_symbol

        # Keep track of gene ids which could be retrieved from CIVIC
        if gene_key not in retrieved_genes:
            retrieved_genes.append(gene_key)

        # NOTE: it seems that only genes having at least 1 variant record available are included in the offline cache (eg. gene ADORA1 has no variants and is found via API but not in the cache)
        # Skip genes that do not have any variants available in CIVIC
        if not gene_variants:
            if gene_key not in no_variants:
                no_variants.append(gene_key)
            continue

        # Iterate variant records associated to the current gene
        # Retrieve all relevant info listed for each variant
        for variant_record in gene_variants:
            # Internal variant id in CIVIC
            variant_id = str(variant_record.id)
            # Keep track of all variants retrieved from CIVIC
            if variant_id not in all_variants:
                all_variants.append(variant_id)
            # Variant name in CIVIC; use uppercase for consistency
            variant_name = variant_record.name.strip().upper()
            hgvs_expressions = variant_record.hgvs_expressions
            # Score to assess the accumulation of evidence for each variant (quantity and quality)
            civic_score = variant_record.civic_actionability_score
            # Sanity check for empty scores
            if civic_score is None:
                civic_score = "NULL"

            # List of evidence records available for the current variant (can be empty)
            evidence_items = variant_record.evidence_items

            # Iterate through the listed evidence items and store relevant information for this variant
            # Variants which do not have any clinical data associated to them will be directly skipped
            for evidence_record in evidence_items:

                # FIXME: Sanity check that all critical elements are present and non-empty (entrez_name, variant name, evidence_type, disease name)

                # Use uppercase for consistency of the tags
                evidence_status = evidence_record.status.strip().upper()
                evidence_type = evidence_record.evidence_type.strip().upper()
                disease = evidence_record.disease.name.strip().upper()

                evidence_level = evidence_record.evidence_level                     # just in case, should never be None
                variant_origin = evidence_record.variant_origin                     # can be None
                evidence_direction = evidence_record.evidence_direction             # can be None
                clinical_significance = evidence_record.clinical_significance       # can be None

                # Skip records that are not accepted evidence
                if (evidence_status != "ACCEPTED"):
                    continue

                # Skip records that correspond to germline variants
                # The variant_origin field might be blank/empty (None)
                if variant_origin:
                    if re.search("GERMLINE", variant_origin):
                        continue

                # Sanity check for empty evidence direction, clinical significance or level
                # 'NULL' is introduced to distinguish from 'N/A' tag
                if evidence_direction is None:
                    evidence_direction = "NULL"
                else:
                    evidence_direction = evidence_direction.strip().upper()
                if clinical_significance is None:
                    clinical_significance = "NULL"
                else:
                    clinical_significance = clinical_significance.strip().upper()
                if evidence_level is None:
                    evidence_level = "NULL"
                else:
                    evidence_level = evidence_level.strip().upper()

                # Combine the direction and significance of the evidence in one term
                evidence = evidence_direction + ':' + clinical_significance

                # At this point, currently evaluate evidence item has passed all the checks/filters
                # Keep track of the current evidence item under the corresponding variant and gene
                if gene_key not in varMap.keys():
                    varMap[gene_key] = {}

                # Variant name should be unique within gene (found some duplicates but all were submitted, not accepted data)
                if variant_name not in varMap[gene_key].keys():
                    varMap[gene_key][variant_name] = {}
                    varMap[gene_key][variant_name]['id'] = variant_id
                    varMap[gene_key][variant_name]['civic_score'] = civic_score

                    ## Generate list of strings that will be used to match our input variants in CIVIC
                    ## Returned list always has at least length=1 (in this case, containing only variant name)
                    ## For CNV, variant matching is not based on HGVS, so matchStrings will only contain the variant name
                    matchStrings = generate_civic_matchStrings(variant_name, hgvs_expressions, dataType)
                    varMap[gene_key][variant_name]['match_hgvs'] = matchStrings
                    # Keep original HGVS annotations (empty list when nothing is available)
                    # Use uppercase to avoid mismatches due to case
                    varMap[gene_key][variant_name]['hgvs'] = [h.strip().upper() for h in hgvs_expressions]

                    # Include associated variant types (sequence ontology terms). There can be multiple terms
                    varMap[gene_key][variant_name]['types'] = []
                    for vartype_record in variant_record.variant_types:
                        varMap[gene_key][variant_name]['types'].append(vartype_record.name.strip().upper())
                    # Account for empty variant types (can happen)
                    # 'NULL' is introduced to distinguish from 'N/A' tag
                    if not varMap[gene_key][variant_name]['types']:
                        varMap[gene_key][variant_name]['types'] = ["NULL"]

                # FIXME: there is no sanity check for detecting possible variant name duplicates
                if evidence_type not in varMap[gene_key][variant_name].keys():
                    varMap[gene_key][variant_name][evidence_type] = {}
                if disease not in varMap[gene_key][variant_name][evidence_type].keys():
                    varMap[gene_key][variant_name][evidence_type][disease] = {}

                drugs = []
                evidence_drugs = evidence_record.drugs
                for evidence_drug in evidence_drugs:
                    drug_name = evidence_drug.name.strip().upper()
                    if drug_name not in drugs:
                        drugs.append(drug_name)

                # When more than 1 drug are listed for the same evidence item, 'drug_interaction_type' is not null and defines the nature of this multiple drug entry
                drug_interaction = evidence_record.drug_interaction_type
                if drug_interaction is not None:
                    drug_interaction = drug_interaction.strip().upper()
                    # 'Substitutes' indicates that drugs can be considered individually
                    if drug_interaction != "SUBSTITUTES":
                        # Remaining terms ('Sequential' and 'Combination') indicate that drugs should be considered together, so join their names into a single tag
                        # Sort drugs alphabetically to ensure that their order in the combination treatment is always the same
                        drugs.sort()
                        drugs = ["+".join(drugs)]

                if not drugs:
                    # Only non-Predictive evidences and Predictive ones without drugs will have this dummy level
                    # Introduced for consistency purposes within the varMap structure
                    drugs = ["NULL"]

                # Iterate through drugs to add evidences associated to them
                #   For non-Predictive evidences or Predictive with empty drugs, drugs=['NULL']
                #   For Predictive and interaction=None, len(drugs) = 1
                #   For Predictive and interaction='Substitutes', len(drugs)>1
                #   For Predictive and interaction!='Substitutes', len(drugs)=1 (combiantion of several using '+')
                for drug in drugs:
                    if drug not in varMap[gene_key][variant_name][evidence_type][disease].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][drug] = {}
                    if evidence not in varMap[gene_key][variant_name][evidence_type][disease][drug].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][drug][evidence] = {}
                    if evidence_level not in varMap[gene_key][variant_name][evidence_type][disease][drug][evidence].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][drug][evidence][evidence_level] = []
                    # Group all publications associated to the same level. Do not check publication status
                    ## On 25.01.2019, source structure was changed to introduce ASCO abstracts as a source type
## TODO: sanity check for empty ID. Check for type of source?
                    varMap[gene_key][variant_name][evidence_type][disease][drug][evidence][evidence_level].append(evidence_record.source.citation_id.strip())

# TODO: iterate through assertions and repeat above process? They are not structured but (highly curated) free text

    return (varMap,retrieved_genes,no_variants,all_variants)


# Given a list of gene identifiers, query CIVIC for known variants and return a structured dictionary with the relevant results
# List of gene ids can be: CIVIC id, entrez id or gene symbol
def query_civic_genes(genes, identifier_type="entrez_symbol"):
    """Given a list of gene identifiers, query CIVIC for known variants and return a structured dictionary with the relevant results

## TODO

    Arguments
    ---------

    Returns
    -------
    """

    # Check that provided argument is a list (even if length = 1)
    if (not isinstance(genes, list)) or (not genes):
        raise TypeError(
            f"'{genes}' is not of type 'list'.\n"
        )

    # Check that id type corresponds to one of the allowed options
    ignore = check_identifier_type(identifier_type)

    ## Offline cache of CIVICdb should already be loaded

    # Querying functionality provided by module `civicpy` can only use internal CIVIC ids
    # NOTE: when using `civic.get_genes_by_ids()` with CIVIC ids, if any of them is not contained in the cache file, then it will try to directly query the CIVICdb, causing this script to crash

    # Workaround: retrieve all gene records available in the offline cache and parse them to match to input genes
    # Avoid at all costs directly querying the CIVIC db: even for CIVIC ids, we will match to available records from the offline cache
    all_results = civic.get_all_genes()

    # Iterate individual records and retrieve only those matching the provided gene ids
    results = []
    for gene_record in all_results:
        toKeep = False

        if (identifier_type == "civic_id"):
            gene_id = str(gene_record.id)                    # expectation is a single number
            if gene_id in genes:
                toKeep = True

        if (identifier_type == "entrez_id"):
            gene_id = str(gene_record.entrez_id)             # expectation is a single number
            if gene_id in genes:
                toKeep = True

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
                    toKeep = True
                    break

        if toKeep:
            results.append(gene_record)

    # At this point, all CIVIC results for queried genes have been retrieved in a list
    # Process gene records into a dictionary with structured format
    # gene -> variants -> evidence_items
    (dict_results,retrieved_genes,no_variants,all_variants) = reformat_results(results, identifier_type)
    return (dict_results,retrieved_genes,no_variants,all_variants)
