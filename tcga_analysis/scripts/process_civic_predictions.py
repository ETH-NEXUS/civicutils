#!/usr/bin/env python

'''
Parse CIViCutils annotations for a set of samples
SNVs and CNVs provided separately

Lourdes Rosano, Feb 2022
'''

import sys
import os
import argparse
import re
import copy


# Define mapping of special cases where known drugs are referred to in CIViC using synonyms or special terms
# CIViC_synonym -> drug_name
# drug_synonyms_mapping = {'DOVITINIB DILACTIC ACID (TKI258 DILACTIC ACID)':'DOVITINIB', '5-FLUOROURACIL':'FLUOROURACIL', '5-FU':'FLUOROURACIL', 'ADO-TRASTUZUMAB EMTANSINE':'TRASTUZUMAB EMTANSINE', 'PD0325901':'PD-0325901', 'PD173074':'PD-173074', 'BGJ-398':'INFIGRATINIB', 'BGJ398':'INFIGRATINIB'}

## Dictionary that allows prioritization of CIViC support categories (in oncoprint) when >1 gene has CIViC info for a given drug+sample
# support_mapping = {'civic_support':1, 'civic_resistance':1, 'civic_conflict':2, 'civic_unknown':3, 'civic_unspecific_vars':4, 'dgidb_only':5}



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


# Process evidence string reported by CIViCutils and retrieve the evidence direction, clinical significance and either the drug name of the cancer-specificity classification of the disease (depending on the given evidence type)
def get_clinical_info(evidence_string, has_drug=False):
    # Sanity check the expected separator character and format
    # Assume direction and clinical significance are always separated by ","
    evidence_split = evidence_string.split(",")

    # NOTE: take into account CIViCutils can report several evidence levels and/or publication ids aggregated using separator character ","
    # E.g.'..(B(PUBMED_17590872:..,PUBMED_19903786:..),C(ASCO_19223544:..))..'
    if len(evidence_split) > 2:
        # Remove entries consisting only of PMIDs to have clinical significance in the last possible position
        # NOTE: assumes that citation ids are always a series of digits and that evidence levels always correspond to one particular letter

# FIXME: review regex below to ensure it will always work for the new format of CIViC evidences
        evidence_split = [x for x in evidence_split if not re.match('([A-Z]\()?[A-Z]+_\d+(\)+)?', x)]

    clinical_signf_split = evidence_split[-1].strip().split("(")
    clinical_signf = clinical_signf_split[0]

    # NOTE: depending on the type of evidence being parsed, the format can be slighly different: either the drug name is available (only 'Predictive' evidence) or the cancer-specificity classification (other evidence types)
    cond_value = None

    if has_drug:
        # NOTE: take into account that drug names in CIViC can contain separator character ","
        # Once publication ids have been removed using get_clinical_significance(), evidences involving a drug name with "," will still have length > 2
        # E.g. "BEZ235 (NVP-BEZ235, DACTOLISIB)"
        if len(evidence_split) > 2:
            evidence_split = [",".join(evidence_split[0:-1])]

        # Split the remaining evidence string further and retrieve the drug name reported from CIViC
        # E.g. "CETUXIMAB(SUPPORTS,RESISTANCE(..."
        interim_direction_split = evidence_split[0].strip().split("(")

        # NOTE: take into account that drug names in CIViC can contain separator characters "(" and ")"
        # E.g. "RAPAMYCIN (SIROLIMUS)"
        interim_drug = interim_direction_split[0].strip()
        if len(interim_direction_split) > 2:
            interim_drug = "(".join(interim_direction_split[0:-1])
        # Use uppercase to avoid mismatches of the drug names due to case
        cond_value = interim_drug.upper()

    else:
        # Split the remaining evidence string further and retrieve the cancer specificity classification reported by CIViCutils
        # E.g. "NCT(SUPPORTS,RESISTANCE(..."
        interim_direction_split = evidence_split[0].strip().split("(")
        # Use lowercase to avoid mismatching of "ct" classes due to case
        ct_type = interim_direction_split[0].strip().lower()
        cond_value = ct_type

    # In both cases (drug available or not), the direction reported from CIViC uses the same format
    direction = interim_direction_split[-1].strip()

    # Sanity check that the term associated with the direction (either drug name or ct type) could be successfully retrieved
    if cond_value is None:
        raise ValueError("Encountered unexpected format while parsing evidence string '%s'" %(evidence_string))

    return (direction, clinical_signf, cond_value)


# Per tier, get total number of features parsed for the sample (only makes sense to compute mean for tiers != 4)
# Feature can be e.g. diseases or drugs
def process_feature_per_tier(input_mapping):
    n_items_tier1 = 0.0
    n_items_tier1b = 0.0
    n_items_tier2 = 0.0
    n_items_tier3 = 0.0
    # tier -> feature -> None
    if "tier_1" in input_mapping.keys():
        n_items_tier1 = float(len(input_mapping["tier_1"].keys()))
    if "tier_1b" in input_mapping.keys():
        n_items_tier1b = float(len(input_mapping["tier_1b"].keys()))
    if "tier_2" in input_mapping.keys():
        n_items_tier2 = float(len(input_mapping["tier_2"].keys()))
    if "tier_3" in input_mapping.keys():
        n_items_tier3 = float(len(input_mapping["tier_3"].keys()))
    return (n_items_tier1, n_items_tier1b, n_items_tier2, n_items_tier3)


# Per "ct" class, get total number of features parsed for the sample
# Feature can be e.g. diseases or drugs
def process_feature_per_ct(input_mapping):
    n_items_ct = 0.0
    n_items_gt = 0.0
    n_items_nct = 0.0
    # ct -> feature -> None
    if "ct" in input_mapping.keys():
        n_items_ct = float(len(input_mapping["ct"].keys()))
    if "gt" in input_mapping.keys():
        n_items_gt = float(len(input_mapping["gt"].keys()))
    if "nct" in input_mapping.keys():
        n_items_nct = float(len(input_mapping["nct"].keys()))
    return (n_items_ct, n_items_gt, n_items_nct)


# Per tier and "ct" class, get total number of features parsed for the sample (only makes sense to compute mean for tiers != 4)
# Feature can be e.g. diseases or drugs
def process_feature_per_tier_and_ct(input_mapping):
    n_items_tier1_ct = 0.0
    n_items_tier1_gt = 0.0
    n_items_tier1_nct = 0.0
    n_items_tier1b_ct = 0.0
    n_items_tier1b_gt = 0.0
    n_items_tier1b_nct = 0.0
    n_items_tier2_ct = 0.0
    n_items_tier2_gt = 0.0
    n_items_tier2_nct = 0.0
    n_items_tier3_ct = 0.0
    n_items_tier3_gt = 0.0
    n_items_tier3_nct = 0.0
    # tier -> ct -> feature -> None
    if "tier_1" in input_mapping.keys():
        (n_items_tier1_ct, n_items_tier1_gt, n_items_tier1_nct) = process_feature_per_ct(input_mapping["tier_1"])
    if "tier_1b" in input_mapping.keys():
        (n_items_tier1b_ct, n_items_tier1b_gt, n_items_tier1b_nct) = process_feature_per_ct(input_mapping["tier_1b"])
    if "tier_2" in input_mapping.keys():
        (n_items_tier2_ct, n_items_tier2_gt, n_items_tier2_nct) = process_feature_per_ct(input_mapping["tier_2"])
    if "tier_3" in input_mapping.keys():
        (n_items_tier3_ct, n_items_tier3_gt, n_items_tier3_nct) = process_feature_per_ct(input_mapping["tier_3"])

    return (n_items_tier1_ct, n_items_tier1_gt, n_items_tier1_nct, n_items_tier1b_ct, n_items_tier1b_gt, n_items_tier1b_nct, n_items_tier2_ct, n_items_tier2_gt, n_items_tier2_nct, n_items_tier3_ct, n_items_tier3_gt, n_items_tier3_nct)


# Per tier, compute mean feature for the sample (only makes sense to compute mean for tiers != 4)
# Mean feature can be e.g. mean number of matched variants, diseases or drugs
def process_mean_feature_per_tier(input_mapping, n_tier1, n_tier1b, n_tier2, n_tier3):
    mean_feature_tier1 = 0.0
    mean_feature_tier1b = 0.0
    mean_feature_tier2 = 0.0
    mean_feature_tier3 = 0.0
    # tier -> # feature
    # Sanity check for divisions by 0
    if ("tier_1" in input_mapping.keys()) and n_tier1:
        mean_feature_tier1 = float(float(input_mapping["tier_1"]) / float(n_tier1))
    if ("tier_1b" in input_mapping.keys()) and n_tier1b:
        mean_feature_tier1b = float(float(input_mapping["tier_1b"]) / float(n_tier1b))
    if ("tier_2" in input_mapping.keys()) and n_tier2:
        mean_feature_tier2 = float(float(input_mapping["tier_2"]) / float(n_tier2))
    if ("tier_3" in input_mapping.keys()) and n_tier3:
        mean_feature_tier3 = float(float(input_mapping["tier_3"]) / float(n_tier3))
    return (mean_feature_tier1, mean_feature_tier1b, mean_feature_tier2, mean_feature_tier3)


# Per "ct" class, compute mean feature for the sample
# Mean feature can be e.g. mean number of matched variants, diseases or drugs
def process_mean_feature_per_ct(input_mapping, n_variants):
    mean_feature_ct = 0.0
    mean_feature_gt = 0.0
    mean_feature_nct = 0.0
    # ct -> # feature
    # Sanity check for divisions by 0
    if ("ct" in input_mapping.keys()) and n_variants:
        mean_feature_ct = float(float(input_mapping["ct"]) / float(n_variants))
    if ("gt" in input_mapping.keys()) and n_variants:
        mean_feature_gt = float(float(input_mapping["gt"]) / float(n_variants))
    if ("nct" in input_mapping.keys()) and n_variants:
        mean_feature_nct = float(float(input_mapping["nct"]) / float(n_variants))
    return (mean_feature_ct, mean_feature_gt, mean_feature_nct)


# Per tier and "ct" class, compute mean feature for the current sample (only makes sense to compute mean for tiers != 4)
# Mean feature can be e.g. mean number of matched variants, diseases or drugs
def process_mean_feature_per_tier_and_ct(input_mapping, n_tier1, n_tier1b, n_tier2, n_tier3):
    mean_feature_tier1_ct = 0.0
    mean_feature_tier1_gt = 0.0
    mean_feature_tier1_nct = 0.0
    mean_feature_tier1b_ct = 0.0
    mean_feature_tier1b_gt = 0.0
    mean_feature_tier1b_nct = 0.0
    mean_feature_tier2_ct = 0.0
    mean_feature_tier2_gt = 0.0
    mean_feature_tier2_nct = 0.0
    mean_feature_tier3_ct = 0.0
    mean_feature_tier3_gt = 0.0
    mean_feature_tier3_nct = 0.0
    # tier -> ct -> # feature
    if "tier_1" in input_mapping.keys():
        (mean_feature_tier1_ct, mean_feature_tier1_gt, mean_feature_tier1_nct) = process_mean_feature_per_ct(input_mapping["tier_1"], n_tier1)
    if "tier_1b" in input_mapping.keys():
        (mean_feature_tier1b_ct, mean_feature_tier1b_gt, mean_feature_tier1b_nct) = process_mean_feature_per_ct(input_mapping["tier_1b"], n_tier1b)
    if "tier_2" in input_mapping.keys():
        (mean_feature_tier2_ct, mean_feature_tier2_gt, mean_feature_tier2_nct) = process_mean_feature_per_ct(input_mapping["tier_2"], n_tier2)
    if "tier_3" in input_mapping.keys():
        (mean_feature_tier3_ct, mean_feature_tier3_gt, mean_feature_tier3_nct) = process_mean_feature_per_ct(input_mapping["tier_3"], n_tier3)

    return (mean_feature_tier1_ct, mean_feature_tier1_gt, mean_feature_tier1_nct, mean_feature_tier1b_ct, mean_feature_tier1b_gt, mean_feature_tier1b_nct, mean_feature_tier2_ct, mean_feature_tier2_gt, mean_feature_tier2_nct, mean_feature_tier3_ct, mean_feature_tier3_gt, mean_feature_tier3_nct)


# Parse and process CIViCutils annotations reported for the variants of a given sample, in the provided input file
def parse_input_file(sample_file, sample_name, civic_info_mapping):
    print("Sample %s. File: %s" %(sample_name, sample_file))
    sorted_cts = ["ct","gt","nct"]              # define order of priority of "ct" classes assigned by CIViCutils

    infile = open(sample_file,'r')
    # Assume header containing a specific format and column names
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

    ## Keep track of several different annotations reported by CIViCutils

    all_variants = 0.0                            # all lines
    all_civic_variants = 0.0                      # all lines having CIViC info available (tier != 4)
    n_tier_1 = 0.0                                # all lines with tier = 1
    n_tier_1b = 0.0                               # all lines with tier = 1b
    n_tier_1_agg = 0.0                            # all lines with tier = 1 or tier = 1b
    n_tier_2 = 0.0                                # all lines with tier = 2
    n_tier_3 = 0.0                                # all lines with tier = 3
    n_tier_4 = 0.0                                # all lines with tier = 4
    n_predictive = 0.0                            # all lines with 'Predictive' CIViC evidence available
    n_diagnostic = 0.0                            # all lines with 'Diagnostic' CIViC evidence available
    n_prognostic = 0.0                            # all lines with 'Prognostic' CIViC evidence available
    n_predisposing = 0.0                          # all lines with 'Predisposing' CIViC evidence available
    n_predictive_no_tier3 = 0.0                   # all lines with 'Predictive' CIViC evidence available (excluding tier3 matches which can introduce biases)
    n_diagnostic_no_tier3 = 0.0                   # all lines with 'Diagnostic' CIViC evidence available (excluding tier3 matches which can introduce biases)
    n_prognostic_no_tier3 = 0.0                   # all lines with 'Prognostic' CIViC evidence available (excluding tier3 matches which can introduce biases)
    n_predisposing_no_tier3 = 0.0                 # all lines with 'Predisposing' CIViC evidence available (excluding tier3 matches which can introduce biases)

    matched_variants = 0.0                        # keep track of the total number of variants matched overall for the current sample
    # tier -> # variants
    matched_variants_mapping = {}                 # keep track of the total number of variants matched per tier for the current sample

    # disease -> ct
    disease_mapping = {}                          # keep track of all disease names parsed across the current sample
    # disease -> ct
    disease_mapping_no_tier3 = {}                 # keep track of all disease names parsed across the current sample (excluding tier3 matches which can introduce biases)
    # tier -> disease -> None
    per_tier_disease_mapping = {}                 # per tier, keep track of all disease names parsed across the current sample
    # ct -> disease -> None
    ct_mapping = {}                               # keep track of all disease names assigned to each "ct" class for the current sample
    # ct -> disease -> None
    ct_mapping_no_tier3 = {}                      # keep track of all disease names assigned to each "ct" class for the current sample (excluding tier3 matches which can introduce biases)
    # tier -> ct -> disease -> None
    per_tier_ct_mapping = {}                      # per tier, keep track of all disease names assigned to each "ct" class for the current sample

    matched_diseases = 0.0                        # keep track of the total number of unique disease names matched overall across all variants for the current sample
    # ct -> # diseases
    matched_diseases_ct_mapping = {}              # keep track of the total number of unique disease names matched per "ct" class across all variants for the current sample
    # ct -> # diseases
    matched_diseases_ct_mapping_no_tier3 = {}     # keep track of the total number of unique disease names matched per "ct" class across all variants for the current sample (excluding tier3 matches which can introduce biases)
    # tier -> # diseases
    per_tier_matched_diseases_mapping = {}        # keep track of the total number of unique disease names matched per tier for the current sample
    # tier -> ct -> # diseases
    per_tier_matched_diseases_ct_mapping = {}     # keep track of the total number of unique disease names matched per tier and "ct" class for the current sample

    n_drug_info = 0.0                             # all lines having consensus drug support info available
    n_drug_info_no_tier3 = 0.0                    # all lines having consensus drug support info available (excluding tier3 matches which can introduce biases)
    # drug -> ct -> consensus_support
    consensus_drug_mapping = {}                   # keep track of the total number of unique drug names parsed across the consensus drug support for the current sample
    # drug -> ct -> consensus_support
    consensus_drug_mapping_no_tier3 = {}          # keep track of the total number of unique drug names parsed across the consensus drug support for the current sample (excluding tier3 matches which can introduce biases)
    # drug -> ct -> consensus_support
    prior_consensus_drug_mapping = {}             # keep track of the total number of unique drug names parsed across the consensus drug support of the highest available "ct" class for the current sample
    # drug -> ct -> consensus_support
    prior_consensus_drug_mapping_no_tier3 = {}    # keep track of the total number of unique drug names parsed across the consensus drug support of the highest available "ct" class for the current sample (excluding tier3 matches which can introduce biases)
    # tier -> drug -> None
    per_tier_consensus_drug_mapping = {}          # per tier, keep track of the total number of unique drug names parsed across the consensus drug support for the current sample
    # tier -> drug -> None
    per_tier_prior_consensus_drug_mapping = {}    # per tier, keep track of the total number of unique drug names parsed across the consensus drug support of the highest available "ct" class for the current sample
    # ct -> drug -> consensus_support
    consensus_ct_mapping = {}                     # keep track of the total number of unique drug names parsed per "ct" class across the consensus drug support for the current sample
    # ct -> drug -> consensus_support
    consensus_ct_mapping_no_tier3 = {}            # keep track of the total number of unique drug names parsed per "ct" class across the consensus drug support for the current sample (excluding tier3 matches which can introduce biases)
    # ct -> drug -> consensus_support
    prior_consensus_ct_mapping = {}               # keep track of the total number of unique drug names parsed per highest available "ct" class across the consensus drug support for the current sample
    # ct -> drug -> consensus_support
    prior_consensus_ct_mapping_no_tier3 = {}      # keep track of the total number of unique drug names parsed per highest available "ct" class across the consensus drug support for the current sample (excluding tier3 matches which can introduce biases)
    # tier -> ct -> drug -> None
    per_tier_consensus_ct_mapping = {}            # per tier, keep track of the total number of unique drug names parsed per "ct" class across the consensus drug support for the current sample
    # tier -> ct -> drug -> None
    per_tier_prior_consensus_ct_mapping = {}      # per tier, keep track of the total number of unique drug names parsed per highest available "ct" class across the consensus drug support for the current sample


    # Each line in the input file corresponds to a single variant in the genome
    for line in infile:
        all_variants += 1.0
        line_split = line.strip().split("\t")

        ## 1) Process tier of the variant match

        tier = str(line_split[tier_pos].strip())
        # Avoid having issues due to using numbers/strings for the tiers
        tier = "tier_" + tier
        # Tier=4 should be skipped as no information was found on CIViC for the current gene
        # (i.e. all associated columns will be empty)
        if tier=="tier_4":
            n_tier_4 += 1.0
            continue

        # Keep track of the number of variants assigned to each tier in the current sample file
        all_civic_variants += 1.0
        if tier=="tier_3":
            n_tier_3 += 1.0
        if tier=="tier_2":
            n_tier_2 += 1.0
        if tier=="tier_1b":
            n_tier_1b += 1.0
        if tier=="tier_1":
            n_tier_1 += 1.0
        if tier=="tier_1" or tier=="tier_1b":
            n_tier_1_agg += 1.0


        ## 2) Process number of CIViC variants matched per line (and associated tier)

        # Parse column 'CIViC_Score', assumed to contain all CIViC variants matched to the current line (and their associated scores in CIViC)
        # NOTE: assume that variant names are listed using ";" as a separator character, and that duplicates are not possible
        civic_scores = str(line_split[civic_score_pos].strip())
        if civic_scores == ".":
            if not (tier=="tier_3" or tier=="tier_4"):
                raise ValueError("Encountered unexpected case of variant with tier!=3 and tier!=4 but no associated variant matches from CIViC in line %s" %(line.strip()))
            n_variants = 0.0
        else:
            # NOTE: for SNVs, all variant matches are ensured to originate from a single variant annotation + gene
            # NOTE: for CNVs, several variant matches arising from different genes (but same associated tier) can happen within the same line
            civic_score_list = civic_scores.split(';')
            n_variants = float(len(civic_score_list))

            # Sanity check there are no duplicated CIViC variant entries reported within the same line
            if (str(float(n_variants)) != str(float(len(set(civic_score_list))))):
                raise ValueError("Encountered duplicated variant matches from CIViC in line %s" %(line.strip()))
            # Sanity check that at least 1 variant should have been matched if column is not empty
            if n_variants == 0.0:
                raise ValueError("Encountered unexpected case of no associated variant matches from CIViC in line %s" %(line.strip()))

        # Keep track of number of matched variants, also keep counts per tier
        # tier -> # CIViC variants matched
        if tier not in matched_variants_mapping.keys():
            matched_variants_mapping[tier] = 0.0
        matched_variants_mapping[tier] += n_variants
        matched_variants += n_variants


        ## 3) Process disease names and associated cancer specificity classifications across all evidence types (columns)

        # disease -> None
        interim_disease_mapping = {}                # keep track of the total number of unique disease names parsed within the current variant line
        # ct -> disease -> None
        interim_ct_mapping = {}                     # keep track of the total number of unique disease names assigned to each "ct" class within the current variant line

        # Iterate across the 4 relevant evidence columns assumed to be present in the file
        # Keep track of all unique instances of disease_name + ct_type parsed across the available CIViC evidences
        for tmp_pos in [drug_pos, diag_pos, prog_pos, pred_pos]:
            # Empty evidence columns should be skipped as no information was found on CIViC
            has_evidence = False
            evidence_infos = str(line_split[tmp_pos].strip())
            if evidence_infos == ".":
                continue

            # NOTE: assume that independent evidence items are listed using ";" as a separator character
            # NOTE: assume that cancer-specificity classifications should always be available for all disease names reported, independently of the given evidence type column
            # E.g. CCND1:AMPLIFICATION:LUNG NON-SMALL CELL CARCINOMA|NCT(SUPPORTS,POOR OUTCOME(B(PUBMED_17070615:ACCEPTED:FULLY CURATED:SOMATIC:3)));..
            evidence_infos_list = evidence_infos.split(";")
            for evidence_info in evidence_infos_list:

                # NOTE: assume that evidence annotations retrieved from CIViC will never contain separator character "|"
                tmp_evidence_split = evidence_info.strip().split("|")

                # Extract relevant information about the gene, variant and disease name reported from CIViC
                variant_and_disease = tmp_evidence_split[0].strip()

                # Sanity check the expected separator character and format
                variant_and_disease_split = variant_and_disease.split(":")
                if len(variant_and_disease_split) < 3:
                    raise ValueError("Encountered unexpected format of CIViC evidence annotations in column %s of line %s" %(str(tmp_pos+1), line.strip()))

                # NOTE: assume gene and disease names from CIViC can never contain separator character ":"
                # Use uppercase to avoid mismatches of the gene and disease name due to case
                gene_name = variant_and_disease_split[0].strip().upper()
                disease_name = variant_and_disease_split[len(variant_and_disease_split)-1].strip().upper()

                # NOTE: take into account that some CIViC variant names can contain separator character ":"
                # E.g. "LMNA::NTRK1 E11-E10:18.5"
                variant_name = variant_and_disease_split[1].strip()
                if len(variant_and_disease_split) > 3:
                    variant_name = [":".join(variant_and_disease_split[1:-1])]

                # a) Special case of drug prediction evidences (column 'CIViC_PREDICTIVE')
                # In this situation, the evidence strings have a slightly different format to include the drug name
                # E.g. CCND1:AMPLIFICATION:BREAST CANCER|NCT|TAMOXIFEN(DOES NOT SUPPORT,RESISTANCE(B(PUBMED_24367492:ACCEPTED:FULLY CURATED:SOMATIC:3)))
                interim_evidence_split = None
                if tmp_pos == drug_pos:
                    # Sanity check the expected separator character and format
                    if len(tmp_evidence_split) != 3:
                        raise ValueError("Encountered unexpected format of CIViC evidence annotations in column %s of line %s" %(str(tmp_pos+1), line.strip()))
                    # Retrieve the ct" class assigned by CIViCutils and use lowercase to avoid mismatchings due to case
                    ct_type = tmp_evidence_split[1].strip().lower()
                    # Split the remaining evidence string further and retrieve the direction, clinical significance and drug name reported from CIViC
                    drug_and_evidence = tmp_evidence_split[2].strip()
                    (direction, clinical_signf, interim_drug) = get_clinical_info(drug_and_evidence, has_drug=True)
                    interim_evidence_split = drug_and_evidence

                # b) All other evidences
                # E.g. CCND1:AMPLIFICATION:LUNG NON-SMALL CELL CARCINOMA|NCT(SUPPORTS,POOR OUTCOME(B(PUBMED_17070615:ACCEPTED:FULLY CURATED:SOMATIC:3)))
                else:
                    # Sanity check the expected separator character and format
                    if len(tmp_evidence_split) != 2:
                        raise ValueError("Encountered unexpected format of CIViC evidence annotations in column %s of line %s" %(str(tmp_pos+1), line.strip()))
                    # Split the remaining evidence string further and retrieve the direction and clinical significance reported from CIViC, as well as the "ct" class assigned by CIViCutils
                    ct_and_evidence =  tmp_evidence_split[1].strip()
                    (direction, clinical_signf, ct_type) = get_clinical_info(ct_and_evidence, has_drug=False)
                    interim_evidence_split = ct_and_evidence

                # Sanity check for expected "ct" classifications provided
                if ct_type not in sorted_cts:
                    raise ValueError("Encountered unexpected cancer-specificity classification '%s' in column %s of line %s" %(ct_type, str(tmp_pos+1), line.strip()))
                # Sanity check that the required split string could be retrieved
                if interim_evidence_split is None:
                    raise ValueError("Encountered unexpected format of CIViC evidence annotations in column %s of line %s" %(str(tmp_pos+1), line.strip()))


                ## At this point, all relevant clinical information for the current variant line has been correctly parsed (tier, gene, variant, disease, ct, and optionally drug if evidence is 'Predictive')

                ## 1) Keep track of disease information across all variants and evidence columns for the current sample

                # Keep track of unique disease names parsed
                # disease -> ct
                if disease_name not in disease_mapping.keys():
                    disease_mapping[disease_name] = ct_type
                else:
                    # Sanity check that the same disease name can only be associated to one "ct" classification across the file
                    parsed_ct_type = disease_mapping[disease_name]
                    if ct_type != parsed_ct_type:
                        raise ValueError("Disease name '%s' was found to be associated to two different 'ct' classifications '%s' and '%s'!" %(disease_name, ct_type, parsed_ct_type))

                # Per tier, keep track of unique disease names parsed across all CIViC predictions
                # tier -> disease -> None
                if tier not in per_tier_disease_mapping.keys():
                    per_tier_disease_mapping[tier] = {}
                if disease_name not in per_tier_disease_mapping[tier].keys():
                    per_tier_disease_mapping[tier][disease_name] = None

                # Keep track of unique disease names associated with each 'ct' classification
                # ct -> disease -> None
                if ct_type not in ct_mapping.keys():
                    ct_mapping[ct_type] = {}
                if disease_name not in ct_mapping[ct_type].keys():
                    ct_mapping[ct_type][disease_name] = None

                # Per tier, keep track of unique disease names associated with each 'ct' classification
                # tier -> ct -> disease -> None
                if tier not in per_tier_ct_mapping.keys():
                    per_tier_ct_mapping[tier] = {}
                if ct_type not in per_tier_ct_mapping[tier].keys():
                    per_tier_ct_mapping[tier][ct_type] = {}
                if disease_name not in per_tier_ct_mapping[tier][ct_type].keys():
                    per_tier_ct_mapping[tier][ct_type][disease_name] = None


                ## Version of stats above excluding tier3 matches which can introduce biases
                # Apply following block only to variants which are tier1, tier1b or tier2 (exclude tier3 and tier4 cannot have disease info available)
                if tier != "tier_3":
                    # Keep track of unique disease names parsed
                    # disease -> ct
                    if disease_name not in disease_mapping_no_tier3.keys():
                        # Already checked above for 1-1 correspondances between disease name and a single "ct" classification
                        disease_mapping_no_tier3[disease_name] = ct_type

                    # Keep track of unique disease names associated with each 'ct' classification
                    # ct -> disease -> None
                    if ct_type not in ct_mapping_no_tier3.keys():
                        ct_mapping_no_tier3[ct_type] = {}
                    if disease_name not in ct_mapping_no_tier3[ct_type].keys():
                        ct_mapping_no_tier3[ct_type][disease_name] = None


                ## 2) Keep track of disease information only within the current variant line (to compute disease stats across all variants and tiers)

                # Keep track of unique disease names parsed
                # disease -> None
                if disease_name not in interim_disease_mapping.keys():
                    interim_disease_mapping[disease_name] = None

                # Keep track of unique disease names associated with each 'ct' classification
                # ct -> disease -> None
                if ct_type not in interim_ct_mapping.keys():
                    interim_ct_mapping[ct_type] = {}
                if disease_name not in interim_ct_mapping[ct_type].keys():
                    interim_ct_mapping[ct_type][disease_name] = None

                # Keep track of whether the currently evaluated column (i.e. evidence type) has at least 1 evidence item available
                # NOTE: at this point, all sanity checks for the expected format of the column and listed evidence items have already been performed (see above)
                has_evidence = True

            ## End of loop that iterates individual evidence items found within the same column

            # Check whether the currently evaluated evidence type had any evidence items listed
            if has_evidence:
                # Keep track of the number of variants having some info available for each evidence type
                if tmp_pos == drug_pos:
                    n_predictive += 1.0
                elif tmp_pos == diag_pos:
                    n_diagnostic += 1.0
                elif tmp_pos == prog_pos:
                    n_prognostic += 1.0
                elif tmp_pos == pred_pos:
                    n_predisposing += 1.0
                else:
                    raise ValueError("Encountered unexpected condition while parsing Evidence type annotations in line %s" %(line.strip()))

                # Also keep track of total number of variants with each evidence type available, excluding tier3 matches which can introduce biases
                if tier != "tier_3":
                    if tmp_pos == drug_pos:
                        n_predictive_no_tier3 += 1.0
                    elif tmp_pos == diag_pos:
                        n_diagnostic_no_tier3 += 1.0
                    elif tmp_pos == prog_pos:
                        n_prognostic_no_tier3 += 1.0
                    elif tmp_pos == pred_pos:
                        n_predisposing_no_tier3 += 1.0

        ## End of loop that iterates across all evidence type columns


# FIXME: Selected info can be used to split the evidence strings in turn
#                 # Use the first part of the string (now already known) to split and retrieve all evidence items supporting the current claim (i.e. different combinations of evidence level + publication id)
#                 if tmp_pos == drug_pos:
#                     # Format of interim_evidence_split (drug_and_evidence): 'DRUG(DIRECTION,CLINICALSIGNF(A(ref1,ref2,..),B(ref1,ref5)..))'
#                     split_pattern = interim_drug + "(" + direction + "," + clinical_signf + "("
#                 else:
#                     # Format of interim_evidence_split (ct_and_evidence): 'CT_TYPE(DIRECTION,CLINICALSIGNF(A(ref1,ref2,..),B(ref1,ref5)..))'
#                     split_pattern = ct_type + "(" + direction + "," + clinical_signf + "("
# 
#                 # Split returns list containing as many elements as evidence items supporting the current claim
#                 # NOTE: this is the case even when the same publication id is reported across several levels
#                 pub_ids_split = interim_evidence_split.split(split_pattern)[1].strip().split(",")

# FIXME: taken from 'https://gitlab.ethz.ch/nexus/thalmann_seiler_tcga_2017/-/blob/master/scripts/combine_all_civic_prediction.py'
#                     # In CIVIC, drugs can form part of a combinatorial treatment, indicated by 'DRUG1+DRUG2+...'
#                     drug_list = interim_drug.strip().split("+")
#                     for drug in drug_list:
# 
#                             # Keep track of how many different evidence items (level+PMIDs) support this direction+clinicalSignificance
#                             # Also, this way we keep track of supporting evidence items for this direction+clinicalSignificance across multiple variants for gene of interest in one sample
#                             for x in range(0, len(pub_ids_split)): # contains as many elements as evidence items support this claim (even when the same reference is used for many levels)
#                                 clinInfoDict[tmp_synonym][sampleName][gene][tier][cancerTag][direction].append(clinSignf)


        ## At this point, all evidence items available across the 4 evidence type columns have been parsed
        ## Process and keep track of disease information parsed for the current variant line

        # Keep track of number of unique matched diseases across all variants, also keep counts per tier
        tmp_n_diseases = float(len(interim_disease_mapping.keys()))
        # tier -> # diseases
        if tier not in per_tier_matched_diseases_mapping.keys():
            per_tier_matched_diseases_mapping[tier] = 0.0
        per_tier_matched_diseases_mapping[tier] += tmp_n_diseases
        matched_diseases += tmp_n_diseases

        # Keep track of number of unique matched diseases per "ct" class across all variants, also keep counts per tier
        # tier -> ct -> # diseases
        if tier not in per_tier_matched_diseases_ct_mapping.keys():
            per_tier_matched_diseases_ct_mapping[tier] = {}

        # Iterate all "ct" classes available for the current variant line and keep track of associated disease info
        # ct -> disease -> None
        for interim_ct_type in interim_ct_mapping.keys():
            interim_n_diseases = float(len(interim_ct_mapping[interim_ct_type].keys()))
            # ct -> # diseases
            if interim_ct_type not in matched_diseases_ct_mapping.keys():
                matched_diseases_ct_mapping[interim_ct_type] = 0.0
            matched_diseases_ct_mapping[interim_ct_type] += interim_n_diseases
             # tier -> ct -> # diseases
            if interim_ct_type not in per_tier_matched_diseases_ct_mapping[tier].keys():
                per_tier_matched_diseases_ct_mapping[tier][interim_ct_type] = 0.0
            per_tier_matched_diseases_ct_mapping[tier][interim_ct_type] += interim_n_diseases

            # Also keep track of number of unique matched diseases per "ct" class across all variants, excluding tier3 matches which can introduce biases
            if tier == "tier_3":
                continue
            # ct -> # diseases
            if interim_ct_type not in matched_diseases_ct_mapping_no_tier3.keys():
                matched_diseases_ct_mapping_no_tier3[interim_ct_type] = 0.0
            matched_diseases_ct_mapping_no_tier3[interim_ct_type] += interim_n_diseases


        ## 4) Process column listing consensus support across available drugs

## TODO: (double check info with the corresponding info parsed from Predictive column?)

        # Only process further variant lines which have consensus drug prediction information available in column 'CIViC_Drug_Support'
        drug_infos = str(line_split[drug_supp_pos].strip())
        if drug_infos == ".":
            continue

        # Assume multiple consensus support strings can be listed separated by ";"
        # E.g.: 'ENTRECTINIB:NCT:CIVIC_SUPPORT;LAROTRECTINIB:NCT:CIVIC_SUPPORT;..'
        drug_infos_list = drug_infos.split(";")
        if not drug_infos_list:
            raise ValueError("Encountered unexpected case of consensus drug support annotations in line %s" %(line.strip()))
        # Keep track of all variant lines having non-empty consensus drug support info
        n_drug_info += 1.0
        # Keep track of all variant lines, excluding those classified as tier3, having non-empty consensus drug support info
        if tier != "tier_3":
            n_drug_info_no_tier3 += 1.0

        # ct -> drug -> [consensus_support_1,..,consensus_support_N]
        interim_consensus_ct_mapping = {}

        # Iterate individual string of consensus drug support
        # Assume one string per combination of drug name + "ct" class (note the same drug can be available for different "ct" classes)
        for drug_info in drug_infos_list:
            # Format: 'DRUG_NAME:CT_TYPE:CIVIC_SUPPORT'
            drug_split = drug_info.strip().split(":")
            # Sanity check the expected separator character and format
            if len(drug_split) != 3:
                raise ValueError("Encountered unexpected format of consensus drug support annotations in line %s" %(line.strip()))
            # Use uppercase for drug names and consensus support strings to avoid mismatches due to case
            drug = drug_split[0].strip().upper()
            consensus_support = drug_split[2].strip()
            # Use lowercase for "ct" classes to avoid mismatches due to case
            ct_type = drug_split[1].strip().lower()

            # Keep track of all disease names, associated "ct" class, and available consensus drug support strings
            # ct -> drug -> [consensus_support_1,..,consensus_support_N]
            if ct_type not in interim_consensus_ct_mapping.keys():
                interim_consensus_ct_mapping[ct_type] = {}
            if drug not in interim_consensus_ct_mapping[ct_type].keys():
                interim_consensus_ct_mapping[ct_type][drug] = []
            # Sanity check for duplicated consensus support strings for the same disease + "ct" class
            if consensus_support in interim_consensus_ct_mapping[ct_type][drug]:
                print("Warning! Skipping duplicated consensus support '%s' for drug '%s' and cancer-specificity classification '%s' encountered in line %s" %(consensus_support, drug, ct_type, line.strip()))
            interim_consensus_ct_mapping[ct_type][drug].append(consensus_support)


        # Sanity check that at this point, all variant lines parsed should have at least one consensus drug prediction associated
        if not interim_consensus_ct_mapping:
            raise ValueError("Encountered unexpected case of consensus drug support annotations in line %s" %(line.strip()))

        # Sort all "ct" classes available for the current sample by the priority order defined at the beginning of this function (i.e. 'sorted_cts')
        sorted_ct_list = sorted(interim_consensus_ct_mapping.keys(), key=lambda x: sorted_cts.index(x))
        # Select the "ct" class with the highest priority for the current line (i.e. at least one will always be available, even if 'nct')
        pick_ct = sorted_ct_list[0]

        # Keep track of the total number of unique (consensus) drug names parsed for the current sample, also keep track per tier and "ct" class available
        # Both at the level of all consensus drug predictions available for the current line, as well as only those drug predictions associated to the "ct" class with highest priority available for the current line
        # ct -> drug -> consensus_support
        for tmp_ct in interim_consensus_ct_mapping.keys():
            for tmp_drug in interim_consensus_ct_mapping[tmp_ct].keys():
                # drug -> ct -> consensus_support
                if tmp_drug not in consensus_drug_mapping.keys():
                    consensus_drug_mapping[tmp_drug] = {}
                if tmp_ct not in consensus_drug_mapping[tmp_drug].keys():
                    consensus_drug_mapping[tmp_drug][tmp_ct] = []
                # tier -> drug -> None
                if tier not in per_tier_consensus_drug_mapping.keys():
                    per_tier_consensus_drug_mapping[tier] = {}
                if tmp_drug not in per_tier_consensus_drug_mapping[tier].keys():
                    per_tier_consensus_drug_mapping[tier][tmp_drug] = None
                # ct -> drug -> consensus_support
                if tmp_ct not in consensus_ct_mapping.keys():
                    consensus_ct_mapping[tmp_ct] = {}
                if tmp_drug not in consensus_ct_mapping[tmp_ct].keys():
                    consensus_ct_mapping[tmp_ct][tmp_drug] = []
                # tier -> ct -> drug -> None
                if tier not in per_tier_consensus_ct_mapping.keys():
                    per_tier_consensus_ct_mapping[tier] = {}
                if tmp_ct not in per_tier_consensus_ct_mapping[tier].keys():
                    per_tier_consensus_ct_mapping[tier][tmp_ct] = {}
                if tmp_drug not in per_tier_consensus_ct_mapping[tier][tmp_ct].keys():
                    per_tier_consensus_ct_mapping[tier][tmp_ct][tmp_drug] = None

                consensus_list = interim_consensus_ct_mapping[tmp_ct][tmp_drug]
                # NOTE: expectation is that each combination of disease name + "ct" class can only have one single consensus support string associated
                if len(consensus_list) > 1:
                    print("Warning! Encountered multiple consensus support strings ('%s') for drug '%s' and cancer-specificity classification '%s' in line %s" %(consensus_list, tmp_drug, tmp_ct, line.strip()))

                # Also keep track of total and per "ct" number of unique (consensus) drug names parsed for the current sample, excluding tier3 matches which can introduce biases
                if tier != "tier_3":
                    # drug -> ct -> consensus_support
                    if tmp_drug not in consensus_drug_mapping_no_tier3.keys():
                        consensus_drug_mapping_no_tier3[tmp_drug] = {}
                    if tmp_ct not in consensus_drug_mapping_no_tier3[tmp_drug].keys():
                        consensus_drug_mapping_no_tier3[tmp_drug][tmp_ct] = []
                    # ct -> drug -> consensus_support
                    if tmp_ct not in consensus_ct_mapping_no_tier3.keys():
                        consensus_ct_mapping_no_tier3[tmp_ct] = {}
                    if tmp_drug not in consensus_ct_mapping_no_tier3[tmp_ct].keys():
                        consensus_ct_mapping_no_tier3[tmp_ct][tmp_drug] = []


                # Keep track of drug information associated with the "ct" class of highest priority for the current variant line
                if tmp_ct != pick_ct:
                    continue

                # drug -> ct -> consensus_support
                if tmp_drug not in prior_consensus_drug_mapping.keys():
                    prior_consensus_drug_mapping[tmp_drug] = {}
                if tmp_ct not in prior_consensus_drug_mapping[tmp_drug].keys():
                    prior_consensus_drug_mapping[tmp_drug][tmp_ct] = []
                # tier -> drug -> None
                if tier not in per_tier_prior_consensus_drug_mapping.keys():
                    per_tier_prior_consensus_drug_mapping[tier] = {}
                if tmp_drug not in per_tier_prior_consensus_drug_mapping[tier].keys():
                    per_tier_prior_consensus_drug_mapping[tier][tmp_drug] = None
                # ct -> drug -> consensus_support
                if tmp_ct not in prior_consensus_ct_mapping.keys():
                    prior_consensus_ct_mapping[tmp_ct] = {}
                if tmp_drug not in prior_consensus_ct_mapping[tmp_ct].keys():
                    prior_consensus_ct_mapping[tmp_ct][tmp_drug] = []
                # tier -> ct -> drug -> None
                if tier not in per_tier_prior_consensus_ct_mapping.keys():
                    per_tier_prior_consensus_ct_mapping[tier] = {}
                if tmp_ct not in per_tier_prior_consensus_ct_mapping[tier].keys():
                    per_tier_prior_consensus_ct_mapping[tier][tmp_ct] = {}
                if tmp_drug not in per_tier_prior_consensus_ct_mapping[tier][tmp_ct].keys():
                    per_tier_prior_consensus_ct_mapping[tier][tmp_ct][tmp_drug] = None

                # Also keep track of drug information associated with the "ct" class of highest priority for the current variant line, excluding tier3 matches which can introduce biases
                if tier != "tier_3":
                    # drug -> ct -> consensus_support
                    if tmp_drug not in prior_consensus_drug_mapping_no_tier3.keys():
                        prior_consensus_drug_mapping_no_tier3[tmp_drug] = {}
                    if tmp_ct not in prior_consensus_drug_mapping_no_tier3[tmp_drug].keys():
                        prior_consensus_drug_mapping_no_tier3[tmp_drug][tmp_ct] = []
                    # ct -> drug -> consensus_support
                    if tmp_ct not in prior_consensus_ct_mapping_no_tier3.keys():
                        prior_consensus_ct_mapping_no_tier3[tmp_ct] = {}
                    if tmp_drug not in prior_consensus_ct_mapping_no_tier3[tmp_ct].keys():
                        prior_consensus_ct_mapping_no_tier3[tmp_ct][tmp_drug] = []


## TODO: check in the consensus drug support column: report, out of all unique drugs, how many (%) have associated: a. civic_support, civic_resistance, etc. Also, take into account situations where the same drug is associated to different classifications depending on the ct, variant, etc. -> how to handle?

    infile.close()


    ## A) Stats on mean number of matched variants

    # Keep track of the total number of variants matched overall for the current sample, excluding tier3 matches which can introduce biases
    matched_variants_no_tier3 = 0.0
    # tier -> # variants
    for tmp_tier in matched_variants_mapping.keys():
        if tmp_tier == "tier_3":
            continue
        matched_variants_no_tier3 += float(matched_variants_mapping[tmp_tier])

    # Compute mean number of matched variants for the sample (only makes sense to compute mean on lines that had CIViC matches available)
    mean_matched_variants = 0.0
    # Sanity check for divisions by 0
    if all_civic_variants:
        mean_matched_variants = float(float(matched_variants) / float(all_civic_variants))

    # Also compute mean number of matched variants excluding tier3 matches (can introduce biases as this is an unspecific match)
    mean_matched_variants_no_tier3 = 0.0
    # Get total number of parsed variants, excluding those classified with a tier3 (unspecific match)
    n_civic_variants_no_tier3 = n_tier_1 + n_tier_1b + n_tier_2
    if n_civic_variants_no_tier3:
        mean_matched_variants_no_tier3 = float(float(matched_variants_no_tier3) / float(n_civic_variants_no_tier3))

    # Per tier, compute mean number of matched variants for the sample
    # tier -> # variants
    (mean_matched_variants_tier1, mean_matched_variants_tier1b, mean_matched_variants_tier2, mean_matched_variants_tier3) = process_mean_feature_per_tier(matched_variants_mapping, n_tier_1, n_tier_1b, n_tier_2, n_tier_3)



    ## B) Stats on mean number of matched diseases

    # Keep track of the total number of unique disease names matched overall across all variants for the current sample (excluding tier3 matches which can introduce biases)
    matched_diseases_no_tier3 = 0.0
    # tier -> # diseases
    for tmp_tier in per_tier_matched_diseases_mapping.keys():
        if tmp_tier == "tier_3":
            continue
        matched_diseases_no_tier3 += float(per_tier_matched_diseases_mapping[tmp_tier])

    # Compute mean number of matched diseases for the sample (only makes sense to compute mean on lines that had CIViC matches available)
    mean_matched_diseases = 0.0
    # Sanity check for divisions by 0
    if all_civic_variants:
        mean_matched_diseases = float(float(matched_diseases) / float(all_civic_variants))

    # Also compute mean number of matched diseases excluding tier3 matches (can introduce biases as this is an unspecific match)
    mean_matched_diseases_no_tier3 = 0.0
    # Base mean on total number of parsed variants with tier1, tier1b or tier2 (computed above)
    if n_civic_variants_no_tier3:
        mean_matched_diseases_no_tier3 = float(float(matched_diseases_no_tier3) / float(n_civic_variants_no_tier3))

    # Per tier, compute mean number of matched diseases for the sample
    # tier -> # diseases
    (mean_matched_diseases_tier1, mean_matched_diseases_tier1b, mean_matched_diseases_tier2, mean_matched_diseases_tier3) = process_mean_feature_per_tier(per_tier_matched_diseases_mapping, n_tier_1, n_tier_1b, n_tier_2, n_tier_3)

    # Per ct, compute mean number of matched diseases across all variants for the current sample
    # ct -> # diseases
    (mean_matched_diseases_ct, mean_matched_diseases_gt, mean_matched_diseases_nct) = process_mean_feature_per_ct(matched_diseases_ct_mapping, all_civic_variants)
    # Also compute version of stats per ct, compute mean number of matched diseases across all variants for the current sample
# excluding tier3 matches which can introduce biases
    # ct -> # diseases
    (mean_matched_diseases_ct_no_tier3, mean_matched_diseases_gt_no_tier3, mean_matched_diseases_nct_no_tier3) = process_mean_feature_per_ct(matched_diseases_ct_mapping_no_tier3, n_civic_variants_no_tier3)

    # Per tier, compute mean number of matched diseases per "ct" class for the current sample
    # tier -> ct -> # diseases
    (mean_matched_diseases_tier1_ct, mean_matched_diseases_tier1_gt, mean_matched_diseases_tier1_nct, mean_matched_diseases_tier1b_ct, mean_matched_diseases_tier1b_gt, mean_matched_diseases_tier1b_nct, mean_matched_diseases_tier2_ct, mean_matched_diseases_tier2_gt, mean_matched_diseases_tier2_nct, mean_matched_diseases_tier3_ct, mean_matched_diseases_tier3_gt, mean_matched_diseases_tier3_nct) = process_mean_feature_per_tier_and_ct(per_tier_matched_diseases_ct_mapping, n_tier_1, n_tier_1b, n_tier_2, n_tier_3)

    # Per tier, compute number of unique disease names per "ct" class for the current sample
    # tier -> ct -> disease -> None
    (n_diseases_tier1_ct, n_diseases_tier1_gt, n_diseases_tier1_nct, n_diseases_tier1b_ct, n_diseases_tier1b_gt, n_diseases_tier1b_nct, n_diseases_tier2_ct, n_diseases_tier2_gt, n_diseases_tier2_nct, n_diseases_tier3_ct, n_diseases_tier3_gt, n_diseases_tier3_nct) = process_feature_per_tier_and_ct(per_tier_ct_mapping)


    # Get total number of unique disease names parsed across the current sample
    # disease -> ct
    n_diseases = len(disease_mapping.keys())

    # Per tier, get total number of unique disease names parsed for the sample
    # tier -> disease -> None
    (n_diseases_tier1, n_diseases_tier1b, n_diseases_tier2, n_diseases_tier3) = process_feature_per_tier(per_tier_disease_mapping)

    # Per ct, get total number of unique disease names parsed for the sample
    # ct -> disease -> None
    (n_diseases_ct, n_diseases_gt, n_diseases_nct) = process_feature_per_ct(ct_mapping)

    # Per tier, get total number of unique disease names parsed for each "ct" class for the current sample (only makes sense to compute mean for tiers != 4)
    # tier -> ct -> disease -> None
    (n_diseases_tier1_ct, n_diseases_tier1_gt, n_diseases_tier1_nct, n_diseases_tier1b_ct, n_diseases_tier1b_gt, n_diseases_tier1b_nct, n_diseases_tier2_ct, n_diseases_tier2_gt, n_diseases_tier2_nct, n_diseases_tier3_ct, n_diseases_tier3_gt, n_diseases_tier3_nct) = process_feature_per_tier_and_ct(per_tier_ct_mapping)


    ## Version of stats above excluding tier3 matches which can introduce biases

    # Get total number of unique disease names parsed across the current sample excluding tier3 matches
    # disease -> ct
    n_diseases_no_tier3 = len(disease_mapping_no_tier3.keys())
    # Per ct, get total number of unique disease names parsed for the sample excluding tier3 matches
    # ct -> disease -> None
    (n_diseases_ct_no_tier3, n_diseases_gt_no_tier3, n_diseases_nct_no_tier3) = process_feature_per_ct(ct_mapping_no_tier3)



    ## a) Version of drug stats based on all consensus drug predictions found for the sample

    # Get total number of unique drug names parsed across the current sample
    # drug -> ct -> consensus_support
    n_drugs_all = len(consensus_drug_mapping.keys())

    # Compute mean number of "ct" classes available per drug parsed across the current sample
    interim_ct_sum = 0.0
    # Per drug, retrieve number of "ct" classes parsed for the current sample (important to use version of dict without prioritization of "ct" class performed!)
    for this_drug in consensus_drug_mapping.keys():
        n_ct_classes_avail = float(len(consensus_drug_mapping[this_drug].keys()))
        interim_ct_sum += n_ct_classes_avail
    mean_ct_classes_avail = 0.0
    # Sanity check for divisions by 0
    if n_drugs_all:
        mean_ct_classes_avail = float(float(interim_ct_sum) / float(n_drugs_all))

    # Per tier, get total number of unique drug names parsed for the sample
    # tier -> drug -> None
    (n_drugs_all_tier1, n_drugs_all_tier1b, n_drugs_all_tier2, n_drugs_all_tier3) = process_feature_per_tier(per_tier_consensus_drug_mapping)

    # Per ct, get total number of unique drug names parsed for the sample
    # ct -> drug -> consensus_support
    (n_drugs_all_ct, n_drugs_all_gt, n_drugs_all_nct) = process_feature_per_ct(consensus_ct_mapping)


    ## Version of stats above excluding tier3 matches which can introduce biases

    # Get total number of unique drug names parsed across the current sample excluding tier3 matches
    # drug -> ct -> consensus_support
    n_drugs_all_no_tier3 = len(consensus_drug_mapping_no_tier3.keys())
    # Per ct, get total number of unique drug names parsed for the sample excluding tier3 matches
    # ct -> drug -> consensus_support
    (n_drugs_all_ct_no_tier3, n_drugs_all_gt_no_tier3, n_drugs_all_nct_no_tier3) = process_feature_per_ct(consensus_ct_mapping_no_tier3)
    # Compute mean number of "ct" classes available per drug parsed across the current sample, excluding tier3 matches which can introduce biases
    interim_ct_sum_no_tier3 = 0.0
    # Per drug, retrieve number of "ct" classes parsed for the current sample (important to use version of dict without prioritization of "ct" class performed, and already excluding info from tier3 matches!)
    for this_drug_no_tier3 in consensus_drug_mapping_no_tier3.keys():
        n_ct_classes_avail_no_tier3 = float(len(consensus_drug_mapping_no_tier3[this_drug_no_tier3].keys()))
        interim_ct_sum_no_tier3 += n_ct_classes_avail_no_tier3
    mean_ct_classes_avail_no_tier3 = 0.0
    # Sanity check for divisions by 0
    if n_drugs_all_no_tier3:
        mean_ct_classes_avail_no_tier3 = float(float(interim_ct_sum_no_tier3) / float(n_drugs_all_no_tier3))



    ## b) Version of drug stats based only on the consensus drug predictions for the "ct" class with highest priority available per variant line

    # Get total number of unique drug names parsed for the highest "ct" class available across the current sample
    # drug -> ct -> consensus_support
    n_drugs_prior = len(prior_consensus_drug_mapping.keys())

    # Per tier, get total number of unique drug names parsed for the highest "ct" class available across the current sample
    # tier -> drug -> None
    (n_drugs_prior_tier1, n_drugs_prior_tier1b, n_drugs_prior_tier2, n_drugs_prior_tier3) = process_feature_per_tier(per_tier_prior_consensus_drug_mapping)

    # Per ct, get total number of unique drug names parsed for the highest "ct" class available across the current sample
    # ct -> drug -> consensus_support
    (n_drugs_prior_ct, n_drugs_prior_gt, n_drugs_prior_nct) = process_feature_per_ct(prior_consensus_ct_mapping)

    # Per tier, get total number of unique drug names parsed for each "ct" class for the current sample
    # tier -> ct -> drug -> None
    (n_drugs_tier1_ct, n_drugs_tier1_gt, n_drugs_tier1_nct, n_drugs_tier1b_ct, n_drugs_tier1b_gt, n_drugs_tier1b_nct, n_drugs_tier2_ct, n_drugs_tier2_gt, n_drugs_tier2_nct, n_drugs_tier3_ct, n_drugs_tier3_gt, n_drugs_tier3_nct) = process_feature_per_tier_and_ct(per_tier_consensus_ct_mapping)

    # Per tier and ct, get total number of unique drug names parsed for the highest "ct" class available across the current sample
    # tier -> ct -> drug -> None
    (n_drugs_tier1_ct_prior, n_drugs_tier1_gt_prior, n_drugs_tier1_nct_prior, n_drugs_tier1b_ct_prior, n_drugs_tier1b_gt_prior, n_drugs_tier1b_nct_prior, n_drugs_tier2_ct_prior, n_drugs_tier2_gt_prior, n_drugs_tier2_nct_prior, n_drugs_tier3_ct_prior, n_drugs_tier3_gt_prior, n_drugs_tier3_nct_prior) = process_feature_per_tier_and_ct(per_tier_prior_consensus_ct_mapping)


    ## Version of stats above excluding tier3 matches which can introduce biases

    # Get total number of unique drug names parsed for the highest "ct" class available across the current sample, excluding tier3 matches
    # drug -> ct -> consensus_support
    n_drugs_prior_no_tier3 = len(prior_consensus_drug_mapping_no_tier3.keys())
    # Per ct, get total number of unique drug names parsed for the highest "ct" class available across the current sample, excluding tier3 matches
    # ct -> drug -> consensus_support
    (n_drugs_prior_ct_no_tier3, n_drugs_prior_gt_no_tier3, n_drugs_prior_nct_no_tier3) = process_feature_per_ct(prior_consensus_ct_mapping_no_tier3)


    ## Keep track of all computed info and stats per sample across the entire patient cohort being processed

    # sample -> [#vars, #civic, #civic_no_tier3, #tier1, #tier1b, #tier1agg, #tier2, #tier3, #tier4, #predictive, #diagnostic, #prognostic, #predisposing, #predictive_no_tier3, #diagnostic_no_tier3, #prognostic_no_tier3, #predisposing_no_tier3, mean_matched_vars, mean_matched_vars_no_tier3, mean_matched_vars_tier1, mean_matched_vars_tier1b, mean_matched_vars_tier2, mean_matched_vars_tier3, mean_matched_diseases, mean_matched_diseases_no_tier3, mean_matched_diseases_tier1, mean_matched_diseases_tier1b, mean_matched_diseases_tier2, mean_matched_diseases_tier3, mean_matched_diseases_ct, mean_matched_diseases_gt, mean_matched_diseases_nct, mean_matched_diseases_ct_no_tier3, mean_matched_diseases_gt_no_tier3, mean_matched_diseases_nct_no_tier3, mean_matched_diseases_tier1_ct, mean_matched_diseases_tier1_gt, mean_matched_diseases_tier1_nct, mean_matched_diseases_tier1b_ct, mean_matched_diseases_tier1b_gt, mean_matched_diseases_tier1b_nct, mean_matched_diseases_tier2_ct, mean_matched_diseases_tier2_gt, mean_matched_diseases_tier2_nct, mean_matched_diseases_tier3_ct, mean_matched_diseases_tier3_gt, mean_matched_diseases_tier3_nct, #diseases, #diseases_no_tier3, #diseases_tier1, #diseases_tier1b, #diseases_tier2, #diseases_tier3, #diseases_ct, #diseases_gt, #diseases_nct, #diseases_ct_no_tier3, #diseases_gt_no_tier3, #diseases_nct_no_tier3, #diseases_tier1_ct, #diseases_tier1_gt, #diseases_tier1_nct, #diseases_tier1b_ct, #diseases_tier1b_gt, #diseases_tier1b_nct, #diseases_tier2_ct, #diseases_tier2_gt, #diseases_tier2_nct, #diseases_tier3_ct, #diseases_tier3_gt, #diseases_tier3_nct, #drugs, #drugs_no_tier3, mean_cts_per_drug, mean_cts_per_drug_no_tier3, #drugs_tier1, #drugs_tier1b, #drugs_tier2, #drugs_tier3, #drugs_ct, #drugs_gt, #drugs_nct, #drugs_ct_no_tier3, #drugs_gt_no_tier3, #drugs_nct_no_tier3, #drugs_prior, #drugs_prior_no_tier3, #drugs_tier1_prior, #drugs_tier1b_prior, #drugs_tier2_prior, #drugs_tier3_prior, #drugs_ct_prior, #drugs_gt_prior, #drugs_nct_prior, #drugs_ct_prior_no_tier3, #drugs_gt_prior_no_tier3, #drugs_nct_prior_no_tier3, #drugs_tier1_ct, #drugs_tier1_gt, #drugs_tier1_nct, #drugs_tier1b_ct, #drugs_tier1b_gt, #drugs_tier1b_nct, #drugs_tier2_ct, #drugs_tier2_gt, #drugs_tier2_nct, #drugs_tier3_ct, #drugs_tier3_gt, #drugs_tier3_nct, drugs_tier1_ct_prior, #drugs_tier1_gt_prior, #drugs_tier1_nct_prior,, #drugs_tier1b_ct_prior, #drugs_tier1b_gt_prior, #drugs_tier1b_nct_prior, #drugs_tier2_ct_prior, #drugs_tier2_gt_prior, #drugs_tier2_nct_prior, #drugs_tier3_ct_prior, #drugs_tier3_gt_prior, #drugs_tier3_nct_prior]
    if sample_name in civic_info_mapping.keys():
        raise ValueError("Sample name '%s' was already parsed!")

    civic_info_mapping[sample_name] = [all_variants, all_civic_variants, n_civic_variants_no_tier3, n_tier_1, n_tier_1b, n_tier_1_agg, n_tier_2, n_tier_3, n_tier_4, n_predictive, n_diagnostic, n_prognostic, n_predisposing, n_predictive_no_tier3, n_diagnostic_no_tier3, n_prognostic_no_tier3, n_predisposing_no_tier3, mean_matched_variants, mean_matched_variants_no_tier3, mean_matched_variants_tier1, mean_matched_variants_tier1b, mean_matched_variants_tier2, mean_matched_variants_tier3, mean_matched_diseases, mean_matched_diseases_no_tier3, mean_matched_diseases_tier1, mean_matched_diseases_tier1b, mean_matched_diseases_tier2, mean_matched_diseases_tier3, mean_matched_diseases_ct, mean_matched_diseases_gt, mean_matched_diseases_nct, mean_matched_diseases_ct_no_tier3, mean_matched_diseases_gt_no_tier3, mean_matched_diseases_nct_no_tier3, mean_matched_diseases_tier1_ct, mean_matched_diseases_tier1_gt, mean_matched_diseases_tier1_nct, mean_matched_diseases_tier1b_ct, mean_matched_diseases_tier1b_gt, mean_matched_diseases_tier1b_nct, mean_matched_diseases_tier2_ct, mean_matched_diseases_tier2_gt, mean_matched_diseases_tier2_nct, mean_matched_diseases_tier3_ct, mean_matched_diseases_tier3_gt, mean_matched_diseases_tier3_nct, n_diseases, n_diseases_no_tier3, n_diseases_tier1, n_diseases_tier1b, n_diseases_tier2, n_diseases_tier3, n_diseases_ct, n_diseases_gt, n_diseases_nct, n_diseases_ct_no_tier3, n_diseases_gt_no_tier3, n_diseases_nct_no_tier3, n_diseases_tier1_ct, n_diseases_tier1_gt, n_diseases_tier1_nct, n_diseases_tier1b_ct, n_diseases_tier1b_gt, n_diseases_tier1b_nct, n_diseases_tier2_ct, n_diseases_tier2_gt, n_diseases_tier2_nct, n_diseases_tier3_ct, n_diseases_tier3_gt, n_diseases_tier3_nct, n_drugs_all, n_drugs_all_no_tier3, mean_ct_classes_avail, mean_ct_classes_avail_no_tier3, n_drugs_all_tier1, n_drugs_all_tier1b, n_drugs_all_tier2, n_drugs_all_tier3, n_drugs_all_ct, n_drugs_all_gt, n_drugs_all_nct, n_drugs_all_ct_no_tier3, n_drugs_all_gt_no_tier3, n_drugs_all_nct_no_tier3, n_drugs_prior, n_drugs_prior_no_tier3, n_drugs_prior_tier1, n_drugs_prior_tier1b, n_drugs_prior_tier2, n_drugs_prior_tier3, n_drugs_prior_ct, n_drugs_prior_gt, n_drugs_prior_nct, n_drugs_prior_ct_no_tier3, n_drugs_prior_gt_no_tier3, n_drugs_prior_nct_no_tier3, n_drugs_tier1_ct, n_drugs_tier1_gt, n_drugs_tier1_nct, n_drugs_tier1b_ct, n_drugs_tier1b_gt, n_drugs_tier1b_nct, n_drugs_tier2_ct, n_drugs_tier2_gt, n_drugs_tier2_nct, n_drugs_tier3_ct, n_drugs_tier3_gt, n_drugs_tier3_nct, n_drugs_tier1_ct_prior, n_drugs_tier1_gt_prior, n_drugs_tier1_nct_prior, n_drugs_tier1b_ct_prior, n_drugs_tier1b_gt_prior, n_drugs_tier1b_nct_prior, n_drugs_tier2_ct_prior, n_drugs_tier2_gt_prior, n_drugs_tier2_nct_prior, n_drugs_tier3_ct_prior, n_drugs_tier3_gt_prior, n_drugs_tier3_nct_prior]
    return civic_info_mapping


# Write to output a series of statistics retrieved from parsing and processing CIViCutils annotations for a set of samples
def write_results_to_output(sample_order, input_mapping, outfile):
    # Iterate sample names to be reported to output in the provided order
    for sample in sample_order:
        # sample -> [infos1,..,infosN]
        if sample not in input_mapping.keys():
            raise ValueError("Provided sample name '%s' was not parsed!")
        # Retrieve and sanity check available CIViCutils info for the current sample
        civic_infos = input_mapping[sample]
        if len(civic_infos) != 121:
            raise ValueError("Expected 121 stat values from processing CIViC annotations for sample '%s'!" %(sample))
        # Reported numeric values must be converted into strings before writing to output
        civic_infos_strings = [str(round(x, 2)) for x in civic_infos]
        # Write each sample in a separate line
        outfile.write("%s\t%s\n" %(sample, "\t".join(civic_infos_strings)))
    return None



'''
Script
'''

parser = argparse.ArgumentParser(description='Parse and process CIViCutils annotations reported for SNVs and CNVs of a set of samples.')
parser.add_argument('--input_dir_civic_snv', dest='input_dir_civic_snv', required=True, help='Input directory with SNV CIViC files of format [sample].civic_snv.txt.')
parser.add_argument('--input_dir_civic_cnv', dest='input_dir_civic_cnv', required=True, help='Input directory with CNV CIViC files of format [sample].civic_cnv.txt.')
parser.add_argument('--file_suffix_snv', dest='file_suffix_snv', required=True, help='To retrieve correct input files, specify the desired file ending ("civic_snv.txt" for SNV data).')
parser.add_argument('--file_suffix_cnv', dest='file_suffix_cnv', required=True, help='To retrieve correct input files, specify the desired file ending ("civic_cnv.txt" for CNV data).')
parser.add_argument('--outfile_tag', dest='outfile_tag', required=True, help='Name prefix of the output files.')

args = parser.parse_args()


## 1) Process input files for SNVs annotated with CIViCutils

civic_info_mapping_snv = {}
seen_samples_snvs = []

for file in os.listdir(args.input_dir_civic_snv):
    sample_file_snv = "%s%s" %(args.input_dir_civic_snv,os.path.basename(file))
    if os.path.isfile(sample_file_snv) and sample_file_snv.endswith(args.file_suffix_snv):
        name_snv = os.path.basename(file).split("_")[0].split("-")[1]
        if name_snv in seen_samples_snvs:
            raise ValueError("Sample name '%s' was already parsed for SNV results!")
            sys.exit(1)
        seen_samples_snvs.append(name_snv)
        civic_info_mapping_snv = parse_input_file(sample_file_snv, name_snv, civic_info_mapping_snv)


## 2) Process input files for CNVs annotated with CIViCutils

civic_info_mapping_cnv = {}
seen_samples_cnvs = []

for file in os.listdir(args.input_dir_civic_cnv):
    sample_file_cnv = "%s%s" %(args.input_dir_civic_cnv,os.path.basename(file))
    if os.path.isfile(sample_file_cnv) and sample_file_cnv.endswith(args.file_suffix_cnv):
        name_cnv = os.path.basename(file).split("_")[0].split("-")[1]
        if name_cnv in seen_samples_cnvs:
            raise ValueError("Sample name '%s' was already parsed for CNV results!")
        seen_samples_cnvs.append(name_cnv)
        civic_info_mapping_cnv = parse_input_file(sample_file_cnv, name_cnv, civic_info_mapping_cnv)


# Sanity check that set of samples provided for SNVs and CNVs are identical (n=412 for TCGA-BLCA)
if (set(seen_samples_snvs) != set(seen_samples_cnvs)):
    raise ValueError("Sample names in provided SNV and CNV input directories do not match!")

# Write stats of CIViCutils annotations for all samples separately for SNVs and CNVs
outfile_snv = open(args.outfile_tag + ".snvs.tsv",'w') # Results from processing SNV annotations from CIViCutils
outfile_cnv = open(args.outfile_tag + ".cnvs.tsv",'w') # Results from processing CNV annotations from CIViCutils

# Header is identical for both output tables
output_header = "sample_name\tall_variants\tall_civic_variants\tall_civic_variants_no_tier3\tn_tier_1\tn_tier_1b\tn_tier_1_agg\tn_tier_2\tn_tier_3\tn_tier_4\tn_predictive_vars\tn_diagnostic_vars\tn_prognostic_vars\tn_predisposing_vars\tn_predictive_vars_no_tier3\tn_diagnostic_vars_no_tier3\tn_prognostic_vars_no_tier3\tn_predisposing_vars_no_tier3\tmean_matched_vars\tmean_matched_vars_no_tier3\tmean_matched_vars_tier1\tmean_matched_vars_tier1b\tmean_matched_vars_tier2\tmean_matched_vars_tier3\tmean_matched_diseases\tmean_matched_diseases_no_tier3\tmean_matched_diseases_tier1\tmean_matched_diseases_tier1b\tmean_matched_diseases_tier2\tmean_matched_diseases_tier3\tmean_matched_diseases_ct\tmean_matched_diseases_gt\tmean_matched_diseases_nct\tmean_matched_diseases_ct_no_tier3\tmean_matched_diseases_gt_no_tier3\tmean_matched_diseases_nct_no_tier3\tmean_matched_diseases_tier1_ct\tmean_matched_diseases_tier1_gt\tmean_matched_diseases_tier1_nct\tmean_matched_diseases_tier1b_ct\tmean_matched_diseases_tier1b_gt\tmean_matched_diseases_tier1b_nct\tmean_matched_diseases_tier2_ct\tmean_matched_diseases_tier2_gt\tmean_matched_diseases_tier2_nct\tmean_matched_diseases_tier3_ct\tmean_matched_diseases_tier3_gt\tmean_matched_diseases_tier3_nct\tn_diseases\tn_diseases_no_tier3\tn_diseases_tier1\tn_diseases_tier1b\tn_diseases_tier2\tn_diseases_tier3\tn_diseases_ct\tn_diseases_gt\tn_diseases_nct\tn_diseases_ct_no_tier3\tn_diseases_gt_no_tier3\tn_diseases_nct_no_tier3\tn_diseases_tier1_ct\tn_diseases_tier1_gt\tn_diseases_tier1_nct\tn_diseases_tier1b_ct\tn_diseases_tier1b_gt\tn_diseases_tier1b_nct\tn_diseases_tier2_ct\tn_diseases_tier2_gt\tn_diseases_tier2_nct\tn_diseases_tier3_ct\tn_diseases_tier3_gt\tn_diseases_tier3_nct\tn_drugs_all\tn_drugs_all_no_tier3\tmean_ct_classes_avail\tmean_ct_classes_avail_no_tier3\tn_drugs_all_tier1\tn_drugs_all_tier1b\tn_drugs_all_tier2\tn_drugs_all_tier3\tn_drugs_all_ct\tn_drugs_all_gt\tn_drugs_all_nct\tn_drugs_all_ct_no_tier3\tn_drugs_all_gt_no_tier3\tn_drugs_all_nct_no_tier3\tn_drugs_prior\tn_drugs_prior_no_tier3\tn_drugs_prior_tier1\tn_drugs_prior_tier1b\tn_drugs_prior_tier2\tn_drugs_prior_tier3\tn_drugs_prior_ct\tn_drugs_prior_gt\tn_drugs_prior_nct\tn_drugs_prior_ct_no_tier3\tn_drugs_prior_gt_no_tier3\tn_drugs_prior_nct_no_tier3\tn_drugs_tier1_ct\tn_drugs_tier1_gt\tn_drugs_tier1_nct\tn_drugs_tier1b_ct\tn_drugs_tier1b_gt\tn_drugs_tier1b_nct\tn_drugs_tier2_ct\tn_drugs_tier2_gt\tn_drugs_tier2_nct\tn_drugs_tier3_ct\tn_drugs_tier3_gt\tn_drugs_tier3_nct\tn_drugs_tier1_ct_prior\tn_drugs_tier1_gt_prior\tn_drugs_tier1_nct_prior\tn_drugs_tier1b_ct_prior\tn_drugs_tier1b_gt_prior\tn_drugs_tier1b_nct_prior\tn_drugs_tier2_ct_prior\tn_drugs_tier2_gt_prior\tn_drugs_tier2_nct_prior\tn_drugs_tier3_ct_prior\tn_drugs_tier3_gt_prior\tn_drugs_tier3_nct_prior"

outfile_snv.write(output_header + "\n")
outfile_cnv.write(output_header + "\n")

# Use same order of the reported samples in both output tables
sorted_samples = sorted(set(seen_samples_snvs))

write_results_to_output(sorted_samples, civic_info_mapping_snv, outfile_snv)
write_results_to_output(sorted_samples, civic_info_mapping_cnv, outfile_cnv)

outfile_snv.close()
outfile_cnv.close()
