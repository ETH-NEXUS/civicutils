#!/usr/bin/env python

'''
Process and combine variant info from CIViCutils
Lourdes Rosano, Feb 2022
'''

import sys
import os
import argparse
import re
import copy


# Define mapping of special cases where known drugs are referred to in CIViC using synonyms or special terms
# CIViC_synonym -> drug_name
drug_synonyms_mapping = {'DOVITINIB DILACTIC ACID (TKI258 DILACTIC ACID)':'DOVITINIB', '5-FLUOROURACIL':'FLUOROURACIL', '5-FU':'FLUOROURACIL', 'ADO-TRASTUZUMAB EMTANSINE':'TRASTUZUMAB EMTANSINE', 'PD0325901':'PD-0325901', 'PD173074':'PD-173074', 'BGJ-398':'INFIGRATINIB', 'BGJ398':'INFIGRATINIB'}

## Dictionary that allows prioritization of CIViC support categories (in oncoprint) when >1 gene has CIViC info for a given drug+sample
supportPrior = {'civic_support':1, 'civic_resistance':1, 'civic_conflict':2, 'civic_unknown':3, 'civic_unspecific_vars':4, 'dgidb_only':5}



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


def get_clinical_info(evidence_string, has_drug=False):
    # Sanity check the expected separator character and format
    # Assume direction and clinical significance are always separated by ","
    evidence_split = evidence_string.split(",")

    # NOTE: take into account CIViCutils can report several evidence levels and/or publication ids aggregated using separator character ","
    # E.g.'..(B(PUBMED_17590872:..,PUBMED_19903786:..),C(PUBMED_19223544:..))..'
    if len(evidence_split) > 2:
        # Remove entries consisting only of PMIDs to have clinical significance in the last possible position
        # NOTE: assumes that citation ids are always a series of digits and that evidence levels always correspond to one particular letter

# FIXME: review regex below to ensure it will always work for the new format of CIViC evidences

        evidence_split = [x for x in evidence_split if not re.match('([A-Z]\()?[A-Z]+_\d+(\)+)?', x)]

#                         new_evidence_split = []
#                         for x in interim_drug_and_evidence_split:
#                             if not re.match('([A-Z]\()?\d+(\)+)?', x):
#                                 new_evidence_split.append(x)
#                         for interim in new_evidence_split:
#                             drug_and_evidence_split.append(interim)
#                     else:
#                         for x in interim_drug_and_evidence_split:
#                             drug_and_evidence_split.append(x)

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

    return (direction, clinical_signf, cond_value, evidence_split)



# FIXME
# Special case Tier='SNV_only' should also be skipped. This only occurs for CNV CIViC data and corresponds to a special tier=3 case
# Also skip line if Tier=1,2 or 3, but still no drug information was available on CIViC for the current gene (eg. prognostic info only)

# Parse and process CIViC annotations available for the variants in the provided input file (corresponding to a given sample name)
def parse_input_file(sample_file, sample_name, civic_info_mapping):
# def parse_input_file(sample_file,sample_name,civic_info_mapping,drugsToVariants):
    print("Sample %s. File: %s" %(sample_name, sample_file))

    sorted_cts = ["ct","gt","nct"]

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

# FIXME
#     civicLines = 0
#     parsedDrugs = []
#     parsedGenes = []


    ## Keep track of several different annotations reported by CIViCutils

    all_variants = 0                            # all lines
    all_civic_variants = 0                      # all lines having CIViC info available (tier != 4)
    n_tier_1 = 0                                # all lines with tier = 1
    n_tier_1b = 0                               # all lines with tier = 1b
    n_tier_1_agg = 0                            # all lines with tier = 1 or tier = 1b
    n_tier_2 = 0                                # all lines with tier = 2
    n_tier_3 = 0                                # all lines with tier = 3
    n_tier_4 = 0                                # all lines with tier = 4

    matched_variants = 0                        # keep track of the total number of variants matched overall for the current sample
    # tier -> # variants
    matched_variants_mapping = {}               # keep track of the total number of variants matched per tier for the current sample

    # disease -> ct
    disease_mapping = {}                        # keep track of all disease names parsed across the current sample
    # tier -> disease -> None
    per_tier_disease_mapping = {}               # per tier, keep track of all disease names parsed across the current sample
    # ct -> disease -> None
    ct_mapping = {}                             # keep track of all disease names assigned to each "ct" class for the current sample
    # tier -> ct -> disease -> None
    per_tier_ct_mapping = {}                    # per tier, keep track of all disease names assigned to each "ct" class for the current sample

    matched_diseases = 0                        # keep track of the total number of unique disease names matched overall across all variants for the current sample
    # ct -> # diseases
    matched_diseases_ct_mapping = {}            # keep track of the total number of unique disease names matched per "ct" class across all variants for the current sample
    # tier -> # diseases
    per_tier_matched_diseases_mapping = {}      # keep track of the total number of unique disease names matched per tier for the current sample
    # tier -> ct -> # diseases
    per_tier_matched_diseases_ct_mapping = {}   # keep track of the total number of unique disease names matched per tier and "ct" class for the current sample

    # drug -> ct -> consensus_support
    consensus_drug_mapping = {}                 # keep track of the total number of unique drug names parsed across the consensus drug support for the current sample
    # drug -> ct -> consensus_support
    prior_consensus_drug_mapping = {}           # keep track of the total number of unique drug names parsed across the consensus drug support of the highest available "ct" class for the current sample
    # tier -> drug -> None
    per_tier_consensus_drug_mapping = {}        # per tier, keep track of the total number of unique drug names parsed across the consensus drug support for the current sample
    # tier -> drug -> None
    per_tier_prior_consensus_drug_mapping = {}  # per tier, keep track of the total number of unique drug names parsed across the consensus drug support of the highest available "ct" class for the current sample
    # ct -> drug -> consensus_support
    consensus_ct_mapping = {}                   # keep track of the total number of unique drug names parsed per "ct" class across the consensus drug support for the current sample
    # ct -> drug -> consensus_support
    prior_consensus_ct_mapping = {}             # keep track of the total number of unique drug names parsed per highest available "ct" class across the consensus drug support for the current sample

    # Each line in the input file corresponds to a single variant in the genome
    for line in infile:
        all_variants += 1
        line_split = line.strip().split("\t")

        ## 1) Process tier of the variant match
        tier = str(line_split[tier_pos].strip())

        # Avoid having issues due to using numbers/strings for the tiers
        tier = "tier_" + tier

        # Tier=4 should be skipped as no information was found on CIViC for the current gene
        # (i.e. all associated columns will be empty)
        if tier=="tier_4":
            n_tier_4 += 1
            continue

        # Keep track of the number of variants assigned to each tier in the current sample file
        all_civic_variants += 1
        if tier=="tier_3":
            n_tier_3 += 1
        if tier=="tier_2":
            n_tier_2 += 1
        if tier=="tier_1b":
            n_tier_1b += 1
        if tier=="tier_1":
            n_tier_1 += 1
        if tier=="tier_1" or tier=="tier_1b":
            n_tier_1_agg += 1


        ## 2) Process number of CIViC variants matched per line (and associated tier)

# FIXME
# Then, to compute mean simply divide the final count of # matches for each tier (e.g. dict["tier_1"]), by the total number of variants tagged as that tier (e.g. n_tier_1)
# we can compute the mean number of matched variants across the whole variant set (patient) by simply aggregating the match count across all tiers, and then dividing by the total number of variants with civic info


        # Parse column 'CIViC_Score', assumed to contain all CIViC variants matched to the current line (and their associated scores in CIViC)
        # NOTE: assume that variant names are listed using ";" as a separator character, and that duplicates are not possible
        civic_scores = str(line_split[civic_score_pos].strip())
        if civic_scores == ".":
            if not (tier=="3" or tier=="4"):
                raise ValueError("Encountered unexpected case of variant with tier!=3 and tier!=4 but no associated variant matches from CIViC in line %s" %(line.strip()))
            n_variants = 0
        else:
            # NOTE: for SNVs, all variant matches are ensured to originate from a single variant annotation + gene
            # NOTE: for CNVs, several variant matches arising from different genes (but same associated tier) can happen within the same line
            civic_score_list = civic_scores.split(';')
            n_variants = len(civic_score_list)

            # Sanity check there are no duplicated CIViC variant entries reported within the same line
            if n_variants != len(set(civic_score_list)):
                raise ValueError("Encountered duplicated variant matches from CIViC in line %s" %(line.strip()))
            # Sanity check that at least 1 variant should have been matched if column is not empty
            if n_variants == 0:
                raise ValueError("Encountered unexpected case of no associated variant matches from CIViC in line %s" %(line.strip()))

        # Keep track of number of matched variants, also keep counts per tier
        # tier -> # CIViC variants matched
        if tier not in matched_variants_mapping.keys():
            matched_variants_mapping[tier] = 0
        matched_variants_mapping[tier] += n_variants
        matched_variants += n_variants



        ## 3) Process disease names and associated cancer specificity classifications across all evidence types (columns)

# TODO: parse all diseases + ct available for each variant line (and associated tier) -> report # diseases for the whole variant set (patient) as well as counts per tier, and also report # diseases per ct classification for the whole variant set (patient) as well as counts per tier
# TODO: do not count unique diseases + ct per variant but rather, look at the whole set and retrieve the unique diseases+ct, and look at the tier subsets and retrieve the unique diseases+ct

        # disease -> None
        interim_disease_mapping = {}                # keep track of the total number of unique disease names parsed within the current variant line
        # ct -> disease -> None
        interim_ct_mapping = {}                     # keep track of the total number of unique disease names assigned to each "ct" class within the current variant line

        # Iterate across the 4 relevant evidence columns assumed to be present in the file
        # Keep track of all unique instances of disease_name + ct_type parsed across the available CIViC evidences
        for tmp_pos in [drug_pos, diag_pos, prog_pos, pred_pos]:
            # Empty evidence columns should be skipped as no information was found on CIViC
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
                # Sanity check that the same disease name can only be associated to one "ct" classification across the file
                else:
# FIXME
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



# FIXME: Selected info can be used to split the evidence strings in turn
                # Use the first part of the string (now already known) to split and retrieve all evidence items supporting the current claim (i.e. different combinations of evidence level + publication id)
                if tmp_pos == drug_pos:
                    # Format of interim_evidence_split (drug_and_evidence): 'DRUG(DIRECTION,CLINICALSIGNF(A(ref1,ref2,..),B(ref1,ref5)..))'
                    split_pattern = interim_drug + "(" + direction + "," + clinical_signf + "("
                else:
                    # Format of interim_evidence_split (ct_and_evidence): 'CT_TYPE(DIRECTION,CLINICALSIGNF(A(ref1,ref2,..),B(ref1,ref5)..))'
                    split_pattern = ct_type + "(" + direction + "," + clinical_signf + "("

                # Split returns list containing as many elements as evidence items supporting the current claim
                # NOTE: this is the case even when the same publication id is reported across several levels
                pub_ids_split = interim_evidence_split.split(split_pattern)[1].strip().split(",")



#                     # In CIVIC, drugs can form part of a combinatorial treatment, indicated by 'DRUG1+DRUG2+...'
#                     drug_list = interim_drug.strip().split("+")
#                     for drug in drug_list:
# FIXME
#                         if drug not in parsedDrugs:
#                             parsedDrugs.append(drug)
# 
#                         # At this step, it is critical to consider drug synonyms and special names that exist in CIViC
#                         drug_synonyms = [drug]

#                         # 1) Handle special cases of the form "SYNONYM 1 (SYNONYM 2)"
#                         # Introduce sanity check for 1 particular case: "BEZ235 (NVP-BEZ235, Dactolisib)"
#                         # Also, careful with "O(6)-BENZYLGUANINE"
#                         if (re.match(r'(\S+)\s+\((.+)\)', drug)) and ("," not in drug):
#                             tmp_drug_list = re.match(r'(\S+)\s+\((.+)\)', drug).groups()
#                             for tmp_drug in tmp_drug_list:
#                                 if tmp_drug not in drug_synonyms:
#                                     drug_synonyms.append(tmp_drug)

#                         # 2) Handle special cases where the drug in CIVIC corresponds to a known synonym (defined at the beginning of the script using drug_synonyms_mapping)
#                         # CIViC_synonym -> drug_name
#                         if drug in drug_synonyms_mapping.keys():
#                             drug_synonym = drug_synonyms_mapping[drug]
#                             if drug_synonym not in drug_synonyms:
#                                 drug_synonyms.append(drug_synonym)

#                         # Now, iterate all potential drug synonyms and associate them to all the available drug information
#                         # Duplication is necessary in this case to avoid missing drugs because of a name mismatch
#                         for tmp_synonym in drug_synonyms:
#                             # 2) Keep track of clinical information
#                             if tmp_synonym not in clinInfoDict.keys():
#                                 clinInfoDict[tmp_synonym] = {}
#                             if sampleName not in clinInfoDict[tmp_synonym].keys():
#                                 clinInfoDict[tmp_synonym][sampleName] = {}
#                             if gene not in clinInfoDict[tmp_synonym][sampleName].keys():
#                                 clinInfoDict[tmp_synonym][sampleName][gene] = {}
#                             if tier not in clinInfoDict[tmp_synonym][sampleName][gene].keys():
#                                 clinInfoDict[tmp_synonym][sampleName][gene][tier] = {}
#                             # Classify cancer into 'ct' (ie. cancer type specific) or 'nct'
#                             if cancerTag not in clinInfoDict[tmp_synonym][sampleName][gene][tier].keys():
#                                 clinInfoDict[tmp_synonym][sampleName][gene][tier][cancerTag] = {}
#                             # clinical significances associated to the same direction
#                             if direction not in clinInfoDict[tmp_synonym][sampleName][gene][tier][cancerTag].keys():
#                                 clinInfoDict[tmp_synonym][sampleName][gene][tier][cancerTag][direction] = []
#                             # Keep track of how many different evidence items (level+PMIDs) support this direction+clinicalSignificance
#                             # Also, this way we keep track of supporting evidence items for this direction+clinicalSignificance across multiple variants for gene of interest in one sample
#                             for x in range(0, len(pub_ids_split)): # contains as many elements as evidence items support this claim (even when the same reference is used for many levels)
#                                 clinInfoDict[tmp_synonym][sampleName][gene][tier][cancerTag][direction].append(clinSignf)


        ## At this point, all evidence items available across the 4 evidence type columns have been parsed
        ## Process and keep track of disease information parsed for the current variant line

        # Keep track of number of unique matched diseases across all variants, also keep counts per tier
        n_diseases = len(interim_disease_mapping.keys())
        # tier -> # diseases
        if tier not in per_tier_matched_diseases_mapping.keys():
            per_tier_matched_diseases_mapping[tier] = 0
        per_tier_matched_diseases_mapping[tier] += n_diseases
        matched_diseases += n_diseases


        # Keep track of number of unique matched diseases per "ct" class across all variants, also keep counts per tier
        # tier -> ct -> # diseases
        if tier not in per_tier_matched_diseases_ct_mapping.keys():
            per_tier_matched_diseases_ct_mapping[tier] = {}

        # Iterate all "ct" classes available for the current variant line and keep track of associated disease info
        # ct -> disease -> None
        for interim_ct_type in interim_ct_mapping.keys():
            interim_n_diseases = len(interim_ct_mapping[interim_ct_type].keys())
            # ct -> # diseases
            if interim_ct_type not in matched_diseases_ct_mapping.keys():
                matched_diseases_ct_mapping[interim_ct_type] = 0
            matched_diseases_ct_mapping[interim_ct_type] += interim_n_diseases
             # tier -> ct -> # diseases
            if interim_ct_type not in per_tier_matched_diseases_ct_mapping[tier].keys():
                per_tier_matched_diseases_ct_mapping[tier][interim_ct_type] = 0
            per_tier_matched_diseases_ct_mapping[tier][interim_ct_type] += interim_n_diseases



        ## 4) Process column listing consensus support across available drugs

        drug_infos = str(line_split[drug_supp_pos].strip())


        ## 3. TODO: (double check info with the corresponding info parsed from Predictive column?)

        # ct -> drug -> consensus_support
        interim_consensus_ct_mapping = {}

        # Process further variant lines which have drug consensus information available
        if drug_infos != ".":
            # Assume multiple consensus support strings can be listed separated by ";"
            # E.g.: 'ENTRECTINIB:NCT:CIVIC_SUPPORT;LAROTRECTINIB:NCT:CIVIC_SUPPORT;..'
            drug_infos_list = drug_infos.split(';')
            if not drug_infos_list:
                raise ValueError("Encountered unexpected case of no associated variant matches in CIViC for line %s" %(line.strip()))

            n_drug_info += 1

            # Iterate individual string of consensus drug support
            # Assume one string per combination of drug name + "ct" class (note the same drug can be available for different "ct" classes)
            for drug_info in drug_infos_list:
                # Format: 'DRUG_NAME:CT_TYPE:CIVIC_SUPPORT'
                drug_split = drug_info.strip().split(":")
                # Sanity check the expected separator character and format
                if len(drug_split) != 3:
                    raise ValueError("Encountered unexpected format of CIViC evidence annotations in column %s of line %s" %(str(tmp_pos+1), line.strip()))
                # Use uppercase for drug names and consensus support strings to avoid mismatches due to case
                drug = drug_split[0].strip().upper()
                consensus_support = drug_split[2].strip()
                # Use lowercase for "ct" classes to avoid mismatches due to case
                ct_type = drug_split[1].strip().lower()

                # Keep track of all disease names, associated "ct" class, and available consensus drug support strings in column 'CIViC_Drug_Support'
                # ct -> drug -> [consensus_support]
                if ct_type not in interim_consensus_ct_mapping.keys():
                    interim_consensus_ct_mapping[ct_type] = {}
                if drug not in interim_consensus_ct_mapping[ct_type].keys():
                    interim_consensus_ct_mapping[ct_type][drug] = []
                # Sanity check for duplicated consensus support strings for the same disease + "ct" class
                if consensus_support in interim_consensus_ct_mapping[ct_type][drug]:
                    print("Warning! Skipping duplicated consensus support '%s' for drug '%s' and cancer-specificity classification '%s' encountered in line %s" %(consensus_support, drug, ct_type, line.strip()))
                interim_consensus_ct_mapping[ct_type][drug].append(consensus_support)

        # TODO: for each variant line, select only drugs associated with the highest available ct classification (if only nct, then all will be taken into account) -> 
        # TODO: count # of unique drugs (already filtered for only highest ct) predicted for the whole variant set (patient), as well as counts per tier
        # TODO: also, per ct classification, report # unique drugs, and also report # unique drugs per tier
        # TODO: taking into account all available drugs + associated ct, for each unique drug, report # of ct classifications associated with it across the file, also, report counts per tier (then do mean across all drugs for the patient, as well as mean across all variants assigned to each tier)

        # Skip variant lines that are not associated to any consensus drug predictions
        if not interim_consensus_ct_mapping:
            continue


        # Sort all "ct" classes available for the current sample by the priority order defined at the beginning of this function (i.e. 'sorted_cts')
        sorted_ct_list = sorted(interim_consensus_ct_mapping.keys(), key=lambda x: sorted_cts.index(x))
        # Select the "ct" class with the highest priority (at least one will always be available)
        pick_ct = sorted_ct_list[0]

        # Keep track of the total number of unique (consensus) drug names parsed for the current sample, also keep track per tier and "ct" class available
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

                consensus_list = interim_consensus_ct_mapping[tmp_ct][tmp_drug]
                # NOTE: expectation is that each combination of disease name + "ct" class can only have one single consensus support string associated
                if len(consensus_list) > 1:
                    print("Warning! Encountered multiple consensus support strings ('%s') for drug '%s' and cancer-specificity classification '%s' in line %s" %(consensus_list, tmp_drug, tmp_ct, line.strip()))

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


        ## 4. TODO: check in the consensus drug support column: report, out of all unique drugs, how many (%) have associated: a. civic_support, civic_resistance, etc. Also, take into account situations where the same drug is associated to different classifications depending on the ct, variant, etc. -> how to handle?


    infile.close()
# FIXME
    print("Parsed %s lines with available CIViC information associated to %s genes and %s drugs" %(civicLines,len(parsedGenes),len(parsedDrugs)))


matched_diseases
matched_diseases_ct_mapping
per_tier_matched_diseases_ct_mapping



# # ct -> # diseases
# matched_diseases_ct_mapping = {}            # keep track of the total number of unique disease names matched per "ct" class across all variants for the current sample

# # tier -> ct -> # diseases
# per_tier_matched_diseases_ct_mapping = {}   # keep track of the total number of unique disease names matched per tier and "ct" class for the current sample



# #diseases, #diseases_tier1, #diseases_tier1b, #diseases_tier2, #diseases_tier3
# #diseases_ct, #diseases_gt, #diseases_nct, #diseases_tier1_ct, #diseases_tier1_gt, #diseases_tier1_nct, #diseases_tier1b_ct, #diseases_tier1b_gt, #diseases_tier1b_nct, #diseases_tier2_ct, #diseases_tier2_gt, #diseases_tier2_nct, #diseases_tier3_ct, #diseases_tier3_gt, #diseases_tier3_nct

# #drugs, #drugs_tier1, #drugs_tier1b, #drugs_tier2, #drugs_tier3
# mean_cts_per_drug
# #drugs_ct, #drugs_gt, #drugs_nct, #drugs_tier1_ct, #drugs_tier1_gt, #drugs_tier1_nct, #drugs_tier1b_ct, #drugs_tier1b_gt, #drugs_tier1b_nct, #drugs_tier2_ct, #drugs_tier2_gt, #drugs_tier2_nct, #drugs_tier3_ct, #drugs_tier3_gt, #drugs_tier3_nct

# #drugs_prior, #drugs_tier1_prior, #drugs_tier1b_prior, #drugs_tier2_prior, #drugs_tier3_prior
# #drugs_ct_prior, #drugs_gt_prior, #drugs_nct_prior, #drugs_tier1_ct_prior, #drugs_tier1_gt_prior, #drugs_tier1_nct_prior, #drugs_tier1b_ct_prior, #drugs_tier1b_gt_prior, #drugs_tier1b_nct_prior, #drugs_tier2_ct_prior, #drugs_tier2_gt_prior, #drugs_tier2_nct_prior, #drugs_tier3_ct_prior, #drugs_tier3_gt_prior, #drugs_tier3_nct_prior


# sample -> [

# #vars, #civic, #tier1, #tier1b, #tier1agg, #tier2, #tier3, #tier4
# mean_matched_vars, mean_matched_vars_tier1, mean_matched_vars_tier1b, mean_matched_vars_tier2, mean_matched_vars_tier3

    if sample_name not in civic_info_mapping.keys():
        civic_info_mapping[sample_name] = [all_variants, all_civic_variants, n_tier_1, n_tier_1b, n_tier_1_agg, n_tier_2, n_tier_3, n_tier_4, matched_variants, matched_diseases]

# TODO: sanity check for divisions by 0 (e.g. no variants in sample or no variants with CIViC info)

    # Compute mean number of matched variants for the sample (only makes sense to compute mean on lines that had CIViC matches available)
    mean_matched_variants = matched_variants / all_civic_variants

    # Per tier, compute mean number of matched variants for the sample (only makes sense to compute mean for tiers != 4)
    mean_matched_variants_tier1 = 0
    mean_matched_variants_tier1b = 0
    mean_matched_variants_tier2 = 0
    mean_matched_variants_tier3 = 0
    # tier -> # variants
    if "tier_1" in matched_variants_mapping.keys():
        mean_matched_variants_tier1 =⋅matched_variants_mapping["tier_1"] / n_tier_1
    if "tier_1b" in matched_variants_mapping.keys():
        mean_matched_variants_tier1b =⋅matched_variants_mapping["tier_1b"] / n_tier_1b
    if "tier_2" in matched_variants_mapping.keys():
        mean_matched_variants_tier2 =⋅matched_variants_mapping["tier_2"] / n_tier_2
    if "tier_3" in matched_variants_mapping.keys():
        mean_matched_variants_tier3 =⋅matched_variants_mapping["tier_3"] / n_tier_3

    # Compute mean number of matched diseases for the sample (only makes sense to compute mean on lines that had CIViC matches available)
    mean_matched_diseases = matched_diseases / all_civic_variants

    # Per tier, compute mean number of matched diseases for the sample (only makes sense to compute mean for tiers != 4)
    mean_matched_diseases_tier1 = 0
    mean_matched_diseases_tier1b = 0
    mean_matched_diseases_tier2 = 0
    mean_matched_diseases_tier3 = 0
    # tier -> # diseases
    if "tier_1" in per_tier_matched_diseases_mapping.keys():
        mean_matched_diseases_tier1 =⋅per_tier_matched_diseases_mapping["tier_1"] / n_tier_1
    if "tier_1b" in per_tier_matched_diseases_mapping.keys():
        mean_matched_diseases_tier1b =⋅per_tier_matched_diseases_mapping["tier_1b"] / n_tier_1b
    if "tier_2" in per_tier_matched_diseases_mapping.keys():
        mean_matched_diseases_tier2 =⋅per_tier_matched_diseases_mapping["tier_2"] / n_tier_2
    if "tier_3" in per_tier_matched_diseases_mapping.keys():
        mean_matched_diseases_tier3 =⋅per_tier_matched_diseases_mapping["tier_3"] / n_tier_3


    # Per ct, compute mean number of matched diseases for the sample
    mean_matched_diseases_ct = 0
    if "ct" in matched_diseases_ct_mapping.keys():
        mean_matched_diseases_ct = matched_diseases_ct_mapping["ct"] /
    mean_matched_diseases_gt = 0
    if "gt" in matched_diseases_ct_mapping.keys():
        mean_matched_diseases_gt = matched_diseases_ct_mapping["gt"] / 
    mean_matched_diseases_nct = 0
    if "nct" in matched_diseases_ct_mapping.keys():
        mean_matched_diseases_nct = matched_diseases_ct_mapping["nct"] /


    # Get total number of unique disease names parsed across the current sample
    n_diseases = len(disease_mapping.keys())


# TODO
    # Get total number of unique disease names parsed across the current sample
    #  len(disease_mapping.keys())
    # Per tier, get total number of unique disease names parsed across the current sample
    #  for tier in per_tier_disease_mapping.keys(): len(per_tier_disease_mapping[tier].keys())
    # Per ct, get total number of unique disease names parsed across the current sample
    #  for ct in ct_mapping.keys(): len(ct_mapping[ct].keys())
    # Per tier and ct, get total number of unique disease names parsed across the current sample
    #  for tier in per_tier_ct_mapping.keys(): for ct in per_tier_ct_mapping[tier].keys(): len(per_tier_ct_mapping[tier][ct].keys())


    # Get total number of unique drug names parsed across the current sample
    #  len(consensus_drug_mapping.keys())
    # Per tier, get total number of unique drug names parsed across the current sample
    #  for tier in per_tier_consensus_drug_mapping.keys(): len(per_tier_consensus_drug_mapping[tier].keys())
    # Per drug, retrieve number of associated "ct" classes across the current sample (no prioritization of "ct" class performed)
    # Compute mean across all drugs available for the sample
    #  sum_n_ct_across_drugs = 0 for drug in consensus_drug_mapping.keys(): n_ct = len(consensus_drug_mapping[drug].keys()) sum_n_ct_across_drugs += n_ct (sum_n_ct_across_drugs/n_drugs)
    # Per ct, get total number of unique drug names parsed across the current sample
    #  for ct in consensus_ct_mapping.keys(): len(consensus_ct_mapping[ct].keys())

    # Get total number of unique drug names parsed for the highest "ct" class available across the current sample
    #  len(prior_consensus_drug_mapping.keys())
    # Per tier, get total number of unique drug names parsed for the highest "ct" class available across the current sample
    #  for tier in per_tier_prior_consensus_drug_mapping.keys(): len(per_tier_prior_consensus_drug_mapping[tier].keys())
    # Per ct, get total number of unique drug names parsed for the highest "ct" class available across the current sample
    #  for ct in prior_consensus_ct_mapping.keys(): len(prior_consensus_ct_mapping[ct].keys())

    return (civic_info_mapping)


'''
Script
'''

parser = argparse.ArgumentParser(description='Combine SNV, CNV, DRS and CIViC drug predictions.')
parser.add_argument('--inputTable', dest='inputTable', required=True, help='Input file with combined SNV+DRS+CNV gene-drug predictions.')
parser.add_argument('--inputDir_civic_snv', dest='inputDir_civic_snv', required=True, help='Input directory with SNV CIViC files of format [sample].civic_snv.txt.')
parser.add_argument('--inputDir_civic_cnv', dest='inputDir_civic_cnv', required=True, help='Input directory with CNV CIViC files of format [sample].civic_cnv.txt.')
parser.add_argument('--fileEnding_snv', dest='fileEnding_snv', required=True, help='To retrieve correct input files, specify the desired file ending ("civic_snv.txt" for SNV data).')
parser.add_argument('--fileEnding_cnv', dest='fileEnding_cnv', required=True, help='To retrieve correct input files, specify the desired file ending ("civic_cnv.txt" for CNV data).')
parser.add_argument('--outFileTag', dest='outFileTag', required=True, help='Name prefix of the output files.')

args = parser.parse_args()


### Process input files for CIViC (pool SNV+CNV data together)

# Global dictionary (pools SNV+CNV CIViC info) to match drugs to the underlying target genes and CIViC support in each sample
## TODO: Beware of inconsistency between tier3 in SNV and tier3 in CNV (all vs. CNV-only CIViC records)
##       Because of this, when pooling and same drug+sample+gene has tier3 for both, CNV entries will be duplicated because they are also in SNV results

civic_info_mapping = {}

## 4a) Process input files for SNV CIViC
seenSamples_snv = []
for file in os.listdir(args.inputDir_civic_snv):
    sampleFile_snv = "%s%s" %(args.inputDir_civic_snv,os.path.basename(file))
    if os.path.isfile(sampleFile_snv) and sampleFile_snv.endswith(args.fileEnding_snv):
        name_snv = os.path.basename(file).split("_")[0].split("-")[1]
        if name_snv not in seenSamples_snv:
            seenSamples_snv.append(name_snv)
        else:
            print("Error! Multiple SNV CIViC files for %s!" %(name_snv))
            sys.exit(1)
        (civic_info_mapping) = parse_input_file(sampleFile_snv,name_snv,civic_info_mapping)

## 4b) Process input files for CNV CIViC
seenSamples_cnv = []
for file in os.listdir(args.inputDir_civic_cnv):
    sampleFile_cnv = "%s%s" %(args.inputDir_civic_cnv,os.path.basename(file))
    if os.path.isfile(sampleFile_cnv) and sampleFile_cnv.endswith(args.fileEnding_cnv):
        name_cnv = os.path.basename(file).split("_")[0].split("-")[1]
        if name_cnv not in seenSamples_cnv:
            seenSamples_cnv.append(name_cnv)
        else:
            print("Error! Multiple CNV CIViC files for %s!" %(name_cnv))
            sys.exit(1)
        (civic_info_mapping) = parse_input_file(sampleFile_cnv,name_cnv,civic_info_mapping)

## SNV and CNV samples must be the same (n=412)
if (set(seenSamples_snv) != set(seenSamples_cnv)):
    print("Error! Sample names in input directories %s and %s do not match." %(args.inputDir_civic_snv,args.inputDir_civic_cnv))
    sys.exit(1)

## TODO: For now, ignore DRS information, ie. only DGIDB+CIViC
infile = open(args.inputTable,'r')
outfile_det = open(args.outFileTag + ".details.tsv",'w') # Details table
outfile_onco = open(args.outFileTag + ".condensed.tsv",'w') # Oncoprint table
outfile_infos = open(args.outFileTag + ".civicInfos.tsv",'w') # CIViC details table

onco_header = "Drug\tNumberSamples_genomic\tNumberSamples_drs\tNumberSamples_civic\tCancer_related\tContained_in_snv\tContained_in_drs\tContained_in_cnv\tContained_in_civic"
infos_header = "GENE\tDRUG\tCONTAINED_IN_CIVIC_RESULTS\tMUTATED_SAMPLES\tCIVIC_SAMPLES\tPOSITIVE_SAMPLES\tNEGATIVE_SAMPLES\tUNKNOWN_CONFLICT_SAMPLES\tUNKNOWN_DNS_SAMPLES\tUNKNOWN_BLANK_SAMPLES\tUNKNOWN_NOVARMATCH_SAMPLES"

# Dictionary to retrieve the corresponding sample name given a position (ie. column) in the line
dict_predPos = {}
headerPred = infile.readline().strip().split('\t')
startSampleCols = 7
for posHead in range(startSampleCols,len(headerPred)):
    ## Add current sample to header of output files
    onco_header += "\t" + headerPred[posHead]
    dict_predPos[posHead] = headerPred[posHead]

## Write complete header to output files
outfile_onco.write(onco_header + "\n")
outfile_det.write(onco_header + "\n")  # header is the same for both the oncoprint and details tables

## Classification of evidences will be counted for each gene and drug (based on structure of supportDict)
## Generate all combinations of direction + clinical significance
countDict = {}
## Use the following dictionary to classify each combination of direction+clinSignf as POSITIVE or NEGATIVE and use this information in the output header column names
colNameDict = {}
for direction in supportDict:
    for clinSignf in supportDict[direction].keys():
        supportType = direction + "_" + clinSignf
        countDict[supportType] = 0
        colNameDict[supportType] = supportDict[direction][clinSignf]
## Add 'UNKNOWN_BLANK' support as well (not present in supportDict)
## Also, this type is not classified as POSITIVE or NEGATIVE for the column names (see below)
countDict['UNKNOWN_BLANK'] = 0

# Dictionary to retrieve the corresponding support count given a position (ie. column) in the line
dict_infosPos ={}
## Generate and write header to output file
for posInfos,supportType in enumerate(countDict.keys()):
    if supportType in colNameDict.keys():
        ## Change UNKNOWN_DNS (due to DOES_NOT_SUPPORT) to simply UNKNOWN
        if colNameDict[supportType] == 'UNKNOWN_DNS':
            infos_header += "\t" + 'UNKNOWN_' + supportType + "_TOTAL"
        else:
            infos_header += "\t" + colNameDict[supportType] + '_' + supportType + "_TOTAL"
    ## For UNKNOWN (due to blanks), there is no classification as POSITIVE OR NEGATIVE, so do not add in column name
    else:
        infos_header += "\t" + supportType + "_TOTAL"
    dict_infosPos[posInfos] = supportType

## Write complete header to output file
outfile_infos.write(infos_header + "\n")


seenDrugs = []
## Iterate drug table (ie. each line corresponds to unique drug)
for line in infile:
    ## Initialize dictionary for counting classification of evidences for current drug
    drugCountMap = {}
    ## Initialize dictionary for counting classification of samples for each drug-gene interaction relating to current drug
    ## array[0] = mutSamples, array[1] = posSamples, array[2] = negSamples, array[3] = unkSamples, array[4] = unkSamples_noMatchVar,
    ## array[5] = civicSamples, array[6] = dnsSamples, array[7] = conflictSamples
    drugCountSampleMap = {}

    line_split = line.strip().split('\t')
    ## Retrieve drug name for current line
    predDrug = line_split[0].upper()
    if predDrug not in seenDrugs:
        seenDrugs.append(predDrug)
    else:
        print("Warning! Drug %s contained twice in input file." %(predDrug))
        continue

    sampleNum_genomic = line_split[1]
    sampleNum_drs = line_split[2]
    cancer_rel = line_split[3]
    matchedSNV = line_split[4]
    matchedDRS = line_split[5]
    matchedCNV = line_split[6]

    matchedCIVIC = 'n'
    ## Pooled CIViC data for SNV and CNV
    if predDrug in civic_info_mapping.keys():
        matchedCIVIC = 'y'

    ## Drug and NumberSamples columns will be added at the end
    outString_det = cancer_rel + '\t' + matchedSNV + '\t' + matchedDRS + '\t' + matchedCNV + '\t' + matchedCIVIC
    outString_onco = cancer_rel + '\t' + matchedSNV + '\t' + matchedDRS + '\t' + matchedCNV + '\t' + matchedCIVIC

    ## Keep track of total no. samples that had support for CIViC (any type)
    samples_civic = 0

    ## Iterate samples and combine predictions for current drug
    for pos in range(startSampleCols,len(line_split)):
        ## Retrieve corresponding sample name for current position
        sample = dict_predPos[pos]
        predDGIDB = line_split[pos]
        predSupport = "no_mutation"
        ctSupport = ''
        tierSupport = '4'

        ## If at least one gene was found, this means there is a DGIDB prediction
        if 'score' in predDGIDB:
            predSupport = "dgidb_only"

        ## Parse gene(s) predicted by DGIDB for currently evaluated sample and drug (if any), and combine CIViC support for each (if available)
        ## The input table in this case corresponds to a combination of SNV+DRS+CNV predictions, in the specified order and separated with '|'
        ## Format eg.: 'PIK3CA (multiple variants, score: 9) | ERBB2;MIR4728 (score: 5) | -0.0183 | TP53 (score: 30) | ERBB2 (score: 7)'
        predSplit = predDGIDB.split(' | ')
        ## There can be several DGIDB gene predictions within one data type (ie. snv or cnv)
        predGenes = []
        genesWithSupport = []
        predStrings = []
        skip_conflict = False
        ## Only parse DGIDB predicionts and skip DRS; anyway only possible to retrieve necessary gene info from DGIDB
        for singlePred in predSplit:
            tempPred = singlePred
            ## Skip accounting for current gene because it was already seen (due to snv+cnv mutation info)
            skip_count = False  
            if 'score' in singlePred:
                gene = singlePred.split('(')[0].strip()
                ## Check for format containing variant annotations
                ## Format: eg. CDKN2A:c.355G>A|p.Glu119Lys;c.365G>A|p.Arg122Gln;c.*245G>A| (multiple variants, score: 4)
                gene_and_var = ''
                if ":" in gene:
                    gene_and_var = singlePred.split('(')[0].strip()
                    gene = gene_and_var.split(":")[0].strip()
                if gene not in predGenes:
                    predGenes.append(gene)
                else:
                    skip_count = True

                ## Combine DGIDB prediction with CIViC information (if available)
                ## Since many genes can be predicted for a single drug+sample, iterate all and report best CIViC support
                tempSupport = "dgidb_only"
                ctGene = ''         # whether CIViC info derives from ct/gt/nct
                tierGene = ''       # whether CIViC info derives from tier 1/2/3
                ## Initialize dictionaries for counting classification of evidences and samples, for current gene (and drug)
                if gene not in drugCountMap.keys():
                    ## One entry per support type (ie. combination of direction+clinical significance), including 'UNKNOWN'
                    ## If gene has no support for CIViC, ie. dgidb_only, then all entries will be 0
                    drugCountMap[gene] = copy.deepcopy(countDict)
                    ## For each gene-drug interaction of current drug, keep track of support across samples
                    ##  array[0] = mutSamples, array[1] = posSamples, array[2] = negSamples, array[3] = unkSamples, array[4] = unkSamples_noMatchVar,
                    ##  array[5] = civicSamples, array[6] = dnsSamples, array[7] = conflictSamples
                    ## If gene has no support for CIViC, dgidb_only, then all entries will be 0 except first one (dgidb)
                    drugCountSampleMap[gene] = [0,0,0,0,0,0,0,0]

                ## Account this gene as predicted by DGIDB for the currently evaluated drug
                if not skip_count:
                    drugCountSampleMap[gene][0] += 1
                geneSupport = []
                ## Check whether there is any CIViC info available for current DGIDB prediciton
                noCIVIC = False
                if predDrug not in civic_info_mapping.keys():
                    noCIVIC = True
                elif sample not in civic_info_mapping[predDrug].keys():
                    noCIVIC = True
                elif gene not in civic_info_mapping[predDrug][sample].keys():
                    noCIVIC = True

                ## If not, then continue iteration of predictions and append (unchanged) prediction details string for reporting in final output
                if noCIVIC:
                    predStrings.append(tempPred)
                    continue

                ## Add sanity check for gene having already been evaluated for CIViC support (ie. when same gene is mutated in SNV and CNV)
                ## Required as we will only add 'MULTIPLE' tag when several **different** genes have been evaluated for CIViC support
                if gene not in genesWithSupport:
                    genesWithSupport.append(gene)
                ## Account this gene as having CIViC info available for the currently evaluated drug
                if not skip_count:
                    drugCountSampleMap[gene][5] += 1

                ## In dict civic_info_mapping, only available tiers are 1,2,3 (4 is excluded since no clinical info in this case)
                allTiers = list(civic_info_mapping[predDrug][sample][gene].keys())
                ## Just report clinical info for a single tier (if many are available). Prioritize available tiers for this (1>2>3)
                ## Only combination of drug+sample+gene that has >1 tier -> AFATINIB + A6B1 + ERBB3 (tier 2, tier 3)
                if (len(allTiers)>1):
                    tier = min(allTiers)
                    tierGene = tier
                else:
                    tier = allTiers[0]
                    tierGene = tier

                ## Classify sample as unknown when only tier 3 is available for drug+sample+gene and no other
                if tier == '3':
                    geneSupport = ['UNKNOWN_TIER3']
                else:
                    ## Do not ignore cancer specificity but generate several 'supports' depending on the ct
                    geneSupport_ct = []
                    geneSupport_gt = []
                    geneSupport_nct = []
                    for ct in civic_info_mapping[predDrug][sample][gene][tier].keys():
                        for direction in civic_info_mapping[predDrug][sample][gene][tier][ct].keys():
                            for clinSig in civic_info_mapping[predDrug][sample][gene][tier][ct][direction]:
                                ## Unknown is assigned when there is at least one 'N/A' or 'NULL for direction or clinical significance
                                if ('NULL' in direction) or ('N/A' in direction) or ('NULL' in clinSig) or ('N/A' in clinSig):
                                    thisSupport = 'UNKNOWN_BLANK'
                                    ## Increase by 1 the counter for unknown significance
                                    if not skip_count:
                                        drugCountMap[gene][thisSupport] += 1
                                else:
                                    if direction not in supportDict.keys():
                                        print("Error! Could not find direction %s in support dictionary." %(direction))
                                        sys.exit(1)
                                    if clinSig not in supportDict[direction].keys():
                                        print("Error! Could not find clinical significance %s in support dictionary." %(clinSig))
                                        sys.exit(1)
                                    thisSupport = supportDict[direction][clinSig]
                                    ## Increase by 1 the counter for current combination of direction+clinical significance
                                    if not skip_count:
                                        drugCountMap[gene][direction + "_" + clinSig] += 1
                                ## Split gene support by cancer-type specificity
                                ## Keep track of all available support evidence
                                if ct == 'ct':
                                    geneSupport_ct.append(thisSupport)
                                elif ct == 'gt':
                                    geneSupport_gt.append(thisSupport)
                                elif ct == 'nct':
                                    geneSupport_nct.append(thisSupport)

                    ## Once all clinical info has been iterated, choose what supporting evidence to report for current gene based on available cancer-specificity
                    ## Apply hierarchy: ct > gt > nct
                    if geneSupport_ct:
                        geneSupport = geneSupport_ct
                        ctGene = 'ct'
                    elif geneSupport_gt:
                        geneSupport = geneSupport_gt
                        ctGene = 'gt'
                    elif geneSupport_nct:
                        geneSupport = geneSupport_nct
                    else:
                        print("Error! No clinical information was found for gene %s." %(gene))
                        sys.exit(1)

                ## Counters for available CIViC evidence items
                count_pos = 0
                count_neg = 0
                count_unk = 0
                count_dns = 0

                ## Classify sample (for given gene) depending on supporting evidence
                ## Tier 3 case: collected evidence is too disperse to be reported
                if ('UNKNOWN_TIER3' in geneSupport):
                    tempSupport = "civic_unspecific_vars"
                    if not skip_count:
                        drugCountSampleMap[gene][4] += 1

                ## Remaining classification includes POSITIVE, NEGATIVE, UNKNOWN_BLANK (ie. due to blank info) and UNKNOWN_DNS (ie. due to DOES_NOT_SUPPORT evidence)
                else:
                    count_pos = geneSupport.count('POSITIVE')
                    count_neg = geneSupport.count('NEGATIVE')
                    count_unk = geneSupport.count('UNKNOWN_BLANK')
                    count_dns = geneSupport.count('UNKNOWN_DNS')
                    ## Pool UNKNOWN_BLANK and UNKNOWN_DNS together (as both result in unknown CIViC support)
                    count_total_unk = count_unk + count_dns
                    ## Sanity check that there is at least some support
                    if (count_pos == 0) and (count_neg == 0) and (count_total_unk == 0):
                        print("Error! Unexpected support case for gene %s." %(gene))
                        sys.exit(1)
                    ## Resolve contradicting evidence (if any) by majority vote
                    ## For this, pool UNKNOWN_BLANK and UNKNOWN_DNS together
                    if (count_total_unk > count_pos) and (count_total_unk > count_neg):
                        ## In oncoprint, simply classify samples as 'civic_unknown'
                        tempSupport = "civic_unknown"
                        ## However, for detailed table, distinguish between UNKNOWN_BLANK and UNKNOWN_DNS
                        ## If count is equal for both types of unknowns, classify as UNKNOWN_DNS
                        if not skip_count:
                            if count_dns >= count_unk:
                                drugCountSampleMap[gene][6] += 1
                            else:
                                drugCountSampleMap[gene][3] += 1
                    elif count_pos == count_neg:
                        tempSupport = "civic_conflict"
                        if not skip_count:
                            drugCountSampleMap[gene][7] += 1
                    elif (count_pos > count_neg) and (count_pos >= count_total_unk):
                        tempSupport = "civic_support"
                        if not skip_count:
                            drugCountSampleMap[gene][1] += 1
                    elif (count_neg > count_pos) and (count_neg >= count_total_unk):
                        tempSupport = "civic_resistance"
                        if not skip_count:
                            drugCountSampleMap[gene][2] += 1
                    else:
                        print("Error! Unexpected support case for drug %s, sample %s and gene %s." %(predDrug,sample,gene))
                        sys.exit(1)

                ## Decide final support to report based on the type of ct and achieved civic support
                ##  - First, decide on cancer specificity: ct > gt > nct (always take the best possible ct no matter the associated CIViC support)
                ##  - Then, decide on civic support (use dictionary 'supportPrior' for this):
                ##      support = resistance > civic_conflict > civic_unknown > civic_unspecific_vars (tier3) > dgidb_only
                if tempSupport != "dgidb_only":
                    ## If available evidence was found to be conflicting across the mutated genes (conflicting signals), skip further classification and report this information
                    ## However, if available evidence was found to be conflicting WITHIN a gene (conflicting evidences), then consider it for prioritization of CIViC support
                    if not skip_conflict:
                        ## Change final civic support if current gene has same or better tier
                        if (tierGene <= tierSupport):
                            ## Case of better tier, always keep it
                            if (tierGene < tierSupport):
                                ctSupport = ctGene
                                tierSupport = tierGene
                                predSupport = tempSupport
                            ## Case of same tier, change final civic support if current gene has same or better ct
                            elif (ctGene == ctSupport) or (ctGene == 'ct' and ctSupport == 'gt') or (ctGene and not ctSupport):
                                ## Case of better ct, always keep it
                                if (ctGene != ctSupport):
                                    ctSupport = ctGene
                                    tierSupport = tierGene
                                    predSupport = tempSupport
                                ## Case of same ct, decide on civic support based on priority of categories
                                else:
                                    predSupport_score = supportPrior[predSupport]
                                    tempSupport_score = supportPrior[tempSupport]
                                    ## Case of both support and resistance: conflicting signals
                                    if (predSupport_score==1 and tempSupport_score==1) and (predSupport != tempSupport):
                                        ctSupport = ctGene
                                        tierSupport = tierGene
                                        predSupport = "civic_conflict"
                                        # Do not prioritize support any further from now on
                                        skip_conflict = True
                                    ## Case of better support, always keep it
                                    elif (predSupport_score > tempSupport_score):
                                        ctSupport = ctGene
                                        tierSupport = tierGene
                                        predSupport = tempSupport
                                    ## Same support will stay the same anyway

                ## Finally, add CIViC support info to current gene prediction details
                ## First, retrieve DGIDB score (really, whole content inside brackets) for currently evaluated DGIDB prediction
                predContent = singlePred.split('(')[1].strip()
                predContent = predContent.split(')')[0].strip()
                ctString = ctGene
                if not ctString:
                    ctString = 'nct'
                ## Format of CIViC details string: tier | ct | no.positives | no.negatives | no.unknown_dns | no.unknown_blank
                civicString = tierGene + '|' + ctString + '|' + str(count_pos) + '|' + str(count_neg) + '|' + str(count_dns) + '|' + str(count_unk)
                ## Reconstruct initial DGIDB prediction details string, now also adding CIViC info
                predContent += ', civic: ' + civicString
                ## Here, tempPred value is changed (if no CIViC details were available, then temPred would correspond to original value in table)
                ## Check for format containing variant annotations
                if gene_and_var:
                    tempPred = gene_and_var + ' (' + predContent + ')'
                else:
                    tempPred = gene + ' (' + predContent + ')'

                #### **End of genes iteration (within single DGIDB prediction)

            ## Append prediction details string for reporting in final output
            ## Note that 'tempPred' may have changed due to addition to CIViC details (if any)
            predStrings.append(tempPred)

        #### **End of DGIDB predictions iteration**

        ## Add support for current sample to output file
        ## Include whether current support was derived from ct/gt information or not
        if ctSupport:
            ctSupport = '_' + ctSupport
        outString_onco += '\t' + str(predSupport + ctSupport)
        ## Add details string (containing CIViC support details as well when available) for current sample to output file
        outString_det += '\t' + ' | '.join(predStrings)

        ## Account for sample having some kind of CIViC support (any) for current drug
        if len(genesWithSupport) > 0:
            samples_civic += 1
        ## Sanity check of whether >1 **different** gene per drug+patient has CIViC support
        ## Cases of same gene being mutated in SNV and CNV do not get MULTIPLE tag
        if len(genesWithSupport) > 1:
            ## Include tag for the oncoprint to indicate that >1 gene had CIViC support for the given drug+sample
            outString_onco += ';MULTIPLE'
            print("Warning! Found >1 gene with CIViC support for drug %s and sample %s." %(predDrug,sample))

    #### **End of sample iteration**

    ## Add drug name and sample count for current drug (line)
    ## Write line for current drug (all samples) to output files
    outfile_onco.write("%s\t%s\t%s\t%s\t%s\n" %(predDrug,sampleNum_genomic,sampleNum_drs,samples_civic,outString_onco))
    outfile_det.write("%s\t%s\t%s\t%s\t%s\n" %(predDrug,sampleNum_genomic,sampleNum_drs,samples_civic,outString_det))

    ## Write to output aggregated info across all samples and associated genes, for current drug
    for gene in drugCountMap.keys():
        geneCounts = drugCountMap[gene]
        ## array[0] = mutSamples, array[1] = posSamples, array[2] = negSamples, array[3] = unkSamples, array[4] = unkSamples_noMatchVar,
        ## array[5] = civicSamples, array[6] = dnsSamples, array[7] = conflictSamples
        geneSampleCounts = drugCountSampleMap[gene]
        ## If there is at least one sample with CIViC information,
        ## that means current gene+drug combination is present in CIViC
        contained_CIVIC = 'n'
        if geneSampleCounts[5]>0:
            contained_CIVIC = 'y'
        outString_infos = gene + '\t' + predDrug + '\t' + contained_CIVIC + '\t' + str(geneSampleCounts[0]) + '\t' + str(geneSampleCounts[5]) + '\t' + str(geneSampleCounts[1]) + '\t' + str(geneSampleCounts[2]) + '\t' + str(geneSampleCounts[7]) + '\t' + str(geneSampleCounts[6]) + '\t' + str(geneSampleCounts[3]) + '\t' +  str(geneSampleCounts[4])
        for i in range(0,len(dict_infosPos)):
            supportType = dict_infosPos[i]
            supportCount = geneCounts[supportType]
            outString_infos += '\t' + str(supportCount)

        ## Add line for support of current gene+drug to output file
        outfile_infos.write(outString_infos + '\n')

outfile_onco.close()
outfile_det.close()
outfile_infos.close()
print("\nParsed %s predicted drugs for %s samples." %(len(seenDrugs),len(dict_predPos.keys())))
