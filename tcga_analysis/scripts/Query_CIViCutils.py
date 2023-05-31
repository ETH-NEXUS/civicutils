#!/usr/bin/env python

"""
Use CIViCutils to query and process CIViC variant information
Requires newer format of variant annotations, i.e.: GENE:variant|variant;GENE:variant;...
Version queries offline cache provided in civicpy module

Lourdes Rosano, Feb 2022
"""

import sys
import argparse

## Load relevant functions from CIViCutils package

sys.path.append("/path/to/civicutils")

from read_and_write import get_dict_support, write_evidences, write_output_line
from query import query_civic
from filtering import filter_civic
from match import match_in_civic, annotate_ct, filter_ct, process_drug_support, reprocess_drug_support_across_selected_variants
from utils import check_match_before_writing, check_keys, check_keys_not, check_data_type, check_dict_entry


"""
Functions
"""

# Given the header, find the input columns containing genes and variants/cnvs values (use provided arguments directly)
# For SNVs, also retrieve the input columns containing variant impacts and exons
def get_column_positions(header, data_type):
    header_split = header.split("\t")
    gene_pos = -1
    var_pos = -1
    impact_pos = -1
    exon_pos = -1
    for pos in range(0,len(header_split)):
        # Avoid mismatches due to case by always using uppercase
        if args.colname_gene.upper() == header_split[pos].upper():
            gene_pos = pos
        if args.colname_data.upper() == header_split[pos].upper():
            var_pos = pos
        if data_type == "SNV":
            if args.colname_impact.upper() == header_split[pos].upper():
                impact_pos = pos
            if args.colname_exon.upper() == header_split[pos].upper():
                exon_pos = pos
    # Check if required columns could be found in the header
    if (gene_pos == -1) or (var_pos == -1):
        print("\nError! Could not match all input columns in header %s." %(header))
        sys.exit(1)
    # Check if optional columns were provided (only for SNVs)
    if data_type == "SNV":
        if (impact_pos == -1) or (exon_pos == -1):
            print("\nError! Could not match all input columns in header %s." %(header))
            sys.exit(1)
    return (gene_pos, var_pos, impact_pos, exon_pos)


# Helper function to parse and process an input file of SNV data into structured dictionaries
# Adapted from CIViCutils to be able to parse lines having multiple gene (and variant) annotations associated
# Assumes header and that relevant column names provided as arguments are found (impact and exon columns are optional)
def read_in_snvs_multiple_genes(input_snv_file):
    # n_line -> [[var, gene, dna, prot, impact, exon], [var, gene, dna, prot, impact, exon], ...]
    raw_mapping = {}
    # gene -> variant (dna|prot|impact|exon|n_line) -> null
    snv_mapping = {}
    infile = open(input_snv_file, "r")
    input_header = infile.readline().strip()
    # Retrieve indices for all relevant columns (i.e. gene and variant annotation; for SNV also variant impact and exon info)
    # If not all are found, script will exit with an error
    # When data_type=CNV, impact_pos=-1 and exon_pos=-1 since they will be disregared
    (gene_pos, var_pos, impact_pos, exon_pos) = get_column_positions(input_header, data_type)
    # Each parsed line corresponds to a different genomic variant (which can have several annotations available)
    for n_line,line in enumerate(infile):
        line_split = line.strip().split("\t")
        all_variants = line_split[var_pos].strip().split(";")
        all_impacts = line_split[impact_pos].strip().split(";")
        all_exons = line_split[exon_pos].strip().split(";")
        # Sanity check the the 3 columns contain the same number of variants annotations, sorted in the same order (1-1 correspondance of positions)
        # Adapted from: https://stackoverflow.com/questions/35791051/better-way-to-check-if-all-lists-in-a-list-are-the-same-length
        check_lists = [all_variants, all_impacts, all_exons]
        it = iter(check_lists)
        n_variants = len(next(it))
        if not all(len(l) == n_variants for l in it):
            raise ValueError("Encountered different number of available annotations in columns '%s', '%s' and '%s'!" %(args.colname_data, args.colname_impact, args.colname_exon))
        # n_line -> [var1,..,varN]}
        raw_mapping[str(n_line)] = []
        # Loop all available items (same number), and retrieve for each variant corresponding info from each list
        # For each element, retrieve single gene (separated by ":"), retrieve single c./g./n. annotation, and retrieve single p. annotation when available
        # For each element, retrieve the corresponding variant impact (separated from gene by ":") and exon info (separated from gene and variant by ":")
        for n_variant in range(0, n_variants):
            variant = all_variants[n_variant]
            variant_split = variant.strip().split(":")
            gene = variant_split[0].strip()
            hgvs_split = variant_split[1].strip().split("|")
            c_variant = hgvs_split[0].strip()
            p_variant = hgvs_split[1].strip()
            # Sanity check required information is never empty
            # NOTE: p_variant is considered optional
            if not gene:
                raise ValueError("Encountered empty gene annotation in column '%s'!" %(args.colname_data))
            if not c_variant:
                raise ValueError("Encountered empty variant annotation in column '%s'!" %(args.colname_data))
            impact_tag = all_impacts[n_variant]
            impact_split = impact_tag.strip().split(":")
            if impact_split[0].strip() != gene:
                raise ValueError("Encountered different genes in positional annotations of columns '%s' and '%s'!" %(args.colname_data, args.colname_impact))
            impact = impact_split[1].strip()
            exon_tag = all_exons[n_variant]
            exon_split = exon_tag.strip().split(":")
            if exon_split[0].strip() != gene:
                raise ValueError("Encountered different genes in positional annotations of columns '%s' and '%s'!" %(args.colname_data, args.colname_exon))
            if exon_split[1].strip() != variant_split[1].strip():
                raise ValueError("Encountered different variants in positional annotations of columns '%s' and '%s'!" %(args.colname_data, args.colname_exon))
            exon = exon_split[2].strip()
            # Sanity check required information is never empty
            # NOTE: exon is considered optional
            if not impact:
                raise ValueError("Encountered empty variant impact annotation in column '%s'!" %(args.colname_impact))
            # NOTE: for now, skip variants which do not have a valid c.HGVS expression available (e.g. n.HGVS)
            # We can do this because we will use alternative custom function "write_match_multiple_annotations" to write results to output, and these variant annotations will not be removed from the affected columns
            if not c_variant.startswith("c."):
                continue
            # Process raw_mapping to have gene-centered dict
            # Returns dict of gene -> [var1,var2,..,varN], where a given var="dna|prot|impact|exon|n_line"
            if gene not in snv_mapping.keys():
                snv_mapping[gene] = {}
            # Collapse variant info separated with "|"
            # Keep track of what line each variant comes from
            variant = c_variant + "|" + p_variant + "|" + impact + "|" + exon + "|" + str(n_line)
        # NOTE: Variants can never be duplicated because of the different row numbers assigned
        # if variant in snv_mapping[gene].keys():
        #     print("Found duplicated variant '%s|%s' for gene '%s' in line '%s'!" %(c_variant, p_variant, gene, str(n_line)))
        #     sys.exit(1)
            snv_mapping[gene][variant] = None
            raw_mapping[str(n_line)].append([variant, gene, c_variant, p_variant, impact, exon])
    infile.close()
    return (raw_mapping, snv_mapping, input_header)


# Helper function to parse and process an input file of CNV data into structured dictionaries
# Adapted from CIViCutils to be able to parse lines having multiple genes annotated for the same CNV
# Assumes header and that relevant column names provided as arguments are found (no optional columns)
def read_in_cnvs_multiple_genes(input_cnv_file):
    # n_line -> {details -> cnv, genes -> [gene1,..,geneN]}
    raw_mapping = {}
    # gene -> variant (cnv|n_line) -> null
    cnv_data = {}
    infile = open(input_cnv_file, "r")
    input_header = infile.readline().strip()
    # Retrieve indices for all relevant columns (i.e. gene and variant annotation; for SNV also variant impact and exon info)
    # If not all are found, script will exit with an error
    # When data_type=CNV, impact_pos=-1 and exon_pos=-1 since they will be disregared
    (gene_pos, var_pos, impact_pos, exon_pos) = get_column_positions(input_header, data_type)
    # Each parsed line corresponds to a different CNV (which can have several genes annotated)
    for n_line,line in enumerate(infile):
        line_split = line.strip().split("\t")
        cnv = line_split[var_pos].strip()
        # n_line -> {details -> cnv, genes -> [gene1,..,geneN]}
        raw_mapping[str(n_line)] = {}
        raw_mapping[str(n_line)]["details"] = cnv
        raw_mapping[str(n_line)]["genes"] = []
        # NOTE: in this case, the same CNV can contain several annotated genes separated with ";" in the Gene column
        genes = line_split[gene_pos].strip().split(";")
        # Collapse variant info separated with "|"
        # Keep track of what line each variant comes from
        variant = cnv + "|" + str(n_line)
        for gene in genes:
            raw_mapping[str(n_line)]["genes"].append(gene)
            # Process raw_mapping to have gene-centered dict
            # Returns dict of gene -> [var1,var2,..,varN], where a given var="cnv|n_line"
            if gene not in cnv_data.keys():
                cnv_data[gene] = {}
        # NOTE: Variants can never be duplicated because of the different row numbers assigned
        # if variant in cnv_data[gene].keys():
        #     print("Found duplicated variant '%s' for gene '%s' in line '%s'!" %(cnv, gene, str(n_line)))
        #     sys.exit(1)
            cnv_data[gene][variant] = None
    infile.close()
    return (raw_mapping, cnv_data, input_header)


# Process and prioritize the tier classifications available for the CIViCutils annotations of a given gene and variant
def check_tier_and_matches(gene, combined_id, match_mapping, has_support=True):
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]
    # Check if match_mapping contains the provided input variants
    if gene not in match_mapping.keys():
        raise ValueError("Provided gene '%s' is not contained in provided 'match_mapping'." %(gene))
    if combined_id not in match_mapping[gene].keys():
        raise ValueError("Provided variant '%s' is not contained in supplied 'match_mapping' for gene '%s'." %(combined_id, gene))
    all_tiers = list(match_mapping[gene][combined_id].keys())
    selected_tier = "tier_4"
    # Pick highest tier available
    for tmp_tier in all_tiers:
        # Only look into tiers that have at least 1 match in CIViC
        skip_tier = False
        if has_support:
            if not match_mapping[gene][combined_id][tmp_tier]["matched"]:
                skip_tier = True
        else:
            if not match_mapping[gene][combined_id][tmp_tier]:
                skip_tier = True
        if skip_tier:
            continue
        # Compare the currently evaluated tier with the highest one encountered (so far)
        # NOTE: use the order of priority already defined in "sorted_tiers"
        tmp_tier_index = sorted_tiers.index(tmp_tier)
        tier_index = sorted_tiers.index(selected_tier)
        if tmp_tier_index < tier_index:
            selected_tier = tmp_tier
    # Return all variants matched for the selected tier and the provided gene+variant
    all_variants = []
    if selected_tier != "tier_4":
        if has_support:
            for tmpVar in match_mapping[gene][combined_id][selected_tier]["matched"]:
                all_variants.append(tmpVar)
        else:
            for tmpVar in match_mapping[gene][combined_id][selected_tier]:
                all_variants.append(tmpVar)
    return (selected_tier, all_variants)


# Write header to output table (add extra column when consensus drug support is provided)
# NOTE: for now, always write column for consensus drug support
def write_header(outfile, write_support=True):
    sorted_evidence_types = ["PREDICTIVE", "DIAGNOSTIC", "PROGNOSTIC", "PREDISPOSING"]
    if write_support:
        out_header = "%s\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\tCIViC_Drug_Support\t%s" %(header, "\t".join(["CIViC_" + x for x in sorted_evidence_types]))
    else:
        out_header = "%s\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\t%s" %(header, "\t".join(["CIViC_" + x for x in sorted_evidence_types]))
    outfile.write(out_header + "\n")
    return None


# Write header to output table (add extra column when consensus drug support is provided)
def write_match_multiple_annotations(match_mapping, variant_mapping, raw_mapping, input_file, outfile, data_type="SNV", has_support=True, has_ct=True, write_ct=False, write_support=True, write_complete=False):
    # NOTE: uppercase is critical for matching!
    sorted_evidence_types = ["PREDICTIVE", "DIAGNOSTIC", "PROGNOSTIC", "PREDISPOSING"]
    # Define the evidence type which contain CIViC drug predictions
    drug_evidence_type = "PREDICTIVE"
    special_cases = ["NON_SNV_MATCH_ONLY", "NON_CNV_MATCH_ONLY", "NON_EXPR_MATCH_ONLY"]
    sorted_cts = ["ct", "gt", "nct"]
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]
    varmap_entries_variant = ["name", "hgvs", "types"]

    # Use CIViCutils functionality to sanity check that all expected data and formats are correct
    check_match_before_writing(match_mapping, variant_mapping, raw_mapping, has_support, has_ct, write_ct, write_support, write_complete)

    if write_support:
        if not has_support:
            raise ValueError("Option 'write_support' cannot be selected when 'has_support'=False!")
    if write_ct:
        if not has_ct:
            raise ValueError("Option 'write_ct' cannot be selected when 'has_ct'=False!")

    # Keep track of all matches and non-matches
    exact_matches = 0           # Tier 1
    syn_matches = 0             # Tier 1b
    pos_matches = 0             # Tier 2
    no_matches = 0              # Tier 3 (includes special cases when there is no match of the same variant type returned)
    gene_not_found = 0          # Tier 4

    # Iterate through the input table once more and simultaneously write results to output table
    # Each parsed line corresponds to a different variant
    infile = open(input_file, "r")
    ignore_header = infile.readline().strip()
    for n_line,line in enumerate(infile):
        n_line = str(n_line)
        if n_line not in raw_mapping.keys():
            raise ValueError("Line %s could not be found in provided 'raw_mapping'!" %(n_line))

        tier_to_write = -1
        # gene -> variant
        genes_and_variants_to_write = {}

        # Process and prioritize CIViC results matched across all variants and genes annotated for the current SNV line
        # NOTE: when reporting results for a variant line, make sure to only report one set of results (i.e. corresponding to one single variant annotation c. + p., if avail). This is to ensure that the same CIViC variant record is not reported multiple times for many variant annotations (as they all refer to the same variant, this would bias the available evidence)
        #   1. if none have CIViC, then variant line is tier4
        #   2. if only one variant annotation has CIViC, then report the corresponding tier (only 1 since we already applied highest priority filter)
        #   3. if >1 variant annotations have CIViC:
        #    a. choose the one with the highest tier
        #    b. if >1 with the same tier, pick the one with the highest number of CIViC variant matches
        #   c. if same number of CIViC variants, choose the first variant annotation encountered (random choice)
        if data_type == "SNV":
            nested_variant_list = raw_mapping[n_line]
            selected_tier = "tier_4"
            n_variants = -1
            selected_variant = []
            for variant_list in nested_variant_list:
                if (len(variant_list) != 6):
                    raise ValueError("Must provide 6 elements to describe a SNV variant (even if some can be empty): variant_id, gene, dna, [prot], impact, [exon], ...")
                combined_id = variant_list[0]
                gene = variant_list[1]
                c_variant = variant_list[2]
                p_variant = variant_list[3]
                impact = variant_list[4]
                exon = variant_list[5]

                # NOTE: for now, skip variants which do not have a valid c.HGVS expression available (e.g. n.HGVS)
                # We can do this because we will use alternative custom function "write_match_multiple_annotations" to write results to output, and these variant annotations will not be removed from the affected columns
                if not c_variant.startswith("c."):
                    continue

                # Process and prioritize available tiers
                (tmp_tier, all_variants) = check_tier_and_matches(gene, combined_id, match_mapping, has_support)
                # Compare the currently evaluated tier with the highest one encountered (so far)
                # NOTE: use the order of priority already defined in "sorted_tiers"
                prioritize_tier = False
                tmp_tier_index = sorted_tiers.index(tmp_tier)
                tier_index = sorted_tiers.index(selected_tier)
                # Keep the highest tier available
                if tmp_tier_index < tier_index:
                    prioritize_tier = True
                # If both tiers are the same, keep the variant annotation with the highest number of CIViC variant ids matched
                elif tmp_tier_index == tier_index:
                    if len(all_variants) > n_variants:
                        prioritize_tier = True
                # If both the tier and the number of matched variants are the same, then pick the first variant annotation encountered (random choice)
                if prioritize_tier:
                    selected_tier = tmp_tier
                    selected_variant = [gene, combined_id]
                    n_variants = len(all_variants)

            if not selected_variant:
                raise ValueError("Could not select a single valid variant annotation for line %s" %(line.strip()))

            # In the end, we only select one single variant annotation per SNV line (even if >1 genes are annotated)
            tier_to_write = selected_tier
            final_gene = selected_variant[0]
            final_variant = selected_variant[1]
            # Keep track of the gene and variant annotation selected for the current line, as well as all the corresponding variant ids matched in CIViC
            genes_and_variants_to_write[final_gene] = {}
            genes_and_variants_to_write[final_gene][final_variant] = []
            # Process and prioritize available tiers
            (tmp_tier, all_variants) = check_tier_and_matches(final_gene, final_variant, match_mapping, has_support)
            for tmp_variant in all_variants:
                genes_and_variants_to_write[final_gene][final_variant].append(tmp_variant)


        # Process and prioritize CIViC results matched across all genes annotated for the current CNV line
        # NOTE: when reporting results for a CNV line, sanity check there are no duplicated genes, and report all matched CIViC information available for the highest tier (as they will always correspond to different genes, hence the retrieved CIViC records will never be duplicated for a given CNV type and tier)
        if data_type == "CNV":
            if ("details" not in raw_mapping[n_line].keys()) or ("genes" not in raw_mapping[n_line].keys()):
                raise ValueError("Must provide 2 elements to describe a CNV variant: gene, cnv, ...")
            cnv = raw_mapping[n_line]["details"]
            combined_id = cnv + "|" + n_line
            all_genes = raw_mapping[n_line]["genes"]
            # Classify all annotated genes for the current CNV according to the corresponding tier from CIViC
            tmp_cnv_results = {}
            for gene in all_genes:
                # Process and prioritize available tiers
                (tmp_tier, all_variants) = check_tier_and_matches(gene, combined_id, match_mapping, has_support)
                if tmp_tier not in tmp_cnv_results.keys():
                    tmp_cnv_results[tmp_tier] = []
                if gene in tmp_cnv_results[tmp_tier]:
                    raise ValueError("Encountered duplicated gene annotation in line %s" %(line.strip()))
                tmp_cnv_results[tmp_tier].append(gene)

            # Select the highest available tier for the current line
            all_available_tiers = list(tmp_cnv_results.keys())
            selected_tier = "tier_4"
            for tmp_avail_tier in all_available_tiers:
                # Compare the currently evaluated tier with the highest one encountered (so far)
                # NOTE: use the order of priority already defined in "sorted_tiers"
                tmp_tier_index = sorted_tiers.index(tmp_avail_tier)
                tier_index = sorted_tiers.index(selected_tier)
                if tmp_tier_index < tier_index:
                    selected_tier = tmp_avail_tier

            # Once the highest (single) tier classification has been selected, keep track of all the annotated genes which have said tier, and keep track of their corresponding variant ids matched in CIViC
            tier_to_write = selected_tier
            for final_gene in tmp_cnv_results[selected_tier]:
                if final_gene not in genes_and_variants_to_write.keys():
                    genes_and_variants_to_write[final_gene] = {}
                    genes_and_variants_to_write[final_gene][combined_id] = []
                # Process and prioritize available tiers
                (tier, all_variants) = check_tier_and_matches(final_gene, combined_id, match_mapping, has_support)
                for tmp_variant in all_variants:
                    genes_and_variants_to_write[final_gene][combined_id].append(tmp_variant)


        ## At this point, the genes and variant annotations which will be reported for the current variant line have already been parsed and selected
        gene_scores = []
        gene_var_types = []
        drug_support = []
        result_mapping = {}
        for gene in genes_and_variants_to_write.keys():
            # Expectation is that there is only 1 variant available per gene
            for variant_annotation in genes_and_variants_to_write[gene].keys():
                all_civic_variants = genes_and_variants_to_write[gene][variant_annotation]
                if tier_to_write == "tier_4":
                    if all_civic_variants:
                        raise ValueError("Unexpectedly found matched variants for a line classified as 'tier_4': %s" %(all_civic_variants))
                elif has_support:
                    this_drug_support = match_mapping[gene][variant_annotation][tier_to_write]["drug_support"]
                    if write_support:
                        for i in this_drug_support:
                            drug_support.append(i.upper())
                for variant_id in all_civic_variants:
                    # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                    # This is a dummy tag and not an actual variant record from CIViC, so skip checking in variant_mapping
                    if variant_id.upper() in special_cases:
                        # In this case, current line will be associated with tier3, but all columns will be empty with "."
                        for evidence_type in sorted_evidence_types:
                            if evidence_type not in result_mapping.keys():
                                result_mapping[evidence_type] = []
                        continue

                    civic_variant = variant_mapping[gene][variant_id]["name"]

                    gene_var_types.append(gene + ":" + civic_variant + ":" + ",".join(variant_mapping[gene][variant_id]["types"]))
                    molecular_profile_ids = set(list(variant_mapping[gene][variant_id].keys())) ^ set(varmap_entries_variant)
                    for molecular_profil_id in molecular_profile_ids:
                        gene_scores.append(gene + ":" + civic_variant + ":" + molecular_profil_id + ":" + str(variant_mapping[gene][variant_id][molecular_profil_id]["civic_score"]))                        
                        for evidence_type in sorted_evidence_types:
                            if evidence_type in variant_mapping[gene][variant_id][molecular_profil_id]["evidence_items"].keys():
                                if evidence_type not in result_mapping.keys():
                                    result_mapping[evidence_type] = []
                                write_drug = False
                                # Check whether current evidence column name corresponds to the one defined at the beginning of the function (i.e. indicates drug prediction evidence)
                                if evidence_type == drug_evidence_type:
                                    write_drug=True
                                if has_ct:
                                    check_keys(list(variant_mapping[gene][variant_id][molecular_profil_id]["evidence_items"][evidence_type].keys()), "variant_mapping", sorted_cts, matches_all=True)
                                    for ct in variant_mapping[gene][variant_id][molecular_profil_id]["evidence_items"][evidence_type].keys():
                                        if write_ct:
                                            results = write_evidences(variant_mapping[gene][variant_id][molecular_profil_id]["evidence_items"][evidence_type][ct], write_drug=write_drug, write_ct=ct, write_complete=write_complete)
                                        else:
                                            results = write_evidences(variant_mapping[gene][variant_id][molecular_profil_id]["evidence_items"][evidence_type][ct], write_drug=write_drug, write_ct=None, write_complete=write_complete)
                                        for x in results:
                                            result_mapping[evidence_type].append(gene + ":" + civic_variant + ":" + x)
                                else:
                                    check_keys_not(list(variant_mapping[gene][variant_id][molecular_profil_id]["evidence_items"][evidence_type].keys()), "variant_mapping", sorted_cts)
                                    if write_ct:
                                        raise ValueError("Option 'write_ct' cannot be selected when 'has_ct'=False!")
                                    results = write_evidences(variant_mapping[gene][variant_id][molecular_profil_id]["evidence_items"][evidence_type], write_drug=write_drug, write_ct=None, write_complete=write_complete)
                                    for x in results:
                                        result_mapping[evidence_type].append(gene + ":" + civic_variant + ":" + x)
        # NOTE: in the above case, if CIViC info is reported across multiple different variant annotations within the same line (e.g. case of the CNVs with >1 genes having CIViC info with the same tier), then the drug support will need to be recomputed to aggregate the information across all reported variants
        # Use custom function to recompute final consensus drug support for the current line
        # Now, support is aggregated across all genes and variants for which CIViC info is to be reported
        # NOTE: currently, for SNV data, only 1 single gene and variant annotation is selected per line
        # NOTE: currently, for CNV data, multiple genes annotated for the same CNV can be reported with CIViC info (only those having the same, highest tier)
        drug_support = []
        drug_support = reprocess_drug_support_across_selected_variants(genes_and_variants_to_write, match_mapping, variant_mapping, support_dict, has_support=True)
        ## Write result to output
        out_line = write_output_line(tier_to_write, line.strip(), gene_scores, gene_var_types, drug_support, result_mapping, write_support)
        # Once here, all genes and variants within the row have been classified into tiers
        if tier_to_write == "tier_4":
            gene_not_found += 1
        elif (tier_to_write == "tier_3") or (tier_to_write in special_cases):
            no_matches += 1
        elif tier_to_write == "tier_2":
            pos_matches += 1
        elif tier_to_write == "tier_1b":
            syn_matches += 1
        elif tier_to_write == "tier_1":
            exact_matches += 1

        outfile.write(out_line + "\n")
    infile.close()

    return (exact_matches, syn_matches, pos_matches, no_matches, gene_not_found)



"""
Script
"""

parser = argparse.ArgumentParser(description="Query CIViC to retrieve drug information for snvs or cnvs.")
parser.add_argument("--infile", dest="infile", required=True, help="Input table with genes and variants, needs to be tab separated.")
parser.add_argument("--outfile", dest="outfile", required=True, help="Name of the output file.")
parser.add_argument("--cancer_type_list", dest="cancer_type_list", required=True, help="Comma-separated list of accepted cancer types. Partial matches will be sought.")
parser.add_argument("--black_list", dest="black_list", required=True, help="Comma-separated list of not accepted cancer types. Partial matches will be sought, and those matched will be excluded from the list of considered evidence.")
parser.add_argument("--high_level_list", dest="high_level_list", required=True, help="Comma-separated list of high level cancer types (e.g. Cancer). Only exact matches will be sought. The user should be aware that results for high level cancer types will only be retrieved when no match for cancer specific types (i.e. --cancer_type_list) is found.")
parser.add_argument("--colname_gene", dest="colname_gene", required=True, help="Name of column containing gene symbols.")
parser.add_argument("--colname_data", dest="colname_data", required=True, help="Name of column containing variant annotations (SNV) or cnv categories (CNV).")
parser.add_argument("--data_type", dest="data_type", required=True, choices=["snv", "cnv"], default="snv", help="Type of data in infile. Possible options are: snv, cnv.")
parser.add_argument("--colname_impact", dest="colname_impact", required=False, help="Name of column containing variant impacts (for SNV only).")
parser.add_argument("--colname_exon", dest="colname_exon", required=False, help="Name of column containing variant exon information (for SNV only).")

args = parser.parse_args()

# Check input parameters
data_type = args.data_type.upper() # SNV or CNV

# Check expected data type using CIViCutils functionality
check_data_type(data_type)

# For SNV, two additional arguments are required (i.e. column names for variant impact and exon information)
if (data_type=="SNV") and (args.colname_impact is None or args.colname_exon is None):
    parser.error("--data_type 'snv' requires --colname_impact and --colname_exon.")

print("\nParameters:\n infile: %s\n outfile: %s\n data_type: %s\n colname_gene: %s\n colname_data: %s\n colname_impact: %s\n colname_exon: %s\n" %(args.infile, args.outfile, args.data_type, args.colname_gene, args.colname_data, args.colname_impact, args.colname_exon))


cancer_type_list = args.cancer_type_list
black_list = args.black_list
high_level_list = args.high_level_list
cancer_type_list = [s.strip().upper() for s in cancer_type_list.split(",")]
black_list = [s.strip().upper() for s in black_list.split(",")]
high_level_list = [s.strip().upper() for s in high_level_list.split(",")]

print("\nWhite listed cancer types: {}".format(",".join(cancer_type_list)))
print("Black listed cancer types: {}".format(",".join(black_list)))
print("High level cancer types: {}".format(",".join(high_level_list)))

# Sanity check for "empty" input lists (i.e. in the form of [""])
# Turn them into "real" empty lists for implementation purposes
if cancer_type_list == [""]:
    cancer_type_list = []
if black_list == [""]:
    black_list = []
if high_level_list == [""]:
    high_level_list = []


# Read in and process file of input SNV variants
if data_type=="SNV":
    (raw_mapping, var_data, header) = read_in_snvs_multiple_genes(args.infile)
# Read in and process file of input CNV variants
if data_type=="CNV":
    (raw_mapping, var_data, header) = read_in_cnvs_multiple_genes(args.infile)


# Already write new header into output file
outfile = open(args.outfile,"w")
write_header(outfile, write_support=True)

# Extract all genes parsed in the provided input file
genes = list(var_data.keys())

# When no genes were found in the input file, write empty files (necessary for snakemake pipeline) and exit without error
# e.g. empty input file containing only header because patient had no variants at all
if not genes:
    print("\nDid not find any genes in column \'{}\' of file {}".format(args.colname_gene, args.infile))
    outfile.close()
    sys.exit(0)


## Query CIViC for the genes of interest

print("\nTotal # genes to query: {}".format(len(genes)))
print("\nRetrieving data from CIViC...")
# Query input genes in CIViC
variant_mapping = query_civic(genes, identifier_type="entrez_symbol")


# (variant_mapping, retrieved_genes, no_variants, all_variants) = query_civic_genes(genes, identifier_type="entrez_symbol", data_type=data_type)
# retrieved = set(retrieved_genes)
# unmatched = list(set(genes) - retrieved)
# print("Found %s/%s genes in CIViC associated to %s variants. Found %s CIViC genes that had no variants available." %(len(retrieved_genes), len(genes), len(all_variants), len(no_variants)))
# print("\nGenes with no CIViC data: {}".format(",".join(unmatched)))

# Filter undesired evidences to avoid matching later on
variant_mapping = filter_civic(variant_mapping, evidence_type_not_in=["FUNCTIONAL", "ONCOGENIC"], evidence_status_in=["ACCEPTED"], var_origin_not_in=["GERMLINE"], output_empty=False)


# Match input SNV variants in CIViC, pick highest tier available per input gene+variant
# Tier hierarchy: 1 > 1b > 2 > 3 > 4
(match_mapping, matched_ids, variant_mapping) = match_in_civic(var_data, data_type=data_type, identifier_type="entrez_symbol", select_tier="highest", var_map=variant_mapping)


# Annotate matched CIViC evidences with cancer specificity of the associated diseases
disease_name_not_in = []
disease_name_in = ["bladder"]
alt_disease_names = ["solid tumor"]
annot_mapping = annotate_ct(variant_mapping, disease_name_not_in, disease_name_in, alt_disease_names)

# Filter CIViC evidences to pick only those for the highest cancer specificity available
# ct hierarchy: ct > gt > nct
annot_mapping = filter_ct(annot_mapping, select_ct="highest")

# Get custom dictionary of support from data.yml (provided within the package)
# This defines how each combination of evidence direction + clinical significance in CIViC is classified in terms of drug support (e.g. sensitivity, resistance, unknown, etc.)
support_dict = get_dict_support()

# Process drug support of the matched variants using the annotated CIViC evidences
annot_match_mapping = process_drug_support(match_mapping, annot_mapping, support_dict)

# Write to output
# Parse input file again, now checking (and prioritizing, if necessary) the available CIViC info per variant line
# Report the CT classification of each disease, and write column with the overall drug support of the match for each available CT class
(exact_matches, syn_matches, pos_matches, no_matches, gene_not_found) = write_match_multiple_annotations(annot_match_mapping, annot_mapping, raw_mapping, args.infile, outfile, data_type=data_type, has_support=True, has_ct=True, write_ct=True, write_support=True, write_complete=True)

outfile.close()

print("\nTotal # exact matches: {}".format(exact_matches))
print("Total # synonym matches: {}".format(syn_matches))
print("Total # positional matches: {}".format(pos_matches))
print("Total # gene-only matches: {}".format(no_matches))
string_not_found = "Total # rows without variant data: {}"
if data_type == "CNV":
    string_not_found = "Total # genes without variant data: {}"
print(string_not_found.format(gene_not_found))
print("---------------------")

