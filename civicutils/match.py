import sys
import os
import re

from civicutils.utils import check_arguments, check_argument, translate_aa, check_is_dict, check_is_list, check_is_str, check_keys, check_keys_not, check_tier_selection, parse_input, check_empty_input, check_dict_entry, check_is_chgvs, check_is_phgvs, check_identifier_type, check_data_type, check_is_bool, uppercase_list, check_is_none, check_logfc


def civic_name_to_hgvs(var_name):
    """
    Given a single CIViC variant name, extract potential HGVS annotations by parsing and modifying the name while using prior knowledge about CIViC naming conventions.
    Function only applicable to SNV data.
    :param var_name:	Name of CIViC variant record.
    :return:		List of corresponding potential HGVS expressions generated based on the provided variant name.
    """
    # Sanity check of provided arguments
    check_argument(var_name, "var_name")
    check_is_str(var_name, "var_name")
    # NOTE: uppercase is critical for the match!
    var_name = var_name.upper()

    hgvs_strings = []
    # 1) HGVS 1-letter protein code, including stop codons and general variants (e.g. V600E, *600E, V600* or V600)
    # NOTE: variant names like T193C (refer to DNA) would be translated to proteins
    if re.match(r'([A-Z*])(\d+)([A-Z*]?)($|\s\()', var_name):
        pos = re.match(r'([A-Z*])(\d+)([A-Z*]?)($|\s\()', var_name).groups()
        aa1 = pos[0]
        npos = pos[1]
        aa2 = pos[2]
        ## Generate the HGVS string by translating 1-letter aa codes into 3-letter

        # Special case of protein extensions where "*" is represented by "Ter" in input table
        if aa1 == "*":
            aa1_new = "TER"
        else:
            # Translate 1-letter aa code to 3-letter code
            aa1_new = translate_aa(aa1)
        # Check for general variants (second aa will be "")
        if aa2:
            # Special cases like p.Ter600Ter happen in input table, so in theory,
            # *600* should also be possible (however does not seem to happen in CIViC)
            if (aa1_new == "TER") and (aa2=="*"):
                aa2_new = "TER"
            else:
                # Translate 1-letter aa code to 3-letter code
                # In cases where aa2="*", it will remain as "*" (using Ter is a special case when aa1=Ter)
                aa2_new = translate_aa(aa2)
        else:
            aa2_new = ""
        # Sanity check that translated strings correspond to valid aa codes
        if (aa1_new is not None) and (aa2_new is not None):
            # Construct new protein annotation
            new_annotation = "P." + aa1_new + npos + aa2_new
            hgvs_strings.append(new_annotation)

    # 2) Extract embedded c. annotation if available i.e. "(c.XXX)"
    # Use re.search to allow for matches beyond the start of the string
    if re.search(r'\((C\..+?)\)', var_name):
        # r = re.findall(r'\((c\..+?)\)', var_name)
        annot = re.search(r'\((C\..+?)\)', var_name).groups()[0]
        hgvs_strings.append(annot)

    # 3) Frameshifts (e.g. T157FS or T157MFS)
    if re.match(r'([A-Z])(\d+)([A-Z]?)FS', var_name):
        pos = re.match(r'([A-Z])(\d+)([A-Z]?)FS', var_name).groups()
        aa = pos[0]
        npos = pos[1]
        # Translate aa in 1st position and if ok, write frameshift in short form
        aa_new = translate_aa(aa)
        if aa_new is not None:
            annot = "P." + aa_new + npos + "FS"
            hgvs_strings.append(annot)

    return hgvs_strings


def civic_hgvs_to_input(civic_hgvs):
    """
    Given a single CIViC HGVS expression, parse and modify it to ensure it complies with the HGVS guidelines expected in the input table.
    Function only applicable to SNV data. NOTE: currently, only modifications of p. HGVS expressions are supported.
    :param civic_hgvs:	One single HGVS expression retrieved from CIViC.
    :return:		Modified HGVS expression complying with user's expectations.
    """
    # Sanity check of provided arguments
    check_argument(civic_hgvs, "civic_hgvs")
    check_is_str(civic_hgvs, "civic_hgvs")
    # NOTE: uppercase is critical for the match!
    civic_hgvs = civic_hgvs.upper()
    ## The following "special hgvs" cases should be mutually exclusive
    ## i.e. a single HGVS expression should match only 1 of the cases
    new_annot = None

    ## Frameshift
    # Input table seems to use short form for frameshift annotations (i.e. p.Glu55fs) while CIViC does not
    if re.match(r'(P\.[A-Z]+[0-9]+)[A-Z]+FS.*', civic_hgvs):
        new_annot = re.sub(r'(P\.[A-Z]+[0-9]+)[A-Z]+FS.*', r'\1FS', civic_hgvs)
        # Only return new HGVS expression if something changed
        if new_annot != civic_hgvs:
            return new_annot

    ## Nonsense mutations (i.e. gain of stop codon)
    # In general, CIViC seems to use "Ter" for stop codons, only 1 exception so far (p.F76Lfs*56)
    # Input table uses "*" for stop codons; exception in protein extensions (e.g. p.Ter370Tyrext*?)
    # "*" is also used for nucleotide numbering in c.,n., but not relevant here (only p.)
    if re.match(r'(P\.[A-Z]+[0-9]+)TER', civic_hgvs):
        new_annot = re.sub(r'(P\.[A-Z]+[0-9]+)TER', r'\1*', civic_hgvs)
        # Only return new HGVS expression if something changed
        if new_annot != civic_hgvs:
            return new_annot

    ## Loss of stop codon
    # In CIViC, a stop codon loss is expressed as p.Ter214Cys
    # In HGVS, this is enconded as an extension (e.g. p.Ter370Tyrext*?)
    # In this case, we need to change annotation from input table (generate_civic_matchStrings)

    ## Silent mutations (i.e. no aa change)
    # CIViC uses a "=" in the 2nd position to represent no change (e.g. p.Pro61=)
    # In the input table, the 2nd aa is specified (e.g. p.Pro61Pro)
    if re.match(r'P\.([A-Z]+)([0-9]+)=', civic_hgvs):
        new_annot = re.sub(r'P\.([A-Z]+)([0-9]+)=', r'P.\1\2\1', civic_hgvs)
        # Only return new HGVS expression if something changed
        if new_annot != civic_hgvs:
            return new_annot

    return new_annot


def extract_p_start(p_hgvs):
    """
    Given a single HGVS p. annotation (of the form p.Pro61Cys...), extract the start of the string composed of the original/affected aminoacid and corresponding position (i.e. p.Pro61).
    Function only applicable to SNV data.
    :param p_hgvs:	One single p.HGVS expression (complete and using 3-letter aminoacid codes). 
    :return:		Modified p.HGVS expression containing only the first aminoacid and the affected position.
    """
    # Sanity check of provided arguments
    check_argument(p_hgvs, "p_hgvs")
    check_is_str(p_hgvs, "p_hgvs")
    # NOTE: uppercase is critical for the match!
    p_hgvs = p_hgvs.upper()

    start_annot = None
    if re.match(r'(P\.[A-Z]+[0-9]+)', p_hgvs):
        start_annot = re.match(r'(P\.[A-Z]+[0-9]+)', p_hgvs).groups()[0]

    return start_annot


def check_general_variant(var_name):
    """
    Check whether a variant name retrieved from CIViC corresponds to a collection of multiple variants (e.g. V600).
    Function only applicable to SNV data.
    :param var_name:	Name of the CIViC variant record.
    :return:		Boolean indicating if the variant name points to a general collection of variants.
    """
    # Sanity check of provided arguments
    check_argument(var_name, "var_name")
    check_is_str(var_name, "var_name")
    # NOTE: uppercase is critical for the match!
    var_name = var_name.upper()

    is_general = False
    if re.match(r'[A-Z]\d+($|\s\()', var_name):
        is_general = True

    return is_general


def cnv_is_exon_string(var_name):
    """
    Given a single CIViC variant name, return whether it corresponds to a CNV variant record related to exons.
    For this, check whether the variant name matches special CNV exon cases known to be contained in CIViC (e.g. EXON 1-2 DELETION, EXON 5 DELETION, 3' EXON DELETION...).
    :param var_name:	Name of the CIViC variant record.
    :return:		Boolean indicating if the variant corresponds to a CNV related to exons.
    """
    # NOTE: uppercase is critical for the match!
    exon_strings = ["^EXON [0-9-]+ DELETION$", "^[35\']+ EXON DELETION$", "^EXON [0-9-]+ SKIPPING MUTATION$"]

    # Sanity check of provided arguments
    check_argument(var_name, "var_name")
    check_is_str(var_name, "var_name")
    var_name = var_name.upper()

    is_exon = False
    for exon_string in exon_strings:
        if re.search(exon_string, var_name):
            is_exon = True

    return is_exon


def expr_is_exon_string(var_name):
    """
    Given a single CIViC record name, return whether it corresponds to an expression record related to exons, and the type of expression change.
    For this, check whether the variant name matches special EXPRESSION exon cases known to be contained in CIViC (e.g. EXON 1-2 EXPRESSION, EXON 5 OVEREXPRESSION...).
    :param var_name:	Name of the CIViC variant record.
    :return:		Tuple of one boolean and one string; first indicates if the record corresponds to a differential expression record, and second indicates what type of expression change the record relates to (i.e. expression, overexpression or underexpression).
    """
    # Sanity check of provided arguments
    check_argument(var_name, "var_name")
    check_is_str(var_name, "var_name")
    var_name = var_name.upper()

    expr_type = ""
    is_exon = False

    # NOTE: uppercase is critical for the match!
    if re.search("^EXON [0-9-]+ EXPRESSION$", var_name):
        is_exon = True
        expr_type = "EXPRESSION"
    elif re.search("^EXON [0-9-]+ OVEREXPRESSION$", var_name):
        is_exon = True
        expr_type = "OVEREXPRESSION"
    elif re.search("^EXON [0-9-]+ UNDEREXPRESSION$", var_name):
        is_exon = True
        expr_type = "UNDEREXPRESSION"

    return (is_exon, expr_type)


def civic_return_all_snvs(gene_data):
    """
    Given all CIViC variant records associated to a given gene, return the set of variants detected to correspond to SNV records.
    For this, remove the most common CNV and EXPRESSION ids existing in CIViC from the complete set of the gene's variants.
    :param gene_data:	Nested dictionary corresponding to a subset of 'var_map' for one single CIViC gene. See README for more details about the structure of 'var_map'.
    :return:		List of variant CIViC ids corresponding to SNV records.
    """
    # Sanity check of provided arguments
    check_argument(gene_data, "gene_data")
    check_is_dict(gene_data, "gene_data")

    # Get all variant names classified as CNV or EXPRESSION (if any)
    # These functions will also attempt matching of variant name to common CNV and EXPRESSION names related to exons
    cnv_ids = civic_return_all_cnvs(gene_data)
    expr_ids = civic_return_all_expr(gene_data)
    # All variant names not matching a CNV or EXPRESSION will be returned
    matches = []
    for var_id in list(gene_data.keys()):
        # Skip variant ids found to match a CNV or EXPRESSION record
        if (var_id in cnv_ids) or (var_id in expr_ids):
            continue
        if var_id not in matches:
            matches.append(var_id)

    return matches


def civic_return_all_cnvs(gene_data):
    """
    Given all CIViC variant records associated to a given gene, return all those variant ids that correspond to CNV records.
    For this, check whether the corresponding gene's variant names match the most common CNV names known to be contained in CIViC.
    Additionally, consider other special CNV cases present in CIViC (e.g. EXON 1-2 DELETION, EXON 5 DELETION, 3' EXON DELETION...).
    :param gene_data:	Nested dictionary corresponding to a subset of 'var_map' for one single CIViC gene. See README for more details about the structure of 'var_map'.
    :return:		List of variant CIViC ids corresponding to CNV records.
    """
    # Common CNV record names in CIViC
    # NOTE: uppercase is critical for the match!
    cnv_names = ["AMPLIFICATION", "DELETION", "LOSS", "COPY NUMBER VARIATION"]

    # Sanity check of provided arguments
    check_argument(gene_data, "gene_data")
    check_is_dict(gene_data, "gene_data")

    # All matched variant ids will be returned
    matches = []
    for var_id in list(gene_data.keys()):
        # Retrieve the name assigned to the current variant id
        check_dict_entry(gene_data[var_id], "gene_data", "name", "name")
        # NOTE: uppercase is critical for the match!
        var_name = gene_data[var_id]["name"].upper()
        # Also attempt matching of variant name to common CNV names related to exons
        is_exon = cnv_is_exon_string(var_name)
        if (var_name in cnv_names) or is_exon:
#             ## Special case for "LOSS": only considered synonym of "DELETION" when a certain variant type is present
#             if var_name == "LOSS" and ("TRANSCRIPT_ABLATION" not in gene_data[var_name]["types"]):
#                 continue
            if var_id not in matches:
                matches.append(var_id)

    return matches


def civic_return_all_expr(gene_data):
    """
    Given all CIViC variant records associated to a given gene, return all those variant ids that correspond to EXPRESSION records.
    For this, check whether the corresponding gene's variant names match the most common EXPRESSION names known to be contained in CIViC.
    Additionally, consider other special EXPRESSION cases present in CIViC (e.g. EXON 1-2 OVEREXPRESSION, EXON 5 UNDEREXPRESSION, ...).
    :param gene_data:	Nested dictionary corresponding to a subset of 'var_map' for one single CIViC gene. See README for more details about the structure of 'var_map'.
    :return:		List of variant CIViC ids corresponding to Expression records.
    """
    # Common EXPRESSION record names in CIViC
    # NOTE: uppercase is critical for the match!
    expr_names = ["OVEREXPRESSION", "UNDEREXPRESSION", "EXPRESSION"]

    # Sanity check of provided arguments
    check_argument(gene_data, "gene_data")
    check_is_dict(gene_data, "gene_data")

    # All matched variant ids will be returned
    matches = []
    for var_id in list(gene_data.keys()):
        # Retrieve the name assigned to the current variant id
        check_dict_entry(gene_data[var_id], "gene_data", "name", "name")
        # NOTE: uppercase is critical for the match!
        var_name = gene_data[var_id]["name"].upper()
        # Also attempt matching of variant name to common EXPRESSION names related to exons
        (is_exon, expr_type) = expr_is_exon_string(var_name)
        if (var_name in expr_names) or is_exon:
            if var_id not in matches:
                matches.append(var_id)

    return matches


def civic_match_strings(var_name, hgvs_expressions, data_type):
    """
    Given a CIViC variant name and its corresponding HGVS expressions (if available), generate a list of potential strings that will be used to match the CIViC record to an input variant.
    Function applicable to all data types (SNV, CNV and EXPR); matching framework is different depending on the given data type.
    For data types 'CNVs' and 'EXPR', matching is not based at all on HGVS expressions, and exclusively the CIViC variant name is returned.
    :param var_name:		Name of CIViC variant record.
    :param hgvs_expressions:	List of HGVS expressions retrieved from CIViC and corresponding to the provided CIViC variant record. Expectation is that this argument will only be non-empty for SNVs.
    :param data_type:		['SNV', 'CNV', 'EXPR']
                        	SNV:   Variant corresponds to genomic single nucleotide mutations and insertions/deletions.
                        	CNV:   Variant corresponds to genomic copy number alterations.
                        	EXPR:  Variant corresponds to differential gene expression data.
				Type of data being queried.
    :return:			List of strings used to match the CIViC variant record to the input variants.
    """
    # Sanity check of provided argument
    check_data_type(data_type)
    check_argument(var_name, "var_name")
    check_is_str(var_name, "var_name")
    # NOTE: uppercase is critical for the match!
    var_name = var_name.upper()
    hgvs_expressions = uppercase_list(hgvs_expressions, "hgvs_expressions")

    match_strings = []
    # For data types "CNV" and "EXPR", only the last step 5) is executed, since variant matching is not done at the HGVS level but only at the level of the record name
    if data_type == "SNV":
        # 1) First, remove reference from annotation (i.e. "transcript_id:")
        # This step will be skipped for CNVs
        for x in hgvs_expressions:
            # CIViC HGVS expressions are of the form reference + single annotation
            # e.g. [NM_007313.2:c.1001C>T, NP_005148.2:p.Thr315Ile]
            new = x.split(":")[-1].upper()
            if new not in match_strings:
                match_strings.append(new)
                # 2) Second, generate modified HGVS strings that comply with input table format
                new_annot = civic_hgvs_to_input(new)
                # Resulting annotation will be None unless something was modified (special cases only)
                if (new_annot is not None) and (new_annot not in match_strings):
                    match_strings.append(new_annot)
        # 3) Third, generate a list of potential HGVS strings based on the actual CIViC variant name
        potential_hgvs = civic_name_to_hgvs(var_name)
        for x in potential_hgvs:
            # Sometimes, duplicated HGVS annotations will get generated (e.g. V600E already has p.Val600Glu)
            if x not in match_strings:
                match_strings.append(x)
        # 4) Last, add potential positional matches for already existing strings (will only affect p. annotations)
        for x in match_strings:
            start = extract_p_start(x)
            # Only p. annotations will return a positional string (e.g. p.Val600Glu -> p.Val600)
            if (start is not None) and (start not in match_strings): 
                match_strings.append(start)
    # 5) Also, add CIViC variant name to allow for matching using input "descriptional" strings (e.g. EXON 15 MUTATION)
    # For data types "CNV" and "EXPR", this string is the only relevant one for matching to CIViC variants
    match_strings.append(var_name)

    return match_strings


def input_match_strings(var_annotations, data_type, impact_annots=[], exon_annots=[]):
    """
    Given a set of one or more genomic alterations (SNV or CNV), generate a list of potential strings that will be used to match the variant in CIViC e.g. EXON 15 MUTATION, AMPLIFICATION.
    Function only appplicable to data types 'SNV' and 'CNV' (the strings generated are different in each case).
    :param var_annotations:	List of genomic alterations (either SNV/InDels or CNVs, depending on the given data type). In the case of type 'SNV', the list can contain HGVS expressions.
    :param data_type:		['SNV', 'CNV', 'EXPR']
                        	SNV:   Variants correspond to genomic single nucleotide mutations and insertions/deletions.
                        	CNV:   Variants correspond to genomic copy number alterations.
                        	EXPR:  Variants correspond to differential gene expression data.
				Type of data being queried.
    :param impact_annots:	Optional (only supported for data_type='SNV'). When provided, variant impacts will be used to generate additional potential record names to match CIViC evidence. See README for more details about these annotations.
    :param exon_annots:		Optional (only supported for data_type='SNV'). When provided, exon information will be used to generate additional potential record names to match CIViC evidence. Note that if exon annotations are provided, then variant impacts must also be provided, and there must be a 1-1 correspondance in the order between the listed variant impacts and the provided exon/intron annotations (as impacts are used to determine if variant is intronic or exonic). See README for more details about these annotations.
    :return:			List of strings that can be used to match input variants to CIViC variant records.
    """
    # Sanity check of provided arguments
    check_data_type(data_type)
    check_argument(var_annotations, "var_annotations")
    # NOTE: uppercase is critical for the match!
    var_annotations = uppercase_list(var_annotations, "var_annotations")
    impact_annots = uppercase_list(impact_annots, "impact_annots")
    exon_annots = uppercase_list(exon_annots, "exon_annots")

    match_strings = []
    # Used for SNV: Keep track of whether a string corresponds to an exact (True) or positional match (False)
    is_exact = []
    # Used for SNV: is_true_exact can be used together with is_exact to distinguish between a true exact match
    # (is_exact=True, is_true_exact=True) corresponding to input HGVS strings, or a more general match
    # (is_exact=True, is_true_exact=False) corresponding to descriptive terms e.g. EXON 1 MUTATION with lower preference
    is_true_exact = []

    ## a) Given a set of HGVS expressions, impacts and exon information for a single variant, generate a list of
    ## potential strings that will be used to match the variant in CIViC e.g. EXON 15 MUTATION, TRUNCATING FRAMESHIFT
    if data_type == "SNV":
        # 1) First, add all (unique) HGVS annotations gathered from input table
        # If matched, they will correspond to exact matches (True in is_exact)
        for var_annot in var_annotations:
            if var_annot not in match_strings:
                match_strings.append(var_annot)
                is_exact.append(True)
                is_true_exact.append(True)
                # Special case for protein extensions (e.g. p.Ter130Tyrext*?) in input table
                # Subset string to match CIViC's format (p.Ter130Tyr)
                if re.match(r'(P\.TER[0-9]+[A-Z]+)EXT', var_annot):
                    new_var_annot = re.match(r'(P\.TER[0-9]+[A-Z]+)EXT', var_annot).groups()[0]
                    if new_var_annot not in match_strings:
                        match_strings.append(new_var_annot)
                        is_exact.append(True)
                        is_true_exact.append(True)
        # 2) Second, add potential positional matches for already existing strings (will only affect p. annotations)
        # If matched, they will correspond to positional matches (False in is_exact)
        for x in match_strings:
            start = extract_p_start(x)
            # Only p. annotations will return a positional string (e.g. p.Val600Glu -> p.Val600)
            if (start is not None) and (start not in match_strings):
                match_strings.append(start)
                is_exact.append(False)
                is_true_exact.append(False)

        # 3) Last, add potential variant synonyms for matching to CIViC's record names (e.g. EXON 15 MUTATION)
        # If matched, they will correspond to exact matches (is_exact=True) but not to true exact matches (is_true_exact=False)
        ## Sanity check that both lists (impact and exon annotations) are the same length
        new_tags = []
        ## Always include "MUTATION" as a potential variant tag
        new_tags.append("MUTATION")

        # Iterate impact information to generate and include additional variant tags
        for impact in impact_annots:
            # Introduce sanity check on variant impact being available, as this field is optional
            if not impact:
                continue
            ## Generate potential additional variant tags based on the variant impact
            ## (can be contained in impact e.g. "frameshift_variant&stop_gained")
            if re.search("3_PRIME_UTR_VARIANT",impact):
                new_tags.append("3' UTR MUTATION")

            if re.search("5_PRIME_UTR_VARIANT",impact):
                new_tags.append("5' UTR MUTATION")

            if re.search("STOP_GAINED",impact):
                new_tags.append("TRUNCATING MUTATION")

            if re.search("FRAMESHIFT_VARIANT",impact):
                new_tags.append("FRAMESHIFT MUTATION")

        # Iterate exon information to generate and include additional variant tags
        for i,exon in enumerate(exon_annots):
            # Introduce sanity check on variant exon being available, as this field is optional
            # Happens with impacts like e.g. downstream_gene_variant, upstream_gene_variant, intergenic_region
            if not exon:
                continue
            # When provided, exon info must be matched with impact information, to indicate if variant is exonic or intronic
            # Sanity check that both lists (impact and exon annotations) are the same length
            if len(impact_annots) != len(exon_annots):
                raise ValueError("Provided 'impact_annots' and 'exon_annots' are not of identical length. Please provide one variant impact annotation for each available exon annotation.")
            ## Based on rank (exon/intron) information, generate an additional tag
            rank = exon.split("/")[0]
            # After checking some examples, it seems that intron information is associated to variant impacts "intron_variant"
            # and "sequence_feature" (can be contained in impact e.g. "splice_donor_variant&intron_variant")
            # NOTE: any variant impact other than the above two will generate a tag "EXON XX MUTATION", even "synonymous" mutations
            if re.search("INTRON_VARIANT", impact_annots[i]) or re.search("SEQUENCE_FEATURE", impact_annots[i]):
                new_tags.append("INTRON " + rank + " MUTATION")
            else:
                new_tags.append("EXON " + rank + " MUTATION")
                ## As of 13/02/2019, only one case of variant name "EXON XX FRAMESHIFT" (gene CALR)
                if re.search("FRAMESHIFT_VARIANT",impact_annots[i]):
                    new_tags.append("EXON " + rank + " FRAMESHIFT")

        ## Add unique tags and assign corresponding information about type of match
        for tag in new_tags:
            if tag not in match_strings:
                match_strings.append(tag)
                is_exact.append(True)
                is_true_exact.append(False)

    # var_annotations can contain one or multiple elements which correspond to CNV categories
    if data_type == "CNV":
        new_tags = []
        # Iterate individual input annotations and attempt to match in CIViC using known CNV record names used in the database
        for var_annot in var_annotations:
            # CIViC seems to consider that GAIN and AMP are the same CNV
            if (var_annot == "AMPLIFICATION") or (var_annot == "AMP") or (var_annot == "GAIN") or (var_annot == "DUPLICATION") or (var_annot == "DUP"):
                new_tags.append("AMPLIFICATION")
            # CIViC seems to consider that DELETION and LOSS are the same CNV
            # Both CNV categories "DELETION" and "LOSS" occur in CIViC
            elif (var_annot == "DELETION") or (var_annot == "DEL") or (var_annot == "LOSS"):
                new_tags.append("DELETION")
                new_tags.append("LOSS")
            new_tags.append("COPY NUMBER VARIATION")
        for tag in new_tags:
            if tag not in match_strings:
                match_strings.append(tag)
                is_exact.append(True)
                is_true_exact.append(True)

    return (match_strings, is_exact, is_true_exact)


def get_expression_strings(gene, logfc):
    """
    Given a differentially expressed gene and its corresponding logFC value, generate a list of potential strings that will be used to match the alteration in CIViC.
    Possible tag combinations are either OVEREXPRESSION and EXPRESSION (when logFC>0), or UNDEREXPRESSION and EXPRESSION (when logFC<0).
    :param gene:	One single gene ID (Entrez symbol, Entrez ID or CIViC ID) that was found to be differentially expressed.
    :param logfc:	Corresponding log fold-change value for the gene, resulting from the differential expression analysis.
    :return:		List of strings used to match input differentially expressed genes and CIViC records.
    """
    # Sanity check of provided arguments
    check_arguments([gene, logfc], ["gene", "logfc"])
    check_is_str(gene, "gene")
    # NOTE: uppercase is critical for the match!
    gene = gene.upper()
    # Sanity checks about non-empty number and sign
    logfc = check_logfc(logfc, gene)

    # Expression tags will be used to match expression records in CIViC
    # expr_change = ""
    match_strings = []
    # Generate expression tags to be able to match record in CIViC
    # As of 29/09/2019, only relevant records in CIViC are: "OVEREXPRESSION", "EXPRESSION", "UNDEREXPRESSION" (ignore remaining)
    if logfc > 0:
        match_strings.append("OVEREXPRESSION")
        # expr_change = "OVEREXPRESSION"
    elif logfc < 0:
        match_strings.append("UNDEREXPRESSION")
        # expr_change = "UNDEREXPRESSION"
    else:
        raise ValueError("Invalid logFC = '%s' for gene '%s'. Only differentially expressed genes are valid." %(logfc, gene))

    # EXPRESSION (directionless) tag will also be taken into account for matching to CIViC records
    match_strings.append("EXPRESSION")

    # Special case for gene "CDKN2A": synonym "p16" is used in CIViC records (e.g. "p16 EXPRESSION")
    # Include these cases as well to be able to match them
    new_tags = []
    if gene == "CDKN2A":
        for this_tag in match_strings:
            new_tag = "P16 " + this_tag
            if new_tag not in match_strings:
                match_strings.append(new_tag)

    return match_strings


def match_variants_in_civic(gene, variants, var_map, data_type, impacts=[], exons=[]):
    """
    Match input variants from a given gene (either SNVs/InDels or CNVs) to provided variant-level CIViC records.
    An exhaustive search of potential matches between input and CIViC is performed by this function, as there could be multiple hits matches in cases of redundancy (e.g. E55FS and E55RFSTER11 from CIViC both translate to p.Glu55fs), as well as multiple positional matches (note that for SNVs, so-called 'general' variants are prioritized over other potential positional hits, e.g. V600). Also, in case of having perfect and non-perfect matches available for SNVs (e.g. V600E and EXON 15 MUTATION), the former are prioritized over the second.
    :param gene:	One single gene identifier (Entrez symbol, Entrez ID or CIViC ID).
    :param variants:	List of variant annotations associated to the provided gene, to be matched to CIViC data (HGVS, synonym descriptive terms e.g. EXON 15 MUTATION and positional strings e.g. p.Val600).
    :param var_map:	Nested dictionary of genes (must use same type of identifier as 'gene') and variants retrieved from CIViC. See README for more details about the specific structure.
    :param data_type:	['SNV', 'CNV']
                        SNV:   Variants correspond to genomic single nucleotide mutations and insertions/deletions.
                        CNV:   Variants correspond to genomic copy number alterations.
			Type of data being queried.
    :param impacts:	Optional (only supported for data_type='SNV'). List of impact annotations associated to the input variant. When provided, variant impacts will be used to generate additional potential record names to match CIViC evidence. See README for more details about these annotations.
    :param exons:	Optional (only supported for data_type='SNV'). List of exon/intron annotations associated to the input variant. When provided, exon/intron information will be used to generate additional potential record names to match CIViC evidence. Note that if exon annotations are provided, then variant impacts must also be provided, and there must be a 1-1 correspondance in the order between the listed variant impacts and the provided exon/intron annotations (as impacts are used to determine if variant is intronic or exonic). See README for more details about these annotations.
    :return:		Nested dictionary with fixed structure containing all tier categories and corresponding list of CIViC variant matches found in each case (if any). See README for more details about the returned dictionary (i.e. 'match_map').
    """
    # Returns dictionary with the following structure
    match = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}
    # Sanity check of provided arguments
    check_data_type(data_type)
    check_arguments([gene, variants], ["gene", "variants"])
    check_is_str(gene, "gene")
    check_is_list(variants, "variants")
    check_is_dict(var_map, "var_map")
    check_is_list(impacts, "impacts")
    check_is_list(exons, "exons")
    # NOTE: uppercase is critical for the match!
    gene = gene.upper()
    variants = uppercase_list(variants, "variants")
    impacts = uppercase_list(impacts, "impacts")
    exons = uppercase_list(exons, "exons")

    # Function for generating additional synonym strings for input annotations (e.g. EXON 15 MUTATION)
    # is_exact is a list of equal length to input_strings, indicating whether a given string corresponds to an exact (True) or positional (False) match
    # variants must refer to the same variant, e.g. "c." and "p." annotations of the same genomic aberration
    # impacts and exons can be empty. When exons is provided, then impacts must be provided as well and elements must have a 1-1 correspondance with the impacts.
    (input_strings, is_exact, is_true_exact) = input_match_strings(variants, data_type, impacts, exons)

    ## If gene is in CIViC, possible tier levels are 1,2,3
    if gene in var_map.keys():
        for var_id in var_map[gene].keys():
            check_dict_entry(var_map[gene][var_id], "var_map", "name", "name")
            variant_name = var_map[gene][var_id]["name"]
            check_dict_entry(var_map[gene][var_id], "var_map", "hgvs", "hgvs")
            hgvs_expressions = var_map[gene][var_id]["hgvs"]

            # Generate list of strings that will be used to match input variants to CIViC database
            # Returned list always has at least length=1 (in this case, containing only variant name)
            # For data_type=CNV, variant matching is not based on HGVS, so match_strings will only contain the variant name
            civic_strings = civic_match_strings(variant_name, hgvs_expressions, data_type)

            # Iterate input annotations strings and attempt match to any CIViC variant string
            # The position of each input string corresponds to the type of match (true exact, synonym exact or positional)
            for indx,input_annot in enumerate(input_strings):
                if input_annot in civic_strings:
                    # Determine type of match using its position and store accordingly
                    # For CNV, all matches will be exact and true exact matches (either variant is in CIViC or not)
                    if is_exact[indx]:
                        # For SNV: give preference to truly exact matches over descriptive synonym terms
                        # e.g. report V600E over EXON 15 MUTATION
                        if is_true_exact[indx]:
                            # For SNV: True exact match (e.g. V600E)
                            # For CNV: all exact matches will also be true exact matches (either CNV is in CIViC or not)
                            # There could be 0,1,>1 exact matches (e.g. DELETION + LOSS + COPY NUMBER for cnvs)
                            if var_id not in match["tier_1"]:
                                match["tier_1"].append(var_id)
                        else:
                            # For SNV: descriptive term match (e.g. EXON 15 MUTATION, FRAMESHIFT MUTATION)
                            # There could be 0,1,>1 synonym matches
                            if var_id not in match["tier_1b"]:
                                match["tier_1b"].append(var_id)
                    else:
                        # Positional match
                        # There could be 0,1,>1 positional matches
                        if var_id not in match["tier_2"]:
                            match["tier_2"].append(var_id)
                # For CNV: When there is no exact match for a "DELETION", also consider special CIViC CNV records related to exons (these will be positional matches)
                else:
                    if (data_type == "CNV") and (input_annot == "DELETION"):
                        # In CNV case, civic_strings corresponds to the CIViC variant name (single string)
                        for temp_string in civic_strings:
                            # Look for special cases like e.g. "EXON 5 DELETION", "EXON 1-2 DELETION" or "3' EXON DELETION"
                            is_exon = cnv_is_exon_string(temp_string)
                            if is_exon:
                                # These special cases are accounted for as positional matches (tier 2)
                                # There could be 0,1,>1 positional matches
                                if var_id not in match["tier_2"]:
                                    match["tier_2"].append(var_id)

        # Once all CIViC variants have been iterated, determine final tier and corresponding matched variants

        # For CNV: either exact or positional matches will occurr due to the implementation design. Positional matches will correspond to EXON records (e.g. EXON 1-2 DELETION, 3' EXON DELETION..)
        # For SNV: if there are positional matches, check for preferential positional matches, i.e. general variants (like V600)
        if match["tier_2"] and (data_type == "SNV"):
            for tmp_id in match["tier_2"]:
                # Variant was already matched, so it was already checked for entry "Name" existing
                tmp_name = var_map[gene][tmp_id]["name"]
                is_general = check_general_variant(tmp_name)
                if is_general:
                    # Stop as soon as a general variant is found and report only this
                    match["tier_2"] = [tmp_id]
                    break
        # If no match was found, then tier case is 3, and return all relevant variants
        if not (match["tier_1"] or match["tier_1b"] or match["tier_2"]):
            # For SNV: when input variant could not be matched in CIViC (tier3), return all CIViC variants associated to the given gene but that do not correspond to a CNV or EXPRESSION related variant
            if data_type == "SNV":
                match["tier_3"] = civic_return_all_snvs(var_map[gene])
                # NOTE: Sanity check for situation where gene does not contain any variant records for the requested data type (e.g. all are CNVs)
                # Use dummy variant name to keep track of these situations
                if not match["tier_3"]:
                    match["tier_3"] = ["NON_SNV_MATCH_ONLY"]

            # For CNV: when input cnv could not be matched in CIViC (tier3), return all CIViC cnvs associated the given gene (if any)
            # i.e. not all CIViC records are returned but only those that are matched to a "CNV" tags
            if data_type == "CNV":
                match["tier_3"] = civic_return_all_cnvs(var_map[gene])
                # NOTE: Sanity check for situation where gene does not contain any variant records for the requested data type (e.g. all are SNVs)
                # Use dummy variant name to keep track of these situations
                if not match["tier_3"]:
                    match["tier_3"] = ["NON_CNV_MATCH_ONLY"]

    # If gene is not in CIViC, tier level is 4 and match will be empty
    else:
        # Sanity check that no other match should have been found
        if (match["tier_1"] or match["tier_1b"] or match["tier_2"] or match["tier_3"]):
            raise ValueError("Encountered unexpected tier case!")
        # No variants will be reported in this case (as the information is not available)
        match["tier_4"] = True
    return match


def match_expression_in_civic(gene, expression_strings, var_map):
    """
    Match input differential expression data for a given gene to provided CIViC records.
    An exhaustive search of potential matches between input and CIViC is performed by this function, as there could be multiple hits (e.g. EXPRESSION and OVEREXPRESSION from CIViC both can be matched to an input gene with logFC>0). Note that tier1b and tier2 are not supported for this kind of data.
    :param gene:		One single gene identifier (Entrez symbol, Entrez ID or CIViC ID).
    :param expression_strings:	List of expression annotations associated to the provided gene, to be matched to CIViC data (i.e. synonym descriptive terms e.g. OVEREXPRESSION, UNDEREXPRESSION or even EXON 15 UNDEREXPRESSION).
    :param var_map:		Nested dictionary of genes (must use same type of identifier as 'gene') and variant-level records retrieved from CIViC. See README for more details about the specific structure.
    :return:			Nested dictionary with fixed structure containing all tier categories and corresponding list of CIViC variant matches found in each case (if any). See README for more details about the returned dictionary (i.e. 'match_map').
    """
    # NOTE: tier1b and tier2 are not supported for this kind of data
    match = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}
    # Sanity check of provided arguments
    check_arguments([gene, expression_strings], ["gene", "expression_strings"])
    check_is_str(gene, "gene")
    check_is_list(expression_strings, "expression_strings")
    check_is_dict(var_map, "var_map")
    # NOTE: uppercase is critical for the match!
    gene = gene.upper()
    expression_strings = uppercase_list(expression_strings, "expression_strings")

    # If gene is in CIViC, only possible tier level is 1
    if gene in var_map.keys():
        for var_id in var_map[gene].keys():
            # Only variant name is relevant in this case
            # Disregard HGVS strings for expression variants
            check_dict_entry(var_map[gene][var_id], "var_map", "name", "name")
            variant_name = var_map[gene][var_id]["name"]

            # For data_type=EXPR (in truth, data_type!=SNV), variant matching is not based on HGVS, so match_strings will only contain the variant name and the list of HGVS is empty
            civic_strings = civic_match_strings(variant_name, hgvs_expressions=[], data_type="EXPR")

            # Iterate available expression tags for the given gene and attempt match to a CIViC record for that gene
            # Opposite to the SNV and CNV case, here there is no tier hierarchy, i.e. gene expression is either matched in CIViC or not
            for expr_tag in expression_strings:
                if expr_tag in civic_strings:
                    # There could be 0,1,>1 exact matches (e.g. OVEREXPRESSION, EXPRESSION)
                    if var_id not in match["tier_1"]:
                        match["tier_1"].append(var_id)
            # Special case for EXPRESSION records related to EXONS: manually add them as matches when necessary
            # Specific expression change must be identical in order to consider them as matched: e.g. "OVEREXPRESSION" and "EXON 18 OVEREXPRESSION"
            for this_string in civic_strings:
                (is_exon, expr_type) = expr_is_exon_string(this_string)
                if is_exon and expr_type:
                    # NOTE: records of the type "P16 EXPRESSION" would be missed
                    if expr_type in expression_strings:
                        if var_id not in match["tier_1"]:
                            match["tier_1"].append(var_id)

        # If no exact match was found, then tier case is 3, and return all relevant EXPRESSION records for the current gene (if any are found in CIViC)
        # If there are no matching records for tier3, then this category will also be empty in the returned dictionary
        if not (match["tier_1"] or match["tier_1b"] or match["tier_2"]):
            # Not all CIViC records are returned but only those that are matched to a "EXPRESSION" tag
            match["tier_3"] = civic_return_all_expr(var_map[gene])
            # NOTE: Sanity check for situation where gene does not contain any variant records for the requested data type (e.g. all are SNVs)
            # Use dummy variant name to keep track of these situations
            if not match["tier_3"]:
                match["tier_3"] = ["NON_EXPR_MATCH_ONLY"]

    # If gene is not in CIViC, tier level is 4 and match will be empty
    else:
        # Sanity check that no other match should have been found
        if (match["tier_1"] or match["tier_1b"] or match["tier_2"] or match["tier_3"]):
            raise ValueError("Encountered unexpected tier case!")
        # No variants will be reported in this case (as the information is not available)
        match["tier_4"] = True
    return match


def add_match(match_map, gene, variant, match):
    """
    :param match_map:	Nested dictionary with fixed structure containing all tier categories and corresponding list of CIViC variant matches found in each case (if any) for the input gene and variant at hand. See README for more details about the specific structure used in 'match_map'.
    :param gene:	Single gene identifier (Entrez symbol, Entrez ID or CIViC ID) to append the given match into in dictionary 'match_map'.
    :param variant:	Single input variant annotation to append the given match into in dictionary 'match_map'.
    :param match:	Nested dictionary with fixed structure containing all tier matches found in CIViC for the given gene and variant annotation (if any).
    :return:		Tuple of one nested dictionary (updated 'match_map' dictionary, now including the new provided 'match' for the given gene and variant) and list (list of variant ids in CIViC which have been matched overall across all tiers for the given gene and variant). 
    """
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]

    # Sanity check of provided arguments
    check_arguments([match_map, gene, variant, match], ["match_map", "gene", "variant", "match"])
    check_is_dict(match_map, "match_map")
    check_is_str(gene, "gene")
    check_is_str(variant, "variant")
    check_is_dict(match, "match")

    # Sanity check that dict contains all expected keys (throws an error if not)
    check_keys(list(match.keys()), "match", sorted_tiers, matches_all=True)
    # Sanity check that provided match dict is not annotated for drug support
    for this_tier in match.keys():
        if this_tier != "tier_4":
            check_is_list(match[this_tier], this_tier)
        else:
            check_is_bool(match[this_tier], this_tier)
    # Sanity check that gene and variant are keys of the dict
    check_dict_entry(match_map, "match_map", gene, "gene")
    check_dict_entry(match_map[gene], "match_map", variant, "variant")
    # Overwritte any existing matches already assigned to the current variant
    match_map[gene][variant] = {}

    # Keep track of all the variant ids matched across tiers for the current gene + variant
    matched_ids = []
    # Add match results to the current gene + variant in match_map
    for x in match.keys():
        if x == "tier_4":
            match_map[gene][variant][x] = match[x]
        else:
            match_map[gene][variant][x] = []
            for var_id in match[x]:
                match_map[gene][variant][x].append(var_id)
                # It is possible that the same CIViC variant id has been matched for several tiers
                if var_id not in matched_ids:
                    matched_ids.append(var_id)

    return (match_map, matched_ids)


def match_in_civic(var_data, data_type, identifier_type, select_tier="all", var_map=None):
    """
    Match input molecular alterations to variant-level CIViC information.
    :param var_data:		Dictionary of input molecular data of a given type (SNV, CNV or EXPR) to be matched in CIViC, with an assumed fixed structure. See README for more details about the specific structure of 'var_data' depending on the data type (i.e. 'snv_data', 'cnv_data', 'expr_data'). 
    :param data_type:		['SNV', 'CNV', 'EXPR']
                        	SNV:   Variants correspond to genomic single nucleotide mutations and insertions/deletions
                        	CNV:   Variants correspond to genomic copy number alterations
                        	EXPR:  Variants correspond to differential gene expression data
				Type of data being queried.
    :param identifier_type:	['entrez_symbol', 'entrez_id', 'civic_id']
                        	entrez_symbol:   Entrez gene symbol
                        	entrez_id: 	 Entrez gene identifier
                        	civic_id: 	 CIViC internal identifier
                        	Type of gene identifier used in your query. Defaults to 'entrez_symbol'.
    :param select_tier:		['highest', 'all', <custom list>]
				highest:	Select the highest encountered tier per variant match case (hierarchy 1>1b>2>3>4)
				all:		Do not apply any filtering and return all available tiers for each variant match.
				<custom_list>:	Alternatively, the user can provide a list of the specific tier categories to select for (if all are provided, then no filtering is applied)
				Type of tier selection to be performed on the returned variant matches. Can be either a list or a string, and defaults to 'all'.
    :param var_map:		Optional. Nested dictionary of genes (must use the provided type of identifier) and variant records retrieved from CIViC, to be used for matching the input molecular alterations. Note: when supplying a custom 'var_map' for perfoming the matching of CIViC data (e.g. when filters need to be applied beforehand to remove undesirable CIViC records), anything that is not provided in this dictionary will be interpreted as not available in the knowledgebase. See README for more details about the specific structure expected for dictionary 'var_map'.
    :return:			Tuple of 3 elements: 'match_map' (nested dictionary with fixed structure containing all tier categories and corresponding matches found in CIViC), 'matched_ids' (list of all variant ids that could be matched in CIViC) and updated 'var_map' (nested dictionary with fixed structure containing all CIViC records associated with the matched variants).
    """
    # Hierarchy of tier classes (descending priority)
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]

    # Sanity check of provided arguments
    check_argument(var_data, "var_data")
    check_is_dict(var_data, "var_data")
    check_data_type(data_type)
    check_identifier_type(identifier_type)
    # Process and sanity check provided select_tier (expected format, valid values, etc.)
    select_tier = check_tier_selection(select_tier, sorted_tiers)

    match_map = {}
    matched_ids = []     # keep track of all matched variant ids across genes and tiers

    # Check if an existing var_map was provided or if CIViC needs to be queried
    if (var_map is None) or (not var_map):
        from query import query_civic
        # gene -> variant -> null
        all_genes = list(var_data.keys())
        var_map = query_civic(all_genes, identifier_type)
    else:
        check_is_dict(var_map, "var_map")

    # gene -> variant -> null
    # where variant -> var="dna|prot|impact|exon|n_line"
    for gene in var_data.keys():
        gene = check_empty_input(gene, "gene", is_required=True)
        if gene not in match_map.keys():
            match_map[gene] = {}
        for variant in var_data[gene].keys():
            variant = check_empty_input(variant, "variant", is_required=True)
            # Overwrite duplicated variant ids (should never happen)
            match_map[gene][variant] = {}
            variants = []
            impact_list = []
            exon_list = []

            variant_list = variant.split("|")
            if data_type == "SNV":
                # Sanity check for expected SNV format (at least 4 fields should exist, even if empty)
                if (len(variant_list) < 4):
                    raise ValueError("Must provide at least 4 fields to describe a SNV variant (even if some can be empty): 'dna|prot|impact|exon|..|'")
                # Format: var="dna|prot|[impact]|[exon]|..|"
                # NOTE: all fields can contain >1 terms separated with "," (no spaces). Fields "dna" and "prot" are required and "impacts" and "exons" are optional
                # NOTE: there might be additional fields after, e.g. "n_line" when data has been parsed from a file
                c_variants = variant_list[0]
                p_variants = variant_list[1]
                impacts = variant_list[2]
                exons = variant_list[3]

                # Sanity check for required and optional fields

                # Field must exist and cannot contain empty values
                c_variants_list = parse_input(c_variants, "c_variants", is_required=True)
                # Field must exist but can contain empty values
                p_variants_list = parse_input(p_variants, "p_variants", is_required=False)
                for c_var in c_variants_list:
                    # Sanity check that c. variant is not empty (as this field is not always required)
                    if not c_var:
                        continue
                    # Sanity check that variant starts with "c."
                    check_is_chgvs(c_var)
                    if c_var not in variants:
                        variants.append(c_var)
                for p_var in p_variants_list:
                    # Sanity check that p. variant is not empty (as this field is not always required)
                    if not p_var:
                        continue
                    # Sanity check that variant starts with "p."
                    check_is_phgvs(p_var)
                    if p_var not in variants:
                        variants.append(p_var)
                # Field is optional and can contain empty values
                impact_list = parse_input(impacts, "impacts", is_required=False)
                # Field is optional and can contain empty values
                exon_list = parse_input(exons, "exons", is_required=False)

                # Since c. and p. annotations can be optionally provided for each variant, sanity check that at least 1 annotation was provided (of either type)
                if not variants:
                    raise ValueError("At least one non-empty variant annotation (either 'c.' or 'p.') must be provided per variant!")

            if data_type == "CNV":
                # Format: var="cnv|.."
                # NOTE: CNV field can contain >1 terms separated with "," (no spaces)
                # NOTE: there might be additional fields after, e.g. "n_line" when data has been parsed from a file
                cnv_variants = variant_list[0]
                variants = parse_input(cnv_variants, "cnv_variants", is_required=True)

            if data_type == "EXPR":
                # Format: expr="logFC|.."
                # NOTE: logFC field should be a single number (sign will determine if gene is overexpressed or underexpressed)
                # NOTE: there might be additional fields after, e.g. "n_line" when data has been parsed from a file
                logfc = variant_list[0]
                variants = get_expression_strings(gene, logfc)

            # tier -> [matched_vars]
            if (data_type == "SNV") or (data_type == "CNV"):
                match = match_variants_in_civic(gene, variants, var_map, data_type, impacts=impact_list, exons=exon_list)

                # Avoid executing unnecessary filter when select_tier="all" (as resulting dict will be identical)
                if select_tier != "all":
                    match = filter_match(match, select_tier)

            # In the case of EXPRESSION data, only tier1 is possible (either CIViC record was matched or not)
            if (data_type == "EXPR"):
                match = match_expression_in_civic(gene, variants, var_map)

            # Add the match to the current entry for gene + variant
            # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
            (match_map, all_ids) = add_match(match_map, gene, variant, match)

            # Keep track of all the variant ids matched across genes and tiers
            for var_id in all_ids:
                if var_id not in matched_ids:
                    matched_ids.append(var_id)

    # At this point, all input gene + variants in var_data have been matched to CIViC
    # Filter var_map used to match CIViC info based on the matched variant ids (to avoid returning unnecessary records)
    from filtering import filter_civic
    var_map = filter_civic(var_map, var_id_in=matched_ids, output_empty=False)

    # Return match_map, list of all matched variant ids, and associated variant records retrieved from CIViC (or provided by user), already filtered for the matched variants
    return (match_map, matched_ids, var_map)


def filter_match(match, select_tier):
    """
    Filter tier matches associated to one single gene and variant from the input, according to the provided tier selection.
    :param match:		Nested dictionary with fixed structure containing all tier matches found in CIViC (if any) for one given gene and variant annotation.
    :param select_tier:		['highest', 'all', <custom list>]
				highest:	Select the highest encountered tier per variant match case (hierarchy 1>1b>2>3>4)
				all:		Do not apply any filtering and return all available tiers for each variant match.
				<custom_list>:	Alternatively, the user can provide a list of the specific tier categories to select for (if all are provided, then no filtering is applied)
				Type of tier selection to be performed on the provided variant matches. Can be either a list or a string.
    :return:			Updated 'match' dictionary after applying the provided tier filter/selection.
    """
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]

    # Sanity check of provided arguments
    check_argument(match, "match")
    check_is_dict(match, "match")
    # Process and sanity check provided select_tier (expected format, valid values, etc.)
    select_tier = check_tier_selection(select_tier, sorted_tiers)
    # Sanity check that dict contains all expected keys (throws an error if not)
    check_keys(list(match.keys()), "match", sorted_tiers, matches_all=True)
    # Sanity check that provided match dict is not annotated for drug support
    for this_tier in match.keys():
        if this_tier != "tier_4":
            check_is_list(match[this_tier], this_tier)
        else:
            check_is_bool(match[this_tier], this_tier)

    new = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}
    keep_tiers = []

    # When select_tier="all", filter is off (keep data for all tiers)
    if isinstance(select_tier, str) and (select_tier == "all"):
        for tmp_tier in sorted_tiers:
            keep_tiers.append(tmp_tier)

    # When select_tier="highest", then keep data only for the highest tier 1>1b>2>3>4
    elif isinstance(select_tier, str) and (select_tier == "highest"):
        for tmp_tier in sorted_tiers:
            if tmp_tier != "tier_4":
                if match[tmp_tier]:
                    keep_tiers.append(tmp_tier)
                    break
            else:
                # If we reached to this point, that means all other tiers were empty (tier_4=True)
                keep_tiers.append(tmp_tier)
                break

    # When select_tier is a list of tiers, then keep data only for those tiers
    elif isinstance(select_tier, list):
        for tmp_tier in sorted_tiers:
            if tmp_tier in select_tier:
                keep_tiers.append(tmp_tier)

    # Sanity check for any duplicated tiers
    keep_tiers = list(set(keep_tiers))

    # Keep variant matches only for the selected tiers (all other data will be excluded from output)
    # Respect the order of tiers while filling the new dictionary
    for tier in sorted_tiers:
        if tier in keep_tiers:
            if tier != "tier_4":
                for var_id in match[tier]:
                    new[tier].append(var_id)

    # Always check if current match corresponds to tier_4 after filtering
    # e.g. if only tier_1 is selected and current match did not have any
    if not (new["tier_1"] or new["tier_1b"] or new["tier_2"] or new["tier_3"]):
        new["tier_4"] = True

    return new


def filter_matches(match_map, select_tier):
    """
    Apply filtering of tier matches (according to the provided tier selection) to the entire nested dictionary of input genes and variants queried in CIViC.
    :param match_map:		Nested dictionary with fixed structure containing all tier matches found in CIViC (if any) for a set of genes and associated variant annotations. See README for more details about the structure of 'match_map'.
    :param select_tier:		['highest', 'all', <custom list>]
				highest:	Select the highest encountered tier per variant match case (hierarchy 1>1b>2>3>4)
				all:		Do not apply any filtering and return all available tiers for each variant match.
				<custom_list>:	Alternatively, the user can provide a list of the specific tier categories to select for (if all are provided, then no filtering is applied)
				Type of tier selection to be performed on the provided variant matches. Can be either a list or a string.
    :return:			Tuple of 2 elements: 'clean_map' (updated 'match_map' dictionary after applying the provided tier filter/selection) and 'matched_ids' (list of all variant ids that could be matched in CIViC).
    """
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]

    # Check provided argument
    check_argument(match_map, "match_map")
    check_is_dict(match_map, "match_map")

    clean_map = {}
    matched_ids = []

    # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
    # where variant -> var="dna|prot|impact|exon|n_line"
    for gene in match_map.keys():
        if gene not in clean_map.keys():
            clean_map[gene] = {}
        for variant in match_map[gene].keys():
            clean_map[gene][variant] = {}
            new = filter_match(match_map[gene][variant], select_tier)
            # Add the filtered match to the current entry for gene + variant
            (clean_map, all_ids) = add_match(clean_map, gene, variant, new)
            # Keep track of all the variant ids matched across genes and tiers
            for var_id in all_ids:
                if var_id not in matched_ids:
                    matched_ids.append(var_id)

    return (clean_map, matched_ids)


def classify_diseases(disease_list, disease_name_not_in, disease_name_in, alt_disease_names):
    """
    Given a list of disease names, first filter out those partially matching the non-allowed terms (if any), and from the remaining set, return either:
       - Subset of diseases names from the list matching white-listed terms - use a partial match. These disease names can be considered as 'ct' (cancer type specific).
       - If previous is not available, subset of disease names matching high-level or general terms - use an exact match. These disease names can be considered as 'gt' (general cancer type).
       - If previous is not available, return all available disease names supplied as input (i.e. thus, original set except those matched to the non-allowed list, if any where provided and foun). These disease names can be considered as 'nct' (non cancer type specific).
    :param disease_list:		List of disease names from CIViC to be classified into 'ct', 'gt' or 'nct' categories.
    :param disease_name_not_in:		List of non-allowed terms to remove undesirable disease names. Partial matching is applied. Can be empty.
    :param disease_name_in:		List of clinically relevant ('ct') terms to select disease names of interest. Partial matching is applied (e.g. returns 'Uveal melanoma' when term 'melanoma' is provided). Can be empty.
    :param alt_disease_names:		List of alternative broad-definition or general ('gt') terms to use as a second-best strategy to select disease names of interest (e.g. 'Cancer', 'Solid tumor'). Exact matching is applied (e.g. returns disease 'cancer' only when 'cancer' is the provided term). Can be empty.
    :return:				Tuple of 3 lists: 'ct_list' (disease names classified as 'ct'), 'gt_list' (disease names classified as 'gt'), 'nct_list' (disease names classified as 'nct').
    """
    # Sanity check and uppercase terms in all the provided lists
    disease_list = uppercase_list(disease_list, "disease_list")
    disease_name_not_in = uppercase_list(disease_name_not_in, "disease_name_not_in")
    disease_name_in = uppercase_list(disease_name_in, "disease_name_in")
    alt_disease_names = uppercase_list(alt_disease_names, "alt_disease_names")

    # Keep track of diseases that did not pass the non-allowed filter
    non_allowed_matched = []
    # Keep track of diseases that passed the non-allowed filter
    clean_set = []
    # Final list of matched diseases (associated to one of 3 categories depending on match type: "ct", "gt" or "nct")
    matched = []
    # Keep track of which type of disease specificy was matched in the end ("ct", "gt" or "nct")
    ct_list = []
    gt_list = []
    nct_list = []

    ## 1) First, remove diseases that partially match non-allowed terms (if any are provided)
    # NOTE: PARTIAL matches to the non-allowed terms are permited! eg:
    #   - including "small" will remove "non-small cell lung cancer" and "lung small cell carcinoma"
    #   - including "non-small" will remove "non-small cell lung cancer" but not "lung small cell carcinoma"
    if disease_name_not_in:
        # Iterate available diseases and keep track of those that partially match to non-allowed terms (there can be several)
        for disease in disease_list:
            # To find partial matches, it is necessary to iterate through both lists (input and civic)
            for non_allowed in disease_name_not_in:
                # Search for partial match of INPUT disease in CIViC disease e.g. "Melanoma" (input list) in "Skin Melanoma" (CIViC) and not opposite
                if non_allowed in disease:
                    if disease not in non_allowed_matched:
                        non_allowed_matched.append(disease)
        # Iterate available diseases once again to retrieve those that passed the non-allowed filter
        for disease in disease_list:
            # Retrieve valid diseases only
            if disease in non_allowed_matched:
                continue
            if disease not in clean_set:
                clean_set.append(disease)

    ## If no terms were provided by the user in the non-allowed list, then all available diseases constitute the clean set
    else:
        clean_set = disease_list

    ## 2) Now, iterate the list of "allowed" diseases (i.e. passing the non-allowed filter) and attempt to match to white-listed terms
    # NOTE: again, PARTIAL matches to the white list are allowed! eg:
    #   - including "melanoma" will match "melanoma", "skin melanoma" and "uveal melanoma", but not "skin cancer" (hypothetical)
    #   - including "uveal melanoma" will only match "uveal melanoma"
    for clean_type in clean_set:
        # To find partial matches, it is necessary to iterate through both lists (input and civic)
        for allowed in disease_name_in:
            # Search for partial match of INPUT disease in CIViC disease e.g. "Melanoma" (input list) in "Skin Melanoma" (CIViC) and not opposite
            if allowed in clean_type:
                # Keep track of diseases that passed the white list filter
                if clean_type not in matched:
                    matched.append(clean_type)
                    if clean_type not in ct_list:
                        ct_list.append(clean_type)

    ## 3) For those not already matched as ct, attempt a second-best match strategy to higher-level diseases
    ## In CIViC, some "general" diseases are included as disease, e.g. "cancer" or "solid tumor". Hence, must be exact matches since they are DB-specific
    # NOTE: here, only PERFECT matches are allowed!
    #   - including "cancer" will only match "cancer" and not "lung cancer"
    for clean_type in clean_set:
        if clean_type in alt_disease_names:
            if clean_type not in matched:
                matched.append(clean_type)
                if clean_type not in gt_list:
                    gt_list.append(clean_type)

    ## 4) For those not already matched as either ct or gt, return all allowed (i.e. not in the non-allowed list) diseases available in CIViC
    ## These will be considered as non-specific diseases (i.e. off label)
    for clean_type in clean_set:
        if clean_type not in matched:
            matched.append(clean_type)
            if clean_type not in nct_list:
                nct_list.append(clean_type)

    ## Return diseases that did not match any non-allowed terms, classified and split into: ct, gt, nct
    return (ct_list, gt_list, nct_list)


def add_ct(diseases, ct, gene, variant, molecular_profile, evidence_type, new_map, var_map, is_annot=False):
    """
    Given a list of disease names classified as one single cancer type specificity (ct, gt or nct), include this disease information in the corresponding CIViC records of the gene, variant and evidence type. 
    :param diseases:		List of disease names to be annotated with the provided cancer type specificity category (i.e. 'ct') in 'var_map'.
    :param ct:			['ct', 'gt', 'nct']
				ct:	Cancer type specific 
				gt:	General cancer type specificity
				nct:	non-specific cancer type
				Category of cancer type specificity, to be annotated on the provided CIViC data.
    :param gene:		Gene identifier to append the given disease specificity annotations into in dictionary 'var_map'.
    :param variant:		Input variant annotation to append the given disease specificity annotations into in dictionary 'var_map'.
    :param evidence_type:	Evidence type to append the given disease specificity annotations into in dictionary 'var_map'.
    :param new_map:		Updated 'var_map' dictionary to append the given disease specificity information into.
    :param var_map:		Nested dictionary of genes and variant-level records retrieved from CIViC. See README for more details about the specific structure.
    :param is_annot:		Boolean indicating whether the provided 'var_map' is already annotated with cancer type specificity information (True) or not (False). 
    :return:			Updated 'var_map' dictionary, including the supplied cancer type specificity category (i.e. 'ct') for all the provided disease names (i.e. 'diseases').
    """
    sorted_cts = ["ct", "gt", "nct"]

    # Sanity check for the expected object types and also expected keys contained in both dicts
    check_is_dict(new_map, "new_map")
    check_is_dict(var_map, "var_map")
    check_is_none(is_annot, "is_annot")
    
    if ct not in sorted_cts:
        raise ValueError("Provided ct '%s' is not valid! Please provide one of: %s" %(ct,sorted_cts))
    check_is_bool(is_annot, "is_annot")
    check_dict_entry(new_map, "new_map", gene, "gene")
    check_dict_entry(new_map[gene], "new_map", variant, "variant")
    check_dict_entry(new_map[gene][variant][molecular_profile], "new_map", "evidence_items", "evidence_items")
    check_dict_entry(new_map[gene][variant][molecular_profile]["evidence_items"], "new_map", evidence_type, "evidence_type")
    check_dict_entry(var_map, "var_map", gene, "gene")
    check_dict_entry(var_map[gene], "var_map", variant, "variant")
    check_dict_entry(var_map[gene][variant][molecular_profile], "var_map", "evidence_items", "evidence_items")
    check_dict_entry(var_map[gene][variant][molecular_profile]["evidence_items"],"var_map", evidence_type, "evidence_type")
    
    # Check whether provided var_map is annotated with disease specificity info (ct/gt/nct) or not
    if is_annot:
        check_keys(list(var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts, matches_all=True)
        # If a non-empty list of diseases was provided, check they exist in var_map
        if diseases:
            diseases = uppercase_list(diseases, "diseases")
            check_keys(list(var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct].keys()), "var_map", diseases, matches_all=False)
    else:
        check_keys_not(list(var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts)
        # If a non-empty list of diseases was provided, check they exist in var_map
        if diseases:
            diseases = uppercase_list(diseases, "diseases")
            check_keys(list(var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type].keys()), "var_map", diseases, matches_all=False)
    
    # Do not check if new_map already contains ct annotations (overwritten in this case)
    new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct] = {}
    
    # Iterate available diseases and append all available info to new_map (overwritte if necessary)
    for disease in diseases:
        new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease] = {}
        if is_annot:
            for drug in var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease].keys():
                new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug] = {}
                for evidence in var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug].keys():
                    new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug][evidence] = {}
                    for evidence_level in var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug][evidence].keys():
                        new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level] = []
                        for thisString in var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level]:
                            new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level].append(thisString)
        else:
            for drug in var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][disease].keys():
                new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug] = {}
                for evidence in var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][disease][drug].keys():
                    new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug][evidence] = {}
                    for evidence_level in var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][disease][drug][evidence].keys():
                        new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level] = []
                        for thisString in var_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][disease][drug][evidence][evidence_level]:
                            new_map[gene][variant][molecular_profile]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level].append(thisString)
    return new_map


def annotate_ct(var_map, disease_name_not_in, disease_name_in, alt_disease_names):
    """
    Annotate a nested dictionary of variant-level CIViC records with cancer type specificity annotations according to the associated disease name in each case.
    :param var_map:			Nested dictionary of genes and variants retrieved from CIViC. See README for more details about the specific structure.
    :param disease_name_not_in:		List of non-allowed terms to remove undesirable disease names. Partial matching is applied. Can be empty.
    :param disease_name_in:		List of clinically relevant ('ct') terms to select disease names of interest. Partial matching is applied (e.g. returns 'Uveal melanoma' when term 'melanoma' is provided). Can be empty.
    :param alt_disease_names:		List of alternative broad-definition or general ('gt') terms to use as a second-best strategy to select disease names of interest (e.g. 'Cancer', 'Solid tumor'). Exact matching is applied (e.g. returns disease 'cancer' only when 'cancer' is the provided term). Can be empty.
    :return:				Updated 'var_map' dictionary after annotating disease specificity according to the provided terms. See README for more details about the specific structure of 'var_map' when disease specificity annotations are included.
    """
    
    sorted_cts = ["ct", "gt", "nct"]
    var_map_entries_variant = ["name", "hgvs", "types"]
    new_map = {}

    # Iterate the complete var_map dict and reorganize it to classify diseases
    for gene in var_map.keys():
        if gene not in new_map.keys():
            new_map[gene] = {}
        for variant in var_map[gene].keys():
            # Overwrite duplicated variant ids (should never happen)
            new_map[gene][variant] = {}
            # Sanity check that some expected fields can be found in the dictionary
            check_keys(list(var_map[gene][variant].keys()), "var_map", var_map_entries_variant, matches_all=False)
            new_map[gene][variant]["name"] = var_map[gene][variant]["name"]
            new_map[gene][variant]["hgvs"] = [a for a in var_map[gene][variant]["hgvs"]]
            new_map[gene][variant]["types"] = [b for b in var_map[gene][variant]["types"]]
            
            molecular_profile_ids = set(list(var_map[gene][variant].keys())) ^ set(var_map_entries_variant)

            for molecular_profile_id in molecular_profile_ids:
                new_map[gene][variant][molecular_profile_id] = {}
                new_map[gene][variant][molecular_profile_id]["civic_score"] = var_map[gene][variant][molecular_profile_id]["civic_score"]
                new_map[gene][variant][molecular_profile_id]["n_evidence_items"] = var_map[gene][variant][molecular_profile_id]["n_evidence_items"]
                new_map[gene][variant][molecular_profile_id]["evidence_items"] = {}
                for evidence_type in var_map[gene][variant][molecular_profile_id]["evidence_items"].keys():
                    new_map[gene][variant][molecular_profile_id]["evidence_items"][evidence_type] = {}
                    # Retrieve all disease names associated witht he current evidence type, and classify them into ct, gt, nct
                    all_diseases = list(var_map[gene][variant][molecular_profile_id]["evidence_items"][evidence_type].keys())
                    # Check that provided var_map is not annotated with disease specificity info (ct/gt/nct)
                    # check_keys_not(all_diseases, "var_map", sorted_cts)
                    #   Classify diseases by provided specificity criteria (ct/gt/nct)
                    (ct_disease_list, gt_disease_list, nct_disease_list) = classify_diseases(all_diseases, disease_name_not_in, disease_name_in, alt_disease_names)
                    # Add new layer of classification within the evidence types: specificity of the disease (ct, gt, nct)
                    for ct in sorted_cts:
                        if ct == "ct":
                            new_map = add_ct(ct_disease_list, ct, gene, variant, molecular_profile_id, evidence_type, new_map, var_map, is_annot=False)
                        if ct == "gt":
                            new_map = add_ct(gt_disease_list, ct, gene, variant, molecular_profile_id, evidence_type, new_map, var_map, is_annot=False)
                        if ct == "nct":
                            new_map = add_ct(nct_disease_list, ct, gene, variant, molecular_profile_id, evidence_type, new_map, var_map, is_annot=False)

    return new_map


def filter_ct(var_map, select_ct):
    """
    Apply filtering of variant-level CIViC records according to their associated disease names, based on the provided cancer type specificity selection.
    :param var_map:		Nested dictionary of genes and variants retrieved from CIViC. Must be annotated with disease specificity information (i.e. 'ct'). See README for more details about the specific structure expected by this function.
    :param select_ct:		['highest', 'all', <custom list>]
				highest:	Select the highest encountered disease specificity category per variant match case (hierarchy ct>gt>nct)
				all:		Do not apply any filtering and return all available disease specificities for each variant match.
				<custom_list>:	Alternatively, the user can provide a list of the specific disease specificity categories to select for (if all are provided, then no filtering is applied)
				Type of cancer type specificity selection to be performed on the provided variant records. Can be either a list or a string.
    :return:			Updated 'var_map' dictionary after filtering disease specificity according to the provided selection.
    """
    sorted_cts = ["ct", "gt", "nct"]
    var_map_entries_variant = ["name", "hgvs", "types"]

    # Sanity check for the expected object types
    check_is_dict(var_map, "var_map")
    # Process and sanity check provided select_ct (expected format, valid values, etc.)
    select_ct = check_tier_selection(select_ct, sorted_cts)
    new_map = {}
    # When select_ct="all", filter is off (keep data for all 3 cts)
    # Do not copy dict again, simply return input dict untouched
    if isinstance(select_ct, str) and (select_ct == "all"):
        # Sanity check that provided var_map has expected format and is annotated for ct
        for gene in var_map.keys():
            for variant in var_map[gene].keys():
                check_keys(list(var_map[gene][variant].keys()), "var_map", var_map_entries_variant, matches_all=False)
                molecular_profile_ids = set(list(var_map[gene][variant].keys())) ^ set(var_map_entries_variant)
                for molecular_profile_id in molecular_profile_ids:
                    for evidence_type in var_map[gene][variant][molecular_profile_id]["evidence_items"].keys():
                        check_keys(list(var_map[gene][variant][molecular_profile_id]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts, matches_all=True)
        return var_map
    
    # Otherwise, some filtering needs to be done, so iterate var_map
    # Iterate the complete var_map dict and reorganize it to classify diseases
    for gene in var_map.keys():
        if gene not in new_map.keys():
            new_map[gene] = {}
        for variant in var_map[gene].keys():
            # Overwrite duplicated variant ids (should never happen)
            new_map[gene][variant] = {}
            # Sanity check that some expected fields can be found in the dictionary
            check_keys(list(var_map[gene][variant].keys()), "var_map", var_map_entries_variant, matches_all=False)
            new_map[gene][variant]["name"] = var_map[gene][variant]["name"]
            new_map[gene][variant]["hgvs"] = [a for a in var_map[gene][variant]["hgvs"]]
            new_map[gene][variant]["types"] = [b for b in var_map[gene][variant]["types"]]
            molecular_profile_ids = set(list(var_map[gene][variant].keys())) ^ set(var_map_entries_variant)
            for molecular_profile_id in molecular_profile_ids:
                new_map[gene][variant][molecular_profile_id] = {}
                new_map[gene][variant][molecular_profile_id]["civic_score"] = var_map[gene][variant][molecular_profile_id]["civic_score"]
                new_map[gene][variant][molecular_profile_id]["n_evidence_items"] = var_map[gene][variant][molecular_profile_id]["n_evidence_items"]
                new_map[gene][variant][molecular_profile_id]["evidence_items"] = {}
                for evidence_type in var_map[gene][variant][molecular_profile_id]["evidence_items"].keys():
                    new_map[gene][variant][molecular_profile_id]["evidence_items"][evidence_type] = {}
                    # Sanity check that input dict is annotated with disease specificity info (ie. has the expected format)
                    check_keys(list(var_map[gene][variant][molecular_profile_id]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts, matches_all=True)
                    skip = False
                    for ct in sorted_cts:
                        new_map[gene][variant][molecular_profile_id]["evidence_items"][evidence_type][ct] = {}
                        ct_list = list(var_map[gene][variant][molecular_profile_id]["evidence_items"][evidence_type][ct].keys())
                        # When select_ct="highest", then keep data only for the highest specificity available ct>gt>nct
                        if isinstance(select_ct, str) and (select_ct == "highest"):
                            # Check if data is available for the currently iterated ct
                            if ct_list and (not skip):
                                new_map = add_ct(ct_list, ct, gene, variant, molecular_profile_id, evidence_type, new_map, var_map, is_annot=True)
                                # Prematurely exit the loop after this successful iteration (to keep only the highest ct)
                                skip = True
                        # When select_ct is a list of cts, then keep data only for those selected
                        elif isinstance(select_ct, list):
                            # Check if the currently iterated ct should be kept or skipped
                            # Iterate all cts and add all those that are select (no premature exit of loop)
                            if (ct in select_ct):
                                new_map = add_ct(ct_list, ct, gene, variant, molecular_profile_id, evidence_type, new_map, var_map, is_annot=True)
    return new_map


def process_drug_support(match_map, var_map, support_dict):
    """
    Given a dictionary of CIViC variant-level records matched to a set of input molecular alterations ('match_map'), and the corresponding set of associated clinical information retrieved from CIViC ('var_map'), compute consensus drug response predictions based on the available 'predictive' CIViC information and a helper dictionary of evidence-to-drug response provided by the user ('support_dict').
    :param match_map:		Nested dictionary with fixed structure containing all tier categories and corresponding list of CIViC variant matches found in each case (if any) for the input gene and molecular alteration at hand. This dictionary will be annotated by the function with consensus drug response information computed based on the provided 'predictive' CIViC evidence. See README for more details about the specific structure expected for this dictionary, before and after the annotation of consensus drug information.
    :param var_map:		Nested dictionary of genes and variant-level evidence records retrieved from CIViC, to be used for computing consensus drug response predictions based on the available 'predictive' evidence of the matched CIViC records (provided in 'match_map'). See README for more details about the specific structure expected for this dictionary.
    :param support_dict:	Dictionary of evidence-to-drug response, provided in the data.yml file, and to be used for computing consensus drug predictions. See README for more details about this topic.
    :return:			Updated 'match_map' dictionary, containing consensus drug response prediction information. Individual consensus predictions use format 'DRUG_NAME:CT_CLASS:CONSENSUS_RESPONSE:#positive|#negative|#unknown|#do_not_support'.
    """
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]
    sorted_cts = ["ct", "gt", "nct"]
    evidence_type = "PREDICTIVE"
    special_cases = ["NON_SNV_MATCH_ONLY", "NON_CNV_MATCH_ONLY", "NON_EXPR_MATCH_ONLY"]
    var_map_entries_variant = ["name", "hgvs", "types"]
    
    # Check provided argument
    check_arguments([match_map, support_dict], ["match_map", "support_dict"])
    check_is_dict(match_map, "match_map")
    check_is_dict(var_map, "var_map")
    check_is_dict(support_dict, "support_dict")
    
    # Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
    new_map = {}
    
    # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
    # where variant -> var="dna|prot|impact|exon|n_line"
    for gene in match_map.keys():
        if gene not in new_map.keys():
            new_map[gene] = {}
        for variant in match_map[gene].keys():
            new_map[gene][variant] = {}
            # Sanity check that provided match_map has expected format (ie. is annotated with drug support)
            check_keys(list(match_map[gene][variant].keys()), "match_map", sorted_tiers, matches_all=True)
            for tier in match_map[gene][variant].keys():
                new_map[gene][variant][tier] = {}
                # Sanity check that provided match_map has expected format
                # NOTE: tier4 has specific format compared to the others (boolean vs. list)
                if tier == "tier_4":
                    check_is_bool(match_map[gene][variant][tier], tier)
                    new_map[gene][variant][tier]["matched"] = False
                else:
                    check_is_list(match_map[gene][variant][tier], tier)
                    new_map[gene][variant][tier]["matched"] = []

                new_map[gene][variant][tier]["drug_support"] = []
                drug_map = {}
                if tier != "tier_4":
                    for var_id in match_map[gene][variant][tier]:
                        new_map[gene][variant][tier]["matched"].append(var_id)
                        # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                        # This is a dummy tag and not an actual variant record from CIViC, so skip checking in var_map
                        if var_id.upper() in special_cases:
                            # Sanity check that no other variant was matched when this special case was matched (length of matches should always be one)
                            if len(match_map[gene][variant][tier]) != 1:
                                raise ValueError("Unexpected: encountered multiple matches in special case of empty tier3 match '%s'!" %(match_map[gene][variant][tier]))
                            continue
                        
                        # For matched variants, they must be contained in the provided var_map
                        check_dict_entry(var_map, "var_map", gene, "gene")
                        check_dict_entry(var_map[gene], "var_map", var_id, "variant")
                        
                        molecular_profile_ids = set(list(var_map[gene][var_id].keys())) ^ set(var_map_entries_variant)
                        for molecular_profile_id in molecular_profile_ids:
                            check_dict_entry(var_map[gene][var_id][molecular_profile_id], "var_map", "evidence_items","evidence_items")
                            if evidence_type in var_map[gene][var_id][molecular_profile_id]["evidence_items"].keys():
                                # Sanity check that provided var_map has expected format (ie. annotated for disease specificity)
                                check_keys(list(var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts, matches_all=True)
                                for ct in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type].keys():
                                    for disease in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct].keys():
                                        for drug in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct][disease].keys():
                                            if drug not in drug_map.keys():
                                                drug_map[drug] = {}
                                            if ct not in drug_map[drug].keys():
                                                drug_map[drug][ct] = []
                                            for evidence in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct][disease][drug].keys():
                                                # Split the evidence direction and clinical significance
                                                evidence_list = evidence.strip().split(":")
                                                if (len(evidence_list) != 2):
                                                    raise ValueError("Unexpected format of evidence '%s'! Please provide string as 'EVIDENCE_DIRECTION:CLINICAL_SIGNIFICANCE'." %(evidence))
                                                direction = evidence_list[0]
                                                clin_signf = evidence_list[1]
                                                
                                                # For each evidence (ie combination of direction+clin_signf), count how many different evidence items support it
                                                # At this stage, we find count evidence items by counting how many different combinations of level+pmids there are for the same drug, disease and evidence
                                                if ("NULL" in direction) or ("N/A" in direction) or ("NULL" in clin_signf) or ("N/A" in clin_signf):
                                                    this_drug_support = "UNKNOWN_BLANK"
                                                else:
                                                    if direction not in support_dict.keys():
                                                        raise ValueError("Could not find evidence direction '%s' in provided 'support_dict'!" %(direction))
                                                    if clin_signf not in support_dict[direction].keys():
                                                        raise ValueError("Could not find clinical significance '%s' in provided 'support_dict'!" %(clin_signf))
                                                    this_drug_support = support_dict[direction][clin_signf]
                                                    
                                                # Keep track of number of occurrences for each support type for the given drug
                                                # Here, take into account the number of supporting PMIDs associated to each evidence item
                                                for evidence_level in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct][disease][drug][evidence].keys():
                                                    for this_evidence_string in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level]:
                                                        drug_map[drug][ct].append(this_drug_support)
                                                        
                    # Process drug support information parsed for the current tier match
                    # Drug support for tier 4 will never be available
                    for this_drug in drug_map.keys():
                        for this_ct in drug_map[this_drug].keys():
                            # Given the selected ct, count number of occurrences for each possible support type (if any)
                            count_pos = drug_map[this_drug][this_ct].count("POSITIVE")
                            count_neg = drug_map[this_drug][this_ct].count("NEGATIVE")
                            count_unk = drug_map[this_drug][this_ct].count("UNKNOWN_BLANK")
                            count_dns = drug_map[this_drug][this_ct].count("UNKNOWN_DNS")
                            # Pool UNKNOWN_BLANK and UNKNOWN_DNS together (as both result in unknown CIViC support)
                            count_total_unk = count_unk + count_dns
                            # Sanity check that there is at least some support
                            if (count_pos == 0) and (count_neg == 0) and (count_total_unk == 0):
                                raise ValueError("Encountered unexpected support case for gene '%s'." %(gene))
                            
                            # Resolve contradicting evidence (if any) by majority vote
                            temp_support = ""
                            # For this, pool UNKNOWN_BLANK and UNKNOWN_DNS together
                            # Whenever there is a tie of "confident" (pos or neg) vs "non-confident" (unk), choose the confident one
                            if (count_total_unk > count_pos) and (count_total_unk > count_neg):
                                temp_support = "CIVIC_UNKNOWN"
                            elif count_pos == count_neg:
                                temp_support = "CIVIC_CONFLICT"
                            elif (count_pos > count_neg) and (count_pos >= count_total_unk):
                                temp_support = "CIVIC_SUPPORT"
                            elif (count_neg > count_pos) and (count_neg >= count_total_unk):
                                temp_support = "CIVIC_RESISTANCE"
                            else:
                                raise ValueError("Encountered unexpected support case for gene '%s'." %(gene))
                            
                            # Build support string for each given combination of drug, ct and matched tier
                            # Format: DRUG:CT:SUPPORT:#pos|#neg|#unk|#dns
                            drug_support = this_drug + ":" + this_ct.upper() + ":" + temp_support
                            new_map[gene][variant][tier]["drug_support"].append(drug_support)

            # Always check if current match corresponds to a tier_4 situation (all other tiers will be empty)
            if not (new_map[gene][variant]["tier_1"]["matched"] or new_map[gene][variant]["tier_1b"]["matched"] or new_map[gene][variant]["tier_2"]["matched"] or new_map[gene][variant]["tier_3"]["matched"]):
                new_map[gene][variant]["tier_4"]["matched"] = True
    return new_map


def reprocess_drug_support_across_selected_variants(input_data, match_map, var_map, support_dict, has_support=True):
    """
    Recompute consensus drug response predictions for a specific set of input genes and associated molecular alterations based on the available predictive evidence matched in CIViC.
    Given a dictionary of input genes and aberrations to be considered for the computation ('input_data'), a dictionary of CIViC variant-level records matched to the input alterations ('match_map'), and the corresponding set of associated clinical information retrieved from CIViC ('var_map'), compute consensus drug response predictions based on the available 'predictive' drug data and a helper dictionary of evidence-to-drug response provided by the user ('support_dict').
    :param input_data:		Nested dictionary with fixed structure identical to that of 'match_map' when consensus drug annotations are not available, used by this function to specify a subset of genes and associated variants that should be considered for the recomputation of consensus drug response predictions. See README for more details about the specific structure expected for this dictionary.
    :param match_map:		Nested dictionary with fixed structure containing all tier categories and corresponding list of CIViC variant matches found in each case (if any) for the input gene and molecular alteration at hand. This dictionary will be annotated by the function with consensus drug response information computed based on the provided 'predictive' CIViC evidence. See README for more details about the specific structure expected for this dictionary, before and after the annotation of consensus drug information.
    :param var_map:		Nested dictionary of genes and variant-level evidence records retrieved from CIViC, to be used for computing consensus drug response predictions based on the available 'predictive' evidence of the matched CIViC records (provided in 'match_map'). See README for more details about the specific structure expected for this dictionary.
    :param support_dict:	Dictionary of evidence-to-drug response, provided in the data.yml file, and to be used for computing consensus drug predictions. See README for more details about this topic.
    :param has_support:		Boolean indicating whether the provided 'match_map' is already annotated with consensus drug response prediction information (True) or not (False). 
    :return:			List of consensus drug response predictions computed based on the aggregation of all available evidences across the supplied input genes and aberrations. Individual consensus predictions use format 'DRUG_NAME:CT_CLASS:CONSENSUS_RESPONSE:#positive|#negative|#unknown|#do_not_support'.
    """
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]
    sorted_cts = ["ct", "gt", "nct"]
    evidence_type = "PREDICTIVE"
    special_cases = ["NON_SNV_MATCH_ONLY", "NON_CNV_MATCH_ONLY", "NON_EXPR_MATCH_ONLY"]
    var_map_entries_variant = ["name", "hgvs", "types"]

    # Check provided argument
    check_arguments([match_map, support_dict], ["match_map", "support_dict"])
    check_is_dict(input_data, "input_data")
    check_is_dict(match_map, "match_map")
    check_is_dict(var_map, "var_map")
    check_is_dict(support_dict, "support_dict")

    # Process and aggregate drug support across all the available evidences for the supplied input_data (genes + vaiants)
    drug_support_strings = []
    # Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
    drug_map = {}

    # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
    # where variant -> var="dna|prot|impact|exon|n_line"
    for gene in input_data.keys():
        if gene not in match_map.keys():
            raise ValueError("Input gene '%s' from 'input_data' could not be found in provided 'match_map'!" %(gene))
        for variant in input_data[gene].keys():
            if variant not in match_map[gene].keys():
                raise ValueError("Input variant '%s' (gene '%s') from 'input_data' could not be found in provided 'match_map'!" %(variant, gene))

            # Sanity check that provided match_map has expected format (i.e. is annotated with drug support)
            check_keys(list(match_map[gene][variant].keys()), "match_map", sorted_tiers, matches_all=True)
            for tier in match_map[gene][variant].keys():
                # Sanity check that provided match_map has expected format
                # NOTE: tier4 has specific format compared to the others (boolean vs. list)
                if tier == "tier_4":
                    if has_support:
                        check_is_bool(match_map[gene][variant][tier]["matched"], tier)
                    else:
                        check_is_bool(match_map[gene][variant][tier], tier)
                else:
                    if has_support:
                        check_is_list(match_map[gene][variant][tier]["matched"], tier)
                    else:
                        check_is_list(match_map[gene][variant][tier], tier)

                if tier != "tier_4":
                    if has_support:
                        variant_list = match_map[gene][variant][tier]["matched"]
                    else:
                        variant_list = match_map[gene][variant][tier]
                    for var_id in variant_list:
                        # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                        # This is a dummy tag and not an actual variant record from CIViC, so skip checking in var_map
                        if var_id.upper() in special_cases:
                            # Sanity check that no other variant was matched when this special case was matched (length of matches should always be one)
                            if len(variant_list) != 1:
                                raise ValueError("Unexpected: encountered multiple matches in special case of empty tier3 match '%s'!" %(match_map[gene][variant][tier]))
                            continue

                        # For matched variants, they must be contained in the provided var_map
                        check_dict_entry(var_map, "var_map", gene, "gene")
                        check_dict_entry(var_map[gene], "var_map", var_id, "variant")
                        
                        molecular_profile_ids = set(list(var_map[gene][var_id].keys())) ^ set(var_map_entries_variant)
                        for molecular_profile_id in molecular_profile_ids:
                            check_dict_entry(var_map[gene][var_id][molecular_profile_id], "var_map", "evidence_items", "evidence_items")
                            if evidence_type in var_map[gene][var_id][molecular_profile_id]["evidence_items"].keys():
                                # Sanity check that provided var_map has expected format (i.e. annotated for disease specificity)
                                check_keys(list(var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts, matches_all=True)
                                for ct in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type].keys():
                                    for disease in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct].keys():
                                        for drug in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct][disease].keys():
                                            if drug not in drug_map.keys():
                                                drug_map[drug] = {}
                                            if ct not in drug_map[drug].keys():
                                                drug_map[drug][ct] = []
                                            for evidence in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct][disease][drug].keys():
                                                # Split the evidence direction and clinical significance
                                                evidence_list = evidence.strip().split(":")
                                                if (len(evidence_list) != 2):
                                                    raise ValueError("Unexpected format of evidence '%s'! Please provide string as 'EVIDENCE_DIRECTION:CLINICAL_SIGNIFICANCE'." %(evidence))
                                                direction = evidence_list[0]
                                                clin_signf = evidence_list[1]

                                                # For each evidence (ie combination of direction+clin_signf), count how many different evidence items support it
                                                # At this stage, we find count evidence items by counting how many different combinations of level+pmids there are for the same drug, disease and evidence
                                                if ("NULL" in direction) or ("N/A" in direction) or ("NULL" in clin_signf) or ("N/A" in clin_signf):
                                                    this_drug_support = "UNKNOWN_BLANK"
                                                else:
                                                    if direction not in support_dict.keys():
                                                        raise ValueError("Could not find evidence direction '%s' in provided 'support_dict'!" %(direction))
                                                    if clin_signf not in support_dict[direction].keys():
                                                        raise ValueError("Could not find clinical significance '%s' in provided 'support_dict'!" %(clin_signf))
                                                    this_drug_support = support_dict[direction][clin_signf]

                                                # Keep track of number of occurrences for each support type for the given drug
                                                # Here, take into account the number of supporting PMIDs associated to each evidence item
                                                for evidence_level in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct][disease][drug][evidence].keys():
                                                    for this_evidence_string in var_map[gene][var_id][molecular_profile_id]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level]:
                                                        drug_map[drug][ct].append(this_drug_support)

    # Process drug support information parsed and aggregated across all provided genes and variants
    # Drug support for tier 4 will never be available
    for this_drug in drug_map.keys():
        for this_ct in drug_map[this_drug].keys():
            # Given the selected ct, count number of occurrences for each possible support type (if any)
            count_pos = drug_map[this_drug][this_ct].count("POSITIVE")
            count_neg = drug_map[this_drug][this_ct].count("NEGATIVE")
            count_unk = drug_map[this_drug][this_ct].count("UNKNOWN_BLANK")
            count_dns = drug_map[this_drug][this_ct].count("UNKNOWN_DNS")

            # Pool UNKNOWN_BLANK and UNKNOWN_DNS together (as both result in unknown CIViC support)
            count_total_unk = count_unk + count_dns
            # Sanity check that there is at least some support
            if (count_pos == 0) and (count_neg == 0) and (count_total_unk == 0):
                raise ValueError("Encountered unexpected support case when aggregating evidences across provided 'input_data'!")

            # Resolve contradicting evidence (if any) by majority vote
            temp_support = ""
            # For this, pool UNKNOWN_BLANK and UNKNOWN_DNS together
            # Whenever there is a tie of "confident" (pos or neg) vs "non-confident" (unk), choose the confident one
            if (count_total_unk > count_pos) and (count_total_unk > count_neg):
                temp_support = "CIVIC_UNKNOWN"
            elif count_pos == count_neg:
                temp_support = "CIVIC_CONFLICT"
            elif (count_pos > count_neg) and (count_pos >= count_total_unk):
                temp_support = "CIVIC_SUPPORT"
            elif (count_neg > count_pos) and (count_neg >= count_total_unk):
                temp_support = "CIVIC_RESISTANCE"
            else:
                raise ValueError("Encountered unexpected support case when aggregating evidences across provided 'input_data'!")

            # Build support string for each given combination of drug, ct and matched tier
            # Format: DRUG:CT:SUPPORT:#pos|#neg|#unk|#dns
            drug_support = this_drug + ":" + this_ct.upper() + ":" + temp_support + ":" + str(count_pos) + "|" + str(count_neg) + "|" + str(count_unk) + "|" + str(count_dns)
            drug_support_strings.append(drug_support)

    return drug_support_strings
