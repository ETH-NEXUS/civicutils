import os
import sys
import yaml
import json

BinPath = os.path.split(os.path.realpath(__file__))[0]


def get_dict_aminoacids():
    """
    Retrieve dictionary from the data.yml file to translate between 1-letter and 3-letter aminoacid codes
    :return:	Dictionary of 1-letter to 3-letter aminoacid codes.
    """
    f = BinPath + "/data/data.yml"
    entry_name = "aminoacids"
    with open(f, "r") as infile:
        try:
            data = yaml.safe_load(infile)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(1)
    # Sanity check that expected entry is contained in the yml file
    if entry_name not in data.keys():
        raise ValueError("Please provide a dictionary of one-letter to three-letter aminoacid codes via the '%s' entry in %s!" %(entry_name,f))
    dict_codes = data[entry_name]

    return dict_codes


def get_dict_support():
    """
    Retrieve dictionary from the data.yml file to translate between CIViC evidence directions and clinical significances and their corresponding drug support (assigned by user)
    :return:	Dictionary of evidence to drug support.
    """
    f = BinPath + "/data/data.yml"
    entry_name = "drug_support"
    with open(f, "r") as infile:
        try:
            data = yaml.safe_load(infile)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(1)
    # Sanity check that expected entry is contained in the yml file
    if entry_name not in data.keys():
        raise ValueError("Please provide a custom dictionary of drug support for CIViC evidences via the '%s' entry in %s!" %(entry_name,f))
    support_dict = data[entry_name]

    return support_dict


def check_header_field(name, header_split, is_required=True):
    """
    Retrieve the position for a given column name from a list of tab-split header fields. Throw error when a required column is not found.
    :param name:		Column name to check in the header to retrieve its position.
    :param header_split:	List of columns from splitting a tab-separated header.
    :param is_required:		Boolean indicating if the given column to check is required. If required, throw an error when column is not found.
    :return:			Either numeric position of the column in the header or None (when non-required column is not found).
    """
    pos = None
    if name in header_split:
        pos = header_split.index(name)
    else:
        if is_required:
            raise ValueError("Required column '%s' could not be found in header '%s'" %(name," ".join(header_split)))

    return pos


def process_snv_header(header_split):
    """
    Retrieve the required column names expected for the header of an input SNV file. Assume column names are: Gene, Variant_dna, Variant_prot, Variant_impact and Variant_impact. Only the last two are not required.
    :param header_split:	List of columns from splitting a tab-separated header.
    :return:			Tuple of column positions for the gene, cHGVS, pHGVS, impact and exon, in that order. Only the last two can be None.
    """
    gene_pos = check_header_field("Gene", header_split, is_required=True)
    c_pos = check_header_field("Variant_dna", header_split, is_required=True)
    p_pos = check_header_field("Variant_prot", header_split, is_required=True)
    impact_pos = check_header_field("Variant_impact", header_split, is_required=False)
    exon_pos = check_header_field("Variant_exon", header_split, is_required=False)

    return (gene_pos, c_pos, p_pos, impact_pos, exon_pos)


def read_in_snvs(infile):
    """
    Read-in input file of SNV data and process it into structured dictionaries. Assumes header and that relevant info is contained in the following columns: Gene, Variant_dna, Variant_prot. Optional columns:  Variant_impact, Variant_exon.
    :param infile:    Path to the input SNV file to read in.
    :return:          Tuple of 3 elements: dictionary containing original rows and fields from input file (raw_data), dictionary of SNV data to be used to match in CIViC (snv_data) and list of additional columns provided in the input file (which are not required but should be written to output).
    """
    # dict n_line -> [gene, dna, prot, (impact, exon)]
    raw_data = {}

    # dict gene -> variant (dna|prot|impact|exon|n_line) -> null
    snv_data = {}

    in_file = open(infile, "r")
    header = in_file.readline().strip()
    header_split = header.strip().split("\t")
    (gene_pos, c_pos, p_pos, impact_pos, exon_pos) = process_snv_header(header_split)

    extra_header = []
    if impact_pos:
        extra_header.append("Variant_impact")
    if exon_pos:
        extra_header.append("Variant_exon")

    extra_pos = []
    for pos,x in enumerate(header_split):
        if (pos == gene_pos) or (pos == c_pos) or (pos == p_pos):
            continue
        if impact_pos:
            if pos == impact_pos:
                continue
        if exon_pos:
            if pos == exon_pos:
                continue
        extra_header.append(x)
        extra_pos.append(pos)

    for n_line,line in enumerate(in_file):
        line_split = line.strip().split("\t")
        c_var = line_split[c_pos].strip()
        p_var = line_split[p_pos].strip()
        gene = line_split[gene_pos].strip()
        impact = ""
        if impact_pos:
            impact = line_split[impact_pos].strip()
        exon = ""
        if exon_pos:
            exon = line_split[exon_pos].strip()
        raw_data[str(n_line)] = [gene, c_var, p_var, impact, exon]
        for p in extra_pos:
            raw_data[str(n_line)].append(line_split[p].strip())

        # Process raw_data to have gene-centered dict
        # Returns dict of gene -> [var1,var2,..,varN], where a given var="dna|prot|impact|exon|n_line"
        if gene not in snv_data.keys():
            snv_data[gene] = {}

        # Collapse variant info separated with "|"
        # Keep track of what line each variant comes from
        variant = c_var + "|" + p_var + "|" + impact + "|" + exon + "|" + str(n_line)
        # NOTE: Variants can never be duplicated because of the different row numbers assigned
        # if variant in snv_data[gene].keys():
        #     print("Found duplicated variant '%s|%s' for gene '%s' in line '%s'!" %(c_var, p_var, gene, str(n_line)))
        #     sys.exit(1)
        snv_data[gene][variant] = None
    in_file.close()

    return (raw_data, snv_data, extra_header)


def process_cnv_header(header_split):
    """
    Retrieve the required column names expected for the header of an input CNV file. Assume column names are: Gene, Variant_cnv (both are required). The type of CNV variant should be one of the following terms: 'GAIN', 'DUPLICATION', 'DUP', 'AMPLIFICATION' or 'AMP' (synonyms for AMPLIFICATION), and 'DELETION', 'DEL' or 'LOSS' (synonyms for DELETION).
    :param header_split:	List of columns from splitting a tab-separated header.
    :return:			Tuple of column positions for the gene and CNV type.
    """
    gene_pos = check_header_field("Gene", header_split, is_required=True)
    cnv_pos = check_header_field("Variant_cnv", header_split, is_required=True)

    return (gene_pos, cnv_pos)


def read_in_cnvs(infile):
    """
    Read-in input file of CNV data and process it into structured dictionaries. Assumes header and that relevant info is contained in the following columns: Gene, Variant_cnv.
    :param infile:    Path to the input CNV file to read in.
    :return:          Tuple of 3 elements: dictionary containing original rows and fields from input file (raw_data), dictionary of CNV data to be used to match in CIViC (cnv_data) and list of additional columns provided in the input file (which are not required but should be written to output).
    """
    # dict n_line -> [gene, cnv]
    raw_data = {}
    # dict gene -> variant (cnv|n_line) -> null
    cnv_data = {}

    in_file = open(infile, "r")
    header = in_file.readline().strip()
    header_split = header.strip().split("\t")
    (gene_pos, cnv_pos) = process_cnv_header(header_split)

    extra_header = []
    extra_pos = []
    for pos,x in enumerate(header_split):
        if (pos == gene_pos) or (pos == cnv_pos):
            continue
        extra_header.append(x)
        extra_pos.append(pos)

    for n_line,line in enumerate(in_file):
        line_split = line.strip().split("\t")
        gene = line_split[gene_pos].strip()
        cnv = line_split[cnv_pos].strip()
        raw_data[str(n_line)] = [gene,cnv]
        for p in extra_pos:
            raw_data[str(n_line)].append(line_split[p].strip())
        # Process raw_data to have gene-centered dict
        # Returns dict of gene -> [var1,var2,..,varN], where a given var="cnv|n_line"
        if gene not in cnv_data.keys():
            cnv_data[gene] = {}

        # Collapse variant info separated with "|"
        # Keep track of what line each variant comes from
        variant = cnv + "|" + str(n_line)
        # NOTE: Variants can never be duplicated because of the different row numbers assigned
        # if variant in cnv_data[gene].keys():
        #     print("Found duplicated variant '%s' for gene '%s' in line '%s'!" %(cnv,gene,str(n_line)))
        #     sys.exit(1)
        cnv_data[gene][variant] = None
    in_file.close()

    return (raw_data, cnv_data, extra_header)


def process_expr_header(header_split):
    """
    Retrieve the required column names expected for the header of an input Expression file. Assume column names are: Gene, logFC (both are required). The log fold-change value of a differentially expressed gene should be numeric and different from zero (to query CIViC, it will be translated into 'OVEREXPRESSION', if positive, and 'UNDEREXPRESSION', if negative).
    :param header_split:	List of columns from splitting a tab-separated header.
    :return:			Tuple of column positions for the gene and log fold-change.
    """
    gene_pos = check_header_field("Gene", header_split, is_required=True)
    logfc_pos = check_header_field("logFC", header_split, is_required=True)

    return (gene_pos, logfc_pos)


def read_in_expr(infile):
    """
    Read-in input file of differentially expressed data and process it into structured dictionaries. Assumes header and that relevant info is contained in the following columns: Gene, logFC.
    :param infile:    Path to the input expression file to read in.
    :return:          Tuple of 3 elements: dictionary containing original rows and fields from input file (raw_data), dictionary of EXPR data to be used to match in CIViC (expr_data) and list of additional columns provided in the input file (which are not required but should be written to output).
    """
    # dict n_line -> [gene, logfc]
    raw_data = {}
    # dict gene -> expression (logfc|n_line) -> null
    expr_data = {}

    in_file = open(infile, "r")
    header = in_file.readline().strip()
    header_split = header.strip().split("\t")
    (gene_pos, logfc_pos) = process_expr_header(header_split)

    extra_header = []
    extra_pos = []
    for pos,x in enumerate(header_split):
        if (pos == gene_pos) or (pos == logfc_pos):
            continue
        extra_header.append(x)
        extra_pos.append(pos)

    for n_line,line in enumerate(in_file):
        line_split = line.strip().split("\t")
        gene = line_split[gene_pos].strip()
        logfc = line_split[logfc_pos].strip()

        raw_data[str(n_line)] = [gene, logfc]
        for p in extra_pos:
            raw_data[str(n_line)].append(line_split[p].strip())
        # Process raw_data to have gene-centered dict
        # Returns dict of gene -> [var1,var2,..,varN], where a given var="logfc|n_line"
        if gene not in expr_data.keys():
            expr_data[gene] = {}

        # Collapse expression info separated with "|"
        # Keep track of what line each variant comes from
        expression = str(logfc) + "|" + str(n_line)
        # NOTE: Expression values can never be duplicated because of the different row numbers assigned
        # if expression in expr_data[gene].keys():
        #     print("Found duplicated expression value '%s' for gene '%s' in line '%s'!" %(logfc,gene,str(n_line)))
        #     sys.exit(1)
        expr_data[gene][expression] = None
    in_file.close()

    return (raw_data, expr_data, extra_header)


def write_to_json(in_dict, outfile, indent=1):
    """
    Write a nested dictionary into a new output file (in JSON format).
    :param in_dict:	Object of type 'dict' (can be nested).
    :param outfile:	Path to the output JSON file to write into.
    :param indent:	Value of indentation to use as margin.
    :return:		None
    """
    with open(outfile, "w") as f:
        json.dump(in_dict, f, ensure_ascii=False, indent=indent)

    return None


def write_to_yaml(in_dict, outfile):
    """
    Write a nested dictionary into a new output file (in YAML format).
    :param in_dict:	Object of type 'dict' (can be nested).
    :param outfile:	Path to the output YAML file to write into.
    :return:		None
    """
    with open(outfile, "w") as f:
        # Preserve the original order of the entries in the manifest template
        yaml.dump(in_dict, f, default_flow_style=False, sort_keys=False)

    return None


def write_header_line(data_type, header, write_support):
    """
    Given a list of sorted column names (from splitting the input header), process it to generate the corresponding output header.
    :param data_type:		['SNV', 'CNV', 'EXPR']
				SNV:   Expects a file of genomic single nucleotide variants and insertions/deletions
				CNV:   Expects a file of genomic copy number alterations
				EXPR:  Expects a file of differential gene expression data
				Data type of the corresponding input file (string).
    :param header:		List of input column names (in order) from splitting the tab-separated header.
    :param write_support:	Boolean indicating if processed drug support from CIViC should be written to output.
    :return:			Tuple of 4 elements (last 3 are only relevant when 'data_type=SNV'): string containing complete header ready to be written to output, list of column names excluding 'Variant_impact' and 'Variant_exon' (if they were present), boolean indicating if column 'Variant_impact' was present, boolean indicating if column '\tVariant_exon' was present.
    """
    sorted_evidence_types = ["PREDICTIVE", "DIAGNOSTIC", "PROGNOSTIC", "PREDISPOSING"]

    # Variables only relevant for data_type="SNV"
    write_impact = False
    write_exon = False

    if data_type == "SNV":
        main_header = "Gene\tVariant_dna\tVariant_prot"
        if "Variant_impact" in header:
            main_header += "\tVariant_impact"
            write_impact = True
        if "Variant_exon" in header:
            main_header += "\tVariant_exon"
            write_exon = True

    if data_type == "CNV":
        main_header = "Gene\tVariant_cnv"

    if data_type == "EXPR":
        main_header = "Gene\tlogFC"

    clean_header = []
    if header:
        for tmp in header:
            if (tmp != "Variant_impact" and tmp != "Variant_exon"):
                clean_header.append(tmp)
    if clean_header:
        main_header += "\t%s" %("\t".join(clean_header))

    if write_support:
        out_header = "%s\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\tCIViC_Drug_Support\t%s" %(main_header,"\t".join(["CIViC_" + x for x in sorted_evidence_types]))
    else:
        out_header = "%s\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\t%s" %(main_header,"\t".join(["CIViC_" + x for x in sorted_evidence_types]))

    return (out_header, clean_header, write_impact, write_exon)


def write_output_line(tier, main_line, gene_scores, gene_var_types, drug_support, result_map, write_support):
    """
    Process CIViC information available for a given input line, and generate the corresponding output line. Empty fields are indicated with '.'.
    :param tier:		String with format 'tier_[N]', where N can be: '1', '1b', '2', '3' or '4'. They indicate different tier cases of the variant matching.
    :param main_line:		String already containing the required columns to be written to output. This will be extended with the available CIViC information to generate the complete output line.
    :param gene_scores:		List containing the CIViC scores of the matched variants for the given 'tier'. Already formatted as required i.e. 'GENE:VARIANT_NAME:SCORE'.
    :param gene_var_types:	List containing the variant types of the matched variants in CIViC. Already formatted as required i.e. 'GENE:VARIANT_NAME:VAR_TYPE1,..'.
    :param drug_support:	List containing the strings of processed drug support for the matched variants. Already formatted as required i.e. 'DRUG:CT_TAG:CIVIC_CONSENSUS_PREDICTION'.
    :param result_map:		Nested dictionary containing the CIViC evidences associated with each evidence type, already formatted for writing to output (assumes a particular structure).
    :param write_support:	Boolean indicating if processed drug support from CIViC should be written to output.
    :return:			String containing complete output line ready to be written to output.
    """
    sorted_evidence_types = ["PREDICTIVE", "DIAGNOSTIC", "PROGNOSTIC", "PREDISPOSING"]

    # Remove the "tier" tag from the assined tier
    if tier.startswith("tier_"):
        tier = tier.replace("tier_", "")

    if gene_scores:
        out_line = main_line + "\t" + tier + "\t" + ";".join(gene_scores)
    else:
        out_line = main_line + "\t" + tier + "\t."

    if gene_var_types:
        out_line += "\t" + ";".join(gene_var_types)
    else:
        out_line += "\t."

    if write_support:
        if drug_support:
            out_line += "\t" + ";".join(drug_support)
        else:
            out_line += "\t."

    for evidence_type in sorted_evidence_types:
        # if evidence_type not in result_map.keys():
        #     raise ValueError("Evidence type '%s' is not found in provided 'result_map'!" %(evidence_type))
        if evidence_type in result_map.keys():
            if result_map[evidence_type]:
                out_line += "\t" + ";".join(result_map[evidence_type])
            else:
                out_line += "\t."
        else:
            out_line += "\t."

    return out_line


def write_evidences(item, write_drug=False, write_ct=None, write_complete=False):
    """
    Process CIViC information available for a given evidence type, and reformat into structured strings for reporting to output.
    Format: DISEASE[|DRUG1,DRUG2..](direction, significance(level(PMID,..,PMID),level(..)));
    :param item:		Nested dictionary containing CIViC evidences associated to a single evidence type (assumes a particular structure).
    :param write_drug:		Boolean indicating if the current type being evaluated corresponds to 'PREDICTIVE' (i.e. these are the only evidences associated to any drug names in CIViC).
    :param write_ct:		Boolean indicating if the available disease specificity annotations should be included in the output table.
    :param write_complete:	Boolean indicating if the complete information string of each CIViC evidence should be written to output. When 'write_complete=False', only the ids of the associated publications will be reported instead, already formatted as 'SOURCE_ID' (e.g. 'PUBMED_12345' or 'ASCO_12345').
    :return:			List of complete evidence strings to be written under the current evidence type (one per column).
    """
    evidences = [] 
    # For each disease found in the provided item
    for disease in item.keys():
        # For each drug associated with the given cancer type
        for drug in item[disease].keys():
            # For each evidence associated with the given drug
            # Evidences are simplified by using the combined form "direction:significance"
            for evidence in item[disease][drug].keys():
                # If drug=True, write drug information, i.e. DISEASE|DRUG(..)
                if write_drug:
                    # Always one drug (single or combination with "+")
                    if write_ct:
                        out_string = disease + "|" + write_ct.upper() + "|" + drug + "("
                    else:
                        out_string = disease + "|" + drug + "("
                else:
                    if write_ct:
                        out_string = disease + "|" + write_ct.upper() + "("
                    else:
                        out_string = disease + "("

                # Split the evidence direction and clinical significance
                evidence_list = evidence.strip().split(":")
                if (len(evidence_list) != 2):
                    raise ValueError("Unexpected format of evidence '%s'! Please provide string as 'EVIDENCE_DIRECTION:CLINICAL_SIGNIFICANCE'." %(evidence))
                direction = evidence_list[0]
                clin_signf = evidence_list[1]
                out_string += direction + "," + clin_signf + "("
                # There may be several levels grouped per evidence
                levels = []
                for level in item[disease][drug][evidence].keys():
                    pmids = []
                    # There may be several publications (i.e. PMIDs) grouped per level
                    for z in item[disease][drug][evidence][level]:
                        if write_complete:
                            pmids.append(z)
                        else:
                            z_split = z.strip().split(":")
                            if (len(z_split) != 5):
                                raise ValueError("Unexpected format of evidence item '%s'! Please provide string as 'EVIDENCE_ID:EVIDENCE_STATUS:SOURCE_STATUS:VARIANT_ORIGIN:RATING'." %(z))
                            pub_id = z_split[0].strip()
                            pmids.append(pub_id)
                    levels.append(level + "(" + ",".join(pmids) + ")")
                out_string += ",".join(levels) + "))"
                evidences.append(out_string)

    return evidences


def write_match(match_map, var_map, raw_map, header, data_type, outfile, has_support=True, has_ct=True, write_ct=False, write_support=True, write_complete=False):
    """
    Process input data and matched CIViC information, and reformat to write into a tab-separated table.
    :param match_map:		Dictionary containing data matched in CIViC (there must be a correspondance of 'match_map' and the variant data in 'var_map'). See README for more details about the specific structure of dictionary 'match_map'.
    :param var_map:		Nested dictionary of results from querying genes in CIViC (there must be a correspondance of 'var_map' and the variants matched in 'match_map'). See README for more details about the specific structure of dictionary 'var_map'.
    :param raw_map:		Dictionary containing the original rows and fields from the processed input file (n_line -> [field1,field2,..]).
    :param header:		List of input column names (in order) from splitting the tab-separated header.
    :param data_type:		['SNV', 'CNV', 'EXPR']
                             	SNV:   Expects a file of genomic single nucleotide variants and insertions/deletions
                             	CNV:   Expects a file of genomic copy number alterations
                             	EXPR:  Expects a file of differential gene expression data
                             	String data type of the corresponding input file.
    :param outfile:          	Path to the output file to write the matched CIViC evidences into (tab-separated table with header).
    :param has_support:       	Boolean indicating if the provided 'match_map' is annotated for drug support.
    :param has_ct:            	Boolean indicating if the provided 'var_map' is annotated for disease specificity.
    :param write_ct:          	Boolean indicating if the available disease specificity annotations should be included in the output table. To use this option, 'has_ct' must be True.
    :param write_support:     	Boolean indicating if the available drug support annotations should be included in the output table. To use this option, 'has_support' must be True.
    :param write_complete:    	Boolean indicating if the complete information string of each CIViC evidence should be written to output. When 'write_complete=False', only the ids of the associated publications will be reported instead, already formatted as 'SOURCE_ID' (e.g. 'PUBMED_12345' or 'ASCO_12345').
    :return:                 	None
    """
    # NOTE: uppercase is critical for matching!
    sorted_evidence_types = ["PREDICTIVE", "DIAGNOSTIC", "PROGNOSTIC", "PREDISPOSING"]
    evidence_type = "PREDICTIVE"
    special_cases = ["NON_SNV_MATCH_ONLY", "NON_CNV_MATCH_ONLY", "NON_EXPR_MATCH_ONLY"]
    sorted_cts = ["ct", "gt", "nct"]
    varmap_entries_variant = ["name", "hgvs", "types"]
    
    from civicutils.utils import check_match_before_writing, check_keys, check_keys_not, check_data_type, check_dict_entry
    check_match_before_writing(match_map, var_map, raw_map, has_support, has_ct, write_ct, write_support, write_complete)
    check_data_type(data_type)
    outfile = open(outfile, "w")
    
    # Retrieve the output header given the argument selection
    (out_header, clean_header, write_impact, write_exon) = write_header_line(data_type, header, write_support)
    outfile.write(out_header + "\n")
    
    for n_line in raw_map.keys():
        line_list = raw_map[n_line]
        extra_line = []
        if data_type == "SNV":
            if (len(line_list) < 5):
                raise ValueError("Must provide at least 5 elements to describe a SNV variant (even if some can be empty): 'gene,dna,[prot],[impact],[exon],..'")
            gene = line_list[0]
            c_var = line_list[1]
            p_var = line_list[2]
            impact = line_list[3]
            exon = line_list[4]
            comb_id = c_var + "|" + p_var + "|" + impact + "|" + exon + "|" + str(n_line)
            # Extract any additional fields that might be present for this line
            for pos in range(5, len(line_list)):
                extra_line.append(line_list[pos])
            # Build line string for writing to output
            main_line = gene + "\t" + c_var + "\t" + p_var
            if write_impact:
                main_line += "\t" + impact
            if write_exon:
                main_line += "\t" + exon

        if data_type == "CNV":
            if (len(line_list) < 2):
                raise ValueError("Must provide at least 2 elements to describe a CNV variant: 'gene,cnv,..'")
            gene = line_list[0]
            cnv = line_list[1]
            comb_id = cnv + "|" + str(n_line)
            # Extract any additional fields that might be present for this line
            for pos in range(2, len(line_list)):
                extra_line.append(line_list[pos])
            # Build line string for writing to output
            main_line = gene + "\t" + cnv

        if data_type == "EXPR":
            if (len(line_list) < 2):
                raise ValueError("Must provide at least 2 elements to describe a EXPR variant: 'gene,logFC,..'")
            gene = line_list[0]
            logfc = line_list[1]
            comb_id = str(logfc) + "|" + str(n_line)
            # Extract any additional fields that might be present for this line
            for pos in range(2, len(line_list)):
                extra_line.append(line_list[pos])
            # Build line string for writing to output
            main_line = gene + "\t" + str(logfc)

        # Sanity check that as many data fields were provided as in the header
        if len(extra_line) != len(clean_header):
            raise ValueError("Number of fields available does not match provided header!")
        # Add extra fields to the current line build
        for extra in extra_line:
            main_line += "\t" + extra

        # Check if match_map contains the provided input variants
        if gene not in match_map.keys():
            raise ValueError("Provided gene '%s' is not contained in 'match_map'." %(gene))
        if comb_id not in match_map[gene].keys():
            raise ValueError("Provided variant '%s' is not contained in 'match_map' of gene '%s'." %(comb_id,gene))
        for tier in match_map[gene][comb_id].keys():
            gene_scores = []
            gene_var_types = []
            drug_support = []
            result_map = {}
            write_line = False
            if tier != "tier_4":
                all_variants = []
                if has_support:
                    this_drug_support = match_map[gene][comb_id][tier]["drug_support"]
                    for tmp_var in match_map[gene][comb_id][tier]["matched"]:
                        all_variants.append(tmp_var)
                    if write_support:
                        for i in this_drug_support:
                            drug_support.append(i.upper())
                else:
                    if write_support:
                        raise ValueError("Option 'write_support' cannot be selected when 'has_support'=False!")
                    for tmp_var in match_map[gene][comb_id][tier]:
                        all_variants.append(tmp_var)

                for var_id in all_variants:
                    # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                    # This is a dummy tag and not an actual variant record from CIViC, so skip checking in var_map
                    if var_id.upper() in special_cases:
                        # In this case, current line will be associated with tier3, but all columns will be empty with "."
                        for evidence_type in sorted_evidence_types:
                            if evidence_type not in result_map.keys():
                                result_map[evidence_type] = []
                        continue

                    variant = var_map[gene][var_id]["name"]
           
                    gene_var_types.append(gene + ":" + variant + ":" + ",".join(var_map[gene][var_id]["types"]))
                    molecular_profile_ids = set(list(var_map[gene][var_id].keys())) ^ set(varmap_entries_variant)
                    for molecular_profil_id in molecular_profile_ids:
                        gene_scores.append(gene + ":" + variant + ":" + molecular_profil_id + ":" + str(var_map[gene][var_id][molecular_profil_id]["civic_score"]))
                        for evidence_type in sorted_evidence_types:
                            if evidence_type in var_map[gene][var_id][molecular_profil_id]["evidence_items"].keys():
                                if evidence_type not in result_map.keys():
                                    result_map[evidence_type] = []
                                write_drug = False
                                if evidence_type == evidence_type:
                                    write_drug=True
                                if has_ct:
                                    check_keys(list(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts, matches_all=True)
                                    for ct in var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys():
                                        if write_ct:
                                            results = write_evidences(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type][ct], write_drug=write_drug, write_ct=ct, write_complete=write_complete)
                                        else:
                                            results = write_evidences(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type][ct], write_drug=write_drug, write_ct=None, write_complete=write_complete)
                                        for x in results:
                                            result_map[evidence_type].append(gene + ":" + variant + ":" + molecular_profil_id + ":" + x)
                                else:
                                    check_keys_not(list(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts)
                                    if write_ct:
                                        raise ValueError("Option 'write_ct' cannot be selected when 'has_ct'=False!")
                                    results = write_evidences(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type], write_drug=write_drug, write_ct=None, write_complete=write_complete)
                                    for x in results:
                                        result_map[evidence_type].append(gene + ":" + variant + ":" + molecular_profil_id + ":" + x)

                # Only write line current tier when there was at least one match for it
                if all_variants:
                    write_line = True

            else:
                if has_support:
                    if match_map[gene][comb_id][tier]["matched"]:
                        write_line = True
                else:
                    if write_support:
                        raise ValueError("Option 'write_support' cannot be selected when 'has_support'=False!")
                    if match_map[gene][comb_id][tier]:
                        write_line = True

            if write_line:
                out_line = write_output_line(tier, main_line, gene_scores, gene_var_types, drug_support, result_map, write_support)
                outfile.write(out_line + "\n")
    outfile.close()

    return None
