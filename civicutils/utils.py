import sys
import os
import re


def check_is_none(argument, arg_name):
    """
    Check that a given argument is not None.
    :param argument:	Object to check.
    :param arg_name:	Name of the argument being checked (for error handling).
    :return:		None
    """
    if argument is None:
        raise ValueError("Argument '%s' must be provided!" %(arg_name))

    return None


def check_argument(argument, arg_name):
    """
    Check that a given argument is not None and is non-empty.
    :param argument:	Object to check.
    :param arg_name:	Name of the argument being checked (for error handling).
    :return:            None
    """
    check_is_none(argument, arg_name)
    if not argument:
        raise ValueError("Argument '%s' must be provided!" %(arg_name))

    return None


def check_arguments(arg_list, name_list):
    """
    Given a list of arguments, check that all exist and are non-empty.
    :param arg_list:	List of objects to check.
    :param name_list:	List of object names (for error handling). There must be a 1-1 correspondance of elements with arg_list.
    :return:		None
    """
    check_argument(arg_list, "arg_list")
    check_argument(name_list, "name_list")
    check_is_list(arg_list, "arg_list")
    check_is_list(name_list, "name_list")
    if (len(arg_list) != len(name_list)):
        raise ValueError("Arguments 'arg_list' and 'name_list' must have the same length!")
    for i,argument in enumerate(arg_list):
        check_argument(argument, name_list[i])

    return None


def check_identifier_type(identifier_type):
    """
    Check that a given identifier type is valid.
    :param identifier_type:	['entrez_symbol', 'entrez_id', 'civic_id']
                        	entrez_symbol:	Entrez gene symbol
                        	entrez_id:	Entrez gene identifier
                        	civic_id:	CIViC internal identifier
                        	String identifier type to check for validity.
    :return:            	None
    """
    check_argument(identifier_type, "identifier_type")
    check_is_str(identifier_type, "identifier_type")
    id_types = ["entrez_id", "entrez_symbol", "civic_id"]
    if identifier_type not in id_types:
        raise ValueError("'%s' is not a valid identifier_type. Please provide one of: %s" %(identifier_type, id_types))

    return None


def check_data_type(data_type):
    """
    Check that a given data type is valid.
    :param data_type:	['SNV', 'CNV', 'EXPR']
			SNV:   Expects genomic single nucleotide variants and insertions/deletions
			CNV:   Expects genomic copy number alterations
			EXPR:  Expects differential gene expression data
			String data type to check for validity.
    :return:            None
    """
    check_argument(data_type, "data_type")
    check_is_str(data_type, "data_type")
    data_types = ["SNV", "CNV", "EXPR"]
    if data_type not in data_types:
        raise ValueError("'%s' is not a valid data_type. Please provide one of: %s" %(data_type, data_types))

    return None


def check_empty_field(in_field):
    """
    Check if the given field is None or empty, and return 'NULL' in this case; otherwise, return the given field.
    :param in_field:	Object to check.
    :return:		Either 'NULL' (if object is empty or None) or the input object.
    """
    if (in_field is None) or (not in_field):
        new_field = "NULL"
    else:
        new_field = in_field

    return new_field


def check_empty_input(in_field, field_name, is_required=True):
    """
    Check if the given field is None, "." or empty, and return a empty string "" in this case. If not empty, return the same given field.
    :param in_field:	Field to check for being empty.
    :param field_name:	Name of the field being checked.
    :param is_required:	Boolean indicating if an error should be thrown when the field is empty.
    :return:		Either an empty string (if object is empty, "." or None) or the input object.
    """
    check_is_str(in_field, field_name)
    if (in_field is None) or (not in_field) or (in_field == "."):
        if is_required:
            raise ValueError("'%s' is required and cannot be empty!" %(field_name))
        new_field = ""
    else:
        new_field = in_field

    return new_field


def parse_input(in_field, field_name, is_required=True):
    """
    Given a list of comma-separated fields, check that all exist and are non-empty and return in a list
    :param in_field:	String to check.
    :param field_name:	Name of the field being checked.
    :param is_required:	Boolean indicating if error should be thrown when a field is empty.
    :return:		List of elements retrieved from parsing and checking 'in_field'.
    """
    out = []
    new_field = check_empty_input(in_field, field_name, is_required)
    split_list = new_field.split(",")
    for tmp_split in split_list:
        new = check_empty_input(tmp_split, field_name, is_required)
        if is_required:
            if new and (new not in out):
                out.append(new)
        else:
            out.append(new)

    return out


def check_logfc(logfc, gene):
    """
    Check that a given logFC value for a gene is valid.
    :param logfc:      Log fold-change value to check (must be numeric or float).
    :param gene:       Name of the gene associated with the log fold-change (for error handling).
    :return:           Float log fold-change.
    """
    check_argument(logfc, "logfc")
    try:
        logfc = float(logfc)
    except ValueError:
        print("Invalid logFC = '%s' for gene '%s'. Please provide a numeric value." %(logfc, gene))
        sys.exit(1)

    return logfc


def check_is_bool(in_bool, bool_name):
    """
    Check an object is of type Boolean.
    :param in_bool:	Object to be checked for Boolean type.
    :param bool_name:	Name of the object being checked (for error handling).
    :return:            None
    """
    if not isinstance(in_bool, bool):
        raise TypeError("'%s' is not of type 'bool'" %(bool_name))

    return None


def check_is_list(in_list, list_name):
    """
    Check an object is of type list (even if empty).
    :param in_list:	Object to be checked for list type.
    :param list_name:	Name of the object being checked (for error handling).
    :return:		None
    """
    if not isinstance(in_list, list):
        raise TypeError("'%s' is not of type 'list'" %(list_name))

    return None


def check_is_dict(in_dict, dict_name):
    """
    Check an object is of type dict (even if empty).
    :param in_dict:	Object to be checked for dict type.
    :param dict_name:	Name of the object being checked (for error handling).
    :return:		None
    """
    if not isinstance(in_dict, dict):
        raise TypeError("'%s' is not of type 'dict'" %(dict_name))

    return None


def check_keys(in_keys, dict_name, key_list, matches_all=True):
    """
    Check if a list of keys matches those provided in a list (retrieved from a dict).
    :param in_keys:		List of keys to check.
    :param dict_name:		Name of the dict being checked for containing the keys.
    :param key_list:		List of keys to compare against.
    :param matches_all:		Boolean indicating if the set of keys in 'in_keys' should be identical to those in 'key_list'.
    :return:			None
    """
    check_arguments([in_keys, key_list], [dict_name, "key_list"])
    check_is_list(in_keys, dict_name)
    check_is_list(key_list, "key_list")

    in1 = set(in_keys)
    in2 = set(key_list)
    if matches_all:
        if (in1 != in2):
            raise ValueError("Dictionary '%s' does not contain all of the following keys: %s" %(dict_name, key_list))
    else:
        not_found = False
        for x in in2:
            if x not in in_keys:
                raise ValueError("Dictionary '%s' does not contain the following key: %s" %(dict_name, x))

    return None


def check_keys_not(in_keys, dict_name, key_list):
    """
    Check a given list of keys are not contained in those provided in a list (retrieved from a dict).
    :param in_keys:	List of keys to check.
    :param dict_name:	Name of the dict being checked for not containing the keys.
    :param key_list:	List of keys to compare against.
    :return:		None
    """
    check_arguments([in_keys, key_list], [dict_name, "key_list"])
    check_is_list(in_keys, dict_name)
    check_is_list(key_list, "key_list")
    for key in key_list:
        if key in in_keys:
            raise ValueError("'%s' cannot contain key '%s'!" %(dict_name, key))

    return None


def check_is_str(in_field, field_name):
    """
    Check an object is of type str (even if empty).
    :param in_field:	Object to be checked for 'str' type.
    :param field_name:	Name of the object being checked (for error handling).
    :return:		None
    """
    if not isinstance(in_field, str):
        raise TypeError("'%s' is not of type 'str'" %(field_name))

    return None


def uppercase_list(in_list, list_name):
    """
    Uppercase all strings contained in a list.
    :param in_list:	List of strings to uppercase.
    :param list_name:	Name of the list being uppercased (for error handling).
    :return:		New list of uppercase elements.
    """
    # in_list should be a list (even if empty)
    check_is_list(in_list, list_name)
    new_list = []
    # Convert individual elements into strings, in case provided list contains numbers (e.g. gene or variant ids)
    for tmp_item in in_list:
        new_item = str(tmp_item)
        new_item = new_item.strip().upper()
        new_list.append(new_item)

    return new_list


def check_string_filter_arguments(in_field, field_name, in_list, list_name):
    """
    Check the arguments provided to a filtering function for having the correct types and format.
    :param in_field:	Field to be checked for being of str type and uppercase.
    :param field_name:	Name of the string field being checked.
    :param in_list:	List to be checked and turn to uppercase.
    :param list_name:	Name of the list being checked.
    :return:		Tuple of uppercase field and uppercase list.
    """
    # Check expected class was provided for each argument
    # in_field should be a non-empty string
    check_is_str(in_field, field_name)
    check_empty_input(in_field, field_name, is_required=True)
    # Use uppercase and remove leading/trailing spaces, for consistency of strings
    new_field = in_field.strip().upper()
    new_list = uppercase_list(in_list, list_name)

    return (new_field, new_list)


def check_cutoff_filter_arguments(value, name):
    """
    Check that a cutoff argument provided to a filtering function is a valid number.
    :param value:	Value to be checked for being a valid number.
    :param name:	Name of the value being checked.
    :return:		Value turned into float type.
    """
    # cutoff can be 0
    if (value is None):
        raise ValueError("Argument '%s' must be provided!" %(name))
    try:
        cutoff_f = float(value)
    except ValueError:
        print("Invalid '%s'! Please provide a numeric value." %(name, value))
        sys.exit(1)

    return cutoff_f


def check_dict_entry(input_dict, dict_name, entry, entry_name):
    """
    Check that a given entry is a key in the provided dictionary.
    :param input_dict:	Dictionary to be checked for containing the key.
    :param dict_name:	Name of the dictionary being checked.
    :param entry:	Key to check in the dictionary.
    :param entry_name:	Name of the entry being checked.
    :return:		None
    """
    check_arguments([input_dict, entry], [dict_name, entry_name])
    check_is_dict(input_dict, dict_name)
    check_is_str(entry, entry_name)
    if entry not in input_dict.keys():
        raise ValueError("Could not find %s '%s' in dict %s!" %(entry_name, entry, dict_name))

    return None


def check_tier_selection(select_tier, all_tiers):
    """
    Check that the provided tier selection for filtering is correct.
    :param select_tier:		Either list or string, indicating the tier selection. Possible values are: 'highest' (select the highest encountered tier, hierarchy 1>1b>2>3>4), 'all' (do not filter) or a list with the individual tiers to select for (if all are provided, then no filtering is applied).
    :param all_tiers:		Sorted list of all available tiers (order must correspond to priority hierarchy).
    :return:			Processed list or string of tier selection.
    """
    # Check provided argument
    check_argument(all_tiers, "all_tiers")
    check_is_list(all_tiers, "all_tiers")

    # Check for empty argument
    # Only valid argument types are str or list
    new_selection = None
    if (select_tier is None) or (not select_tier):
        raise ValueError("Please provide a non-empty tier selection!")
    elif isinstance(select_tier, str):
        # Check a valid value was provided for select_tier
        if select_tier not in ["all", "highest"]:
            raise ValueError("Unknown tier option provided: '%s'. Possible options are: 'all' (to return all tiers) or 'highest' (to return only match for the highest tier)." %(select_tier))
        new_selection = select_tier
    elif isinstance(select_tier, list):
        unique_tiers = set(select_tier)
        # Check that valid values were provided in select_tier
        for tmp_tier in list(unique_tiers):
            if tmp_tier not in all_tiers:
                raise ValueError("Provided list contains invalid tier: '%s'. Possible options are: %s" %(tmp_tier, all_tiers))
        new_selection = []
        # Sanity check for cases when all tiers are provided in the list
        if unique_tiers == set(all_tiers):
            # Use correct select_tier value for this situation
            new_selection = "all"
        else:
            for tmp_tier in all_tiers:
                if tmp_tier in unique_tiers: 
                    new_selection.append(tmp_tier)
    else:
        raise TypeError("Please provide either a str ('all', 'highest') or a list of tiers.")

    return new_selection


def check_is_chgvs(hgvs):
    """
    Check that a given HGVS string corresponds to a coding DNA reference sequence i.e. starts with 'c.'
    :param hgvs:     HGVS string to be checked.
    :return:         None
    """
    if not hgvs.startswith("c."):
        raise ValueError("HGVS string '%s' does not start with 'c.'!" %(hgvs))

    return None


def check_is_phgvs(hgvs):
    """
    Check that a given HGVS string corresponds to a protein reference sequence i.e. starts with 'p.'
    :param hgvs:     HGVS string to be checked.
    :return:         None
    """
    if not hgvs.startswith("p."):
        raise ValueError("HGVS string '%s' does not start with 'p.'!" %(hgvs))

    return None


def translate_aa(aminoacid):
    """
    Given a 1-letter aminoacid code, check if it is a valid one and translate into a 3-letter code.
    :param aminoacid:     Aminoacid code in 1-letter format.
    :return:              Translated 3-letter aminoacid code.
    """
    from read_and_write import get_dict_aminoacids
    # Import the dictionary of aminoacid codes provided in the config "data.yml" file
    dict_codes = get_dict_aminoacids()
    aa_new = None
    check_empty_input(aminoacid, "aminoacid", is_required=True)
    if aminoacid.upper() in dict_codes.keys():
        aa_new = dict_codes[aminoacid.upper()]

    return aa_new


def check_match_before_writing(match_map, var_map, raw_map, has_support=True, has_ct=True, write_ct=False, write_support=True, write_complete=False):
    """
    Check that the relevant arguments provided to write_match() are valid and have the expected structure.
    :param match_map:        Dictionary containing data matched in CIViC. See more details about the specific structure of 'match_map' in the README.
    :param var_map:          Dictionary containing CIViC records from querying genes in the database. See more details about the specific structure of 'var_map' in the README.
    :param raw_map:          Dictionary containing input data from parsing the input file.
    :param has_support:      Boolean indicating if the provided 'match_map' is annotated for drug support.
    :param has_ct:           Boolean indicating if the provided 'var_map' is annotated for disease specificity.
    :param write_ct:         Boolean indicating if the available disease specificity annotations should be included in the output table. To use this option, 'has_ct' must be True.
    :param write_support:    Boolean indicating if the available drug support annotations should be included in the output table. To use this option, 'has_support' must be True.
    :param write_complete:   Boolean indicating if the complete evidence annotations should be included in the output table, instead of the reduced version containing only the publication ids.
    :return:                 None
    """
    sorted_tiers = ["tier_1", "tier_1b", "tier_2", "tier_3", "tier_4"]
    var_map_entries_variant = ["name", "hgvs", "types"]
    special_cases = ["NON_SNV_MATCH_ONLY", "NON_CNV_MATCH_ONLY", "NON_EXPR_MATCH_ONLY"]

    check_arguments([match_map, raw_map], ["match_map", "raw_map"])
    check_is_none(has_support, "has_support")
    check_is_none(has_ct, "has_ct")
    check_is_none(write_ct, "write_ct")
    check_is_none(write_support, "write_support")
    check_is_none(write_complete, "write_complete")
    check_is_dict(match_map, "match_map")
    check_is_dict(var_map, "var_map")
    check_is_dict(raw_map, "raw_map")
    check_is_bool(has_support, "has_support")
    check_is_bool(has_ct, "has_ct")
    check_is_bool(write_ct, "write_ct")
    check_is_bool(write_support, "write_support")
    check_is_bool(write_complete, "write_complete")

    for gene in match_map.keys():
        matched = []
        for variant in match_map[gene].keys():
            check_keys(list(match_map[gene][variant].keys()), "match_map", sorted_tiers, matches_all=True)
            for tier in match_map[gene][variant].keys():
                if has_support:
                    check_keys(list(match_map[gene][variant][tier].keys()), "match_map", ["matched", "drug_support"], matches_all=True)
                    check_is_list(match_map[gene][variant][tier]["drug_support"], tier)
                    if tier == "tier_4":
                        check_is_bool(match_map[gene][variant][tier]["matched"], tier)
                    else:
                        check_is_list(match_map[gene][variant][tier]["matched"], tier)
                        for tmp_var in match_map[gene][variant][tier]["matched"]:
                            # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                            if tmp_var.upper() in special_cases:
                                # Check that no other variant was matched when this special case was matched (length of matches should always be one)
                                if len(match_map[gene][variant][tier]["matched"]) != 1:
                                    raise ValueError("Unexpected: encountered multiple matches in special case of empty tier3 match '%s'!" %(match_map[gene][variant][tier]["matched"]))
                                continue
                            if tmp_var not in matched:
                                matched.append(tmp_var)
                else:
                    if tier == "tier_4":
                        check_is_bool(match_map[gene][variant][tier], tier)
                    else:
                        check_is_list(match_map[gene][variant][tier], tier)
                        for tmp_var in match_map[gene][variant][tier]:
                            # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                            if tmp_var.upper() in special_cases:
                                # Check that no other variant was matched when this special case was matched (length of matches should always be one)
                                if len(match_map[gene][variant][tier]) != 1:
                                    raise ValueError("Unexpected: encountered multiple matches in special case of empty tier3 match '%s'!" %(match_map[gene][variant][tier]))
                                continue
                            if tmp_var not in matched:
                                matched.append(tmp_var)
            if matched:
                check_dict_entry(var_map, "var_map", gene, "gene")
            for var_id in matched:
                check_dict_entry(var_map[gene], "var_map", var_id, "variant")
                check_keys(list(var_map[gene][var_id].keys()), "var_map", var_map_entries_variant, matches_all=False)

    return None
