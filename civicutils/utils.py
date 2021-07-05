import sys
import os
import re

def check_is_none(argument,argName):
    """
    Check that a given argument is not None.
    :param argument:    An object to check.
    :param argName:     The name of the argument being checked (for error handling).
    :return:            None
    """
    if argument is None:
        raise ValueError("Argument '%s' must be provided!" %(argName))
    return None


def check_argument(argument,argName):
    """
    Check that a given argument is not None and is non-empty.
    :param argument:    An object to check.
    :param argName:     The name of the argument being checked (for error handling).
    :return:            None
    """
    check_is_none(argument,argName)
    if not argument:
        raise ValueError("Argument '%s' must be provided!" %(argName))
    return None


def check_arguments(argList,nameList):
    """
    Given a list of arguments, check that all exist and are non-empty.
    :param argList:     A list of objects to check.
    :param nameList:    A list of object names (for error handling). There must be a 1-1 correspondance of elements with argList.
    :return:            None
    """
    check_argument(argList,"argList")
    check_argument(nameList,"nameList")
    check_is_list(argList,"argList")
    check_is_list(nameList,"nameList")
    if (len(argList) != len(nameList)):
        raise ValueError("Arguments 'argList' and 'nameList' must have the same length!")
    for i,argument in enumerate(argList):
        check_argument(argument,nameList[i])
    return None


def check_identifier_type(identifier_type):
    """
    Check that a given identifier type is valid.
    :param identifier_type:    ['entrez_symbol', 'entrez_id', 'civic_id']
                        entrez_symbol:   Entrez gene symbol.
                        entrez_id: Entrez gene identifier.
                        civic_id: CIViCdb internal identifier.
                        String identifier type to check for validity.
    :return:            None
    """
    check_argument(identifier_type,"identifier_type")
    check_is_str(identifier_type,"identifier_type")
    idTypes = ["entrez_id", "entrez_symbol", "civic_id"]
    if identifier_type not in idTypes:
        raise ValueError("'%s' is not a valid identifier_type. Please provide one of: %s" %(identifier_type,idTypes))
    return None


def check_data_type(dataType):
    """
    Check that a given data type is valid.
    :param dataType:    ['SNV', 'CNV', 'EXPR']
                        SNV:   Expect genomic single nucleotide variants and insertions/deletions.
                        CNV:   Expect genomic copy number variation data.
                        EXPR:  Expect differential expression gene data.
                        String data type to check for validity.
    :return:            None
    """
    check_argument(dataType,"dataType")
    check_is_str(dataType,"dataType")
    dataTypes = ["SNV", "CNV", "EXPR"]
    if dataType not in dataTypes:
        raise ValueError("'%s' is not a valid dataType. Please provide one of: %s" %(dataType,dataTypes))
    return None


def check_empty_field(inField):
    """
    Check if the given field is None or empty, and return 'NULL' in this case; otherwise, return the given field.
    :param inField:    Object to check.
    :return:           Either 'NULL' (if object is empty or None) or the input object.
    """
    if (inField is None) or (not inField):
        newField = "NULL"
    else:
        newField = inField
    return newField


def check_empty_input(inField, fieldName, isRequired=True):
    """
    Check if the given field is None, "." or empty, and return a empty string "" in this case. If not empty, return the same given field.
    :param inField:     Field to check for being empty.
    :param fieldName:   Name of the field being checked.
    :param isRequired:  Boolean indicating if an error should be thrown when the field is empty.
    :return:            Either an empty string (if object is empty, "." or None) or the input object.
    """
    check_is_str(inField,fieldName)
    if (inField is None) or (not inField) or (inField == "."):
        if isRequired:
            raise ValueError("'%s' is required and cannot be empty!" %(fieldName))
        newField = ""
    else:
        newField = inField
    return newField


def parse_input(inField, fieldName, isRequired=True):
    """
    Given a list of comma-separated fields, check that all exist and are non-empty and return in a list
    :param inField:      A string to check.
    :param fieldName:    Name of the field being checked.
    :param isRequired:   A boolean indicating if error should be thrown when a field is empty.
    :return:             List of elements retrieved from parsing and checking 'inField'.
    """
    out = []
    newField = check_empty_input(inField, fieldName, isRequired)
    splitArr = newField.split(",")
    for tmpSplit in splitArr:
        new = check_empty_input(tmpSplit, fieldName, isRequired)
        if isRequired:
            if new and (new not in out):
                out.append(new)
        else:
            out.append(new)
    return out


def check_logFC(logfc,gene):
    """
    Check that a given logFC value for a gene is valid.
    :param logfc:      A log fold-change value to check (must be numeric or float).
    :param gene:       Name of the gene associated with the log fold-change (for error handling).
    :return:           Float log fold-change.
    """
    check_argument(logfc,"logFC")
    try:
        logfc = float(logfc)
    except ValueError:
        print("Invalid logFC = '%s' for gene '%s'. Please provide a numeric value." %(logFC,gene))
        sys.exit(1)
    return logfc


def check_is_bool(inBool,boolName):
    """
    Check an object is of type Boolean.
    :param inBool:      Object to be checked for Boolean type.
    :param boolName:    Name of the object being checked (for error handling).
    :return:            None
    """
    if not isinstance(inBool, bool):
        raise TypeError("'%s' is not of type 'bool'" %(boolName))
    return None


def check_is_list(inList,listName):
    """
    Check an object is of type list (even if empty).
    :param inList:      Object to be checked for list type.
    :param listName:    Name of the object being checked (for error handling).
    :return:            None
    """
    if not isinstance(inList, list):
        raise TypeError("'%s' is not of type 'list'" %(listName))
    return None


def check_is_dict(inDict,dictName):
    """
    Check an object is of type dict (even if empty).
    :param inDict:      Object to be checked for dict type.
    :param dictName:    Name of the object being checked (for error handling).
    :return:            None
    """
    if not isinstance(inDict, dict):
        raise TypeError("'%s' is not of type 'dict'" %(dictName))
    return None


def check_keys(inKeys,dictName,keyList,matches_all=True):
    """
    Check if a list of keys matches those provided in a list (retrieved from a dict).
    :param inKeys:         List of keys to check.
    :param dictName:       Name of the dict being checked for containing the keys.
    :param keyList:        List of keys to compare against.
    :param matches_all:    Boolean indicating if the set of keys in 'inKeys' should be identical to those in 'keyList'.
    :return:               None
    """
    check_arguments([inKeys,keyList],[dictName,"keyList"])
    check_is_list(inKeys,dictName)
    check_is_list(keyList,"keyList")

    in1 = set(inKeys)
    in2 = set(keyList)
    if matches_all:
        if (in1 != in2):
            raise ValueError("Dictionary '%s' does not contain all of the following keys: %s" %(dictName,keyList))
    else:
        notFound = False
        for x in in2:
            if x not in inKeys:
                raise ValueError("Dictionary '%s' does not contain the following key: %s" %(dictName,x))
    return None


def check_keys_not(inKeys,dictName,keyList):
    """
    Check a given list of keys are not contained in those provided in a list (retrieved from a dict).
    :param inKeys:         List of keys to check.
    :param dictName:       Name of the dict being checked for not containing the keys.
    :param keyList:        List of keys to compare against.
    :return:               None
    """
    check_arguments([inKeys,keyList],[dictName,"keyList"])
    check_is_list(inKeys,dictName)
    check_is_list(keyList,"keyList")
    for key in keyList:
        if key in inKeys:
            raise ValueError("'%s' cannot contain key '%s'!" %(dictName,key))
    return None


def check_is_str(inField,fieldName):
    """
    Check an object is of type str (even if empty).
    :param inField:      Object to be checked for str type.
    :param fieldName:    Name of the object being checked (for error handling).
    :return:             None
    """
    if not isinstance(inField, str):
        raise TypeError("'%s' is not of type 'str'" %(fieldName))
    return None


def uppercase_list(inList,listName):
    """
    Uppercase all strings contained in a list.
    :param inList:         List of strings to uppercase.
    :param listName:       Name of the list being uppercased (for error handling).
    :return:               New list of uppercase elements.
    """
    # inList should be a list (even if empty)
    check_is_list(inList,listName)
    newList = []
    # Convert individual elements into strings, in case provided list contains numbers (eg. gene or variant ids)
    for tmpItem in inList:
        newItem = str(tmpItem)
        newItem = newItem.strip().upper()
        newList.append(newItem)
    return newList


def check_string_filter_arguments(inField, fieldName, inList, listName):
    """
    Check the arguments provided to a filtering function for having the correct types and format.
    :param inField:        Field to be checked for being of str type and uppercase.
    :param fieldName:      Name of the string field being checked.
    :param inList:         List to be checked and turn to uppercase.
    :param listName:       Name of the list being checked.
    :return:               Tuple of uppercase field and uppercase list.
    """
    # Check expected class was provided for each argument
    # inField should be a non-empty string
    check_is_str(inField,fieldName)
    check_empty_input(inField,fieldName,isRequired=True)
    # Use uppercase and remove leading/trailing spaces, for consistency of strings
    newField = inField.strip().upper()
    newList = uppercase_list(inList,listName)
    return (newField,newList)


def check_cutoff_filter_arguments(value,name):
    """
    Check that a cutoff argument provided to a filtering function is a valid number.
    :param value:     Value to be checked for being a valid number.
    :param name:      Name of the value being checked.
    :return:          Value turned into float type.
    """
    # cutoff can be 0
    if (value is None):
        raise ValueError("Argument '%s' must be provided!" %(name))
    try:
        cutoff_f = float(value)
    except ValueError:
        print("Invalid '%s'! Please provide a numeric value." %(name,value))
        sys.exit(1)
    return cutoff_f


def check_dict_entry(myDict,dictName,entry,entryName):
    """
    Check that a given entry is a key in the provided dictionary.
    :param myDict:     Dictionary to be checked for containing the key.
    :param dictName:   Name of the dictionary being checked.
    :param entry:      Key to check in the dictionary.
    :param entryName:  Name of the entry being checked.
    :return:           None
    """
    check_arguments([myDict,entry],[dictName,entryName])
    check_is_dict(myDict,dictName)
    check_is_str(entry,entryName)
    if entry not in myDict.keys():
        raise ValueError("Could not find %s '%s' in dict %s!" %(entryName,entry,dictName))
    return None


def check_tier_selection(select_tier,all_tiers):
    """
    Check that the provided tier selection for filtering is correct.
    :param select_tier:     Either list or string, indicating the tier selection. Possible values are: "highest" (select the highest encountered tier, hierarchy 1>1b>2>3>4), "all" (do not filter) or a list with the individual tiers to select for (if all are provided, then no filtering is applied).
    :param all_tiers:       Sorted list of all available tiers (order must correspond to priority hierarchy).
    :return:                Processed list or string of tier selection.
    """
    # Check provided argument 'all_tiers'
    check_argument(all_tiers,"all_tiers")
    check_is_list(all_tiers,"all_tiers")

    # Check for empty argument
    # Only valid argument types are str or list
    new_selection = None
    if (select_tier is None) or (not select_tier):
        raise ValueError("Please provide a non-empty tier selection!")
    elif isinstance(select_tier, str):
        # Check a valid value was provided for select_tier
        if select_tier not in ["all","highest"]:
            raise ValueError("Unknown tier option provided: '%s'. Possible options are: 'all' (to return all tiers) or 'highest' (to return only match for the highest tier)." %(select_tier))
        new_selection = select_tier
    elif isinstance(select_tier, list):
        unique_tiers = set(select_tier)
        # Check that valid values were provided in select_tier
        for tmp_tier in list(unique_tiers):
            if tmp_tier not in all_tiers:
                raise ValueError("Provided list contains invalid tier: '%s'. Possible options are: %s" %(tmp_tier,all_tiers))
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
        raise TypeError("Please provide either a str ('all','highest') or a list of tiers.")
    return new_selection


def check_is_cHGVS(hgvs):
    """
    Check that a given HGVS string corresponds to a coding DNA reference sequence ie. starts with "c."
    :param hgvs:     HGVS string to be checked.
    :return:         None
    """
    if not hgvs.startswith("c."):
        raise ValueError("HGVS string '%s' does not start with 'c.'!" %(hgvs))
    return None


def check_is_pHGVS(hgvs):
    """
    Check that a given HGVS string corresponds to a protein reference sequence ie. starts with "p."
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
    # Import the dictionary of aminoacid codes provided in the data.yml file
    dict_codes = get_dict_aminoacids()
    aaNew = None
    check_empty_input(aminoacid,"aminoacid",isRequired=True)
    if aminoacid.upper() in dict_codes.keys():
        aaNew = dict_codes[aminoacid.upper()]
    return aaNew


def check_match_before_writing(matchMap, varMap, rawMap, hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False):
    """
    Check that the relevant arguments provided to write_match() are valid and have the expected structure.
    :param matchMap:        Dictionary containing data matched in CIVIC.
    :param varMap:          Dictionary containing CIVIC records from querying genes in the database.
    :param rawMap:          Dictionary containing input data from parsing file.
    :param hasSupport:      Boolean indicating if the provided 'matchMap' is annotated for drug support.
    :param hasCt:           Boolean indicating if the provided 'varMap' is annotated for disease specificity.
    :param writeCt:         Boolean indicating if the available disease specificity annotations should be included in the output table. To use this option, 'hasCt' must be True.
    :param writeSupport:    Boolean indicating if the available drug support annotations should be included in the output table. To use this option, 'hasSupport' must be True.
    :param writeComplete:   Boolean indicating if the complete evidence annotations should be included in the output table, instead of the reduced version containing only the publication ids.
    :return:                None
    """
    sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]
    varmap_entries = ['name','civic_score','hgvs','types','n_evidence_items','evidence_items']
    special_cases = ["NON_SNV_MATCH_ONLY","NON_CNV_MATCH_ONLY","NON_EXPR_MATCH_ONLY"]

    check_arguments([matchMap,varMap,rawMap],["matchMap","varMap","rawMap"])
    check_is_none(hasSupport,"hasSupport")
    check_is_none(hasCt,"hasCt")
    check_is_none(writeCt,"writeCt")
    check_is_none(writeSupport,"writeSupport")
    check_is_none(writeComplete,"writeComplete")
    check_is_dict(matchMap,"matchMap")
    check_is_dict(varMap,"varMap")
    check_is_dict(rawMap,"rawMap")
    check_is_bool(hasSupport,"hasSupport")
    check_is_bool(hasCt,"hasCt")
    check_is_bool(writeCt,"writeCt")
    check_is_bool(writeSupport,"writeSupport")
    check_is_bool(writeComplete,"writeComplete")

    for gene in matchMap.keys():
        matched = []
        for variant in matchMap[gene].keys():
            check_keys(list(matchMap[gene][variant].keys()),"matchMap",sorted_tiers,matches_all=True)
            for tier in matchMap[gene][variant].keys():
                if hasSupport:
                    check_keys(list(matchMap[gene][variant][tier].keys()),"matchMap",["matched","drug_support"],matches_all=True)
                    check_is_list(matchMap[gene][variant][tier]["drug_support"],tier)
                    if tier == "tier_4":
                        check_is_bool(matchMap[gene][variant][tier]["matched"],tier)
                    else:
                        check_is_list(matchMap[gene][variant][tier]["matched"],tier)
                        for tmpVar in matchMap[gene][variant][tier]["matched"]:
                            # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                            if tmpVar.upper() in special_cases:
                                # Check that no other variant was matched when this special case was matched (length of matches should always be one)
                                if len(matchMap[gene][variant][tier]["matched"]) != 1:
                                    raise ValueError("Unexpected: encountered multiple matches in special case of empty tier3 match '%s'!" %(matchMap[gene][variant][tier]["matched"]))
                                continue
                            if tmpVar not in matched:
                                matched.append(tmpVar)
                else:
                    if tier == "tier_4":
                        check_is_bool(matchMap[gene][variant][tier],tier)
                    else:
                        check_is_list(matchMap[gene][variant][tier],tier)
                        for tmpVar in matchMap[gene][variant][tier]:
                            # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                            if tmpVar.upper() in special_cases:
                                # Check that no other variant was matched when this special case was matched (length of matches should always be one)
                                if len(matchMap[gene][variant][tier]) != 1:
                                    raise ValueError("Unexpected: encountered multiple matches in special case of empty tier3 match '%s'!" %(matchMap[gene][variant][tier]))
                                continue
                            if tmpVar not in matched:
                                matched.append(tmpVar)
            if matched:
                check_dict_entry(varMap,"varMap",gene,"gene")
            for varId in matched:
                check_dict_entry(varMap[gene],"varMap",varId,"variant")
                check_keys(list(varMap[gene][varId].keys()),"varMap",varmap_entries,matches_all=True)
# TODO: check evidence types are correct? Since new categories were added recently in CIVICdb

    return None
