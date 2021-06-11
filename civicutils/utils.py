import sys
import os
import re

def check_identifier_type(identifier_type):
    """Check that a given identifier type is valid.

## TODO

    """

    if identifier_type not in ["entrez_id", "entrez_symbol", "civic_id"]:
# FIXME
        raise ValueError("")
#             f"'{identifier_type}' is not a valid identifier type.\n"
#             f"Please, provide one of the following identifier types: 'entrez_id','entrez_symbol','civic_id'.\n"


def check_empty_field(inField):
    """If a given inField is None or empty string, return 'NULL'; otherwise return the same inField string.
## TODO

    """

    if (inField is None) or (not inField):
        newField = "NULL"
    else:
        newField = inField
    return newField


def check_empty_input(inField, fieldName, isRequired=True):
    """If a given inField is None, '.' or empty string, return empty ''; otherwise return the same inField string. For required fields, throw a warning
    """
    if (inField is None) or (not inField) or (inField == "."):
        if isRequired:
# TODO
            raise ValueError("Field '%s' is required and cannot be empty!" %(fieldName))
        newField = ""
    else:
        newField = inField
    return newField


def parse_input(inField, fieldName, isRequired=True):
    """TODO
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


# Check object is a list (even if empty)
def check_is_list(inList):
    if not isinstance(inList, list):
        raise TypeError("")
        #     f"Provided 'inList' is not of type 'list': {inList}.\n"

# Check object is a dict (even if empty)
def check_is_dict(inDict):
    if not isinstance(inDict, dict):
        raise TypeError("")
#             f"Provided 'inDict' is not of type 'dict': {inDict}.\n"

def check_keys(inDict,keyList):
    check_is_dict(inDict)
    check_is_list(keyList)
    if set(inDict.keys()) != set(keyList):
        raise ValueError("")
#             f"Provided 'inDict' does not containg all of the following keys: {keyList}.\n"

# Check object is a str (even if empty)
def check_is_str(inField):
    if not isinstance(inField, str):
        raise TypeError("")
#             f"Provided 'inField' is not of type 'str': {inField}.\n"

def check_string_filter_arguments(inField, inList):
    # Check expected class was provided for each argument
    # inList should be a list (even if empty)
    check_is_list(inList)
    # inField should be a non-empty string
    check_is_str(inField)
    if not inField:
        raise ValueError("")
#             f"Please provided a non-empty 'inField'.\n"
    # Use uppercase and remove leading/trailing spaces, for consistency of strings
    newField = inField.strip().upper()
    newList = []
    # Convert individual elements into strings, in case provided list contains numbers (eg. gene or variant ids)
    for tmpItem in inList:
        newItem = str(tmpItem)
        newItem = newItem.strip().upper()
        newList.append(newItem)
    return (newField,newList)


# FIXME generalize and include sanity check in all functions handling dicts
def check_variant_level_entry(subDict,entry):
# TODO: sanity check that arguments are of the correct type

    if entry not in subDict.keys():
        raise ValueError("")
#             f"Variant-level entry '{entry}' could not be found!\n"


# Possible values for select_tier can be a string or a list of tiers
#       - "highest" = highest encountered tier for each record, order 1>1b>2>3>4
#       - list of tiers eg. ["tier_1","tier_1b"]
#       - "all" = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]
#       - single tier eg. "tier1"
def check_tier_selection(select_tier,all_tiers):
    # Check for empty argument
    # Only valid argument types are str or list
    if select_tier is None:
        raise ValueError("")
#             f"Argument 'select_tier' must be provided: {select_tier}.\n"
#             f"Please, provide either a str ('all','highest') or a list of tiers.\n"
    elif isinstance(select_tier, str):
        # Check a valid value was provided for select_tier
        if select_tier not in ["all","highest"]:
            raise ValueError("")
#                 f"'Unknown str option provided to argument 'select_tier': {select_tier}'.\n"
#                 f"Please, select either 'all' (to return all tiers) or 'highest' (to return only match for the highest tier).\n"
    elif isinstance(select_tier, list):
        # Check that valid values were provided in select_tier
        unique_tiers = set(select_tier)
        if not unique_tiers:
            raise ValueError("")
#                 f"Empty list provided to argument 'select_tier': {select_tier}.\n"
#                 f"Please, provide a non-empty list of tiers.\n"
        for tmp_tier in list(unique_tiers):
            if tmp_tier not in all_tiers:
                raise ValueError("")
#                     f"Invalid tier provided to argument 'select_tier': {tmp_tier}.\n"
#                     f"Please, provide a selection of valid tiers from: {all_tiers}.\n"
        # Sanity check for cases when all tiers are provided in the list
        if unique_tiers == set(all_tiers):
            # Use correct select_tier value for this situation
            select_tier = "all"
        else:
            select_tier = list(unique_tiers)
    else:
        raise TypeError("")
#             f"Invalid type provided to argument 'select_tier': {select_tier}.\n"
#             f"Please, provide either a str ('all','highest') or a list of tiers.\n"
    return select_tier


# Sanity check that a given HGVS string corresponds to a coding DNA reference sequence ie. starts with "c."
def check_is_cHGVS(hgvs):
    if not hgvs.startswith("c."):
        raise ValueError("HGVS string '%s' does not start with 'c.'!" %(hgvs))
    return None


# Sanity check that a given HGVS string corresponds to a protein reference sequence ie. starts with "p."
def check_is_pHGVS(hgvs):
    if not hgvs.startswith("p."):
        raise ValueError("HGVS string '%s' does not start with 'p.'!" %(hgvs))
    return None


# Translate a 1-letter aminoacid code (if it exists) into a 3-letter code
def translate_aa(oneLetter):
    aaNew = None
    if oneLetter.upper() in dict_codes.keys():
        aaNew = dict_codes[oneLetter.upper()]
    return aaNew
