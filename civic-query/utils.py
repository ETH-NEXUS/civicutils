import sys
import os
import re

# import read_and_write
# dict_coodes = get_dict_aminoacids()


def check_identifier_type(identifier_type):
    """Check that a given identifier type is valid.

## TODO

    """

    if identifier_type not in ["entrez_id", "entrez_symbol", "civic_id"]:
        raise ValueError(
            f"'{identifier_type}' is not a valid identifier type.\n"
            f"Please, provide one of the following identifier types: 'entrez_id','entrez_symbol','civic_id'.\n"
        )


def check_empty_field(inField):
    """If a given inField is None or empty string, return 'NULL'; otherwise return the same inField string.
## TODO

    """

    if (inField is None) or (not inField):
        newField = "NULL"
    else:
        newField = inField
    return newField


# Check object is a list (even if empty)
def check_is_list(inList):
    if not isinstance(inList, list):
        raise TypeError(
            f"Provided 'inList' is not of type 'list': {inList}.\n"
        )

# Check object is a dict (even if empty)
def check_is_dict(inDict):
    if not isinstance(inDict, dict):
        raise TypeError(
            f"Provided 'inDict' is not of type 'dict': {inDict}.\n"
        )

def check_keys(inDict,keyList):
    check_is_dict(inDict)
    check_is_list(keyList)
    if set(inDict.keys()) != set(keyList):
        raise ValueError(
            f"Provided 'inDict' does not containg all of the following keys: {keyList}.\n"
        )

# Check object is a str (even if empty)
def check_is_str(inField):
    if not isinstance(inField, str):
        raise TypeError(
            f"Provided 'inField' is not of type 'str': {inField}.\n"
        )

def check_string_filter_arguments(inField, inList):
    # Check expected class was provided for each argument
    # inList should be a list (even if empty)
    check_is_list(inList)
    # inField should be a non-empty string
    check_is_str(inField)
    if not inField:
        raise ValueError(
            f"Please provided a non-empty 'inField'.\n"
        )
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
        raise ValueError(
            f"Variant-level entry '{entry}' could not be found!\n"
        )


# Possible values for select_tier can be a string or a list of tiers
#       - "highest" = highest encountered tier for each record, order 1>1b>2>3>4
#       - list of tiers eg. ["tier_1","tier_1b"]
#       - "all" = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]
#       - single tier eg. "tier1"
def check_tier_selection(select_tier,all_tiers):
    # Check for empty argument
    # Only valid argument types are str or list
    if select_tier is None:
        raise ValueError(
            f"Argument 'select_tier' must be provided: {select_tier}.\n"
            f"Please, provide either a str ('all','highest') or a list of tiers.\n"
        )
    elif isinstance(select_tier, str):
        # Check a valid value was provided for select_tier
        if select_tier not in ["all","highest"]:
            raise ValueError(
                f"'Unknown str option provided to argument 'select_tier': {select_tier}'.\n"
                f"Please, select either 'all' (to return all tiers) or 'highest' (to return only match for the highest tier).\n"
            )
    elif isinstance(select_tier, list):
        # Check that valid values were provided in select_tier
        unique_tiers = set(select_tier)
        if not unique_tiers:
            raise ValueError(
                f"Empty list provided to argument 'select_tier': {select_tier}.\n"
                f"Please, provide a non-empty list of tiers.\n"
            )
        for tmp_tier in list(unique_tiers):
            if tmp_tier not in all_tiers:
                raise ValueError(
                    f"Invalid tier provided to argument 'select_tier': {tmp_tier}.\n"
                    f"Please, provide a selection of valid tiers from: {all_tiers}.\n"
                )
        # Sanity check for cases when all tiers are provided in the list
        if unique_tiers == set(all_tiers):
            # Use correct select_tier value for this situation
            select_tier = "all"
        else:
            select_tier = list(unique_tiers)
    else:
        raise TypeError(
            f"Invalid type provided to argument 'select_tier': {select_tier}.\n"
            f"Please, provide either a str ('all','highest') or a list of tiers.\n"
        )
    return select_tier


# Translate a 1-letter aminoacid code (if it exists) into a 3-letter code
def translate_aa(oneLetter):
    aaNew = None
    if oneLetter.upper() in dict_codes,dict_codes:
        aaNew = dict_codes[oneLetter.upper()]
    return aaNew

# Given a single CIVIC variant name, extract potential HGVS annotations by
# parsing and modifying this string using knowledge on CIVIC naming conventions
# Function only applicable to SNV data
def civicName_to_hgvs(varName):
    # TODO: variant names like T193C (refer to DNA) would be translated to proteins
    hgvsStrings = []
    # 1) HGVS 1-letter protein code, including stop codons and general variants (eg. V600E, *600E, V600* or V600)
    # TODO: multiple cases, what to do? pos = re.findall
    if re.match(r'([A-Z*])(\d+)([A-Z*]?)($|\s\()', varName):
        pos = re.match(r'([A-Z*])(\d+)([A-Z*]?)($|\s\()', varName).groups()
        aa1 = pos[0]
        npos = pos[1]
        aa2 = pos[2]

        ## Generate the HGVS string by translating 1-letter aa codes into 3-letter

        # Special case of protein extensions where '*' is represented by 'Ter' in input table
        if aa1 == '*':
            aa1New = 'TER'
        else:
            # Translate 1-letter aa code to 3-letter code
            aa1New = translate_aa(aa1)
        # Check for general variants (second aa will be '')
        if aa2:
            # Special cases like p.Ter600Ter happen in input table, so in theory,
            # *600* should also be possible (however does not seem to happen in CIVIC)
            if (aa1New == 'TER') and (aa2=='*'):
                aa2New = 'TER'
            else:
                # Translate 1-letter aa code to 3-letter code
                # In cases where aa2='*', it will remain as '*' (using Ter is a special case when aa1=Ter)
                aa2New = translate_aa(aa2)
        else:
            aa2New = ''
        # Sanity check that translated strings correspond to valid aa codes
        if (aa1New is not None) and (aa2New is not None):
            # Construct new protein annotation
            newAnnotation = 'P.' + aa1New + npos + aa2New
            hgvsStrings.append(newAnnotation)

    # 2) Extract embedded c. annotation if available ie. '(c.XXX)'
    # Use re.search to allow for matches not in the string start
    if re.search(r'\((C\..+?)\)', varName):
        # TODO: at least two multiple cases, what to do?
        # 'A50A (c.150C>G); Splicing alteration (c.463-1G>T)'
        # r = re.findall(r'\((c\..+?)\)', varName)
        annot = re.search(r'\((C\..+?)\)', varName).groups()[0]
        hgvsStrings.append(annot)

    # 3) Frameshifts (eg. T157FS or T157MFS)
    if re.match(r'([A-Z])(\d+)([A-Z]?)FS', varName):
        pos = re.match(r'([A-Z])(\d+)([A-Z]?)FS', varName).groups()
        aa = pos[0]
        npos = pos[1]
        # Translate aa in 1st position and if ok, write frameshift in short form
        aaNew = translate_aa(aa)
        if aaNew is not None:
            annot = 'P.' + aaNew + npos + 'FS'
            hgvsStrings.append(annot)

        # TODO: extend cases?

    return hgvsStrings


# Given a single CIVIC hgvs expression, parse and modify it to ensure
# it complies with the HGVS format followed by the input table
# Function only applicable to SNV data
# TODO: Currently, only modification of p. annotations are supported
def civic_hgvs_to_input(civic_hgvs):
    ## The following "special hgvs" cases should be mutually exclusive
    ## ie. a single HGVS expression should match only 1 of the cases
    newAnnot = None

    ## Frameshift
    # Input table seems to use short form for frameshift annotations (ie. p.Glu55fs) while CIVIC does not
    if re.match(r'(P\.[A-Z]+[0-9]+)[A-Z]+FS.*', civic_hgvs):
        newAnnot = re.sub(r'(P\.[A-Z]+[0-9]+)[A-Z]+FS.*', r'\1FS', civic_hgvs)
        # Only return new HGVS expression if something changed
        if newAnnot != civic_hgvs:
            return newAnnot

    ## Nonsense mutations (gain of stop codon)
    # In general, CIViC seems to use 'Ter' for stop codons, only 1 exception so far (p.F76Lfs*56)
    # Input table uses '*' for stop codons; exception in protein extensions (eg. p.Ter370Tyrext*?)
    # '*' is also used for nucleotide numbering in c.,n., but not relevant here (only p.)
    if re.match(r'(P\.[A-Z]+[0-9]+)TER', civic_hgvs):
        newAnnot = re.sub(r'(P\.[A-Z]+[0-9]+)TER', r'\1*', civic_hgvs)
        # Only return new HGVS expression if something changed
        if newAnnot != civic_hgvs:
            return newAnnot

    ## Loss of stop codon
    # In CIVIC, a stop codon loss is expressed as p.Ter214Cys
    # In HGVS, this is enconded as an extension (eg. p.Ter370Tyrext*?)
    # In this case, we need to change annotation from input table (generate_civic_matchStrings)

    ## Silent mutations (no aa change)
    # CIVIC uses a '=' in the 2nd position to represent no change (eg. p.Pro61=)
    # In the input table, the 2nd aa is specified (eg. p.Pro61Pro)
    if re.match(r'P\.([A-Z]+)([0-9]+)=', civic_hgvs):
        newAnnot = re.sub(r'P\.([A-Z]+)([0-9]+)=', r'P.\1\2\1', civic_hgvs)
        # Only return new HGVS expression if something changed
        if newAnnot != civic_hgvs:
            return newAnnot

    return newAnnot


# Given a single HGVS p. annotation (of the form p.Pro61Cys...), extract the start of the string (ie. p.Pro61)
# Function only applicable to SNV data
def extract_p_start(pAnnotation):
    startAnnot = None
    if re.match(r'(P\.[A-Z]+[0-9]+)', pAnnotation):
        startAnnot = re.match(r'(P\.[A-Z]+[0-9]+)', pAnnotation).groups()[0]

    return startAnnot


# Given a single CIVIC variant name, return whether it corresponds to a CNV variant record related to exons
# For this, attemp to match the variant name to special CNV exon cases present in CIVIC (eg. EXON 1-2 DELETION, EXON 5 DELETION, 3' EXON DELETION...)
def cnv_is_exon_string(varName):
    exonStrings = ['^EXON [0-9-]+ DELETION$', '^[35\']+ EXON DELETION$', '^EXON [0-9-]+ SKIPPING MUTATION$']
    isExon = False
    for exonString in exonStrings:
        if re.search(exonString, varName):
            isExon = True

    return isExon

# Given a single CIVIC record name, return whether it corresponds to a EXPRESSION record related to exons
# For this, attemp to match the variant name to special EXPRESSION exon cases present in CIVIC (eg. EXON 1-2 EXPRESSION, EXON 5 OVEREXPRESSION...)
def expr_is_exon_string(varName):
    exonStrings = ['^EXON [0-9-]+ EXPRESSION$', '^EXON [0-9-]+ OVEREXPRESSION$', '^EXON [0-9-]+ UNDEREXPRESSION$']
    isExon = False
    for exonString in exonStrings:
        if re.search(exonString, varName):
            isExon = True

    return isExon



# Given all CIVIC variant records associated to a given gene, return all those variant names that should correspond to SNV records
# For this, remove the most common CNV and EXPRESSION names existing in CIVIC from the complete set of the gene's variant names
# Additionally, consider other special cases present in CIVIC (eg. EXON 1-2 DELETION/EXPRESSION, EXON 5 DELETION/OVEREXPRESSION, 3' EXON DELETION/UNDEREXPRESSION...)
def civic_return_all_snvs(geneData):
    # Common CNV record names in CIVIC
    cnvNames = ['AMPLIFICATION','DELETION','LOSS','COPY NUMBER VARIATION']
    # Common EXPRESSION record names in CIVIC
    exprNames = ['OVEREXPRESSION','UNDEREXPRESSION','EXPRESSION']

    # All variant names not matching a CNV or EXPRESSION will be returned
    matches = []
    for varName in list(geneData.keys()):
        ## Also attempt matching of variant name to common CNV and EXPRESSION names related to exons
        isExon_cnv = cnv_is_exon_string(varName)
        isExon_expr = expr_is_exon_string(varName)
        ## Skip variant names matching a CNV or EXPRESSION
        if (varName in cnvNames) or isExon_cnv:
            continue
        if (varName in exprNames) or isExon_expr:
            continue
        if varName not in matches:
            matches.append(varName)

    return matches


# Given all CIVIC variant records associated to a given gene, return all those variant names that correspond to CNV records
# For this, attemp to match the gene's variant names to the most common CNV names existing in CIVIC
# Additionally, consider other special CNV cases present in CIVIC (eg. EXON 1-2 DELETION, EXON 5 DELETION, 3' EXON DELETION...)
def civic_return_all_cnvs(geneData):
    # Common CNV record names in CIVIC
    cnvNames = ['AMPLIFICATION','DELETION','LOSS','COPY NUMBER VARIATION']

    # All matched variant names will be returned
    matches = []
    for varName in list(geneData.keys()):
        ## Also attempt matching of variant name to common CNV names related to exons
        isExon = cnv_is_exon_string(varName)
        if (varName in cnvNames) or isExon:
#             ## Special case for 'LOSS': only considered synonym of 'DELETION' when a certain variant type is present
#             if varName == 'LOSS' and ('TRANSCRIPT_ABLATION' not in geneData[varName]['types']):
#                 continue
            if varName not in matches:
                matches.append(varName)

    return matches


# Check whether a variant name in CIVIC refers to a group of variants (eg. V600)
# Function only applicable to SNV data
def check_general_variant(varName):
    isGeneral = False
    if re.match(r'[A-Z]\d+($|\s\()', varName):
        isGeneral = True
    return isGeneral



