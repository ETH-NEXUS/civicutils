#!/usr/bin/env python

'''
Utility functions
Lourdes Rosano, November 2020
'''


'''
Required Python modules
'''

import sys
import re

'''
Required internal functions
'''

import read_and_write


'''
Required data
'''

# TODO define dict_codes as global variable that translate_aa can access
dict_coodes = get_dict_aminoacids()


'''
Functions
'''


# Given a list of gene identifiers, remove those that are empty or contain any dots "." (as they cause the query to fail). Also, remove any duplicated identifiers.
def validate_genes(genes):
    out = []    # keep track of unique gene symbols to query
    skip = []   # keep track of unique gene symbols to skip
    # Iterate individual genes
    for gene in genes:
        # Check that gene is not empty
        if gene:
            # Check if current gene should be skipped
            if "." not in gene:
                if gene not in out:
                    out.append(gene)
                else:
                    print("Warning! Removed duplicated gene '%s'!" %(gene))
            else:
                if gene not in skip:
                    skip.append(gene)
        else:
            print("Warning! Removed empty gene '%s'!" %(gene))
    if skip:
        print("Omitted gene symbols: %s" %(",".join(skip)))
    return(out)


def get_field_from_record(record, field_name):
    # Check that provided field_name can be found in record
    if field_name not in record.keys():
        print("Error! Provided field '%s' could not be found in the following record: %s." %(field_name,record))
        sys.exit(1)
    else:
        field = record[field_name]
    return field

# Empty list of records will return empty dict
def get_field_from_records(records, id_field, field_name):
    dictFields = {}
    for record in records:
        thisId = get_field_from_record(record,id_field)
        field = get_field_from_record(record,field_name)
        if thisId not in dictFields.keys():
            dictFields[thisId] = field
        else:
            print("Warning! Skipped record for duplicated id '%s'!" %(thisId))
    return dictFields


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


# Given a CIVIC variant name and its corresponding hgvs expressions (if available),
# generate a list of potential strings that will be used to match the variant to our input
# For CNVs, variant matching is not based on HGVS, since input data is different
def generate_civic_matchStrings(varName, hgvsExpressions, dataType):
    matchStrings = []
    ## For CNVs, only the last step 5) is executed, since variant matching is not done at the HGVS level
    if dataType == 'SNV':
        # 1) First, remove reference from annotation (ie. 'transcriptID:')
        # This step will be skipped for CNVs
        for x in hgvsExpressions:
            # CIViC HGVS expressions are of the form reference + single annotation
            # eg. [NM_007313.2:c.1001C>T, NP_005148.2:p.Thr315Ile]
            new = x.split(':')[-1].upper()
            if new not in matchStrings:
                matchStrings.append(new)
                # 2) Second, generate modified HGVS strings that comply with input table format
                newAnnot = civic_hgvs_to_input(new)
                # Resulting annotation will be None unless something was modified (special cases only)
                if (newAnnot is not None) and (newAnnot not in matchStrings):
                    matchStrings.append(newAnnot)
        # 3) Third, generate a list of potential HGVS strings based on the actual CIVIC variant name
        potentialHGVS = civicName_to_hgvs(varName)
        for x in potentialHGVS:
            # Sometimes, duplicated HGVS annotations will get generated (eg. V600E already has p.Val600Glu)
            if x not in matchStrings:
                matchStrings.append(x)
        # 4) Last, add potential positional matches for already existing strings (will only affect p. annotations)
        for x in matchStrings:
            start = extract_p_start(x)
            # Only p. annotations will return a positional string (eg. p.Val600Glu -> p.Val600)
            if (start is not None) and (start not in matchStrings): 
                matchStrings.append(start)
    # 5) Also, add CIVIC variant name to allow for matching using input 'descriptional' strings (eg. EXON 15 MUTATION)
    ## Only this step is executed for CNV
    matchStrings.append(varName)

    return matchStrings


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


# Given a set of one or more annotations (SNV or CNV), generate a list of potential strings
# that will be used to match the variant in CIVIC eg. EXON 15 MUTATION, AMPLIFICATION
def generate_input_matchStrings(varAnnotations, dataType, impactAnnots, exonAnnots):
    matchStrings = []
    # Used for SNV: Keep track of whether a string corresponds to an exact (True) or positional match (False)
    isExact = []
    # Used for SNV: isTrueExact can be used together with isExact to distinguish between a true exact match
    # (isExact=True, isTrueExact=True) corresponding to input HGVS strings, or a more general match
    # (isExact=True, isTrueExact=False) corresponding to descriptive terms eg. EXON 1 MUTATION with lower preference
    isTrueExact = []

    ## a) Given a set of HGVS expressions, impacts and exon information for a single variant, generate a list of
    ## potential strings that will be used to match the variant in CIVIC eg. EXON 15 MUTATION, TRUNCATING FRAMESHIFT
    if dataType == 'SNV':
        # 1) First, add all (unique) HGVS annotations gathered from input table
        # If matched, they will correspond to exact matches (True in isExact)
        for varAnnot in varAnnotations:
            if varAnnot not in matchStrings:
                matchStrings.append(varAnnot)
                isExact.append(True)
                isTrueExact.append(True)
                # Special case for protein extensions (eg. p.Ter130Tyrext*?) in input table
                # Subset string to match CIVIC's format (p.Ter130Tyr)
                if re.match(r'(P\.TER[0-9]+[A-Z]+)EXT', varAnnot):
                    newVarAnnot = re.match(r'(P\.TER[0-9]+[A-Z]+)EXT', varAnnot).groups()[0]
                    if newVarAnnot not in matchStrings:
                        matchStrings.append(newVarAnnot)
                        isExact.append(True)
                        isTrueExact.append(True)
        # 2) Second, add potential positional matches for already existing strings (will only affect p. annotations)
        # If matched, they will correspond to positional matches (False in isExact)
        for x in matchStrings:
            start = extract_p_start(x)
            # Only p. annotations will return a positional string (eg. p.Val600Glu -> p.Val600)
            if (start is not None) and (start not in matchStrings): 
                matchStrings.append(start)
                isExact.append(False)
                isTrueExact.append(False)
        # 3) Last, add potential variant synonyms for matching to CIVIC's record names (eg. EXON 15 MUTATION)
        # If matched, they will correspond to exact matches (isExact=True) but not to true exact matches (isTrueExact=False)
        ## Sanity check that both lists (impact and exon annotations) are the same length
        if len(impactAnnots) != len(exonAnnots):
            print("\nError! Number of variant impacts and exon annotations are not equal")
            sys.exit(1)
        newTags = []
        ## Always include 'MUTATION' as a potential variant tag
        newTags.append('MUTATION')

        ## Iterate exon and impact information to generate and include additional variant tags
        for i,impact in enumerate(impactAnnots):
            ## Introduce sanity check on variant impact being available
            if not impact:
                continue
            ## First, generate potential additional variant tags based on the variant impact
            ## (can be contained in impact eg. 'frameshift_variant&stop_gained')
            if re.search('3_PRIME_UTR_VARIANT',impact):
                newTags.append("3' UTR MUTATION")

            if re.search('5_PRIME_UTR_VARIANT',impact):
                newTags.append("5' UTR MUTATION")

            if re.search('STOP_GAINED',impact):
                newTags.append('TRUNCATING MUTATION')

            if re.search('FRAMESHIFT_VARIANT',impact):
                newTags.append('FRAMESHIFT MUTATION')

            ## Second, get rank (exon/intron) information to generate an additional tag
            rank = exonAnnots[i]
            ## Skip cases where exon information is not available
            ## eg. with impacts downstream_gene_variant, upstream_gene_variant, intergenic_region
            if not rank:
                continue
            ## Get the actual rank (eg. 2/30)
            rank = rank.split('/')[0]
            ## After checking some examples, it seems that intron information is associated to variant impacts 'intron_variant'
            ## and 'sequence_feature' (can be contained in impact eg. 'splice_donor_variant&intron_variant')
            ## NOTE: any variant impact other than the above two will generate a tag 'EXON XX MUTATION',
            ## even 'synonymous' mutations
            if re.search('INTRON_VARIANT',impact) or re.search('SEQUENCE_FEATURE',impact):
                newTags.append('INTRON ' + rank + ' MUTATION')
            else:
                newTags.append('EXON ' + rank + ' MUTATION')
                ## As of 13/02/2019, only one case of variant name 'EXON XX FRAMESHIFT' (gene CALR)
                if re.search('FRAMESHIFT_VARIANT',impact):
                    newTags.append('EXON ' + rank + ' FRAMESHIFT')

        ## Add unique tags and assign corresponding information about type of match
        for tag in newTags:
            if tag not in matchStrings:
                matchStrings.append(tag)
                isExact.append(True)
                isTrueExact.append(False)

    ## For data type CNV, follow a different procedure
    else:
        newTags = []

        ## b) Given a single CNV annotation (AMP, GAIN, DEL, LOSS), generate a list of potential strings that will be used to match CNV in CIVIC
        if dataType == 'CNV':
            ## Opposite to the SNV case, for CNVs, varAnnotations should correspond to a single element (ie. CNV category)
            ## CIVIC seems to consider that GAIN and AMP are the same CNV
            if (varAnnotations == 'AMP') or (varAnnotations == 'GAIN'):
                newTags.append('AMPLIFICATION')
            ## CIVIC seems to consider that DELETION and LOSS are the same CNV
            ## Both CNV categories 'DELETION' and 'LOSS' occur in CIVIC
            elif (varAnnotations == 'DEL') or (varAnnotations == 'LOSS'):
                newTags.append('DELETION')
                newTags.append('LOSS')
            newTags.append('COPY NUMBER VARIATION')

        for tag in newTags:
            if tag not in matchStrings:
                matchStrings.append(tag)
                isExact.append(True)
                isTrueExact.append(True)

    return matchStrings,isExact,isTrueExact


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


# Given a list of variant annotations for a gene in the input table (HGVS but also synonym descriptive terms eg. EXON 15 MUTATION and positional
# strings eg. p.Val600), and all CIVIC variant records for said gene, attempt mapping of any input annotations to one or more CIVIC records.
# Return tier of the match (exact, positional, no match) and corresponding CIVIC variant records that were matched

# Do exhaustive search as there could be >1 perfect match in cases of redundancy (eg. E55FS and E55RFSTER11
# both translate to p.Glu55fs) and there can be >1 positional matches (for SNV case, 'general' variants
# eg. V600 have preference over the rest). Also, for SNVs in case of having two types of exact matches (truly exact
# eg. V600E and descriptive term eg. EXON 15 MUTATION), give preference to the former as descriptive terms should only
# be used when no true exact match was found.

# For CNV, inputStrings will be limited (usually only 1 or 2), and all will always correspond to exact matches (no positional)
def matchVariants(inputStrings, isExact, isTrueExact, geneData, dataType):
    matchMap = {}
    exactMatches = []
    synMatches = []
    posMatches = []
    # Iterate CIVIC variants for given gene and attempt match of any input string to any CIVIC string
    for varName in geneData.keys():
        # Retrieve all available CIVIC strings (HGVS, manually generated, name, etc.) for the given CIVIC variant
        # For CNV, civicStrings will contain one single element corresponding to the CIVIC variant name
        civicStrings = geneData[varName]['match_hgvs']
        # Iterate input annotations strings and attempt match to any
        # The position of each input string corresponds to the type of match (true exact, synonym exact or positional)
        for indx,inputAnnot in enumerate(inputStrings):
            if inputAnnot in civicStrings:
#                 ## For CNV: especial case for cnv 'LOSS', depending on variant type we should disregard the match or not
#                 if (dataType == 'CNV') and (inputAnnot == 'LOSS'):
#                     ## CIVIC CNVs of type 'LOSS' will only be considered when they are associated to variant type 'transcript_ablation'
#                     if 'TRANSCRIPT_ABLATION' not in geneData[varName]['types']:
#                         ## Disregard incorrect CNV match
#                         continue

                # Determine type of match using its position and store accordingly
                # For CNV, all matches will be exact and true exact matches (either variant is in CIVIC or not)
                if isExact[indx]:
                    # For SNV: give preference to truly exact matches over descriptive synonym terms
                    # eg. report V600E over EXON 15 MUTATION
                    if isTrueExact[indx]:
                        # For SNV: True exact match (eg. V600E)
                        # For CNV: all matches will be exact and true exact matches (either variant is in CIVIC or not)
                        if varName not in exactMatches:
                            exactMatches.append(varName)
                    else:
                        # For SNV: descriptive term match (eg. EXON 15 MUTATION, FRAMESHIFT MUTATION)
                        if varName not in synMatches:
                            synMatches.append(varName)
                else:
                    # Positional match
                    if varName not in posMatches:
                        posMatches.append(varName)
            ## For CNV: When there is no exact match for a 'DELETION', also consider special CIVIC CNV records related to exons (these will be positional matches)
            else:
                ## TODO CNV: Generalize for DELETION and AMPLIFICATION as well
                if (dataType == 'CNV') and (inputAnnot == 'DELETION'):
                    ## In CNV case, civicStrings corresponds to the CIVIC variant name (single string)
                    for tempString in civicStrings:
                        ## Look for special cases like eg. 'EXON 5 DELETION', 'EXON 1-2 DELETION' or '3' EXON DELETION'
                        isExon = cnv_is_exon_string(tempString)
                        if isExon:
                            ## These special cases are accounted for as positional matches (tier 2)
                            if varName not in posMatches:
                                posMatches.append(varName)

    ## Once all CIVIC variants have been iterated, determine final tier and corresponding matched variants
    ## For CNV: either exactMatches or posMatches will occurr due to the implementation design

    # If there is at least one truly exact match, return this (highest preference)
    # For CNV: all exact matches will also be true exact matches (either CNV is in CIVIC or not). Multiple matches could occur for a single CNV (eg. DELETION + LOSS + COPY NUMBER)
    if exactMatches:
        tier = 1
        # There could be >1 perfect matches
        matchedVars = exactMatches
    # If there are no true exact matches but at least one synonym exact match, return this (still exact match)
    elif synMatches:
        tier = 1
        # There could be >1 perfect matches
        matchedVars = synMatches
    # If there are no exact matches of any kind but at least one positional, return this
    elif posMatches:
        tier = 2
        # For CNV: positional matches will correspond to EXON records (eg. EXON 1-2 DELETION, 3' EXON DELETION..)
        # There can be >1 positional matches
        matchedVars = posMatches
        # For SNV: check for presence of preferential positional matches, ie. general variants (like V600)
        if dataType == 'SNV':
            for var in posMatches:
                isGeneral = check_general_variant(var)
                if isGeneral:
                    # Stop as soon as a general variant is found and report only this
                    matchedVars = [var]
                    break

    # If no match was found, then tier case is 3, and return all relevant variants
    else:
        tier = 3
        ## For SNV: when input variant could not be matched in CIVIC (tier3), return all CIVIC variants associated
        ## to the given gene but that do not correspond to a CNV or EXPRESSION related variant
        if dataType == 'SNV':
#             matchedVars = list(geneData.keys())
            matchedVars = civic_return_all_snvs(geneData)
        ## For CNV: when input cnv could not be matched in CIVIC (tier3), return all CIVIC cnvs associated the given gene (if any)
        ## ie. not all CIVIC records are returned but only those that are matched to a 'CNV' tags
        elif dataType == 'CNV':
            matchedVars = civic_return_all_cnvs(geneData)

    return tier,matchedVars


# Check whether a variant name in CIVIC refers to a group of variants (eg. V600)
# Function only applicable to SNV data
def check_general_variant(varName):
    isGeneral = False
    if re.match(r'[A-Z]\d+($|\s\()', varName):
        isGeneral = True
    return isGeneral

