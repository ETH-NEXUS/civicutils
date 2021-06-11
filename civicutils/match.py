import sys
import os
import re

from utils import translate_aa,check_is_dict,check_is_list,check_is_str,check_keys,check_tier_selection,parse_input,check_empty_input,check_dict_entry,check_is_cHGVS,check_is_pHGVS,check_identifier_type,check_data_type
# import utils
from query import query_civic
from filtering import filter_civic


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


# Check whether a variant name in CIVIC refers to a group of variants (eg. V600)
# Function only applicable to SNV data
def check_general_variant(varName):
    isGeneral = False
    if re.match(r'[A-Z]\d+($|\s\()', varName):
        isGeneral = True
    return isGeneral


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
    # Get all variant names classified as CNV or EXPRESSION (if any)
    # These functions will also attempt matching of variant name to common CNV and EXPRESSION names related to exons
    cnvNames = civic_return_all_cnvs(geneData)
    exprNames = civic_return_all_expr(geneData)
    # All variant names not matching a CNV or EXPRESSION will be returned
    matches = []
    for varName in list(geneData.keys()):
        # Skip variant names matching a CNV or EXPRESSION
        if (varName in cnvNames) or (varName in exprNames):
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
        # Also attempt matching of variant name to common CNV names related to exons
        isExon = cnv_is_exon_string(varName)
        if (varName in cnvNames) or isExon:
#             ## Special case for 'LOSS': only considered synonym of 'DELETION' when a certain variant type is present
#             if varName == 'LOSS' and ('TRANSCRIPT_ABLATION' not in geneData[varName]['types']):
#                 continue
            if varName not in matches:
                matches.append(varName)
    return matches


# Given all CIVIC variant records associated to a given gene, return all those variant names that correspond to EXPRESSION records
# For this, attemp to match the gene's variant names to the most common EXPRESSION names existing in CIVIC
# Additionally, consider other special EXPRESSION cases present in CIVIC (eg. EXON 1-2 OVEREXPRESSION, EXON 5 UNDEREXPRESSION, ...)
def civic_return_all_expr(geneData):
    # Common EXPRESSION record names in CIVIC
    exprNames = ['OVEREXPRESSION','UNDEREXPRESSION','EXPRESSION']
    # All matched variant names will be returned
    matches = []
    for varName in list(geneData.keys()):
        # Also attempt matching of variant name to common EXPRESSION names related to exons
        isExon = expr_is_exon_string(varName)
        if (varName in exprNames) or isExon:
            if varName not in matches:
                matches.append(varName)
    return matches


# TODO: expand framework to provide variant type (SNV or CNV) along with each individual variant, always generate matching strings for both types, and match to one depending on retrieved type -> would allow mixing SNV and CNV variants in the same file/input object
# Note that rowMap will have a different structure depending on dataType:
#  - For SNV: gene -> annotationType -> value (either SNV annotation, variant impact or variant exon)
#  - For CNV: gene -> value (CNV category)


# Given a CIVIC variant name and its corresponding hgvs expressions (if available),
# generate a list of potential strings that will be used to match the variant to our input
# For CNVs, variant matching is not based on HGVS, since input data is different
def civic_matchStrings(varName, hgvsExpressions, dataType):
    check_data_type(dataType)
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
    ## For CNVs and EXPRESSION data, this string is the only relevant one for matching to CIVIC variants
    matchStrings.append(varName)

    return matchStrings



# Given a set of one or more genomic alterations (SNV or CNV), generate a list of potential strings that will be used to match the variant in CIVIC eg. EXON 15 MUTATION, AMPLIFICATION
# Function only appplicable to data types "SNV","CNV" (strings generated are different in each case)
# Arguments impactAnnots and exonAnnots are optional for "SNV" argument. When provided, variant impacts will be used to generate additional potential record names. If exon annotations are provided, then impacts must also be provided, and there must be a 1-1 correspondance with this list (to determine if variant is intronic or exonic, from the impact)
def input_matchStrings(varAnnotations, dataType, impactAnnots=[], exonAnnots=[]):
    check_data_type(dataType)
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
        newTags = []
        ## Always include 'MUTATION' as a potential variant tag
        newTags.append('MUTATION')

        # Iterate impact information to generate and include additional variant tags
        for impact in impactAnnots:
            # Introduce sanity check on variant impact being available, as this field is optional
            if not impact:
                continue
            ## Generate potential additional variant tags based on the variant impact
            ## (can be contained in impact eg. 'frameshift_variant&stop_gained')
            if re.search('3_PRIME_UTR_VARIANT',impact):
                newTags.append("3' UTR MUTATION")

            if re.search('5_PRIME_UTR_VARIANT',impact):
                newTags.append("5' UTR MUTATION")

            if re.search('STOP_GAINED',impact):
                newTags.append('TRUNCATING MUTATION')

            if re.search('FRAMESHIFT_VARIANT',impact):
                newTags.append('FRAMESHIFT MUTATION')

        ## When provided, exon info must be matched with impact information, to indicate if variant is exonic or intronic
        if exonAnnots:
            ## Sanity check that both lists (impact and exon annotations) are the same length
            if len(impactAnnots) != len(exonAnnots):
                print("\nError! A variant impact must be provided for each of the provided exon annotations (to guess whether exonic or intronic variant)")
                sys.exit(1)

        ## Iterate exon information to generate and include additional variant tags
        for i,exon in enumerate(exonAnnots):
            # Introduce sanity check on variant exon being available, as this field is optional
            # Happens with impacts like eg. downstream_gene_variant, upstream_gene_variant, intergenic_region
            if not exon:
                continue
            ## Based on rank (exon/intron) information, generate an additional tag
            rank = exon.split('/')[0]
            # After checking some examples, it seems that intron information is associated to variant impacts 'intron_variant'
            # and 'sequence_feature' (can be contained in impact eg. 'splice_donor_variant&intron_variant')
            # NOTE: any variant impact other than the above two will generate a tag 'EXON XX MUTATION',
            # even 'synonymous' mutations
            if re.search('INTRON_VARIANT',impactAnnots[i]) or re.search('SEQUENCE_FEATURE',impactAnnots[i]):
                newTags.append('INTRON ' + rank + ' MUTATION')
            else:
                newTags.append('EXON ' + rank + ' MUTATION')
                ## As of 13/02/2019, only one case of variant name 'EXON XX FRAMESHIFT' (gene CALR)
                if re.search('FRAMESHIFT_VARIANT',impactAnnots[i]):
                    newTags.append('EXON ' + rank + ' FRAMESHIFT')

        ## Add unique tags and assign corresponding information about type of match
        for tag in newTags:
            if tag not in matchStrings:
                matchStrings.append(tag)
                isExact.append(True)
                isTrueExact.append(False)

    ## Opposite to the SNV case, for CNVs, varAnnotations should correspond to a single element (ie. CNV category)
    if dataType == 'CNV':
        newTags = []
        ## CIVIC seems to consider that GAIN and AMP are the same CNV
        if (varAnnotations == 'AMPLIFICATION') or (varAnnotations == 'AMP') or (varAnnotations == 'GAIN') or (varAnnotations == 'DUPLICATION') or (varAnnotations == 'DUP'):
            newTags.append('AMPLIFICATION')
        ## CIVIC seems to consider that DELETION and LOSS are the same CNV
        ## Both CNV categories 'DELETION' and 'LOSS' occur in CIVIC
        elif (varAnnotations == 'DELETION') or (varAnnotations == 'DEL') or (varAnnotations == 'LOSS'):
            newTags.append('DELETION')
            newTags.append('LOSS')
        newTags.append('COPY NUMBER VARIATION')

        for tag in newTags:
            if tag not in matchStrings:
                matchStrings.append(tag)
                isExact.append(True)
                isTrueExact.append(True)

    return (matchStrings,isExact,isTrueExact)


# 
# Possible tag combinations are either OVEREXPRESSION,EXPRESSION (when logFC>0) or UNDEREXPRESSION,EXPRESSION (when logFC<0)
def getExpressionStrings(gene,logfc):
    # TODO: sanity checks about non-empty number and sign
    logfc = check_logFC(logfc,gene)

    # Expression tags will be used to match expression records in CIVIC
    # expr_change = ""
    matchStrings = []
    # Generate expression tags to be able to match record in CIVIC
    # As of 29/09/2019, only relevant records in CIVIC are: 'OVEREXPRESSION', 'EXPRESSION', 'UNDEREXPRESSION' (ignore remaining)
    if logfc > 0:
        matchStrings.append('OVEREXPRESSION')
        # expr_change = 'OVEREXPRESSION'
    elif logfc < 0:
        matchStrings.append('UNDEREXPRESSION')
        # expr_change = 'UNDEREXPRESSION'
    else:
        raise ValueError("Invalid logFC = '%s' for gene '%s'. Only differentially expressed genes are valid." %(logfc,gene))

    # EXPRESSION (directionless) tag will also be taken into account for matching to CIVIC records
    matchStrings.append('EXPRESSION')

    # Special case for gene 'CDKN2A': synonym 'p16' is used in CIVIC records (eg. 'p16 EXPRESSION')
    # Include these cases as well to be able to match them
    new_tags = []
    if gene == 'CDKN2A':
        for this_tag in matchStrings:
            new_tag = 'P16 ' + this_tag
            if new_tag not in matchStrings:
                matchStrings.append(new_tag)

    return matchStrings


# FIXME
# Attempt match to a CIViC record using all available variant annotations for the given gene
# A tier and match will always be retrieved (either exact, positional or all gene variants when no match)
# FIXME
# Given a list of variant annotations for a gene in the input table (HGVS but also synonym descriptive terms eg. EXON 15 MUTATION and positional
# strings eg. p.Val600), and all CIVIC variant records for said gene, attempt mapping of any input annotations to one or more CIVIC records.
# Return tier of the match (exact, positional, no match) and corresponding CIVIC variant records that were matched

# Do exhaustive search as there could be >1 perfect match in cases of redundancy (eg. E55FS and E55RFSTER11
# both translate to p.Glu55fs) and there can be >1 positional matches (for SNV case, 'general' variants
# eg. V600 have preference over the rest). Also, for SNVs in case of having two types of exact matches (truly exact
# eg. V600E and descriptive term eg. EXON 15 MUTATION), give preference to the former as descriptive terms should only
# be used when no true exact match was found.

# For CNV, inputStrings will be limited (usually only 1 or 2), and all will always correspond to exact matches (no positional)
def match_variants_in_civic(gene, variants, varMap, dataType, impacts=[], exons=[]):
    match = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}

# TODO: sanity check that arguments are of the correct type
    check_is_str(gene,"gene")
    check_is_list(variants,"variants")
    check_is_dict(varMap,"varMap")
    check_data_type(dataType)

    # Function for generating additional synonym strings for input annotations (eg. EXON 15 MUTATION)
    # isExact is a list of equal length to inputStrings, indicating whether a given string corresponds
    # to an exact (True) or positional (False) match
    # allVariants must refer to the same variant, eg. "c." and "p." annotations of the same variant
    # FIXME: allImpacts and allExons can be empty
    (inputStrings,isExact,isTrueExact) = input_matchStrings(variants, dataType, impacts, exons)

    ## If gene is in CIVIC, possible tier levels are 1,2,3
    if gene in varMap.keys():
        for var_id in varMap[gene].keys():
            check_dict_entry(varMap[gene][var_id],"varMap","name","name")
            variant_name = varMap[gene][var_id]["name"]
            check_dict_entry(varMap[gene][var_id],"varMap","hgvs","hgvs")
            hgvs_expressions = varMap[gene][var_id]["hgvs"]

            # Generate list of strings that will be used to match input variants to CIVIC database
            # Returned list always has at least length=1 (in this case, containing only variant name)
            # For dataType=CNV, variant matching is not based on HGVS, so matchStrings will only contain the variant name
            civicStrings = civic_matchStrings(variant_name, hgvs_expressions, dataType)

## TODO: in the future, change this if block to situation of 1) having impact and exon 2) or not

            # Iterate input annotations strings and attempt match to any CIVIC variant string
            # The position of each input string corresponds to the type of match (true exact, synonym exact or positional)
            for indx,inputAnnot in enumerate(inputStrings):
                if inputAnnot in civicStrings:
                    # Determine type of match using its position and store accordingly
                    # For CNV, all matches will be exact and true exact matches (either variant is in CIVIC or not)
                    if isExact[indx]:
                        # For SNV: give preference to truly exact matches over descriptive synonym terms
                        # eg. report V600E over EXON 15 MUTATION
                        if isTrueExact[indx]:
                            # For SNV: True exact match (eg. V600E)
                            # For CNV: all exact matches will also be true exact matches (either CNV is in CIVIC or not)
                            # There could be 0,1,>1 exact matches (eg. DELETION + LOSS + COPY NUMBER for cnvs)
                            if var_id not in match["tier_1"]:
                                match["tier_1"].append(var_id)
                        else:
                            # For SNV: descriptive term match (eg. EXON 15 MUTATION, FRAMESHIFT MUTATION)
                            # There could be 0,1,>1 synonym matches
                            if var_id not in match["tier_1b"]:
                                match["tier_1b"].append(var_id)
                    else:
                        # Positional match
                        # There could be 0,1,>1 positional matches
                        if var_id not in match["tier_2"]:
                            match["tier_2"].append(var_id)
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
                                # There could be 0,1,>1 positional matches
                                if var_id not in match["tier_2"]:
                                    match["tier_2"].append(var_id)

        ## Once all CIVIC variants have been iterated, determine final tier and corresponding matched variants
        ## For CNV: either exact or positional matches will occurr due to the implementation design

        # For CNV: positional matches will correspond to EXON records (eg. EXON 1-2 DELETION, 3' EXON DELETION..)
        # For SNV: if there are positional matches, check for preferential positional matches, ie. general variants (like V600)
        if match["tier_2"] and (dataType == 'SNV'):
            for var in match["tier_2"]:
                isGeneral = check_general_variant(var)
                if isGeneral:
                    # Stop as soon as a general variant is found and report only this
                    match["tier_2"] = [var]
                    break

        # If no match was found, then tier case is 3, and return all relevant variants
        if not (match["tier_1"] or match["tier_1b"] or match["tier_2"]):
            ## For SNV: when input variant could not be matched in CIVIC (tier3), return all CIVIC variants associated
            ## to the given gene but that do not correspond to a CNV or EXPRESSION related variant
            if dataType == 'SNV':
                match["tier_3"] = civic_return_all_snvs(varMap[gene])
            ## For CNV: when input cnv could not be matched in CIVIC (tier3), return all CIVIC cnvs associated the given gene (if any)
            ## ie. not all CIVIC records are returned but only those that are matched to a 'CNV' tags
            elif dataType == 'CNV':
                match["tier_3"] = civic_return_all_cnvs(varMap[gene])

    ## If gene is not in CIVIC, tier level is 4 and match will be empty
    else:
        # Sanity check that no other match should have been found
        if (match["tier_1"] or match["tier_1b"] or match["tier_2"] or match["tier_3"]):
            raise ValueError("")
#                 f"Unexpected variant match condition was met.\n"
        # No variants will be reported in this case (as the information is not available)
        match["tier_4"] = True

    return match


# For EXPRESSION data
# TODO: either tier1 match or tier4 (expression record is available in CIVIC for the given gene or not)
def match_expression_in_civic(gene, expression_strings, varMap):
    match = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}

# TODO: sanity check that arguments are of the correct type
    check_is_dict(varMap)

    # If gene is in CIVIC, only possible tier level is 1
    if gene in varMap.keys():
        for var_id in varMap[gene].keys():
            # Only variant name is relevant in this case
            # Disregard HGVS strings for expression variants
            check_dict_entry(varMap[gene][var_id],"varMap","name","name")
            variant_name = varMap[gene][var_id]["name"]

            # For dataType=EXPR (in truth, dataType!=SNV), variant matching is not based on HGVS, so matchStrings will only contain the variant name
            civicStrings = civic_matchStrings(variant_name, hgvs_expressions, dataType="EXPR")

            # Iterate available expression tags for the given gene and attempt match to a CIVIC record for that gene
            # Opposite to the SNV and CNV case, here there is no tier hierarchy, ie. gene expression is either matched in CIVIC or not
            for expr_tag in expression_strings:
                if expr_tag in civicStrings:
                    # There could be 0,1,>1 exact matches (eg. OVEREXPRESSION, EXPRESSION)
                    if var_id not in match["tier_1"]:
                        match["tier_1"].append(var_id)

            # Special case for EXPRESSION records related to EXONS: manually add them as matches when necessary
            # Specific expression change must be identical in order to consider them as matched: eg. 'OVEREXPRESSION' and 'EXON 18 OVEREXPRESSION'
            for this_string in civicStrings:
                (isExon,exprType) = expr_is_exon_string(this_string)
                if isExon and exprType:
                    # FIXME: Records of the type 'P16 EXPRESSION' would be missed
                    if exprType in expression_strings:
                        if var_id not in match["tier_1"]:
                            match["tier_1"].append(var_id)


        # If no exact match was found, then tier case is 3, and return all relevant EXPRESSION records for the current gene (if any are found in CIVIC)
        # If there are no matching records for tier3, then this category will also be empty in the returned dictionary
        if not (match["tier_1"] or match["tier_1b"] or match["tier_2"]):
            # Not all CIVIC records are returned but only those that are matched to a 'EXPRESSION' tag
            match["tier_3"] = civic_return_all_expr(varMap[gene])

    # If gene is not in CIVIC, tier level is 4 and match will be empty
    else:
        # Sanity check that no other match should have been found
        if (match["tier_1"] or match["tier_1b"] or match["tier_2"] or match["tier_3"]):
            raise ValueError("")
#                 f"Unexpected expression match condition was met.\n"
        # No variants will be reported in this case (as the information is not available)
        match["tier_4"] = True

    return match


def add_match(matchMap,gene,variant,match):
    sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]
    # Sanity check that arguments are of the correct type
    check_is_dict(matchMap,"matchMap")
    check_is_dict(match,"match")
    check_is_str(gene,"gene")
    check_is_str(variant,"variant")
    # Sanity check that dict contains all expected keys (throws an error if not)
    check_keys(match,"match",sorted_tiers,matches_all=True)
    # Sanity check that gene and variant are keys of the dict
    check_dict_entry(matchMap,"matchMap",gene,"gene")
    check_dict_entry(matchMap[gene],"matchMap",variant,"variant")

    # Keep track of all the variant ids matched across tiers for the current gene + variant
    matchedIds = []
    # Add match results to the current gene + variant in matchMap
    for x in match.keys():
        if x == "tier_4":
            matchMap[gene][variant][x] = match[x]
        else:
            matchMap[gene][variant][x] = []
            for varId in match[x]:
                matchMap[gene][variant][x].append(varId)
                # It is possible that the same CIVIC variant id has been matched for several tiers
                if varId not in matchedIds:
                    matchedIds.append(varId)

    return (matchMap,matchedIds)


# FIXME
def match_in_civic(varData, dataType, identifier_type, select_tier="all", varMap=None):
    sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]

    # Sanity check that arguments are of the correct type
    check_is_dict(varData,"varData")
    check_is_str(dataType,"dataType")
    check_identifier_type(identifier_type)
    check_data_type(dataType)

    # Process and sanity check provided select_tier (expected format, valid values, etc.)
    select_tier = check_tier_selection(select_tier,sorted_tiers)

    matchMap = {}
    matchedIds = []     # keep track of all matched variant ids across genes and tiers

# TODO: add option to provide user-specified varMap (eg. when filters need to be applied) -> add NOTE or warning about anything not being provided will be interpreted as not available in CIVIC

    if (varMap is None) or (not varMap):
        # gene -> variant -> null
        all_genes = list(varData.keys())
        varMap = query_civic(all_genes, identifier_type)
    # TODO: Check correct structure of varMap?
    else:
        check_is_dict(varMap,"varMap")

    # gene -> variant -> null
    # where variant -> var="dna|prot|impact|exon|lineNumber"
    for gene in varData.keys():
        gene = check_empty_input(gene, "Gene", isRequired=True)
        if gene not in matchMap.keys():
            matchMap[gene] = {}
        for variant in varData[gene].keys():
            variant = check_empty_input(variant, "Variant", isRequired=True)
            # Overwrite duplicated variant ids (should never happen)
            matchMap[gene][variant] = {}
            variants = []
            impactArr = []
            exonArr = []

            varArr = variant.split("|")
            if dataType == "SNV":
                # Sanity check for expected SNV format (at least 4 fields should exist, even if empty)
                if (len(varArr) < 4):
# TODO
                    raise ValueError("Must provide at least 4 fields to describe a SNV variant (even if some can be empty): 'dna|[prot]|[impact]|[exon]|..|'")
                # Format: var="dna|prot|[impact]|[exon]|..|"
                # NOTE: all fields can contain >1 terms separated with ',' (no spaces). Fields 'dna' and 'prot' are required and 'impacts' and 'exons' are optional
                # NOTE: there might be additional fields after, eg. 'lineNumber' when data has been parsed from a file
                cVars = varArr[0]
                pVars = varArr[1]
                impacts = varArr[2]
                exons = varArr[3]

                # Sanity check for required and optional fields

                # Field 'Variant_dna' must exist but cannot contain empty values
                cVarArr = parse_input(cVars, "Variant_dna", isRequired=True)
                # Field 'Variant_prot' must exist but can contain empty values
                pVarArr = parse_input(pVars, "Variant_prot", isRequired=False)
# TODO: use a function for this parsing
                for cVar in cVarArr:
                    # Sanity check that variant starts with "c."
                    check_is_cHGVS(cVar)
                    if cVar not in variants:
                        variants.append(cVar)
                for pVar in pVarArr:
                    # Sanity check that variant starts with "p."
                    # Sanity check that variant starts with "p."
                    if not pVar:
                        continue
                    check_is_pHGVS(pVar)
                    if pVar not in variants:
                        variants.append(pVar)
                # Field 'Variant_impact' is optional and can contain empty values
                impactArr = parse_input(impacts, "Variant_impact", isRequired=False)
                # Field 'Variant_exon' is optional and can contain empty values
                exonArr = parse_input(exons, "Variant_exon", isRequired=False)

            if dataType == "CNV":
                # Format: var="cnv|.."
                # NOTE: CNV field can contain >1 terms separated with ',' (no spaces)
                # NOTE: there might be additional fields after, eg. 'lineNumber' when data has been parsed from a file
                cnvVars = varArr[0]
                variants = parse_input(cnvVars, "Variant_cnv", isRequired=True)

            if dataType == "EXPR":
                # Format: expr="logFC|.."
                # NOTE: logFC field should be a single number (sign will determine if gene is overexpressed or underexpressed)
                # NOTE: there might be additional fields after, eg. 'lineNumber' when data has been parsed from a file
                logFC = varArr[0]
                variants = getExpressionStrings(gene,logFC)

# TODO: check that data type is valid and one of the 3 options
            # tier -> [matched_vars]
            if (dataType == "SNV") or (dataType == "CNV"):
                print("Matching variant of type '%s'..." %(dataType))
                print("Gene: %s" %(gene))
                print("Variants: %s" %(variants))
                print("Impacts: %s" %(impactArr))
                print("Exons: %s" %(exonArr))
                match = match_variants_in_civic(gene, variants, varMap, dataType, impacts=impactArr, exons=exonArr)

                # Avoid executing unnecessary filter when select_tier='all' (as resulting dict will be identical)
                if select_tier != "all":
                    match = filter_match(match, select_tier)

            # In the case of EXPRESSION data, only tier1 is possible (either CIVIC record was matched or not)
            if (dataType == "EXPR"):
                match = match_expression_in_civic(gene, variants, varMap)

            # Add the match to the current entry for gene + variant
            # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
            (matchMap,allIds) = add_match(matchMap,gene,variant,match)

            # Keep track of all the variant ids matched across genes and tiers
            for varId in allIds:
                if varId not in matchedIds:
                    matchedIds.append(varId)

    # At this point, all input gene + variants in varData have been matched to CIVIC
    # Filter varMap used to match CIVIC info based on the matched variant ids (to avoid returning unnecessary records)
    varMap = filter_civic(varMap, var_id_in=matchedIds, output_empty=False)

    # Return matchMap, list of all matched variant ids, and associated variant records retrieved from CIVIC (or provided by user), already filtered for the matched variants
    return (matchMap,matchedIds,varMap)


# TODO: make match of tier case-insensitive?
def filter_match(match, select_tier):
    sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]
    # Sanity check that dict contains all expected keys (throws an error if not)
    check_keys(match,"match",sorted_tiers,matches_all=True)
    # Process and sanity check provided select_tier (expected format, valid values, etc.)
    select_tier = check_tier_selection(select_tier,sorted_tiers)

    new = {"tier_1":[], "tier_1b":[], "tier_2":[], "tier_3":[], "tier_4":False}
    keepTiers = []

    # When select_tier="all", filter is off (keep data for all tiers)
    if isinstance(select_tier, str) and (select_tier == "all"):
        for tmpTier in sorted_tiers:
            keepTiers.append(tmpTier)

    # When select_tier="highest", then keep data only for the highest tier 1>1b>2>3>4
    elif isinstance(select_tier, str) and (select_tier == "highest"):
        for tmpTier in sorted_tiers:
            if tmpTier != "tier_4":
                check_is_list(match[tmpTier])
                if match[tmpTier]:
                    keepTiers.append(tmpTier)
                    break
            else:
                # If we reached to this point, that means all other tiers were empty (tier_4=True)
                keepTiers.append(tmpTier)
                break

    # When select_tier is a list of tiers, then keep data only for those tiers
    elif isinstance(select_tier, list):
        for tmpTier in sorted_tiers:
            if tmpTier in select_tier:
                keepTiers.append(tmpTier)

    # Sanity check for any duplicated tiers
    keepTiers = list(set(keepTiers))

    # Keep variant matches only for the selected tiers (all other data will be excluded from output)
    # Respect the order of tiers while filling the new dictionary
    for tier in sorted_tiers:
        if tier in keepTiers:
            if tier != "tier_4":
                for varId in match[tier]:
                    new[tier].append(varId)

    # Always check if current match corresponds to tier_4 after filtering
    # eg. if only tier_1 is selected and current match did not have any
    if not (new["tier_1"] or new["tier_1b"] or new["tier_2"] or new["tier_3"]):
        new["tier_4"] = True

    return new


def filter_matches(matchMap, select_tier):
    cleanMap = {}
    matchedIds = []

    # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
    # where variant -> var="dna|prot|impact|exon|lineNumber"
    for gene in matchMap.keys():
        if gene not in cleanMap.keys():
            cleanMap[gene] = {}
        for variant in matchMap[gene].keys():
            cleanMap[gene][variant] = {}
            new = filter_match(matchMap[gene][variant], select_tier)
            # Add the filtered match to the current entry for gene + variant
            (cleanMap,allIds) = add_match(cleanMap,gene,variant,new)
            # Keep track of all the variant ids matched across genes and tiers
            for varId in allIds:
                if varId not in matchedIds:
                    matchedIds.append(varId)

    return (cleanMap,matchedIds)


# TODO
# For 'predictive' evidence (writeDrug=True), keep dictionary of drug support for the current variant match (one support per tier)
# Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
def process_drug_support(matchMap, varMap, supportDict):
    evidenceType = 'PREDICTIVE'
    newMap = {}

    # gene -> variant -> {tier1,tier1b..} -> [matched_vars]
    # where variant -> var="dna|prot|impact|exon|lineNumber"
    for gene in matchMap.keys():
        if gene not in varMap.keys():
            raise ValueError("")
            # TODO
        if gene not in newMap.keys():
            newMap[gene] = {}
        for variant in matchMap[gene].keys():
            newMap[gene][variant] = {}
            for tier in matchMap[gene][variant].keys():
                newMap[gene][variant][tier] = {}
                # NOTE: tier4 has specific format compared to the others (boolean vs. list)
                if tier == "tier_4":
                    newMap[gene][variant][tier]['matched'] = False
                else:
                    newMap[gene][variant][tier]['matched'] = []

# TODO: drug support for tier 1, 1b and 2 is clear
# TODO: drug support for tier 3 will not make sense (averaging across variants) -> leave it in case the user is interested
# TODO: drug support for tier 4 should directly be empty
                newMap[gene][variant][tier]['drug_support'] = []
                drugMap = {}
                if tier != "tier_4":
                    check_is_list(matchMap[gene][variant][tier])
                    for varId in matchMap[gene][variant][tier]:
                        if varId not in varMap[gene].keys():
                            raise ValueError("")
                            # TODO
                        newMap[gene][variant][tier]['matched'].append(varId)
# TODO: ensure that filtering returnEmpty=False does affect other functions
                        if evidenceType in varMap[gene][varId]['evidence_items'].keys():
                            for ct in varMap[gene][varId]['evidence_items'][evidenceType].keys():
                                for disease in varMap[gene][varId]['evidence_items'][evidenceType][ct].keys():
                                    for drug in varMap[gene][varId]['evidence_items'][evidenceType][ct][disease].keys():
                                        if drug not in drugMap.keys():
                                            drugMap[drug] = {}
                                        if ct not in drugMap[drug].keys():
                                            drugMap[drug][ct] = []
                                        for evidence in varMap[gene][varId]['evidence_items'][evidenceType][ct][disease][drug].keys():
                                            # Split the evidence direction and clinical significance
                                            if ":" not in evidence:
                                                raise ValueError("")
                                                # TODO
                                            (direction, clin_signf) = evidence.strip().split(':')

                                            # For each evidence (ie combination of direction+clin_signf), count how many different evidence items support it
                                            # At this stage, we find count evidence items by counting how many different combinations of level+pmids there are for the same drug, disease and evidence
                                            if ('NULL' in direction) or ('N/A' in direction) or ('NULL' in clin_signf) or ('N/A' in clin_signf):
# FIXME: have this as part of input config
                                                thisSupport = 'UNKNOWN_BLANK'
                                            else:
                                                if direction not in supportDict.keys():
                                                    # TODO
                                                    print("\nError! Could not find direction %s in support dictionary." %(direction))
                                                    sys.exit(1)
                                                if clin_signf not in supportDict[direction].keys():
                                                    # TODO
                                                    print("\nError! Could not find clinical significance %s in support dictionary." %(clin_signf))
                                                    sys.exit(1)
                                                thisSupport = supportDict[direction][clin_signf]

                                            # Keep track of number of occurrences for each support type for the given drug
                                            # Here, take into account the number of supporting PMIDs associated to each evidence item
                                            for evidence_level in varMap[gene][varId]['evidence_items'][evidenceType][ct][disease][drug][evidence].keys():
                                                for thisString in varMap[gene][varId]['evidence_items'][evidenceType][ct][disease][drug][evidence][evidence_level].keys():
                                                    drugMap[drug][ct].append(thisSupport)

                # Process drug support information parsed for the current tier match
                # Drug support for tier 4 will never be available
                for thisDrug in drugMap.keys():
                    for thisCt in drugMap[thisDrug].keys():
                        # Given the selected ct, count number of occurrences for each possible support type (if any)
# FIXME: have this as part of input config
                        count_pos = drugMap[thisDrug][thisCt].count('POSITIVE')
                        count_neg = drugMap[thisDrug][thisCt].count('NEGATIVE')
                        count_unk = drugMap[thisDrug][thisCt].count('UNKNOWN_BLANK')
                        count_dns = drugMap[thisDrug][thisCt].count('UNKNOWN_DNS')

                        # Pool UNKNOWN_BLANK and UNKNOWN_DNS together (as both result in unknown CIVIC support)
                        count_total_unk = count_unk + count_dns
                        # Sanity check that there is at least some support
                        if (count_pos == 0) and (count_neg == 0) and (count_total_unk == 0):
                            print("Error! Unexpected support case for gene %s." %(gene))
                            sys.exit(1)

                        # Resolve contradicting evidence (if any) by majority vote
                        tempSupport = ''
                        # For this, pool UNKNOWN_BLANK and UNKNOWN_DNS together
                        # Whenever there is a tie of "confident" (pos or neg) vs "non-confident" (unk), choose the confident one
                        if (count_total_unk > count_pos) and (count_total_unk > count_neg):
                            tempSupport = "CIVIC_UNKNOWN"
                        elif count_pos == count_neg:
                            tempSupport = "CIVIC_CONFLICT"
                        elif (count_pos > count_neg) and (count_pos >= count_total_unk):
                            tempSupport = "CIVIC_SUPPORT"
                        elif (count_neg > count_pos) and (count_neg >= count_total_unk):
                            tempSupport = "CIVIC_RESISTANCE"
                        else:
                            print("Error! Unexpected support case for gene %s." %(gene))
                            sys.exit(1)

                        # Build support string for each given combination of drug, ct and matched tier
                        # Format: DRUG:CT:SUPPORT
                        drugSupport = thisDrug + ':' + thisCt.upper() + ':' + tempSupport
                        newMap[gene][variant][tier]['drug_support'].append(drugSupport)

# TODO: Loop for all available tiers that are not tier 4 (could be indicated in the config) -> to make it robust to changes in the tier categories
            # Always check if current match corresponds to a tier_4 situation (all other tiers will be empty)
            if not (newMap[gene][variant]["tier_1"]["matched"] or newMap[gene][variant]["tier_1b"]["matched"] or newMap[gene][variant]["tier_2"]["matched"] or newMap[gene][variant]["tier_3"]["matched"]):
                newMap[gene][variant]["tier_4"]["matched"] = True

    return newMap


# Given a list of disease names, first filter out those matching the black-listed terms (if any), and from the remaining set, return either:
#   - Subset of diseases matching white-listed terms - partial match (eg. 'Melanoma')
#   - If previous is not available, subset of matching high-level diseases - exact match (eg. 'Cancer', 'Solid tumor')
#   - If previous is not available, return all available diseases (ie. so, original set except those matched to the black list, if any)
def classify_diseases(diseaseArr, disease_name_not_in, disease_name_in, alt_disease_names):

    ## Keep track of diseases that did not pass the black-list filter
    blackMatched = []
    ## Keep track of diseases that passed the black-list filter
    cleanSet = []
    ## Final list of matched diseases (associated to one of 3 categories depending on match type: 'ct', 'gt' or 'nct')
    matched = []
    ## Keep track of which type of disease specificy was matched in the end ('ct', 'gt' or 'nct')
    ctArr = []
    gtArr = []
    nctArr = []

    ## 1) First, remove diseases that partially match black-listed terms (if any are provided)
    # NOTE: PARTIAL matches to the black list are allowed! eg:
    #   - including 'small' will remove 'non-small cell lung cancer' and 'lung small cell carcinoma'
    #   - including 'non-small' will remove 'non-small cell lung cancer' but not 'lung small cell carcinoma'
    if disease_name_not_in:
        # Iterate available diseases and keep track of those that partially match to black list (there can be several)
        for disease in diseaseArr:
            # To find partial matches, it is necessary to iterate through both lists (input and civic)
            for blackListed in disease_name_not_in:
                # Search for partial match of INPUT disease in CIVIC disease e.g. 'Melanoma' (input list) in 'Skin Melanoma' (CIVIC) and not opposite
                if blackListed in disease:
                    if disease not in blackMatched:
                        blackMatched.append(disease)
        # Iterate available diseases once again to retrieve those that passed the black list filter
        for disease in diseaseArr:
            # Retrieve valid diseases only
            if disease in blackMatched:
                continue
            if disease not in cleanSet:
                cleanSet.append(disease)

    ## If no black list was provided, then all available diseases constitute the clean set
    else:
        cleanSet = diseaseArr

    ## 2) Now, iterate the list of "allowed" diseases (ie. passing the black list filter) and attempt to match to white-listed terms
    # NOTE: again, PARTIAL matches to the white list are allowed! eg:
    #   - including 'melanoma' will match 'melanoma', 'skin melanoma' and 'uveal melanoma', but not 'skin cancer' (hypothetical)
    #   - including 'uveal melanoma' will only match 'uveal melanoma'
    for cleanType in cleanSet:
        # To find partial matches, it is necessary to iterate through both lists (input and civic)
        for whiteListed in disease_name_in:
            # Search for partial match of INPUT disease in CIVIC disease e.g. 'Melanoma' (input list) in 'Skin Melanoma' (CIVIC) and not opposite
            if whiteListed in cleanType:
                # Keep track of diseases that passed the white list filter
                if cleanType not in matched:
                    matched.append(cleanType)
                    if cleanType not in ctArr:
                        ctArr.append(cleanType)

    ## 3) For those not already matched as ct, attempt a second-best match strategy to higher-level diseases
    ## In CIVIC, some 'general' diseases are included as disease, eg. 'cancer' or 'solid tumor'. Hence, must be exact matches since they are DB-specific
    # NOTE: here, only PERFECT matches are allowed!
    #   - including 'cancer' will only match 'cancer' and not 'lung cancer'
    for cleanType in cleanSet:
        if cleanType in alt_disease_names:
            if cleanType not in matched:
                matched.append(cleanType)
                if cleanType not in gtArr:
                    gtArr.append(cleanType)

    ## 4) For those not already matched as either ct or gt, return all 'allowed' (ie. not black-listed) diseases available in CIVIC
    ## These will be considered as non-specific diseases (ie. off label)
    for cleanType in cleanSet:
        if cleanType not in matched:
            matched.append(cleanType)
            if cleanType not in nctArr:
                nctArr.append(cleanType)

    ## Return diseases that did not match any blacklisted terms, classified and split into: ct, gt, nct
    return (ctArr,gtArr,nctArr)


def add_ct(diseases,ct,gene,variant,evidence_type,newMap,varMap):
    check_is_dict(newMap)
    check_is_dict(varMap)

    # Sanity check that gene and variant are keys of the dict
    check_dict_entry(newMap,"newMap",gene,"gene")
    check_dict_entry(newMap[gene],"newMap",variant,"variant")
    check_dict_entry(newMap[gene][variant],"newMap","evidence_items","key")
    check_dict_entry(newMap[gene][variant]["evidence_items"],"newMap",evidence_type,"evidence type")

    check_dict_entry(varMap,"varMap",gene,"gene")
    check_dict_entry(varMap[gene],"varMap",variant,"variant")
    check_dict_entry(varMap[gene][variant],"varMap","evidence_items","key")
    check_dict_entry(varMap[gene][variant]["evidence_items"],"varMap",evidence_type,"evidence type")

    if ct not in newMap[gene][variant]["evidence_items"][evidence_type].keys():
        newMap[gene][variant]["evidence_items"][evidence_type][ct] = {}
    else:
# FIXME
        raise ValueError("TODO")
    for disease in diseases:
        if disease not in varMap[gene][variant]["evidence_items"][evidence_type].keys():
            # TODO
            sys.exit(1)
        newMap[gene][variant]["evidence_items"][evidence_type][ct][disease] = {}
        for drug in varMap[gene][variant]["evidence_items"][evidence_type][disease].keys():
            newMap[gene][variant]["evidence_items"][evidence_type][ct][disease][drug] = {}
            for evidence in varMap[gene][variant]["evidence_items"][evidence_type][disease][drug].keys():
                newMap[gene][variant]["evidence_items"][evidence_type][ct][disease][drug][evidence] = {}
                for evidence_level in varMap[gene][variant]["evidence_items"][evidence_type][disease][drug][evidence].keys():
                    newMap[gene][variant]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level] = []
                    for thisString in varMap[gene][variant]["evidence_items"][evidence_type][disease][drug][evidence][evidence_level]:
                        newMap[gene][variant]["evidence_items"][evidence_type][ct][disease][drug][evidence][evidence_level].append(thisString)

    return newMap


def annotate_ct(varMap, disease_name_not_in, disease_name_in, alt_disease_names):
    sorted_cts = ["ct","gt","nct"]
    varmap_entries = ['name','civic_score','hgvs','types','n_evidence_items','evidence_items']
    newMap = {}

    # Iterate the complete varMap dict and reorganize it to classify diseases
    for gene in varMap.keys():
        if gene not in newMap.keys():
            newMap[gene] = {}
        for variant in varMap[gene].keys():
            # Overwrite duplicated variant ids (should never happen)
            newMap[gene][variant] = {}
            # Sanity check that some expected fields can be found in the dictionary
            check_keys(varMap[gene][variant],"varMap",varmap_entries,matches_all=False)
            newMap[gene][variant]['name'] = varMap[gene][variant]['name']
            newMap[gene][variant]['civic_score'] = varMap[gene][variant]['civic_score']
            newMap[gene][variant]['hgvs'] = [a for a in varMap[gene][variant]['hgvs']]
            newMap[gene][variant]['types'] = [b for b in varMap[gene][variant]['types']]
            newMap[gene][variant]['n_evidence_items'] = varMap[gene][variant]['n_evidence_items']
            newMap[gene][variant]['evidence_items'] = {}
            for evidence_type in varMap[gene][variant]['evidence_items'].keys():
                newMap[gene][variant]['evidence_items'][evidence_type] = {}
                # Retrieve all disease names associated witht he current evidence type, and classify them into ct, gt, nct
# TODO what happens when list of diseases is empty?
                allDiseases = list(varMap[gene][variant]['evidence_items'][evidence_type].keys())
                (ctDis,gtDis,nctDis) = classify_diseases(allDiseases, disease_name_not_in, disease_name_in, alt_disease_names)
                # Add new layer of classification within the evidence types: specificity of the disease (ct, gt, nct)
                for ct in sorted_cts:
# TODO what happens when list of diseases to add is empty?
                    if ct == "ct":
                        newMap = add_ct(ctDis,ct,gene,variant,evidence_type,newMap,varMap)
                    if ct == "gt":
                        newMap = add_ct(gtDis,ct,gene,variant,evidence_type,newMap,varMap)
                    if ct == "nct":
                        newMap = add_ct(nctDis,ct,gene,variant,evidence_type,newMap,varMap)

    return newMap


def filter_ct(varMap, select_ct):
    sorted_cts = ["ct","gt","nct"]
    # Process and sanity check provided select_ct (expected format, valid values, etc.)
    select_ct = check_tier_selection(select_ct,sorted_cts)

    newMap = {}

    # When select_ct="all", filter is off (keep data for all 3 cts)
    # Do not copy dict again, simply return input dict untouched
    if isinstance(select_ct, str) and (select_ct == "all"):
        return varMap

    # Otherwise, some filtering needs to be done, so iterate varMap
    # Iterate the complete varMap dict and reorganize it to classify diseases
    for gene in varMap.keys():
        if gene not in newMap.keys():
            newMap[gene] = {}
        for variant in varMap[gene].keys():
            # Overwrite duplicated variant ids (should never happen)
            newMap[gene][variant] = {}
# FIXME: This would currently fail: adapt to have strict option
            check_keys(varMap[gene][variant], ['name','civic_score','hgvs','types','n_evidence_items','evidence_items'])
            newMap[gene][variant]['name'] = varMap[gene][variant]['name']
            newMap[gene][variant]['civic_score'] = varMap[gene][variant]['civic_score']
            newMap[gene][variant]['hgvs'] = [a for a in varMap[gene][variant]['hgvs']]
            newMap[gene][variant]['types'] = [b for b in varMap[gene][variant]['types']]
            newMap[gene][variant]['n_evidence_items'] = varMap[gene][variant]['n_evidence_items']
            newMap[gene][variant]['evidence_items'] = {}
            for evidence_type in varMap[gene][variant]['evidence_items'].keys():
                newMap[gene][variant]['evidence_items'][evidence_type] = {}
# TODO: check that input dict has the expected format
                check_keys(varMap[gene][variant]['evidence_items'][evidence_type], ['ct','gt','nct'])
                for ct in sorted_cts:
                    newMap[gene][variant]['evidence_items'][evidence_type][ct] = {}
                    ctArr = list(varMap[gene][variant]['evidence_items'][evidence_type][ct].keys())
                    # When select_ct="highest", then keep data only for the highest specificity available ct>gt>nct
                    if isinstance(select_ct, str) and (select_ct == "highest"):
                        # Check if data is available for the currently iterated ct
                        if ctArr:
                            newMap = add_ct(ctArr,ct,gene,variant,evidence_type,newMap,varMap)
                            # Prematurely exit the loop after this successful iteration (to keep only the highest ct)
                            break

                    # When select_ct is a list of cts, then keep data only for those selected
                    elif isinstance(select_ct, list):
                        # Check if the currently iterated ct should be kept or skipped
                        # Iterate all cts and add all those that are select (no premature exit of loop)
                        if (ct in select_ct):
                            newMap = add_ct(ctArr,ct,gene,variant,evidence_type,newMap,varMap)

    return newMap

