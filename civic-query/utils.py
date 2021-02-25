import sys
import os
import re

# import read_and_write
# dict_coodes = get_dict_aminoacids()

def check_is_list(myList):
    # Check object is a list (even if empty)
    if not isinstance(myList, list):
        raise TypeError(
            f"'{myList}' is not of type 'list'.\n"
        )


def check_string_filter_arguments(field, myList):
    # Check expected class was provided for each argument
    # myList should be a list (even if empty)
    check_is_list(myList)
    # field should be a non-empty string
    if not isinstance(field, str):
        raise TypeError(
            f"'{field}' is not of type 'str'.\n"
        )
    if not field:
        raise TypeError(
            f"'{field}' cannot be an empty string.\n"
        )
    # Use uppercase and remove leading/trailing spaces, for consistency of strings
    newField = field.strip().upper()
# TODO: convert every element into string (in case numbers were provided)
    newList = [x.strip().upper for x in myList]
    return (newField,newList)


def apply_in_filter(field, inList, matchType="exact"):
    (field, inList) = check_string_filter_arguments(field, inList)
    match = False
    # Only perform check when inList is not empty
    if inList:
        # If a list was provided and field is empty, then filter is not passed
        if field != "NULL":
            if (matchType == "exact"):
                if field in inList:
                    match = True
            if (matchType == "partial"):
                for tmp in inList:
                    if field in tmp:
                        match = True
    else:
        match = True
    return match


def apply_not_in_filter(field, outList, matchType="exact"):
    (field, outList) = check_string_filter_arguments(field, outList)
    match = False
    # If no list was provided, even if field is empty, filter is always "passed" (ie. False to not remove variant)
    # This allows to filter record if a particular field is "NULL"
    if outList:
        if (matchType == "exact"):
            if field in outList:
                match = True
        if (matchType == "partial"):
            for tmp in outList:
                if field in tmp:
                    match = True
    return match


# Only keep if cutoff or more instances
# Ignore filter when cutoff=0
def apply_cutoff_filter(field, cutoff):

## TODO check for correct classes provided
# field should be numeric or float
# cutoff should be numeric or float

    match = True
    cutoff_f = float(cutoff_f)
    field_f = float(field)
    if (cutoff_f != float(0)):
        if (field_f < cutoff_f):
            match = False
    return match



## TODO check types of sources are even a thing?

def filter_civic_results(varMap, gene_id_in=[], gene_id_not_in=[], variants=0, var_id_in=[], var_id_not_in=[], var_name_in=[], var_name_not_in=[], civic_score=0, var_origin_in=[], var_origin_not_in=[], var_type_in=[], var_type_not_in=[], evidence_items=0, evidence_status_in=[], evidence_status_not_in=[], source_in=[], source_not_in=[], disease_in=[], disease_not_in=[], evidence_rating=0, evidence_type_in=[], evidence_type_not_in=[], evidence_dir_in=[], evidence_dir_not_in=[], evidence_clinsig_in=[], evidence_clinsig_not_in=[], evidence_level_in=[], evidence_level_not_in=[], drugs=0, drug_interaction=[], drug_name_in=[], drug_name_not_in=[]):

## filters are applied in the same order as their corresponding arguments
## so, logic is always AND for all selected filters?

## if desired filter logic is not possible, then the function would need to be run several times, applyin the filters subsequently

    cleanMap = {}


    ## Process selected filters

## Assume varMap has the following structure:
# TODO
# gene -> variant_name -> id,civic_score,match_hgvs,hgvs,types,origin,

    # gene id should always be available
    for gene_id in varMap.keys():
        gene_id_str = str(gene_id)
        keepGene = apply_in_filter(gene_id_str, gene_id_in, matchType="exact")
        if not keepGene:
            continue
        removeGene = apply_not_in_filter(gene_id_str, gene_id_not_in, matchType="exact")
        if removeGene:
            continue

        # allow number of variants to be 0
        n_variants = len(varMap[gene_id].keys())
        keepGene = apply_cutoff_filter(n_variants, variants)
        if not keepGene:
            continue

        # current gene has passed all gene-level filters
        if gene_id_str not in cleanMap.keys():
            cleanMap[gene_id_str] = {}

        # variant id should always be available
        for var_id in varMap[gene_id].keys():
            var_id_str = str(var_id)
            keepVar = apply_in_filter(var_id_str, var_id_in, matchType="exact")
            if not keepVar:
                continue
            removeVar = apply_not_in_filter(var_id_str, var_id_not_in, matchType="exact")
            if removeVar:
                continue

            # variant name should always be available (never "NULL")
            variant = varMap[gene_id][var_id]["name"]
            keepVar = apply_in_filter(variant, var_name_in, matchType="partial")
            if not keepVar:
                continue
            removeVar = apply_not_in_filter(variant, var_name_not_in, matchType="partial")
            if removeVar:
                continue

            # civic score is always a number (can be 0)
            var_score = varMap[gene_id][var_id]["civic_score"]
            keepVar = apply_cutoff_filter(var_score, civic_score)
            if not keepVar:
                continue

            # variant types will never be empty list
            # use ["NULL"] when not available
            variant_types = varMap[gene_id][var_id]["types"]
            nKeep = 0
            nRemove = 0
            for var_type in variant_types:
                keepType = apply_in_filter(var_type, var_type_in, matchType="partial")
                if keepType:
                    nKeep += 0
                removeType = apply_not_in_filter(var_type, var_type_not_in, matchType="partial")
                if removeType:
                    nRemove += 0
            if (nKeep == 0):
                continue
            if (nRemove > 0):
                continue

            # allow number of evidence items to be 0
            n_evidence_items = varMap[gene_id][var_id]["n_evidence_items"]
            keepVar = apply_cutoff_filter(n_evidence_items, evidence_items)
            if not keepGene:
                continue

            evidence_types = varMap[gene_id][var_id]["evidence_items"].keys()
# 



            # current variant has passed all variant-level filters
            if var_id_str not in cleanMap[gene_id_str].keys():
                cleanMap[gene_id_str][var_id_str] = {}

            for evidence_type in varMap[gene_id][var_id]["evidence_items"].keys():

                # variant origin can be "NULL" when not available
                variant_origin = varMap[gene_id][var_id]["origin"]
                keepVar = apply_in_filter(variant_origin, var_origin_in, matchType="partial")
                if not keepVar:
                    continue
                removeVar = apply_not_in_filter(variant_origin, var_origin_not_in, matchType="partial")
                if removeVar:
                    continue









            # List of evidence records available for the current variant (can be empty)
            evidence_items = variant_record.evidence_items

            # Iterate through the listed evidence items and store relevant information for this variant
            # Variants which do not have any clinical data associated to them will be directly skipped
            for evidence_record in evidence_items:

                # FIXME: Sanity check that all critical elements are present and non-empty (entrez_name, variant name, evidence_type, disease name)

                # Use uppercase for consistency of the tags
                evidence_status = evidence_record.status.strip().upper()
                evidence_type = evidence_record.evidence_type.strip().upper()
                disease = evidence_record.disease.name.strip().upper()

                evidence_level = evidence_record.evidence_level                     # just in case, should never be None
                variant_origin = evidence_record.variant_origin                     # can be None
                evidence_direction = evidence_record.evidence_direction             # can be None
                clinical_significance = evidence_record.clinical_significance       # can be None

                # Skip records that are not accepted evidence
                if (evidence_status != "ACCEPTED"):
                    continue

                # Skip records that correspond to germline variants
                # The variant_origin field might be blank/empty (None)
                if variant_origin:
                    if re.search("GERMLINE", variant_origin):
                        continue

                # Sanity check for empty evidence direction, clinical significance or level
                # 'NULL' is introduced to distinguish from 'N/A' tag
                if evidence_direction is None:
                    evidence_direction = "NULL"
                else:
                    evidence_direction = evidence_direction.strip().upper()
                if clinical_significance is None:
                    clinical_significance = "NULL"
                else:
                    clinical_significance = clinical_significance.strip().upper()
                if evidence_level is None:
                    evidence_level = "NULL"
                else:
                    evidence_level = evidence_level.strip().upper()

                # Combine the direction and significance of the evidence in one term
                evidence = evidence_direction + ':' + clinical_significance

                # At this point, currently evaluate evidence item has passed all the checks/filters
                # Keep track of the current evidence item under the corresponding variant and gene
                if gene_key not in varMap.keys():
                    varMap[gene_key] = {}

                # Variant name should be unique within gene (found some duplicates but all were submitted, not accepted data)
                if variant_name not in varMap[gene_key].keys():
                    varMap[gene_key][variant_name] = {}
                    varMap[gene_key][variant_name]['id'] = variant_id
                    varMap[gene_key][variant_name]['civic_score'] = civic_score

                    ## Generate list of strings that will be used to match our input variants in CIVIC
                    ## Returned list always has at least length=1 (in this case, containing only variant name)
                    ## For CNV, variant matching is not based on HGVS, so matchStrings will only contain the variant name
                    matchStrings = generate_civic_matchStrings(variant_name, hgvs_expressions, dataType)
                    varMap[gene_key][variant_name]['match_hgvs'] = matchStrings
                    # Keep original HGVS annotations (empty list when nothing is available)
                    # Use uppercase to avoid mismatches due to case
                    varMap[gene_key][variant_name]['hgvs'] = [h.strip().upper() for h in hgvs_expressions]

                    # Include associated variant types (sequence ontology terms). There can be multiple terms
                    varMap[gene_key][variant_name]['types'] = []
                    for vartype_record in variant_record.variant_types:
                        varMap[gene_key][variant_name]['types'].append(vartype_record.name.strip().upper())
                    # Account for empty variant types (can happen)
                    # 'NULL' is introduced to distinguish from 'N/A' tag
                    if not varMap[gene_key][variant_name]['types']:
                        varMap[gene_key][variant_name]['types'] = ["NULL"]

                # FIXME: there is no sanity check for detecting possible variant name duplicates
                if evidence_type not in varMap[gene_key][variant_name].keys():
                    varMap[gene_key][variant_name][evidence_type] = {}
                if disease not in varMap[gene_key][variant_name][evidence_type].keys():
                    varMap[gene_key][variant_name][evidence_type][disease] = {}

                drugs = []
                evidence_drugs = evidence_record.drugs
                for evidence_drug in evidence_drugs:
                    drug_name = evidence_drug.name.strip().upper()
                    if drug_name not in drugs:
                        drugs.append(drug_name)

                # When more than 1 drug are listed for the same evidence item, 'drug_interaction_type' is not null and defines the nature of this multiple drug entry
                drug_interaction = evidence_record.drug_interaction_type
                if drug_interaction is not None:
                    drug_interaction = drug_interaction.strip().upper()
                    # 'Substitutes' indicates that drugs can be considered individually
                    if drug_interaction != "SUBSTITUTES":
                        # Remaining terms ('Sequential' and 'Combination') indicate that drugs should be considered together, so join their names into a single tag
                        # Sort drugs alphabetically to ensure that their order in the combination treatment is always the same
                        drugs.sort()
                        drugs = ["+".join(drugs)]

                if not drugs:
                    # Only non-Predictive evidences and Predictive ones without drugs will have this dummy level
                    # Introduced for consistency purposes within the varMap structure
                    drugs = ["NULL"]

                # Iterate through drugs to add evidences associated to them
                #   For non-Predictive evidences or Predictive with empty drugs, drugs=['NULL']
                #   For Predictive and interaction=None, len(drugs) = 1
                #   For Predictive and interaction='Substitutes', len(drugs)>1
                #   For Predictive and interaction!='Substitutes', len(drugs)=1 (combiantion of several using '+')
                for drug in drugs:
                    if drug not in varMap[gene_key][variant_name][evidence_type][disease].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][drug] = {}
                    if evidence not in varMap[gene_key][variant_name][evidence_type][disease][drug].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][drug][evidence] = {}
                    if evidence_level not in varMap[gene_key][variant_name][evidence_type][disease][drug][evidence].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][drug][evidence][evidence_level] = []
                    # Group all publications associated to the same level. Do not check publication status
                    ## On 25.01.2019, source structure was changed to introduce ASCO abstracts as a source type
## TODO: sanity check for empty ID. Check for type of source?
                    varMap[gene_key][variant_name][evidence_type][disease][drug][evidence][evidence_level].append(evidence_record.source.citation_id.strip())



    return


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

