#!/usr/bin/env python

'''
Query CIViC DB
Requires newer format of variant annotations, i.e.: GENE:variant|variant;GENE:variant;...
Version queries offline cache provided in civicpy module (updating cache is managed by Vipin)
Lourdes Rosano, Nov 2018
'''

import sys
import argparse
from civicpy import civic
import re


### Define dictionary of 1-letter to 3-letter codes for aminoacids. Includes custom aminoacids defined only within CIVIC.
## In CIViC variant names, 'X' is used to indicate a stop codon
aaCodes = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS', 'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 'G': 'GLY', 'H': 'HIS',
           'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET', '*': '*', 'X': '*'}

### Define what combination of directions and clinical significances define 'POSITIVE' and 'NEGATIVE' support
supportDict = {'SUPPORTS': {'SENSITIVITY/RESPONSE':'POSITIVE', 'RESISTANCE':'NEGATIVE', 'REDUCED SENSITIVITY':'NEGATIVE', 'ADVERSE RESPONSE':'NEGATIVE'}, 'DOES NOT SUPPORT': {'RESISTANCE':'UNKNOWN_DNS', 'SENSITIVITY/RESPONSE':'UNKNOWN_DNS', 'REDUCED SENSITIVITY':'UNKNOWN_DNS', 'ADVERSE RESPONSE':'UNKNOWN_DNS'}}

### Define global variable to ensure cache file is only loaded once even if several queries are performed
global isLoad
isLoad = False


'''
Functions
'''

### Functions to query CIVIC and process result into structured dictionary

# Check that a given identifier type is contained in the list of allowed values
def sanity_check_identifier_type(identifier_type):
    if identifier_type not in ["entrez_id", "entrez_symbol", "civic_id"]:
        print("\nError! '%s' is not a valid identifier_type. Please provide one of 'entrez_id', 'entrez_symbol' or 'civic_id'!" %(identifier_type))
        sys.exit(1)
    return None


# Given a list of CIVIC gene records, parse and reformat into an structured dictionary
# Assume gene records already contain all information about associated variants and corresponding evidence items
def reformat_results(results, identifier_type):

# TODO: for now, keep the same filters as we had in the previous version
# TODO: to be done soon; do not filter anything but retrieve all the gene-variant records unchanged. we will apply filtering using functions at a later step

    varMap = {}
    retrieved_genes = []    # keep track of genes that could be retrieved from CIVIC
    no_variants = []        # keep track of genes retrieved from CIVIC but with no variants available
    all_variants = []       # keep track of all variants retrieved from CIVIC

    # Check that id type corresponds to one of the allowed options
    ignore = sanity_check_identifier_type(identifier_type)

    # Iterate individual gene records to retrieve associated variants and evidence information
    for gene_record in results:
        # Retrieve all ids associated to this gene (CIVIC, entrez, symbol)
        gene_civic_id = str(gene_record.id)
        gene_id = str(gene_record.entrez_id)
        # Use uppercase for consistency of the gene symbols
        gene_symbol = gene_record.name.strip().upper()
        # Retrieve variants associated to this gene (can be empty)
        gene_variants = gene_record.variants

        # Use the provided gene id (civic, entrez or symbol) to uniquely identify genes
        if identifier_type == "civic_id":
            gene_key = gene_civic_id
        if identifier_type == "entrez_id":
            gene_key = gene_id
        if identifier_type == "entrez_symbol":
            gene_key = gene_symbol

        # Keep track of gene ids which could be retrieved from CIVIC
        if gene_key not in retrieved_genes:
            retrieved_genes.append(gene_key)

        # NOTE: it seems that only genes having at least 1 variant record available are included in the offline cache (eg. gene ADORA1 has no variants and is found via API but not in the cache)
        # Skip genes that do not have any variants available in CIVIC
        if not gene_variants:
            if gene_key not in no_variants:
                no_variants.append(gene_key)
            continue

        # Iterate variant records associated to the current gene
        # Retrieve all relevant info listed for each variant
        for variant_record in gene_variants:
            # Internal variant id in CIVIC
            variant_id = str(variant_record.id)
            # Keep track of all variants retrieved from CIVIC
            if variant_id not in all_variants:
                all_variants.append(variant_id)
            # Variant name in CIVIC; use uppercase for consistency
            variant_name = variant_record.name.strip().upper()
            hgvs_expressions = variant_record.hgvs_expressions
            # Score to assess the accumulation of evidence for each variant (quantity and quality)
            civic_score = variant_record.civic_actionability_score
            # Sanity check for empty scores
            if civic_score is None:
                civic_score = "NULL"

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

# TODO: iterate through assertions and repeat above process? They are not structured but (highly curated) free text

    return (varMap,retrieved_genes,no_variants,all_variants)


# Given a list of gene identifiers, query CIVIC for known variants and return a structured dictionary with the relevant results
# List of gene ids can be: CIVIC id, entrez id or gene symbol
def query_civic_genes(genes, identifier_type="entrez_symbol"):

    # Check that provided argument is a list (even if length = 1)
    if (not isinstance(genes, list)) or (not genes):
        print("\nError! Please provide a list of gene ids to query in CIVIC.")
        sys.exit(1)

    # Check that id type corresponds to one of the allowed options
    ignore = sanity_check_identifier_type(identifier_type)

    ## Load offline cache of CIVICdb
    # Ensure only loaded once
    # NOTE: it is important that on_stale='ignore' to avoid attempts to update the cache via internet access
    global isLoad
    if not isLoad:
        print("Loading offline CIVIC cache...")
        is_successful = civic.load_cache(on_stale='ignore')
        if is_successful:
            isLoad = True
        else:
            print("\nError! Local cache file could not be successfully loaded!")
            sys.exit(1)

    # Querying functionality provided by module `civicpy` can only use internal CIVIC ids
    # NOTE: when using `civic.get_genes_by_ids()` with CIVIC ids, if any of them is not contained in the cache file, then it will try to directly query the CIVICdb, causing this script to crash

    # Workaround: retrieve all gene records available in the offline cache and parse them to match to input genes
    # Avoid at all costs directly querying the CIVIC db: even for CIVIC ids, we will match to available records from the offline cache
    all_results = civic.get_all_genes()

    # Iterate individual records and retrieve only those matching the provided gene ids
    results = []
    for gene_record in all_results:
        toKeep = False

        if (identifier_type == "civic_id"):
            gene_id = str(gene_record.id)                    # expectation is a single number
            if gene_id in genes:
                toKeep = True

        if (identifier_type == "entrez_id"):
            gene_id = str(gene_record.entrez_id)             # expectation is a single number
            if gene_id in genes:
                toKeep = True

        if (identifier_type == "entrez_symbol"):
            gene_id = gene_record.name.strip()          # expectation is single string
            aliases = gene_record.aliases               # expectation is list of strings
            # Perform union of all gene symbols available for the current record
            tmp_symbols = list(set([gene_id]) | set(aliases))
            # Try to match current gene record (any alias) to the provided gene symbols
            for tmp_symbol in tmp_symbols:
                # Use uppercase to ensure consistency of gene symbols
                this_symbol = tmp_symbol.strip().upper()
                if this_symbol in genes:
                    toKeep = True
                    break

        if toKeep:
            results.append(gene_record)

    # At this point, all CIVIC results for queried genes have been retrieved in a list
    # Process gene records into a dictionary with structured format
    # gene -> variants -> evidence_items
    (dict_results,retrieved_genes,no_variants,all_variants) = reformat_results(results, identifier_type)
    return (dict_results,retrieved_genes,no_variants,all_variants)


### Other functions

## Given the header, find the input columns containing genes and variants/cnvs values (use provided arguments directly)
## For SNVs, also retrieve the input columns containing variant impacts and exons
def getColumnIndices(firstInputLine, dataType):
    firstLineSplit = firstInputLine.split('\t')
    index_geneCol = -1
    index_varCol = -1
    index_impCol = -1
    index_exonCol = -1
    for pos in range(0,len(firstLineSplit)):
        # Avoid mismatches due to case by always using uppercase
        if args.colName_gene.upper() == firstLineSplit[pos].upper():
            index_geneCol = pos
        if args.colName_data.upper() == firstLineSplit[pos].upper():
            index_varCol = pos
        if dataType == 'SNV':
            if args.colName_varImpact.upper() == firstLineSplit[pos].upper():
                index_impCol = pos
            if args.colName_exon.upper() == firstLineSplit[pos].upper():
                index_exonCol = pos
    if (index_geneCol == -1) or (index_varCol == -1):
        print("\nError! Could not match all input columns in header %s." %(firstInputLine))
        sys.exit(1)
    if dataType == 'SNV':
        if (index_impCol == -1) or (index_exonCol == -1):
            print("\nError! Could not match all input columns in header %s." %(firstInputLine))
            sys.exit(1)

    return (index_geneCol,index_varCol,index_impCol,index_exonCol)


## Translate a 1-letter aminoacid code (if it exists) into a 3-letter code
def translate_aa(oneLetter):
    aaNew = None
    if oneLetter.upper() in aaCodes:
        aaNew = aaCodes[oneLetter.upper()]
    return aaNew


## Given a single CIVIC variant name, extract potential HGVS annotations by
## parsing and modifying this string using knowledge on CIVIC naming conventions
## Function only applicable to SNV data
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


## Given a single CIVIC hgvs expression, parse and modify it to ensure
## it complies with the HGVS format followed by the input table
## Function only applicable to SNV data
## TODO: Currently, only modification of p. annotations are supported
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


## Given a single HGVS p. annotation (of the form p.Pro61Cys...),
## extract the start of the string (ie. p.Pro61)
## Function only applicable to SNV data
def extract_p_start(pAnnotation):
    startAnnot = None
    if re.match(r'(P\.[A-Z]+[0-9]+)', pAnnotation):
        startAnnot = re.match(r'(P\.[A-Z]+[0-9]+)', pAnnotation).groups()[0]

    return startAnnot


## Given a CIVIC variant name and its corresponding hgvs expressions (if available),
## generate a list of potential strings that will be used to match the variant to our input
## For CNVs, variant matching is not based on HGVS, since input data is different
def generate_civic_matchStrings(varName, hgvsExpressions, dataType):
    matchStrings = []
    ## For CNVs, only the last step 5) is executed, since variant matching is not done at the HGVS level
    if dataType == 'SNV':
        # 1) First, remove reference from annotation (ie. 'transcriptID:')
        ## This step will be skipped for CNVs
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


## Given a set of SNV annotations (either variant annotations, impacts or exons) in the form 'GENE1:ANNOTATION;GENE2:ANNOTATION...',
## create a dictionary of all existing genes and classify associated annotations under their corresponding type
## Function only applicable to SNV data
def process_lineInfo_snv(rowMap,fieldName,fieldAnnotations):
    # For SNV, there can be both multiple genes and variants per line (both separated with ;). When multiple genes, always multiple variants
    # This is also applicable to variant impact and exon information (correspond to single variants)
    splitAnnotations = fieldAnnotations.split(';')
    for annotation in splitAnnotations:
        ## TODO: Skip empty annotations '.'? Can it happen for variant impacts and exons?
        elements = annotation.strip().split(':')
        # Sanity check for annotation format 'GENE1:annotation;..'
        # Always expect at least 2 elements (variants and impacts). For exon information, len(elements) will be 3
        if len(elements) < 2:
            print("\nError! Unexpected annotation format \'{}\'".format(annotation))
            sys.exit(1)
        # Retrieve gene name from field annotation (always present in required format)
        # Use uppercase to avoid mismatches due to case
        gene = elements[0].upper()
        # Retrieve the last splitted string (always contains the relevant information)
        singleAnnots = elements[-1]
        if gene not in rowMap.keys():
            rowMap[gene] = {}
        if fieldName not in rowMap[gene].keys():
            rowMap[gene][fieldName] = []
        # Now, depending on which field (ie. column) we are processing, we need to process differently
        if (fieldName == "variants"):
            # Split annotation into separate annotation levels (e.g. c.503A>C|p.Val200Val)
            singleAnnots = singleAnnots.strip().split('|')
            # Add single variant annotations to their corresponding gene
            for x in singleAnnots:
                # Remove possible blanks caused by splitting and avoid duplication of variant annotations
                # Use uppercase to avoid mismatches due to case
                if x and (x not in rowMap[gene][fieldName]):
                    rowMap[gene][fieldName].append(x.upper())
        else:
            # Here, no need to split further as we already have desired information; Exon eg. 2/3 or Impact eg. missense_variant&stop_gained
            # Do not remove possible blanks as we need this information as well (ie. if exon information
            # not available, will need to remove corresponding impact as well; and viceversa)
            # Use uppercase to avoid mismatches due to case (only applicable to impacts as exon info is numeric)
            rowMap[gene][fieldName].append(singleAnnots.upper())

    return rowMap


## Given multiple (CNV) genes and one associated annotation (CNV category), create a dictionary of all existing genes and associated CNV annotations
## Function only applicable to CNV data
def process_lineInfo_cnv(rowMap, genes, data):
    # For CNV, every line corresponds to a single CNV category. Usually there are multiple genes per line (separated with ;).
    # Possible CNV categories are: AMP, DEL, GAIN and LOSS
    geneSplit = genes.split(';')
    for gene in geneSplit:
        gene = gene.strip().upper()
        if gene not in rowMap.keys():
            rowMap[gene] = data

    return rowMap


## Given a single CIVIC variant name, return whether it corresponds to a CNV variant record related to exons
## For this, attemp to match the variant name to special CNV exon cases present in CIVIC (eg. EXON 1-2 DELETION, EXON 5 DELETION, 3' EXON DELETION...)
def cnv_is_exon_string(varName):
    exonStrings = ['^EXON [0-9-]+ DELETION$', '^[35\']+ EXON DELETION$', '^EXON [0-9-]+ SKIPPING MUTATION$']
    isExon = False
    for exonString in exonStrings:
        if re.search(exonString, varName):
            isExon = True

    return isExon

## Given a single CIVIC record name, return whether it corresponds to a EXPRESSION record related to exons
## For this, attemp to match the variant name to special EXPRESSION exon cases present in CIVIC (eg. EXON 1-2 EXPRESSION, EXON 5 OVEREXPRESSION...)
def expr_is_exon_string(varName):
    exonStrings = ['^EXON [0-9-]+ EXPRESSION$', '^EXON [0-9-]+ OVEREXPRESSION$', '^EXON [0-9-]+ UNDEREXPRESSION$']
    isExon = False
    for exonString in exonStrings:
        if re.search(exonString, varName):
            isExon = True

    return isExon


## Given a set of one or more annotations (SNV or CNV), generate a list of potential strings
## that will be used to match the variant in CIVIC eg. EXON 15 MUTATION, AMPLIFICATION
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

    return (matchStrings,isExact,isTrueExact)


## Given all CIVIC variant records associated to a given gene, return all those variant names that should correspond to SNV records
## For this, remove the most common CNV and EXPRESSION names existing in CIVIC from the complete set of the gene's variant names
## Additionally, consider other special cases present in CIVIC (eg. EXON 1-2 DELETION/EXPRESSION, EXON 5 DELETION/OVEREXPRESSION, 3' EXON DELETION/UNDEREXPRESSION...)
def civic_return_all_snvs(geneData):
    ## Common CNV record names in CIVIC
    cnvNames = ['AMPLIFICATION','DELETION','LOSS','COPY NUMBER VARIATION']
    ## Common EXPRESSION record names in CIVIC
    exprNames = ['OVEREXPRESSION','UNDEREXPRESSION','EXPRESSION']

    ## All variant names not matching a CNV or EXPRESSION will be returned
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


## Given all CIVIC variant records associated to a given gene, return all those variant names that correspond to CNV records
## For this, attemp to match the gene's variant names to the most common CNV names existing in CIVIC
## Additionally, consider other special CNV cases present in CIVIC (eg. EXON 1-2 DELETION, EXON 5 DELETION, 3' EXON DELETION...)
def civic_return_all_cnvs(geneData):
    ## Common CNV record names in CIVIC
    cnvNames = ['AMPLIFICATION','DELETION','LOSS','COPY NUMBER VARIATION']

    ## All matched variant names will be returned
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


## Given a list of variant annotations for a gene in the input table (HGVS but also synonym descriptive terms eg. EXON 15 MUTATION and positional
## strings eg. p.Val600), and all CIVIC variant records for said gene, attempt mapping of any input annotations to one or more CIVIC records.
## Return tier of the match (exact, positional, no match) and corresponding CIVIC variant records that were matched

## Do exhaustive search as there could be >1 perfect match in cases of redundancy (eg. E55FS and E55RFSTER11
## both translate to p.Glu55fs) and there can be >1 positional matches (for SNV case, 'general' variants
## eg. V600 have preference over the rest). Also, for SNVs in case of having two types of exact matches (truly exact
## eg. V600E and descriptive term eg. EXON 15 MUTATION), give preference to the former as descriptive terms should only
## be used when no true exact match was found.

## For CNV, inputStrings will be limited (usually only 1 or 2), and all will always correspond to exact matches (no positional)
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

    return (tier,matchedVars)


## Check whether a variant name in CIVIC refers to a group of variants (eg. V600)
## Function only applicable to SNV data
def check_general_variant(varName):
    isGeneral = False
    if re.match(r'[A-Z]\d+($|\s\()', varName):
        isGeneral = True
    return isGeneral


def cnv_group_lineInfo_by_gene(lineSplit, index_column):
    ## The first half of the output line (up until index_column) remains the same across all new split lines
    lineRoot = "\t".join(lineSplit[0:index_column])
    ## The second half of the output line (from index_column on) will be split into new lines by the different genes
    ## lineMap is a dictionary containing all line info grouped by gene
    lineMap = {}
    ## Keep track of all drugs associated to the genes in the line. Necessary to generate lineMap
    drugDict = {}
    ## Iterate fields in the second half of the line and create dictionary lineMap
    for indx in range(index_column+1,len(lineSplit)):
        lineMap[indx] = {}
        fieldInfo = lineSplit[indx]
        ## All fields contain annotations separated with ';'
        splitAnnotations = fieldInfo.split(';')
        for annotation in splitAnnotations:
            ## Empty annotations are displayed as '.'
            if annotation == '.':
                continue
            ## All annotations have the structure 'ELEMENT1:ELEMENT2', ie. 'GENE1:annotation;..' or 'DRUG1:annotation;..'
            elements = annotation.strip().split(':')
            # Sanity check for annotation format 'GENE1:annotation;..' or 'DRUG1:annotation;..'
            # Always expect at least 2 elements (gene/drug and info).
            if len(elements) < 2:
                print("\nError! Unexpected annotation format \'{}\'".format(annotation))
                sys.exit(1)

            ## Special annotation case for drugs ie. 'DRUG1:annotation;..'
            if (indx == index_column+5) or (indx == index_column+6):
                # Retrieve drug name from field annotation
                drug = elements[0].upper()
                # Retrieve the second element (contains the relevant ClinTrial information)
                singleAnnot = elements[1]
                if drug not in drugDict.keys():
                    ## Drug should be present in drugDict already, as we should have parsed drug column previously
                    print("\nError! Encountered unknown drug in {} from ClinTrial columns".format(annotation))
                    sys.exit(1)
                ## Assign drug annotation to all associated genes (info in drugDict, from previously parsed columns)
                for tempGene in drugDict[drug]:
                    if tempGene not in lineMap[indx].keys():
                        lineMap[indx][tempGene] = []
                    ## Store the complete annotation including the drug name
                    ## TODO: is the following uniqueness ensurement necessary?
                    if annotation not in lineMap[indx][tempGene]:
                        lineMap[indx][tempGene].append(annotation)

            ## Typical annotation case ie. 'GENE1:annotation;..'
            else:
                # Retrieve gene name from field annotation
                # Use uppercase to avoid mismatches due to case
                gene = elements[0].upper()
                # Retrieve the second element (always contains the relevant information)
                singleAnnot = elements[1]
                if gene not in lineMap[indx].keys():
                    lineMap[indx][gene] = []
                ## Store the annotation
                ## TODO: is the following uniqueness ensurement necessary?
                if singleAnnot not in lineMap[indx][gene]:
                    lineMap[indx][gene].append(singleAnnot)
                ## Special annotation case for DGIDB-drug column, need to keep track of all drugs associated to each gene
                if indx == index_column+2:
                    ## eg. '.;.;RXRA:ALITRETINOIN(11,agonist);RXRA:ADAPALENE( 5,.);.;.'
                    tempSplit = singleAnnot.strip().split('(')
                    # Use uppercase to avoid mismatches due to case
                    predDrug = tempSplit[0].upper()
                    # Means there are more brackets included, so drug name containes brackets
                    if len(tempSplit) > 2:
                        predDrug = "(".join(tempSplit[0:-1]).upper()
                    if predDrug not in drugDict.keys():
                        drugDict[predDrug] = []
                    if gene not in drugDict[predDrug]:
                        drugDict[predDrug].append(gene)

        ## When all annotations in a field are empty (ie. '.;.;.;.'), then the corresponding column index in lineMap is empty as well

    return (lineRoot, lineMap)


# Write information about a single cancer type item into one or more structured strings
# i.e. DISEASE[|DRUG1,DRUG2..](direction, significance(level(PMID,..,PMID),level(..)));
## For 'predictive' evidence (writeDrug=True), keep dictionary of drug support for the current gene/line
## Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
def write_evidences(item, cancer, ct, drugMap, writeDrug=False):
    evidences = []
    # For each drug associated with the given cancer type
    for drug in item.keys():
        # For each evidence associated with the given drug
        # Evidences are simplified by using the combined form 'direction:significance'
        for evidence in item[drug].keys():
            ## For each evidence (ie combination of direction+clin_signf), count how many different evidence items support it
            ## At this stage, we find count evidence items by counting how many different combinations of level+pmids there are for the same drug, disease and evidence
            pmids = []
            # If drug=True, write drug information, i.e. DISEASE|DRUG(..)
            if writeDrug:
                # Always one drug (single or combination with '+')
                outString = cancer + '|' + drug + '('
                ## For 'predictive' evidence, keep dictionary of drug support for the current gene/line
                if drug not in drugMap.keys():
                    drugMap[drug] = {}
                ## Possible values for cancer specificity: 'ct' (specific), 'gt' (general), 'nct' (non-specific)
                if ct not in drugMap[drug].keys():
                    drugMap[drug][ct] = []
            else:
                outString = cancer + '('
            # Split the evidence direction and clinical significance
            direction, clin_signf = evidence.split(':')
            outString += direction + ',' + clin_signf + '('
            # There may be several levels grouped per evidence
            levels = []
            for level in item[drug][evidence].keys():
                # There may be several publications (i.e. PMIDs) grouped per level
                levels.append(level + '(' + ','.join(item[drug][evidence][level]) + ')')
                # Count how many different evidence items support this particular evidence item
                for z in item[drug][evidence][level]:
                    # Distinguish cases where the same publication is used to support different and identical evidence levels (they nonetheless count as separate evidence items)
                    pmids.append(z)
#                     new_z = level + '_' + z
#                     if new_z not in pmids:
#                         pmids.append(new_z)
            outString += ','.join(levels) + '))'
            evidences.append(outString)

            ## For 'predictive' evidence, keep dictionary of drug support for the current gene/line
            ## Each combination of direction + clinSignf has an associated support: POSITIVE, NEGATIVE, UNKNOWN_DNS or UNKNOWN_BLANK
            if writeDrug:
                if ('NULL' in direction) or ('N/A' in direction) or ('NULL' in clin_signf) or ('N/A' in clin_signf):
                    thisSupport = 'UNKNOWN_BLANK'
                else:
                    if direction not in supportDict.keys():
                        print("\nError! Could not find direction %s in support dictionary." %(direction))
                        sys.exit(1)
                    if clin_signf not in supportDict[direction].keys():
                        print("\nError! Could not find clinical significance %s in support dictionary." %(clin_signf))
                        sys.exit(1)
                    thisSupport = supportDict[direction][clin_signf]

                ## Keep track of number of occurrences for each support type for the given drug
                ## Here, take into account the number of supporting PMIDs associated to each evidence item
                for z in pmids:
                    drugMap[drug][ct].append(thisSupport)

    return (evidences, drugMap)


# Return in a list all evidence items (in written form) for a given evidence type 
# i.e. DISEASE1[|DRUG1,DRUG2..](direction1,significance1(level1(PMID,..,PMID),level2(..)));
#      DISEASE1[|DRUG1,DRUG2..](direction2,significance2(level1(PMID,..,PMID),level2(..)));
# For each evidence type, return either:
#   - Info on white listed cancer types (eg. 'Melanoma')
#   - If previous is not available, info on high level cancer types (eg. 'Cancer', 'Solid tumor')
#   - If previous is not available, info on all available cancer types for the given variant (except those included in the black list)
## For 'predictive' evidence (writeDrug=True), keep dictionary of drug support for the current gene/line
## Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
def match_cancer_and_return_evidences(vardata, evidence_type, cancerTypes, cancerTypes_notSpec, highLevelTypes, drugMap, writeDrug=False):
    ## Evidences returned for the matched cancer specificity
    evidences = []

    ## Keep track of cancer types that did not pass the black-list filter
    blackMatched = []
    ## Keep track of cancer types that passed the black-list filter
    cleanSet = []
    ## Final list of matched cancer types (associated to either 'ct', 'gt' or 'nct')
    matched = []
    ## Keep track of which type of cancer specificy was matched in the end ('ct', 'gt' or 'nct')
    ct = ''

    # If the given evidence type is not present for the variant, list 'evidences' will be empty
    if evidence_type in vardata.keys():
        ## 1) First, remove cancer types that partially match black-listed terms (if any are provided)
        # NOTE: PARTIAL matches to the black list are allowed! eg:
        #   - including 'small' will remove 'non-small cell lung cancer' and 'lung small cell carcinoma'
        #   - including 'non-small' will remove 'non-small cell lung cancer' but not 'lung small cell carcinoma'
        if cancerTypes_notSpec:
            # Iterate available cancer types and keep track of those that partially match to black list (there can be several)
            for cancerType in vardata[evidence_type].keys():
                # To find partial matches, it is necessary to iterate through both lists (input and civic)
                for blackListed in cancerTypes_notSpec:
                    # Search for partial match of INPUT cancer type in CIVIC cancer type e.g. 'Melanoma' (input list) in 'Skin Melanoma' (CIVIC) and not opposite
                    if blackListed in cancerType:
                        if cancerType not in blackMatched:
                            blackMatched.append(cancerType)
            # Iterate available cancer types once again to retrieve those that passed the black list filter
            for cancerType in vardata[evidence_type].keys():
                # Retrieve valid cancer types only
                if cancerType in blackMatched:
                    continue
                if cancerType not in cleanSet:
                    cleanSet.append(cancerType)

        ## If no black list was provided, then all available cancer types constitute the clean set
        else:
            cleanSet = list(vardata[evidence_type].keys())

        ## 2) Now, iterate the list of "allowed" cancer types (ie. passing the black list filter) and attempt to match to white-listed terms
        # NOTE: again, PARTIAL matches to the white list are allowed! eg:
        #   - including 'melanoma' will match 'melanoma', 'skin melanoma' and 'uveal melanoma', but not 'skin cancer' (hypothetical)
        #   - including 'uveal melanoma' will only match 'uveal melanoma'
        for cleanType in cleanSet:
            # To find partial matches, it is necessary to iterate through both lists (input and civic)
            for whiteListed in cancerTypes:
                # Search for partial match of INPUT cancer type in CIVIC cancer type e.g. 'Melanoma' (input list) in 'Skin Melanoma' (CIVIC) and not opposite
                if whiteListed in cleanType:
                    # Keep track of cancer types that passed the white list filter
                    if cleanType not in matched:
                        ct = 'ct'
                        matched.append(cleanType)

        ## 3) When nothing could be matched to the white-list terms, attempt a second-best match strategy to higher-level cancer types
        ## In CIVIC, some 'general' cancer types are included as disease, eg. 'cancer' or 'solid tumor'. Hence, must be exact matches since they are DB-specific
        # NOTE: here, only PERFECT matches are allowed!
        #   - including 'cancer' will only match 'cancer' and not 'lung cancer'
        if not matched:
            for cleanType in cleanSet:
                if cleanType in highLevelTypes:
                    if cleanType not in matched:
                        ct = 'gt'
                        matched.append(cleanType)

        ## 4) If nothing could be matched to either the white-list or higher-level cancer types, return all 'allowed' (ie. not black-listed) cancer types available in CIVIC
        ## These will be considered as non-specific cancer types (ie. off label)
        if not matched:
            for cleanType in cleanSet:
                if cleanType not in matched:
                    ct = 'nct'
                    matched.append(cleanType)

        ## Now, list 'matched' contains all cancer types for which their evidences items will be returned
        ## They can correspond to either 'ct' (from white-list), 'gt' (from high-level list) or 'nct' (if nothing could be matched, return all that is available)
        for cancerType in matched:
            ## Return evidence items in already formatted strings
            ## Also, return dictionary of drug support for the current gene/line
            ## Format: drug -> ct -> [support1,support2,..,support1] (keep track of total number of evidence items)
            (strings,drugMap) = write_evidences(vardata[evidence_type][cancerType], cancerType, ct, drugMap, writeDrug)
            for s in strings:
                evidences.append(s)

    return (evidences, drugMap)



## Parse and process dictionary of drug support containing available drug prediction support evidence (if any)
## Format: drug -> ct -> [support1,support2,..,support1] (keep track of total number of evidence items)
## Prioritize based on ct and apply majority vote to agree on a CIVIC support decision; one support decision per available drug
def process_drug_support(drugMap):

    supportStrings = []
    ## If no predictive evidence is available, dictionary will be empty (as well as returned list)
    for drug in drugMap.keys():
        ## Prioritize cancer specificity: ct > gt > nct
        thisCT = ''
        if 'ct' in drugMap[drug].keys():
            thisCT = 'ct'
        elif 'gt' in drugMap[drug].keys():
            thisCT = 'gt'
        elif 'nct' in drugMap[drug].keys():
            thisCT = 'nct'
        else:
            print("\nError! Unexpected ct case for gene %s." %(gene))
            sys.exit(1)

        ## Given the selected ct, count number of occurrences for each possible support type (if any)
        count_pos = drugMap[drug][thisCT].count('POSITIVE')
        count_neg = drugMap[drug][thisCT].count('NEGATIVE')
        count_unk = drugMap[drug][thisCT].count('UNKNOWN_BLANK')
        count_dns = drugMap[drug][thisCT].count('UNKNOWN_DNS')

        ## Pool UNKNOWN_BLANK and UNKNOWN_DNS together (as both result in unknown CIVIC support)
        count_total_unk = count_unk + count_dns
        ## Sanity check that there is at least some support
        if (count_pos == 0) and (count_neg == 0) and (count_total_unk == 0):
            print("\nError! Unexpected support case for gene %s." %(gene))
            sys.exit(1)

        ## Resolve contradicting evidence (if any) by majority vote
        tempSupport = ''
        ## For this, pool UNKNOWN_BLANK and UNKNOWN_DNS together
        ## Whenever there is a tie of "confident" (pos or neg) vs "non-confident" (unk), choose the confident one
        if (count_total_unk > count_pos) and (count_total_unk > count_neg):
            tempSupport = "CIVIC_UNKNOWN"
        elif count_pos == count_neg:
            tempSupport = "CIVIC_CONFLICT"
        elif (count_pos > count_neg) and (count_pos >= count_total_unk):
            tempSupport = "CIVIC_SUPPORT"
        elif (count_neg > count_pos) and (count_neg >= count_total_unk):
            tempSupport = "CIVIC_RESISTANCE"
        else:
            print("\nError! Unexpected support case for gene %s." %(gene))
            sys.exit(1)

        ## Build support string for current drug(+gene in line)
        ## Format: DRUG:CT:SUPPORT
        drugSupport = drug + ':' + thisCT.upper() + ':' + tempSupport
        supportStrings.append(drugSupport)

    return supportStrings


'''
Script
'''

parser = argparse.ArgumentParser(description='Query CIViC to retrieve drug information for snvs or cnvs.')
parser.add_argument('--inputTable', dest='inputFile', required=True, help='Input table with genes and variants, needs to be tab separated.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of the output file.')
parser.add_argument('--cancerTypeList', dest='cancerTypeList', required=True, help='Comma-separated list of accepted cancer types. Partial matches will be sought.')
parser.add_argument('--blackList', dest='blackList', required=True, help='Comma-separated list of not accepted cancer types. Partial matches will be sought, and those matched will be excluded from the list of considered evidence.')
parser.add_argument('--highLevelList', dest='highLevelList', required=True, help='Comma-separated list of high level cancer types (e.g. Cancer). Only exact matches will be sought. The user should be aware that results for high level cancer types will only be retrieved when no match for cancer specific types (ie. --cancerTypeList) is found.')
parser.add_argument('--colName_gene', dest='colName_gene', required=True, help='Name of column containing gene symbols.')
parser.add_argument('--colName_data', dest='colName_data', required=True, help='Name of column containing variant annotations (SNV) or cnv categories (CNV).')
parser.add_argument('--dataType', dest='dataType', required=True, choices=['snv','cnv'], default='snv', help='Type of data in inputTable. Possible options are: snv, cnv.')
parser.add_argument('--colName_varImpact', dest='colName_varImpact', required=False, help='Name of column containing variant impacts (for SNV only).')
parser.add_argument('--colName_exon', dest='colName_exon', required=False, help='Name of column containing variant exon information (for SNV only).')

args = parser.parse_args()

## Check input parameters
dataType = args.dataType.upper() # SNV or CNV

## For SNV, two additional arguments are required (ie. column names for variant impact and exon information)
if (dataType=='SNV') and (args.colName_varImpact is None or args.colName_exon is None):
    parser.error("--dataType 'snv' requires --colName_varImpact and --colName_exon.")

print("\nParameters:\n inputTable: %s\n outFile: %s\n dataType: %s\n colName_gene: %s\n colName_data: %s\n colName_varImpact: %s\n colName_exon: %s\n" %(args.inputFile,args.outFile,args.dataType,args.colName_gene,args.colName_data,args.colName_varImpact,args.colName_exon))

infile = open(args.inputFile,'r')
firstInputLine = infile.readline().strip()

## Retrieve indices for all relevant columns (ie. gene and variant annotation; for SNV also variant impact and exon info)
## If not all are found, script will exit with an error
## When dataType is CNV, index_impCol=-1 and index_exonCol=-1 since they will be disregared
(index_geneCol,index_varCol,index_impCol,index_exonCol) = getColumnIndices(firstInputLine, dataType)

## Already write new header into output file
outfile = open(args.outFile,'w')
outHeader = firstInputLine
outHeader += "\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\tCIViC_Drug_Support\tCIViC_Predictive\tCIViC_Diagnostic\tCIViC_Prognostic\tCIViC_Predisposing"
outfile.write(outHeader + "\n")

cancerTypeList = args.cancerTypeList
blackList = args.blackList
highLevelList = args.highLevelList
cancerTypeList = [s.strip().upper() for s in cancerTypeList.split(',')]
blackList = [s.strip().upper() for s in blackList.split(',')]
highLevelList = [s.strip().upper() for s in highLevelList.split(',')]

print('\nWhite listed cancer types: {}'.format(','.join(cancerTypeList)))
print('Black listed cancer types: {}'.format(','.join(blackList)))
print('High level cancer types: {}'.format(','.join(highLevelList)))

## Sanity check for "empty" input lists (ie. in the form of [''])
## Turn them into "real" empty lists for implementation purposes
if cancerTypeList == ['']:
    cancerTypeList = []
if blackList == ['']:
    blackList = []
if highLevelList == ['']:
    highLevelList = []

### Retrieve all input genes to query CIViC for variant info
genes = []
for line in infile:
    lineSplit = line.strip().split("\t")
    # Account for several genes per line (separated with ';')
    # Uppercase as a sanity check
    multiple = lineSplit[index_geneCol].strip().upper().split(";")
    for gene in multiple:
        # Check gene is not empty as a sanity check
        if gene:
            # Uppercase as a sanity check
            gene = gene.strip().upper()
            if gene not in genes:
                genes.append(gene)
infile.close()

# When no genes were found in the input file, write empty files (necessary for snakemake pipeline) and exit without error
# eg. empty input file containing only header because patient had no variants at all
if not genes:
    print("\nDid not find any genes in column \'{}\' of file {}".format(args.colName_gene, args.inputFile))
    outfile.close()
    sys.exit(0)


### Query CIViC for the genes of interest

print('\nTotal # genes to query: {}'.format(len(genes)))
print('\nRetrieving data from CIViC...')

(varMap,retrieved_genes,no_variants,all_variants) = query_civic_genes(genes, identifier_type="entrez_symbol")

retrieved = set(retrieved_genes)
unmatched = list(set(genes) - retrieved)
print("Found %s/%s genes in CIVIC associated to %s variants. Found %s CIVIC genes that had no variants available." %(len(retrieved_genes),len(genes),len(all_variants),len(no_variants)))
print('\nGenes with no CIVIC data: {}'.format(','.join(unmatched)))


### Iterate through input table once more and create output table
### Only new format of input table is allowed!

# Keep track of all matches and non-matches
exactMatches = 0           # Tier 1
posMatches = 0             # Tier 2
noMatches = 0              # Tier 3
notFound = 0               # Tier 4
unmatchedVars = []         # Rows classified as Tier 3 (whole annotation field)
genesNotFound = []         # Genes classified as Tier 4

print('\nMatching input variants...')
infile = open(args.inputFile,'r')
next(infile)
for line in infile:
    lineSplit = line.strip().split("\t")
    outfileLine = line.strip()
    ## For CNV: format of output file will be slightly different in this case, having just one gene per row
    ## Already process current line to return the root (ie. part of the line that will stay unchanged in the output)
    ## and a dictionary containing all the remaining info (ie. part of the line that will be split by row) grouped by gene
    if dataType == 'CNV':
        (outfileLine_cnv_root,lineMapCnv) = cnv_group_lineInfo_by_gene(lineSplit, index_geneCol)

    ## Process current line and create a dictionary of all existing data annotations (SNV, CNV) grouped by gene
    rowMap = {}
    # Note that rowMap will have a different structure depending on dataType:
    #  - For SNV: gene -> annotationType -> value (either SNV annotation, variant impact or variant exon)
    #  - For CNV: gene -> value (CNV category)
    if dataType == 'SNV':
        # For SNV: Sequentially process available annotations within line to create a dictionary of all present genes,variants,impacts and exons in current line
        rowMap = process_lineInfo_snv(rowMap, "variants", lineSplit[index_varCol])
        rowMap = process_lineInfo_snv(rowMap, "impacts", lineSplit[index_impCol])
        rowMap = process_lineInfo_snv(rowMap, "exons", lineSplit[index_exonCol])
    else:
        # For CNV: Create a dictionary of all genes and associated annotations (CNV category) in current line
        rowMap = process_lineInfo_cnv(rowMap, lineSplit[index_geneCol], lineSplit[index_varCol])

    # Iterate genes and variants within the row to determine their 'tier' and matches in CIVIC
    # For SNV: a tier hierarchy will be applied when several genes are present in the same line, so only the highest tier of the row is reported
    tierToWrite = None
    toWrite = {}
    for gene in rowMap.keys():
        ## If gene is in CIVIC, possible tier levels are 1,2,3
        if gene in varMap.keys():
            # For SNV: gene -> annotationType -> value (either SNV annotation, variant impact or variant exon)
            if dataType == 'SNV':
                    # Get all variant HGVS annotations under the same gene symbol
                    # From e.g. GENE1:c.503A>C|p.Val200Val;GENE1:c.502A>C|p.Val200Val, we would get: [c.503A>C, c.502A>C, p.Val200Val]
                    allVariants = rowMap[gene]["variants"]
                    # Also get variant impacts and exons, which will be used as well to match the variants in CIVIC
                    allImpacts = rowMap[gene]["impacts"]
                    allExons = rowMap[gene]["exons"]
            # For CNV: gene -> value (CNV category)
            else:
                    ## Opposite to the SNV case, for CNVs there is only one available annotation ie. string for CNV (cnv category)
                    allVariants = rowMap[gene]
                    allImpacts = []
                    allExons = []

            # Function for generating additional synonym strings for input annotations (eg. EXON 15 MUTATION)
            # isExact is a list of equal length to inputStrings, indicating whether a given string corresponds
            # to an exact (True) or positional (False) match
            inputStrings,isExact,isTrueExact = generate_input_matchStrings(varAnnotations=allVariants, dataType=dataType, impactAnnots=allImpacts, exonAnnots=allExons)

            # Attempt match to a CIViC record using all available variant annotations for the given gene
            # A tier and match will always be retrieved (either exact, positional or all gene variants when no match)
            tier,matched = matchVariants(inputStrings, isExact, isTrueExact, varMap[gene], dataType)

        ## If gene is not in CIVIC, tier level is 4
        else:
            # Keep track of genes with no variant data in CIViC
            if gene not in genesNotFound:
                genesNotFound.append(gene)
            tier = 4
            # No variants will be reported in this case (as the information is not available)
            matched = []

        ## For SNV: tier hierarchy will be applied within the line, and only genes associated to the highest tier will be reported
        if dataType == 'SNV':
            # When first match of the row, just take it
            if tierToWrite is None:
                tierToWrite = tier
                # Add results under current gene (allows for same-tier results with different genes)
                toWrite[gene] = {'matched':[], 'tier':None}
                toWrite[gene]['matched'] = matched
                toWrite[gene]['tier'] = tierToWrite
            # When already had a match previously in the row, apply hierarchy of tiers
            else:
                ## When returned tier is the same as current match, add current gene to results as well
                if tier == tierToWrite:
                    # Add results under current gene (previous matches with the same tier will be under a different gene symbol)
                    # TODO: we are iterating different genes, so should not exist in the dictionary
                    toWrite[gene] = {'matched':[], 'tier':None}
                    toWrite[gene]['matched'] = matched
                    toWrite[gene]['tier'] = tierToWrite
                ## When returned tier is better than current match, exchange it
                # TODO: is it dangerous to treat tiers (classification) as numbers?
                elif tier < tierToWrite:
                    tierToWrite = tier
                    # Stored results are not valid anymore since current tier has higher priority,
                    # so initialize dictionary and store new results under current gene
                    toWrite = {}
                    toWrite[gene] = {'matched':[], 'tier':None}
                    toWrite[gene]['matched'] = matched
                    toWrite[gene]['tier'] = tierToWrite

        ## For CNV: in this case there is no tier hierarchy. All genes correspond to different effects of the same CNV, hence all genes should be reported.
        ## Format of output file will be slightly different in this case (ie. info of one gene per row)
        elif dataType == 'CNV':
            toWrite[gene] = {'matched':[], 'tier':None}
            toWrite[gene]['matched'] = matched
            toWrite[gene]['tier'] = tier


    ### Write result to output
    ### Once here, all genes and variants within the row have been classified into tiers

    ## For SNV: the final genes and variants to report in output have been selected through tier hierarchy (1>2>3>4)
    if dataType == 'SNV':
        ## Only one tier to report per line
        outfileLine += '\t' + str(tierToWrite)
        ## Keep track of the number and types of matches
        if tierToWrite == 4:
            notFound += 1
            ## If tier=4, this means no gene from current line was found in CIVIC, so output line should be empty
            outfileLine += '\t.\t.\t.\t.\t.\t.\t.'
        else:
            if tierToWrite == 1:
                exactMatches += 1
            elif tierToWrite == 2:
                posMatches += 1
            elif tierToWrite == 3:
                noMatches += 1
                # Keep track of unmatched variants for development purposes (whole row)
                if lineSplit[index_varCol] not in unmatchedVars:
                    unmatchedVars.append(lineSplit[index_varCol])

    ## Iterate all genes to be reported
    geneScores = []     # to keep track of the CIVIC score for each matched variant record
    geneVarTypes = []   # to keep track of the CIVIC variant type for each matched variant record
    resultMap = {}      # to keep track of all evidence items (for each evidence type) retrieved from each matched CIVIC record
    drugMap = {}        # to keep track of support for CIVIC drug predictions for the given line (across all matched records)

    ## For SNV: only those with highest tier
    ## For CNV: all genes in the line
    for geneToWrite in toWrite.keys():
        ## For CNV: format of output file is one line per single gene, and all genes from each line will be reported
        ## Output line will be different than SNV format: split current line into several lines, one per gene (group line info by gene)
        ## Also, there will be a different tier reported per gene/line
        if dataType == 'CNV':
            geneScores = []
            geneVarTypes = []
            resultMap = {}
            drugMap = {}

            ## Split current line into several lines, one per gene
            ## The first half of the line (up until 'Genes' column) remains the same across all new split lines (root)
            outfileLine_cnv = outfileLine_cnv_root + '\t' + geneToWrite
            ## The second half of the output line (from 'Genes' column on, ie. DGIDB, ClinicalTrials, etc.) will be split into new lines by the different genes
            ## Write all the info from original line grouped by the current gene
            for colIndx in sorted(lineMapCnv.keys()):
                if geneToWrite in lineMapCnv[colIndx].keys():
                    geneInfo = []
                    ## TODO CNV: introduce exception for clinical trial columns? Drugs were not associated to genes in original input file
                    for annotation in lineMapCnv[colIndx][geneToWrite]:
                        geneInfo.append(geneToWrite + ':' + annotation)
                    if geneInfo:
                        outfileLine_cnv += '\t' + ';'.join(geneInfo)
                    ## TODO: following is internal sanity check, geneInfo should always contain info
                    else:
                        print('\nWarning! No annotations available for column {} and gene {}'.format(colIndx,geneToWrite))
                        sys.exit(1)
                else:
                    outfileLine_cnv += '\t.'

            ## Now that original line has been modified to oblige with new output format, add CIVIC annotations
            ## Report CIVIC tier for the current gene, as in the CNV case there is no tier hierarchy
            tierToWrite = toWrite[geneToWrite]['tier']
            if not toWrite[geneToWrite]['matched'] and (tierToWrite == 3):
                outfileLine_cnv += '\tSNV_only'
            else:
                outfileLine_cnv += '\t' + str(tierToWrite)
            ## Keep track of the number and types of matches
            if tierToWrite == 4:
                notFound += 1
                ## If tier=4, this means current gene was not found in CIVIC, so output line should be empty
                outfileLine_cnv += '\t.\t.\t.\t.\t.\t.\t.'
            else:
                if tierToWrite == 1:
                    exactMatches += 1
                elif tierToWrite == 2:
                    posMatches += 1
                elif tierToWrite == 3:
                    noMatches += 1
                    # Keep track of unmatched CNVs for development purposes (gene + cnv_category)
                    unmatchedPair = geneToWrite + ':' + rowMap[geneToWrite]
                    if unmatchedPair not in unmatchedVars:
                        unmatchedVars.append(unmatchedPair)

        ## Iterate all matched variants for the current reported gene
        for varName in toWrite[geneToWrite]['matched']:
            vardata = varMap[geneToWrite][varName]
            # Store the score for each individual variant for the current gene
            geneScores.append(geneToWrite + ':' + varName + ':' + str(vardata['civic_score']))
            # Store the types for each individual variant (there may be more than 1) for the current gene
            geneVarTypes.append(geneToWrite + ':' + varName + ':' + ','.join(vardata['types']))
            # Process CIVIC evidence associated to the current variant, and prepare results to be written to output
            # Important to be consistent with the column order in the outFile header
            for evidence_type in ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']:
                evidence_items = []
                if evidence_type not in resultMap.keys():
                    resultMap[evidence_type] = []
                if evidence_type == 'PREDICTIVE':
                    writeDrug = True
                else:
                    writeDrug = False
                # Returns a list of evidence items for the given CIVIC record and evidence type, already in written form
                # Only evidences coming from the matched cancer specificity are returned (either 'ct' if in white list, 'gt' in in high level list or 'nct' if unspecific)
                # Also, returns a dictionary of drug support. Format: drug -> ct -> [support1,support2,..,support1] (keep track of total number of evidence items)
                (evidence_items,drugMap) = match_cancer_and_return_evidences(vardata, evidence_type, cancerTypeList, blackList, highLevelList, drugMap, writeDrug)
                # Add resulting evidence strings for the current variant under the appropriate evidence type
                for i in evidence_items:
                    # Add matched CIVIC variant name to each evidence item string
                    resultMap[evidence_type].append(geneToWrite + ':' + varName + ':' + i)


        ## For CNV: reporting of results is done at the gene level, since all genes are reported and each gene placed in a different line
        if dataType == 'CNV':
            ## After parsing all matched variants for the current gene, parse and process available drug prediction support evidence (if any)
            ## Prioritize based on ct and apply majority vote to agree on a CIVIC support decision; one support decision per available drug
            ## If no predictive evidence is available, returned list will be empty
            support_strings = process_drug_support(drugMap)

            ### Write to output results for current gene/line
            if tierToWrite < 4:
                ## Report variant scores for the current gene
                if geneScores:
                    outfileLine_cnv += '\t' + ';'.join(geneScores)
                else:
                    outfileLine_cnv += '\t.'
                ## Report variant types for the current gene
                if geneVarTypes:
                    outfileLine_cnv += '\t' + ';'.join(geneVarTypes)
                else:
                    outfileLine_cnv += '\t.'

                ## Report CIVIC support decision for all drugs predicted in the current line (if any)
                ## Format: DRUG1:CT:SUPPORT;DRUG2:CT:SUPPORT;..;DRUGN:CT:SUPPORT
                if support_strings:
                    outfileLine_cnv += '\t' + ';'.join(support_strings)
                else:
                    outfileLine_cnv += '\t.'

                ## Report CIVIC evidence associated to variants for the current gene
                ## Important to be consistent with the column order in the outFile header
                for evidence_type in ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']:
                    isEmpty = True
                    if evidence_type in resultMap.keys():
                        if resultMap[evidence_type]:
                            isEmpty = False
                            outfileLine_cnv += '\t' + ';'.join(resultMap[evidence_type])
                    # Check for cases when a field is empty for all variants
                    if isEmpty:
                        outfileLine_cnv += '\t.'

            ## Write CNV line (for one single gene) to output
            outfile.write(outfileLine_cnv + '\n')


    ## For SNV: reporting of results is done at the line level, and only genes associated to the line's tier level will be reported, all in the same line
    ## Now that all genes and variants to write for current line have been iterated, write results to output
    if dataType == 'SNV':
        ## After parsing all matched records for the current variant/line, parse and process available drug prediction support evidence (if any)
        ## Prioritize based on ct and apply majority vote to agree on a CIVIC support decision; one support decision per available drug
        ## If no predictive evidence is available, returned list will be empty
        support_strings = process_drug_support(drugMap)

        ### Write to output results for current variant/line
        if tierToWrite < 4:
            ## Report variant scores for all relevant genes from the line
            if geneScores:
                outfileLine += '\t' + ';'.join(geneScores)
            else:
                outfileLine += '\t.'
            ## Report variant types for all relevant genes from the line
            if geneVarTypes:
                outfileLine += '\t' + ';'.join(geneVarTypes)
            else:
                outfileLine += '\t.'

            ## Report CIVIC support decision for all drugs predicted in the current line (if any)
            ## Format: DRUG1:CT:SUPPORT;DRUG2:CT:SUPPORT;..;DRUGN:CT:SUPPORT
            if support_strings:
                outfileLine += '\t' + ';'.join(support_strings)
            else:
                outfileLine += '\t.'

            ## Report CIVIC evidence associated to variants for all relevant genes from the line
            ## Important to be consistent with the column order in the outFile header
            for evidence_type in ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']:
                isEmpty = True
                if evidence_type in resultMap.keys():
                    if resultMap[evidence_type]:
                        isEmpty = False
                        outfileLine += '\t' + ';'.join(resultMap[evidence_type])
                # Check for cases when a field is empty for all variants
                if isEmpty:
                    outfileLine += '\t.'

        ## Write SNV line (for all relevant genes from the line)
        outfile.write(outfileLine + '\n')

infile.close()
outfile.close()

print('\nTotal # exact matches: {}'.format(exactMatches))
print('Total # positional matches: {}'.format(posMatches))
print('Total # gene-only matches: {}'.format(noMatches))
stringNotFound = 'Total # rows without variant data: {}'
if dataType == 'CNV':
    stringNotFound = 'Total # genes without variant data: {}'
print(stringNotFound.format(notFound))
print('---------------------')
# TODO: only report # of genes instead of whole list?
print('Genes with no CIViC variant data:\n{}'.format(','.join(genesNotFound)))
print('\nUnmatched variants:\n {}'.format('\n '.join(unmatchedVars)))
