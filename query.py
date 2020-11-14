#!/usr/bin/env python

'''
Functions for asynchronous querying to CIVIC database
Lourdes Rosano, November 2020
'''

'''
Required Python modules
'''

import sys
import asyncio
import aiohttp

'''
Required internal functions
'''

import read_and_write


'''
Functions
'''

async def fetch(url, session):
    async with session.get(url) as response:
        # TODO: allow for more status codes (only 4XX and 5xx raise error)
        if response.status != 200:
            print("Error! Query {} failed (Status code: {})".format(url, response.status))
            sys.exit(1)
        return await response.json()

async def bound_fetch(sem, url, session):
    # Getter function with semaphore.
    async with sem:
        return await fetch(url, session)

# Query a list of elements (ids or id batches, ie. several ids per call allowed) using urlbase
# batch indicates whether elements correspond to single ids or batches of ids
async def run(elements, urlbase, batch=False, semaphore=100):
    tasks = []

    # Create instance of Semaphore to limit number of concurrent requests
    sem = asyncio.Semaphore(semaphore)

    # Fetch all responses within one Client session, keep connection alive for all requests
    async with aiohttp.ClientSession() as session:
        for e in elements:
            if batch:
                # Join all IDs from the same batch for the corresponding call
                idString = ','.join(e)
            else:
                idString = e
            task = asyncio.ensure_future(bound_fetch(sem, urlbase.format(idString), session))
            tasks.append(task)

        responses = await asyncio.gather(*tasks, return_exceptions=True)
        return responses

# Execute asynchronous API calls to query CIVIC (can be used for both genes and variants)
# Included into a function in order to make error-handling work
def query_civic(elements, urlbase, batch=False):
    old_loop = asyncio.get_event_loop()
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    future = asyncio.ensure_future(run(elements, urlbase, batch))
    results = loop.run_until_complete(future)
    loop.close()
    asyncio.set_event_loop(old_loop)
    return results

# Ensure that given URL contains trailing "/"
def check_url_end(url):
    new_url = url
    if not new_url.endswith("/"):
        new_url += "/"
    return new_url


# Build API URL for the genes endpoint
def get_api_genes():
    (url,end_genes,end_variants) = get_api_access_info()
    # Ensure that both URLs contain trailing "/"
    url = check_url_end(url)
    end_genes = check_url_end(end_genes)
    url_genes = url + end_genes
    return url_genes

# Build API URL for the variants endpoint
def get_api_variants():
    (url,end_genes,end_variants) = get_api_access_info()
    # Ensure that both URLs contain trailing "/"
    url = check_url_end(url)
    end_variants = check_url_end(end_variants)
    url_variants = url + end_variants
    return url_variants

# Error handling function
# When an exception occurs for a given query, it will be returned within the results
# TODO: check that error in check_query would be propagated
# TODO: check that None is a valid return
def check_query(results):
    for r in results:
        if isinstance(r,Exception):
            print("Something went wrong! Exception: %s" %(r))
            sys.exit(1)
    return None

# Query a list of genes
# For identifier_type="entrez_symbol", only perfect matches are allowed (ie. even known aliases will not give a match)
# TODO: check if there is a limitation on the number of batches provided to the gene batch query
def query_civic_genes(genes, identifier_type="entrez_symbol", batch=True, batch_size=200):
    # Sanity check that provided identifier_type is valid
    # Valid identifiers are supplied via the civic.yml file
    dict_ids = get_identifier_type_to_civic()
    if identifier_type not in dict_ids.keys():
        print("Error! Provided identifier_type '%s' is unknown. Please use one of the following: %s." %(identifier_type, list(dict_ids.keys())))
        sys.exit(1)

    base_url = get_api_genes()
    url = base_url + "{}?identifier_type=" + identifier_type

    # Validate gene identifiers before querying
    # Symbols containing a dot "." cause the query to fail, so remove them from the list (if any)
    # Remove empty or duplicated identifiers from the list (if any)
    print("Validating input genes...")
    genes = validate_genes(genes)

    # Check that at least one valid gene was provided
    # NOTE: for gene symbols, check that not all were removed due to the presence of dots
    if not genes:
        print("Error! Please provide at least one valid gene identifier!")
        sys.exit(1)

    # If batch=True, generate batches of genes according to batch_size
    if batch:
        geneBatches = [genes[x:x+batch_size] for x in range(0, len(genes), batch_size)]
    else:
# TODO: check that geneBatches is identical to genes in this case
        geneBatches = genes
    # Returns nested dict of gene records
    # NOTE: query will crash when only a single gene was provided and it could not be found in CIVIC
    print("Querying CIViC...")
    tmp = query_civic(elements=geneBatches, urlbase=url, batch=batch)
    print("Done")
    check_query(tmp)

    # Flatten the nested list to have plain list of gene records
    results = []
    # NOTE: results will be empty list when >1 genes were queried and none could be found in CIVIC
    for recordList in tmp:
        # NOTE: queries containing one single gene will contain a single list (ie. gene) instead of a nested list (ie. batches and genes)
        if isinstance(recordList, dict):
            results.append(recordList)
        else if isinstance(recordList, list):
            for record in recordList:
                if isinstance(record, dict):
                    results.append(record)
    return results

# Query a list of variants CIVIC ids
def query_civic_variants(variants):
    url = get_api_variants()
    print("Querying CIViC...")
    results = query_civic(elements=variants, urlbase=url, batch=False)
    print("Done")
    check_query(results)
    return results


# NOTE: there can be gene records without any variant associated
# Given a single gene record, retrieve its associated variant records and return dict of variant -> value
# NOTE: genes without any associated variants will return an empty dict
def get_variants_from_gene(gene_record, var_id="id", var_value="name"):
    entry_name = "variants"
    if entry_name not in gene_record.keys():
        print("Error! Entry '%s' not found in the following gene record: %s" %(entry_name,gene_record))
        sys.exit(1)
    # Gene records without any variants will have an empty list
    variant_records = gene_record[entry_name]
    dictVars = get_field_from_records(variant_records, id_field=var_id, field_name=var_value)
    return(dictVars)

# Given a list of gene records, return nested dict of gene -> {variants} (dict of variant -> value)
# NOTE: genes without any associated variants will return an empty dict
def get_variants_from_genes(gene_records, gene_id="entrez_symbol", var_id="id", var_value="name"):
    # Sanity check that provided gene_id is valid
    # Valid identifiers are supplied via the civic.yml file
    dict_ids = get_identifier_type_to_civic()
    if gene_id not in dict_ids.keys():
        print("Error! Provided gene_id '%s' is unknown. Please use one of the following: %s." %(gene_id, list(dict_ids.keys())))
        sys.exit(1)
    # dict_ids contains mapping of identifier_type to name of corresponding CIVIC gene record entry
    civic_name = dict_ids[gene_id]
    dictGenes = {}
    # Iterate individual gene records and retrieve the associated variants
    for gene_record in gene_records:
        if civic_name not in gene_record.keys():
            print("Error! Provided gene_id '%s' could not be found in the following gene record: %s" %(gene_id,gene_record))
            sys.exit(1)
        # Use provided gene_id as keys of the output dict
        geneId = gene_record[civic_name]
        if geneId not in dictGenes.keys():
            # NOTE: genes without any associated variants will return an empty dict
            dictGenes[geneId] = get_variants_from_gene(gene_record, var_id=var_id, var_value=var_value)
        else:
            print("Warning! Skipped record for duplicated gene '%s'!" %(geneId))
    return(dictGenes)


def get_all_variant_ids(gene_records):
    # Extract nested dict of gene -> {variants} (where keys are variant ids associated to the gene)
    # NOTE: in this case, the identifier_type used for the gene keys is irrelevant
    dictGenes = get_variants_from_genes(gene_records=gene_records, var_id="id")

    allvariants = []        # keep track of unique variant civic ids
    emptyGenes = []         # keep track of gene ids that did not have any variants associated
    # Iterate individual gene entries in the dict
    for gene in dictGenes.keys():
        # Extract the variants associated to the current gene
        variants = dictGenes[gene]
        # Keep track of gene ids without any associated variant 
        if not variants:
            # Gene keys in the dict will be unique by definition
            emptyGenes.append(gene)
        # Keep track of unique variant civic ids across the complete dict
        for varId in variants:
            if varId not in allvariants:
                allvariants.append(varId)
    print("Parsed a total of %s genes and found %s associated variant ids." %(len(dictGenes.keys()),len(allvariants)))
    if emptyGenes:
        print("Found %s genes that did not have any variants available: %s" %(len(emptyGenes),",".join(emptyGenes)))
    return(allvariants)


def compare_query_and_return(genes, gene_records, identifier_type="entrez_symbol"):
    print("Comparing query and results...")
    print("Validating input genes...")
    validated = validate_genes(genes)
    print("Parsing results...")
    # Return dict of gene -> {variants} (for all genes contained in the provided records)
    # Genes without any associated variants will also be returned
    dictGenes = get_variants_from_genes(gene_records=gene_records, gene_id=identifier_type)
    # Extract all genes which could be matched and retrieved by the query
    retrieved = set(dictGenes.keys())
    # Compare set of input genes to set of retrieved genes
    unmatched = list(set(validated) - retrieved)
    print("Number of input genes (validated):\t%s" %(len(validated)))
    print("Number of retrieved genes:\t%s" %(len(retrieved)))
    if unmatched:
        print("Unmatched genes: %s" %(",".join(unmatched)))
    return None



varMap = {}
for dv in results:
    # Skip variants which do not have any clinical data associated to them
    if not dv['evidence_items']:
#     if not (dv['evidence_items'] or dv['assertions']):
        continue
    # Iterate through the evidence items and store relevant information
    for item in dv['evidence_items']:
        # Sanity check that all critical elements are present and non-empty
        if not (dv['entrez_name'] and dv['name'] and item['evidence_type'] and item['disease']['name']):
            continue
        # Skip records that are not accepted evidence
        if item['status'].upper() != 'ACCEPTED':
            continue
        # Skip records that correspond to germline variants
        # The variant_origin field might be blank/empty
        if item['variant_origin']:
            # TODO: should use re.search instead? To allow anywhere in string
            if re.match('GERMLINE', item['variant_origin'].upper()):
                continue

        # Gene names in CIVIC are HUGO symbols (uppercase) but do a sanity check nevertheless
        gene = dv['entrez_name'].upper()
        variant = dv['name'].upper()
        evidenceType = item['evidence_type'].upper()
        cancerType = item['disease']['name'].upper()
        # Sanity check for empty evidence direction, clinical significance or level
        # 'NULL' is introduced to distinguish from 'N/A' tag
        if item['evidence_direction'] is None:
            item['evidence_direction'] = 'NULL'
        if item['clinical_significance'] is None:
            item['clinical_significance'] = 'NULL'
        if item['evidence_level'] is None:
            item['evidence_level'] = 'NULL'
        # Combine the direction and significance of the evidence in one term
        evidence = item['evidence_direction'].upper() + ':' + item['clinical_significance'].upper()
        level = item['evidence_level'].upper()
        if gene not in varMap.keys():
            varMap[gene] = {}
        # Variant name should be unique within gene
        # (found some duplicates but all were submitted, not accepted data)
        if variant not in varMap[gene].keys():
            varMap[gene][variant] = {}
            # Internal CIViC ID
            varMap[gene][variant]['id'] = dv['id']
            # Score to assess the accumulation of evidence for each variant (quantity and quality)
            # Sanity check for empty scores
            if dv['civic_actionability_score'] is not None:
                varMap[gene][variant]['civic_score'] = dv['civic_actionability_score']
            else:
                varMap[gene][variant]['civic_score'] = 'NULL'

            ## Generate list of strings that will be used to match our input variants in CIVIC
            ## Returned list always has at least length=1 (in this case, containing only variant name)
            ## For CNV, variant matching is not based on HGVS, so matchStrings will only contain the variant name
            matchStrings = generate_civic_matchStrings(variant, dv['hgvs_expressions'], dataType)
            varMap[gene][variant]['match_hgvs'] = matchStrings
            # Keep original HGVS annotations (empty list when nothing is available)
            # Use uppercase to avoid mismatches due to case
            varMap[gene][variant]['hgvs'] = [h.upper() for h in dv['hgvs_expressions']]

            # Include associated variant types (sequence ontology terms). There can be multiple terms
            varMap[gene][variant]['types'] = []
            for vartype in dv['variant_types']:
                varMap[gene][variant]['types'].append(vartype['name'].upper())
            # Account for empty variant types (can happen)
            # 'NULL' is introduced to distinguish from 'N/A' tag
            if not varMap[gene][variant]['types']:
                varMap[gene][variant]['types'] = ['NULL']
        # TODO: there is no sanity check for detecting possible variant name duplicates
        if evidenceType not in varMap[gene][variant].keys():
            varMap[gene][variant][evidenceType] = {}
        if cancerType not in varMap[gene][variant][evidenceType].keys():
            varMap[gene][variant][evidenceType][cancerType] = {}
        if item['drugs']:
            drugs = [d['name'].upper() for d in item['drugs']]
            # When more than 1 drug are listed for the same evidence item, 'drug_interaction_type' is not null and defines the nature of this multiple drug entry
            if item['drug_interaction_type'] is not None:
                # 'Substitutes' indicates that drugs can be considered individually
                if item['drug_interaction_type'].upper() != 'SUBSTITUTES':
                    # Remaining terms ('Sequential' and 'Combination') indicate that drugs should be considered together, so join their names into a single tag
                    # Sort drugs alphabetically to ensure that their order in the combination treatment is always the same
                    drugs.sort()
                    drugs = ['+'.join(drugs)]
#                     # Consider all possible permutations of the drug list
#                     drugs = ['+'.join(per) for per in permutations(drugs)]
#                     drugMatch = None
#                     for drug in drugs:
#                         if drug in varMap[gene][variant][evidenceType][cancerType].keys():
#                             drugMatch = drug
#                             break
#                     # If drug combination is new, get the first permutation
#                     if drugMatch is None:
#                         drugs = [drugs[0]]
#                     else:
#                         drugs = [drugMatch]
        else:
            # Only non-Predictive evidences and Predictive ones without drugs will have this dummy level
            # Introduced for consistency purposes within the varMap structure
            drugs = ['NULL']

        # Iterate through drugs to add evidences associated to them
        #   For non-Predictive evidences or Predictive with empty drugs, drugs=['NULL']
        #   For Predictive and interaction=None, len(drugs) = 1
        #   For Predictive and interaction='Substitutes', len(drugs)>1
        #   For Predictive and interaction!='Substitutes', len(drugs)=1 (combiantion of several using '+')
        for drug in drugs:
            if drug not in varMap[gene][variant][evidenceType][cancerType].keys():
                varMap[gene][variant][evidenceType][cancerType][drug] = {}
            if evidence not in varMap[gene][variant][evidenceType][cancerType][drug].keys():
                varMap[gene][variant][evidenceType][cancerType][drug][evidence] = {}
            if level not in varMap[gene][variant][evidenceType][cancerType][drug][evidence].keys():
                varMap[gene][variant][evidenceType][cancerType][drug][evidence][level] = []
            # Group all publications associated to the same level. Do not check publication status
            ## On 25.01.2019, source structure was changed to introduce ASCO abstracts as a source type
            ## TODO: sanity check for empty ID. Check for type of source?
            varMap[gene][variant][evidenceType][cancerType][drug][evidence][level].append(item['source']['citation_id'])
#             varMap[gene][variant][evidenceType][cancerType][drug][evidence][level].append(item['source']['pubmed_id'])

    # TODO: iterate through assertions and repeat above process
    # for item in dv['assertions']:


### Iterate through input table once more and create output table
### Only new format of input table is allowed!

# Keep track of all matches and non-matches
exactMatches = 0           # Tier 1
posMatches = 0             # Tier 2
noMatches = 0              # Tier 3
notFound = 0               # Tier 4
unmatchedVars = []         # Rows classified as Tier 3 (whole annotation field)
genesNotFound = []         # Genes classified as Tier 4


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
                        print('Warning! No annotations available for column {} and gene {}'.format(colIndx,geneToWrite))
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

outfile.close()

