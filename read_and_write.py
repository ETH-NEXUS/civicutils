#!/usr/bin/env python

'''
Functions for I/O operations
Lourdes Rosano, November 2020
'''

'''
Required Python modules
'''

import sys
import yaml
import json

'''
Required files
'''

# TODO: how to specify? or can they be imported as modules are instead?
# data.yml

'''
Functions
'''

def get_api_access_info():
    f = "data/civic.yml"
    url_entry = "url"
    end_entry = "endpoints"
    gene_entry = "genes"
    var_entry = "variants"
    data = yaml.safe_load(f)

    # Sanity checks that expected entries are contained in the yml file
    if url_entry not in data.keys():
        print("Error! The URL for the query needs to be provided via the '%s' entry in %s!" %(url_entry,f))
        sys.exit(1)
    if end_entry not in data.keys():
        print("Error! The available endpoints for the query need to be provided via the '%s' entry in %s!" %(end_entry,f))
        sys.exit(1)

    url = data[url_entry]
    endpoints = data[end_entry]

    if not url:
        print("Error! Please provide a URL for the query via the '%s' entry in %s!" %(url_rentry,f))
        sys.exit(1)
    if not endpoints:
        print("Error! Please provide the available endpoints for the query via the '%s' entry in %s!" %(end_entry,f))
        sys.exit(1)

    if gene_entry not in endpoints.keys():
        print("Error! Endpoint for querying genes needs to be provided via the '%s' entry of '%s' in %s!" %(gene_entry,end_entry,f))
        sys.exit(1)
    if var_entry not in endpoints.keys():
        print("Error! Endpoint for querying variants needs to be provided via the '%s' entry of '%s' in %s!" %(var_entry,end_entry,f))
        sys.exit(1)

    end_genes = endpoints[gene_entry]
    end_variants = endpoints[var_entry]

    if not end_genes:
        print("Error! Please provide endpoint for querying genes via the '%s' entry of '%s' in %s!" %(gene_entry,end_entry,f))
        sys.exit(1)
    if not end_variants:
        print("Error! Please provide endpoint for querying variants via the '%s' entry of '%s' in %s!" %(var_entry,end_entry,f))
        sys.exit(1)

    return(url,end_genes,end_variants)


def get_identifier_type_to_civic():
    f = "data/civic.yml"
    entry_name = "identifier_types"
    data = yaml.safe_load(f)

    # Sanity check that expected entry is contained in the yml file
    if entry_name not in data.keys():
        print("Error! The set of available identifier types in CIVIC need to be provided via the '%s' entry in %s!" %(entry_name,f))
        sys.exit(1)

    dict_ids = data[entry_name]
    # Sanity check that entry is not empty
    if not dict_ids:
        print("Error! Entry '%s' provided in %s is empty! Please list the set of available identifier types in CIVIC as separate entries under '%s'!" %(entry_name,f,entry_name))
        sys.exit(1)
    # Sanity check that provided identifier_type entries have non-empty values
    for k,v in dict_ids.items():
        if not v:
            print("Error! Entry '%s' of '%s' in %s is empty! Please provide a valid mapping from identifier_type to corresponding field name in CIVIC gene record." %(k,entry_name,f))
            sys.exit(1)

    return(dict_ids)


def get_dict_aminoacids():
    f = "data/data.yml"
    entry_name = "aminoacids"
    data = yaml.safe_load(f)

    # Sanity check that expected entry is contained in the yml file
    if entry_name not in data.keys():
        print("Error! Please provide a dictionary of one-letter to three-letter aminoacid codes via the '%s' entry in %s!" %(entry_name,f))
        sys.exit(1)

    # TODO add sanity checks for all aa being present?
    # TODO add sanity check for contained values eg. special cases like "*"
    dict_codes = data[entry_name]
    return(dict_codes)


# Retrieve a given column name from list of header fields
# Throw a warning when required column is not found
def checkHeaderField(name,headerSplit,isRequired=True):
    pos = None
    if name in headerSplit:
        pos = headerSplit.index(name)
    else:
        if isRequired:
            print("Error! Required column '%s' could not be found in header %s" %(name,"\t".join(headerSplit)))
            sys.exit(1)
    return pos

# Retrieve the following column names from the given header:
# - Variant_dna: HGVS c. annotation for the variant (one per row)
# - Variant_prot: HGVS p. annotation (if available) for the variant (one per row)
# - Gene: gene affected by the variant (one gene per row, so variants affecting multiple genes will be reported as separate lines)
# - Variant_impact: optional. Impact of the variant.
# - Variant_exon: optional. Exon or intron of the variant.
def processSnvHeader(header):
    headerSplit = header.strip().split('\t')
    cPos = checkHeaderField("Variant_dna",headerSplit,isRequired=True)
    pPos = checkHeaderField("Variant_prot",headerSplit,isRequired=True)
    genePos = checkHeaderField("Gene",headerSplit,isRequired=True)
    impactPos = checkHeaderField("Variant_impact",headerSplit,isRequired=False)
    exonPos = checkHeaderField("Variant_exon",headerSplit,isRequired=False)
    return (cPos,pPos,genePos,impactPos,exonPos)



# Assumes header and that relevant info is contained in the following columns: Variant_dna, Variant_prot, Gene. Optional columns: Variant_impact, Variant_exon
# For further info, see docs of processSnvHeader()
def readInSnvs(infile):
    # dict lineNumber -> [dna,prot,gene,impact,exon]
    rawData = {}
    inFile = open(infile,'r')
    header = inFile.readline().strip()
    (cPos,pPos,genePos,impactPos,exonPos) = processSnvHeader(header)
    for nLine,line in enumerate(inFile):
        lineSplit = line.strip().split("\t")
# TODO: sanity check of strings beginning with "c." or "p."
        cVar = lineSplit[cPos].strip()
        pVar = lineSplit[pPos].strip()
        gene = lineSplit[genePos].strip()
        impact = ""
        if impactPos:
            impact = lineSplit[impactPos].strip()
        exon = ""
        if exonPos:
            exon = lineSplit[exonPos].strip()
        rawData[nLine] = [cVar,pVar,gene,impact,exon]
    infile.close()
    return(rawData)

# Process rawData to have gene-centered dict
# Returns dict of gene -> [var1,var2,..,varN], where a given var="lineNumber|dna|prot|impact|exon"
def processSnvData(rawData):
    snvData = {}
    for nLine in rawData.keys():
        # lineNumber -> [dna,prot,gene,impact,exon]
        data = rawData[nLine]
        # Retrieve and remove gene from data list
        gene = data.pop(2)
        if gene not in snvData.keys():
            snvData[gene] = []
        # Collapse variant info separated with "|"
        variant = "|".join(data)
        # Keep track of what line each variant comes from
        mapped_var = nLine + "|" + variant
        # This way, we ensure all variants will be unique in dict snvData
        snvData[gene].append(mapped_var)
    return(snvData)


def match_snv_variants_in_civic(snvData, identifier_type="entrez_symbol", batch=True, batch_size=200):
    matchData = {}
# TODO: Sanity check that snvData has expected format (also, this serves as check that input data corresponds to snv)
    if not snvData:
        print("Error! Dictionary provided as 'snvData' is empty!")
        sys.exit(1)
    genes = list(snvData.keys())
    gene_results = query_civic_genes(genes, identifier_type)
    # Compare query to results (ie. informative messages to user)
    compare_query_and_return(genes, gene_results, identifier_type=identifier_type, batch=batch, batch_size=batch_size)
    # Additional sanity check when no results were returned by the genes query
    # Return empty dict in this case
    if not gene_results:
        print("Warning! No gene records were returned by the query!")
        return (matchData)
    allvariants = get_all_variant_ids(gene_results)
    var_results = query_civic_variants(allvariants)
    # Sanity check when no results were returned by the variants query
    # Return empty dict in this case
    if not var_results:
        print("Warning! No variant records were returned by the query!")
        return (matchData)

# TODO

    return matchData





# TODO: add more options from json.dump?
def write_query_to_json(query, outfile, indent=4):
    with open(outfile, 'w') as f:
        json.dump(query, f, ensure_ascii=False, indent=indent)
    return None

def write_civic_results(data, header, outfile):
#     outFile = open(args.outfile,'w')
#     outHeader = header
#     outHeader += "\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\tCIViC_Drug_Support\tCIViC_Predictive\tCIViC_Diagnostic\tCIViC_Prognostic\tCIViC_Predisposing"
#     outFile.write(outHeader + "\n")

# TODO

    return None
