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
    f = "data/api_access.yml"
    api = yaml.safe_load(f)

    # Sanity checks that expected entries are contained in the yml file
    if "url" not in api.keys():
        print("Error! The URL for the query needs to be provided via the 'url' entry in api_access.yml!")
        sys.exit(1)
    if "endpoints" not in api.keys():
        print("Error! The available endpoints for the query need to be provided via the 'endpoints' entry in api_access.yml!")
        sys.exit(1)

    url = api["url"]
    endpoints = api["endpoints"]

    if not url:
        print("Error! Please provide a URL for the query via the 'url' entry in api_access.yml!")
        sys.exit(1)
    if not endpoints:
        print("Error! Please provide the available endpoints for the query via the 'endpoints' entry in api_access.yml!")
        sys.exit(1)

    if "genes" not in endpoints.keys():
        print("Error! Endpoint for querying genes needs to be provided via the 'genes' entry of 'endpoints' in api_access.yml!")
        sys.exit(1)
    if "variants" not in endpoints.keys():
        print("Error! Endpoint for querying variants needs to be provided via the 'variants' entry of 'endpoints' in api_access.yml!")
        sys.exit(1)

    end_genes = endpoints["genes"]
    end_variants = endpoints["variants"]

    if not end_genes:
        print("Error! Please provide endpoint for querying genes via the 'genes' entry of 'endpoints' in api_access.yml!")
        sys.exit(1)
    if not end_variants:
        print("Error! Please provide endpoint for querying variants via the 'variants' entry of 'endpoints' in api_access.yml!")
        sys.exit(1)

    return(url,end_genes,end_variants)


def get_dict_aminoacids():
    f = "data/data.yml"
    data = yaml.safe_load(f)

    # Sanity check that expected entry is contained in the yml file
    if "aminoacids" not in data.keys():
        print("Error! Please provide a dictionary of one-letter to three-letter aminoacid codes via the 'aminoacids' entry in api_access.yml!")
        sys.exit(1)

    # TODO add sanity checks for all aa being present?
    # TODO add sanity check for contained values eg. special cases like "*"
    dict_codes = data["aminoacids"]
    return(dict_codes)


# TODO: add more options from json.dump?
def write_query(query, outfile, indent=4):
    with open(outfile, 'w') as f:
        json.dump(query, f, ensure_ascii=False, indent=indent)
    return None

def write_civic_results(data, header, outfile):
    outFile = open(args.outfile,'w')
    outHeader = header
    outHeader += "\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\tCIViC_Drug_Support\tCIViC_Predictive\tCIViC_Diagnostic\tCIViC_Prognostic\tCIViC_Predisposing"
    outFile.write(outHeader + "\n")

    # FIXME
