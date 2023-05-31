#!/usr/bin/env python

'''
Retrieve list of available disease names in CIViC cache file
Lourdes Rosano, Feb 2021
'''

import sys
import argparse
from civicpy import civic


'''
Script
'''

parser = argparse.ArgumentParser(description="Parse CIViC cache file and extract list of available disease names.")
parser.add_argument("--outfile", dest="outfile", required=True, help="Output file listing disease names available in CIViC cache file.")

args = parser.parse_args()

# Load offline CIViC cache file
civic.load_cache(on_stale="ignore")

# Retrieve all gene records available in the cache file
all_results = civic.get_all_genes()

diseases = []
# Iterate all available gene records and retrieve all available disease names
for gene_record in all_results:
    gene_variants = gene_record.variants
    for variant_record in gene_variants:
        molecular_profiles = variant_record.molecular_profiles
        for molecular_profile_record in molecular_profiles:
            evidence_items = molecular_profile_record.evidence_items
            for evidence_record in evidence_items:
                if isinstance(evidence_record.disease, civic.Disease):
                    # Use uppercase for consistency of disease names
                    disease_name = evidence_record.disease.name.strip().upper()
                    if disease_name not in diseases:
                        diseases.append(disease_name)

print("Total # CIViC diseases: %s" %(len(diseases)))
# Sort alphabetically
diseases.sort()

# Write output list
outfile = open(args.outfile, "w")
for disease in diseases:
    outfile.write(disease + "\n")
outfile.close()

