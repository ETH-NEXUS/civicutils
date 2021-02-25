#!/usr/bin/env python

'''
Parse complete CIVIC cache file
Retrieve available drugs and diseases
Lourdes Rosano, Feb 2021
'''

import sys
import argparse
from civicpy import civic


'''
Script
'''

parser = argparse.ArgumentParser(description='Parse CIVIC cache file and extract info about drugs and diseases available.')
# parser.add_argument('--outTag', dest='outTag', required=True, help='Name tag of the output files.')

args = parser.parse_args()

# Load CIVIC cache file
civic.load_cache(on_stale='ignore')

# Retrieve all gene records available in the cache file
all_results = civic.get_all_genes()

drugs = []
diseases = []
all_status = []
nEmpty = 0
for gene_record in all_results:
    gene_variants = gene_record.variants
    for variant_record in gene_variants:
        evidence_items = variant_record.evidence_items
        for evidence_record in evidence_items:
            disease_name = evidence_record.disease.name.strip().upper()
            source_status = evidence_record.source.status
            if (source_status is None) or (not source_status):
                nEmpty += 1
            else:
                source_status = source_status.strip().upper()
                if source_status not in all_status:
                    all_status.append(source_status)
            if disease_name not in diseases:
                diseases.append(disease_name)
            evidence_drugs = evidence_record.drugs
            for evidence_drug in evidence_drugs:
                drug_name = evidence_drug.name.strip().upper()
                if drug_name not in drugs:
                    drugs.append(drug_name)

print("Total # CIVIC drugs: %s" %(len(drugs)))
print("Total # CIVIC diseases: %s" %(len(diseases)))
print("All source status:")
print(all_status)
print("Total # empty source status: %s" %(str(nEmpty)))

drugs.sort()
diseases.sort()

# outfile_drugs = open(args.outTag + '.drugs.txt', 'w')
# for drug in drugs:
#     outfile_drugs.write(drug + "\n")
# outfile_drugs.close()

# outfile_diseases = open(args.outTag + '.diseases.txt', 'w')
# for disease in diseases:
#     outfile_diseases.write(disease + "\n")
# outfile_diseases.close()
