srun --cpus-per-task 2 --mem-per-cpu 4000 --time 04:00:00 --pty bash
export CIVICPY_CACHE_FILE=/cluster/customapps/biomed/nexus/sharedutils/db_cache/civicDB/cache.pkl
# run an interactive node
python
 # activate python 3.8.10

from civicpy import civic

civic.load_cache("/cluster/customapps/biomed/nexus/sharedutils/db_cache/civicDB/cache.pkl")

import sys
import os

sys.path.append("/cluster/work/nexus/antoine/Projects/2023_05_Molecular_Profil_CIViCutils/civicutils/civicutils")
# sys.path.append("/Users/hanns/Documents/Work/2023_05_CIViCutils_packages/civicutils/civicutils")

# Import relevant functions
from read_and_write import read_in_snvs,get_dict_support,write_match, write_header_line, write_evidences
from query import query_civic
from filtering import filter_civic
from match import match_in_civic,annotate_ct,filter_ct,process_therapy_support


# Read in file of input SNV variants
(raw_data, snv_data, extra_header) = read_in_snvs("/cluster/work/nexus/antoine/Projects/2023_05_Molecular_Profil_CIViCutils/civicutils/civicutils/data/example_snv.txt")

(rawData,snvData,extraHeader) = read_in_snvs("/Users/hanns/Documents/Work/2023_05_CIViCutils_packages/civicutils/civicutils/data/example_snv.txt")
(rawData,snvData,extraHeader) = read_in_snvs("/Users/hanns/Documents/Work/seiler_civicquery_Nov2021/seiler_civic_input.tsv")


# Query input genes in CIVIC
identifier_type="entrez_symbol"
genes = list(snvData.keys())
all_results = civic.get_all_genes()

# Iterate individual records and retrieve only those matching the provided gene ids
results = []

for gene_record in all_results:
    toKeep = False
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

varMap = {}

gene_record = results[4]

gene_civic_id = str(gene_record.id)
gene_id = str(gene_record.entrez_id)
gene_symbol = gene_record.name.strip().upper()
gene_variants = gene_record.variants

gene_key = gene_symbol
varMap[gene_key] = {}

variant_record = gene_variants[160]

variant_id = str(variant_record.id)
variant_name = variant_record.name.strip().upper()
hgvs_expressions = variant_record.hgvs_expressions
molecular_profiles = variant_record.molecular_profiles


varMap[gene_key][variant_id] = {}
varMap[gene_key][variant_id]['name'] = variant_name
varMap[gene_key][variant_id]['hgvs'] = [h.strip().upper() for h in hgvs_expressions]
varMap[gene_key][variant_id]['types'] = []

for vartype_record in variant_record.variant_types:
        varMap[gene_key][variant_id]['types'].append(vartype_record.name.strip().upper())

molecular_profile = molecular_profiles[2]

molecular_profile_name = molecular_profile.name.strip().upper()
molecular_profile_id = str(molecular_profile.id)
civic_score = molecular_profile.molecular_profile_score
if (civic_score is None) or (not civic_score):
        civic_score = 0.0

evidence_items = molecular_profile.evidence_items
n_items = len(evidence_items)

varMap[gene_key][variant_id][molecular_profile_id] = {}
varMap[gene_key][variant_id][molecular_profile_id]['name'] = molecular_profile_name
varMap[gene_key][variant_id][molecular_profile_id]['civic_score'] = civic_score

varMap[gene_key][variant_id][molecular_profile_id]['n_evidence_items'] = n_items
varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'] = {}

evidence_record = evidence_items[0]

evidence_type = evidence_record.evidence_type.strip().upper()
if hasattr(evidence_record.disease, 'name'):
        disease = evidence_record.disease.name.strip().upper()
else:
        disease = "NULL"

evidence_status = evidence_record.status.strip().upper()
source_type = evidence_record.source.source_type.strip().upper()
source_status = evidence_record.status.strip().upper()
evidence_id = str(evidence_record.source.citation_id).strip()        # expected to be entirely numeric

if evidence_type not in varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'].keys():
        varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type] = {}

if disease not in varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type].keys():
        varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type][disease] = {}

def check_empty_field(inField):
    """
    Check if the given field is None or empty, and return 'NULL' in this case; otherwise, return the given field.
    :param inField:    Object to check.
    :return:           Either 'NULL' (if object is empty or None) or the input object.
    """
    if (inField is None) or (not inField):
        newField = "NULL"
    else:
        newField = inField
    return newField


evidence_level = check_empty_field(evidence_record.evidence_level).strip().upper()                     # string (just in case, should never be None)
variant_origin = check_empty_field(evidence_record.variant_origin).strip().upper()                     # string, can be None
evidence_direction = check_empty_field(evidence_record.evidence_direction).strip().upper()
clinical_significance = check_empty_field(evidence_record.significance).strip().upper()       # string, can be None
evidence_rating = check_empty_field(evidence_record.rating)      

evidence = evidence_direction + ':' + clinical_significance

therapies = []
evidence_therapies = evidence_record.therapies
for evidence_therapy in evidence_therapies:
        therapy_name = evidence_therapy.name.strip().upper()
        if therapy_name not in therapies:
                therapies.append(therapy_name)

therapy_interaction = evidence_record.therapy_interaction_type
if therapy_interaction is not None:
        therapy_interaction = therapy_interaction.strip().upper()

if therapy_interaction != "SUBSTITUTES":
        therapies.sort()
        therapies = ["+".join(therapies)]

if not therapies:
        therapies = ["NULL"]

if (evidence_type != "PREDICTIVE") and (therapies != ["NULL"]):
        raise ValueError("Only evidences of type 'PREDICTIVE' can have therapies associated!")

for therapy in therapies:
        if therapy not in varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type][disease].keys():
                varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type][disease][therapy] = {}
        if evidence not in varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type][disease][therapy].keys():
                varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type][disease][therapy][evidence] = {}
        if evidence_level not in varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type][disease][therapy][evidence].keys():
                varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type][disease][therapy][evidence][evidence_level] = []
        varMap[gene_key][variant_id][molecular_profile_id]['evidence_items'][evidence_type][disease][therapy][evidence][evidence_level].append(source_type + "_" + str(evidence_id) + ":" + evidence_status + ":" + source_status + ":" + variant_origin + ":" + str(evidence_rating))


# at this point, theoricaly query_civic function should work now 
##### # Filter undesired evidences to avoid matching later on 

gene_id_in=[]
gene_id_not_in=[]
min_variants=0
var_id_in=[]
var_id_not_in=[]
var_name_in=[]
var_name_not_in=[]
min_civic_score=0
var_type_in=[]
var_type_not_in=[]
min_evidence_items=0
evidence_type_in=[]
evidence_type_not_in=[]
disease_in=[]
disease_not_in=[]
therapy_name_in=[]
therapy_name_not_in=[]
evidence_dir_in=[]
evidence_dir_not_in=[]
evidence_clinsig_in=[]
evidence_clinsig_not_in=[]
evidence_level_in=[]
evidence_level_not_in=[]
evidence_status_in = ['ACCEPTED']
evidence_status_not_in=[]
source_status_in=[]
source_status_not_in=[]
var_origin_in=[]
var_origin_not_in = ['GERMLINE']
source_type_in=[]
source_type_not_in=[]
min_evidence_rating=0
output_empty=False



varmap_entries_variant = ['name','hgvs','types']
varmap_entries_molecular_profile = ['name','civic_score','n_evidence_items','evidence_items']
sorted_cts = ["ct","gt","nct"]

# copy past all the function from utils.py + function filter_in , filter_not_in and filter_cutoff
from utils import check_string_filter_arguments,check_cutoff_filter_arguments,check_is_bool,check_keys,check_keys_not,check_is_none


check_is_none(output_empty,"output_empty")
check_is_bool(output_empty,"output_empty")

cleanMap = {}
gene_id = "TP53"
gene_id_str = str(gene_id)

keepGene = filter_in(gene_id_str, "gene_id", gene_id_in, "gene_id_in", matchType="exact")
removeGene = filter_not_in(gene_id_str, "gene_id", gene_id_not_in, "gene_id_not_in", matchType="exact")

n_variants = len(varMap[gene_id].keys())
keepGene = filter_cutoff(n_variants, "n_variants", min_variants, "min_variants")

var_id = "918"
var_id_str = str(var_id)
keepVar = filter_in(var_id_str, "var_id", var_id_in, "var_id_in", matchType="exact")
removeVar = filter_not_in(var_id_str, "var_id", var_id_not_in, "var_id_not_in", matchType="exact")


check_keys(list(varMap[gene_id][var_id].keys()),"varMap",varmap_entries_variant,matches_all=False)

variant = varMap[gene_id][var_id]["name"]
keepVar = filter_in(variant, "variant", var_name_in, "var_name_in", matchType="partial")
removeVar = filter_not_in(variant, "variant", var_name_not_in, "var_name_not_in", matchType="partial")


variant_types = varMap[gene_id][var_id]["types"]
nKeep = 0
nRemove = 0
for var_type in variant_types:
        keepType = filter_in(var_type, "var_type", var_type_in, "var_type_in", matchType="partial")
        if keepType:
                nKeep += 1
        removeType = filter_not_in(var_type, "var_type", var_type_not_in, "var_type_not_in", matchType="partial")
        if removeType:
                nRemove += 1
        if (nKeep == 0):
                continue
        if (nRemove > 0):
                continue
        if output_empty:
                if var_id not in cleanMap[gene_id].keys():
                        cleanMap[gene_id][var_id] = {}
                        cleanMap[gene_id][var_id]["name"] = variant
                        #cleanMap[gene_id][var_id]["civic_score"] = var_score
                        cleanMap[gene_id][var_id]["hgvs"] = varMap[gene_id][var_id]['hgvs']
                        cleanMap[gene_id][var_id]["types"] = variant_types
                        #cleanMap[gene_id][var_id]["n_evidence_items"] = 0
                        cleanMap[gene_id][var_id]["evidence_items"] = {}

molecular_profiles_ids = set(list(varMap[gene_id][var_id].keys())) ^ set(varmap_entries_variant)
molecular_profil_id = "892"

n_evidence_items = varMap[gene_id][var_id][molecular_profil_id]["n_evidence_items"]
keepMP = filter_cutoff(n_evidence_items, "n_evidence_items", min_evidence_items, "min_evidence_items")

var_score = varMap[gene_id][var_id][molecular_profil_id]["civic_score"]
keepMP = filter_cutoff(var_score, "var_score", min_civic_score, "min_civic_score")

n_evidence_items_after = 0
evidence_types = list(varMap[gene_id][var_id][molecular_profil_id]["evidence_items"].keys())
allowed_evidence_types = []
for evidence_type in evidence_types:
        keepType = filter_in(evidence_type, "evidence type", evidence_type_in, "evidence_type_in", matchType="exact")
        if not keepType:
                continue
        removeType = filter_not_in(evidence_type, "evidence type", evidence_type_not_in, "evidence_type_not_in", matchType="exact")
        if removeType:
                continue
        allowed_evidence_types.append(evidence_type)

evidence_type = allowed_evidence_types[0]
type_diseases = list(varMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys())
check_keys_not(type_diseases,"varMap",sorted_cts)

allowed_diseases = []
type_disease = type_diseases[0]
keepDisease = filter_in(type_disease, "disease type", disease_in, "disease_in", matchType="partial")
removeDisease = filter_not_in(type_disease, "disease type", disease_not_in, "disease_not_in", matchType="partial")
allowed_diseases.append(type_disease)

disease = allowed_diseases[0]
disease_therapies = list(varMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease].keys())

allowed_therapies = []
disease_therapy = disease_therapies[0]
keeptherapy = filter_in(disease_therapy, "therapy", therapy_name_in, "therapy_name_in", matchType="partial")
removetherapy = filter_not_in(disease_therapy, "therapy", therapy_name_not_in, "therapy_name_not_in", matchType="partial")
allowed_therapies.append(disease_therapy)

therapy = allowed_therapies[0]
evidences = list(varMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy].keys())
allowed_evidences = []
evidence = evidences[0]
evidenceArr = evidence.strip().split(':')

direction = evidenceArr[0]
keepDir = filter_in(direction, "evidence direction", evidence_dir_in, "evidence_dir_in", matchType="exact")
removeDir = filter_not_in(direction, "evidence direction", evidence_dir_not_in, "evidence_dir_not_in", matchType="exact")

clin_signf = evidenceArr[1]
keepClin = filter_in(clin_signf, "clinical significance", evidence_clinsig_in, "evidence_clinsig_in", matchType="exact")
removeClin = filter_not_in(clin_signf, "clinical significance", evidence_clinsig_not_in, "evidence_clinsig_not_in", matchType="exact")
allowed_evidences.append(evidence)

this_evidence = allowed_evidences[0]
evidence_levels = list(varMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy][this_evidence].keys())
allowed_levels = []
evidence_level = evidence_levels[0]
keepLevel = filter_in(evidence_level, "evidence level", evidence_level_in, "evidence_level_in",  matchType="exact")
removeLevel = filter_not_in(evidence_level, "evidence level", evidence_level_not_in, "evidence_level_not_in", matchType="exact")
allowed_levels.append(evidence_level)


level = allowed_levels[0]
all_evidence_items = varMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy][this_evidence][level]

evidence_item = all_evidence_items[0]
itemArr = evidence_item.strip().split(":")
tmpId = itemArr[0]
evidence_status = itemArr[1]
source_status = itemArr[2]
var_origin = itemArr[3]
rating = itemArr[4]

tmpIdArr = tmpId.split("_")
idType = tmpIdArr[0]
this_id = tmpIdArr[1]

keepStatus = filter_in(evidence_status, "evidence status", evidence_status_in, "evidence_status_in", matchType="exact")
removeStatus = filter_not_in(evidence_status, "evidence status", evidence_status_not_in, "evidence_status_not_in", matchType="exact")
keepSource = filter_in(source_status, "source status", source_status_in, "source_status_in", matchType="partial")
removeSource = filter_not_in(source_status, "source status", source_status_not_in, "source_status_not_in", matchType="partial")
keepOrigin = filter_in(var_origin, "variant origin", var_origin_in, "var_origin_in", matchType="partial")
removeOrigin = filter_not_in(var_origin, "variant origin", var_origin_not_in, "var_origin_not_in", matchType="partial")
keepType = filter_in(idType, "source id type", source_type_in, "source_type_in", matchType="exact")
removeType = filter_not_in(idType, "source id type", source_type_not_in, "source_type_not_in", matchType="exact")

keepRating = filter_cutoff(rating, "rating", min_evidence_rating, "min_evidence_rating")


cleanMap[gene_id] = {}
if var_id not in cleanMap[gene_id].keys():
    cleanMap[gene_id][var_id] = {}
    cleanMap[gene_id][var_id]["name"] = variant
    cleanMap[gene_id][var_id]["hgvs"] = varMap[gene_id][var_id]['hgvs']
    cleanMap[gene_id][var_id]["types"] = variant_types


if molecular_profil_id not in cleanMap[gene_id][var_id].keys():
    cleanMap[gene_id][var_id][molecular_profil_id] = {}
    cleanMap[gene_id][var_id][molecular_profil_id]["civic_score"] = var_score
    cleanMap[gene_id][var_id][molecular_profil_id]["n_evidence_items"] = 0
    cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"] = {}    


                     
if evidence_type not in cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"].keys():
    cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type] = {}


if disease not in cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys():
    cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease] = {}


if therapy not in cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease].keys():
    cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy] = {}


if this_evidence not in cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy].keys():
    cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy][this_evidence] = {}

if level not in cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy][this_evidence].keys():
    cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy][this_evidence][level] = []

cleanMap[gene_id][var_id][molecular_profil_id]["evidence_items"][evidence_type][disease][therapy][this_evidence][level].append(evidence_item)
n_evidence_items_after += 1

# at this point, theoricaly filter_civic function should work now 
# Match input SNV variants in CIVIC, pick highest tier available per input gene+variant
# Tier hierarchy: 1 > 1b > 2 > 3 > 4


# load function from match.py

varMap = cleanMap
dataType="SNV"
identifier_type="entrez_symbol"
select_tier="all"
sorted_tiers = ["tier_1","tier_1b","tier_2","tier_3","tier_4"]

varData= snvData

check_argument(varData,"varData")
check_is_dict(varData,"varData")
check_data_type(dataType)
check_identifier_type(identifier_type)
select_tier = check_tier_selection(select_tier,sorted_tiers)

matchMap = {}
matchedIds = []

gene = "TP53"
matchMap[gene] = {}
variant = "c.818G>T|p.Arg273Leu|Missense||0"

matchMap[gene][variant] = {}
variants = []
impactArr = []
exonArr = []

varArr = variant.split("|")
cVars = varArr[0]
pVars = varArr[1]
impacts = varArr[2]
exons = varArr[3]

cVarArr = parse_input(cVars, "Variant_dna", isRequired=False)
pVarArr = parse_input(pVars, "Variant_prot", isRequired=False)

cVar = "c.818G>T"
pVar = "p.Arg273Leu"
variants.append(cVar)
variants.append(pVar)

impactArr = parse_input(impacts, "Variant_impact", isRequired=False)
exonArr = parse_input(exons, "Variant_exon", isRequired=False)
match = match_variants_in_civic(gene, variants, varMap, dataType, impacts=impactArr, exons=exonArr)

(matchMap,allIds) = add_match(matchMap,gene,variant,match)
varId = "12"
matchedIds.append(varId)


# Annotate matched CIVIC evidences with cancer specificity of the associated diseases
disease_name_not_in = []
disease_name_in = ['bladder']
alt_disease_names = ['solid tumor']


sorted_cts = ["ct","gt","nct"]
varmap_entries_variant = ['name','hgvs','types']
newMap = {}
gene = "BRAF"
newMap[gene] = {}

newMap[gene][variant] = {}
variant = "12" 
newMap[gene][variant] = {}

newMap[gene][variant]['name'] = varMap[gene][variant]['name']
newMap[gene][variant]['hgvs'] = [a for a in varMap[gene][variant]['hgvs']]
newMap[gene][variant]['types'] = [b for b in varMap[gene][variant]['types']]

molecular_profile_ids = set(list(varMap[gene][variant].keys())) ^ set(varmap_entries_variant)
molecular_profile_id = "4173"

newMap[gene][variant][molecular_profil_id] = {}
newMap[gene][variant][molecular_profil_id]['civic_score'] = varMap[gene][variant][molecular_profile_id]['civic_score']
newMap[gene][variant][molecular_profil_id]['n_evidence_items'] = varMap[gene][variant][molecular_profile_id]['n_evidence_items']
newMap[gene][variant][molecular_profil_id]['evidence_items'] = {}

evidence_type = "PREDICTIVE"
newMap[gene][variant][molecular_profile_id]['evidence_items'][evidence_type] = {}
allDiseases = list(varMap[gene][variant][molecular_profile_id]['evidence_items'][evidence_type].keys())

(ctDis,gtDis,nctDis) = classify_diseases(allDiseases, disease_name_not_in, disease_name_in, alt_disease_names)

for ct in sorted_cts:
    if ct == "ct":
        newMap = add_ct(ctDis,ct,gene,variant,molecular_profile_id,evidence_type,newMap,varMap,isAnnot=False)
    if ct == "gt":
        newMap = add_ct(gtDis,ct,gene,variant,molecular_profile_id,evidence_type,newMap,varMap,isAnnot=False)
    if ct == "nct":
        newMap = add_ct(nctDis,ct,gene,variant,molecular_profile_id,evidence_type,newMap,varMap,isAnnot=False)



annotMap = newMap
##### filter_ct 
annotMap = filter_ct(annotMap,select_ct="highest")

supportDict = get_dict_support()

annotMatch = process_therapy_support(matchMap,annotMap,supportDict)

# write output
matchMap=annotMatch
varMap=annotMap
rawMap=rawData
header=extraHeader
dataType="SNV"
outfile="/cluster/work/nexus/antoine/Projects/2023_05_Molecular_Profil_CIViCutils/Test.tsv"
hasSupport=True
hasCt=True
writeCt=False
writeSupport=True
writeComplete=False

sorted_evidence_types = ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']
evidenceType = "PREDICTIVE"
special_cases = ["NON_SNV_MATCH_ONLY","NON_CNV_MATCH_ONLY","NON_EXPR_MATCH_ONLY"]
sorted_cts = ["ct","gt","nct"]



# check check_match_before_writing
varmap_entries = ['name','hgvs','types'] 
special_cases = ["NON_SNV_MATCH_ONLY","NON_CNV_MATCH_ONLY","NON_EXPR_MATCH_ONLY"]
check_arguments([matchMap,rawMap],["matchMap","rawMap"])
check_is_none(hasSupport,"hasSupport")
check_is_none(hasCt,"hasCt")
check_is_none(writeCt,"writeCt")
check_is_none(writeSupport,"writeSupport")
check_is_none(writeComplete,"writeComplete")
check_is_dict(matchMap,"matchMap")
check_is_dict(varMap,"varMap")
check_is_dict(rawMap,"rawMap")
check_is_bool(hasSupport,"hasSupport")
check_is_bool(hasCt,"hasCt")
check_is_bool(writeCt,"writeCt")
check_is_bool(writeSupport,"writeSupport")
check_is_bool(writeComplete,"writeComplete") 


write_match(matchMap=annotMatch, varMap=annotMap, rawMap=rawData, header=extraHeader, outfile="/cluster/work/nexus/antoine/Projects/2023_05_Molecular_Profil_CIViCutils/Test.tsv", dataType="SNV", hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False)



write_match(annotMatch, annotMap, rawData, extraHeader, dataType="SNV", outfile="/Users/hanns/Documents/Work/2023_05_CIViCutils_packages/try_seiler_output_snv.tsv", hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False)

















##### Read CNV
python

import sys
import os

sys.path.append("/Users/hanns/Documents/Work/2023_05_CIViCutils_packages/civicutils/civicutils")

# Import relevant functions
from read_and_write import readInCnvs,get_dict_support,write_match
from query import query_civic
from filtering import filter_civic
from match import match_in_civic,annotate_ct,filter_ct,process_therapy_support

# Read in file of input SNV variants
(rawData,snvData,extraHeader) = readInCnvs("/Users/hanns/Documents/Work/2023_05_CIViCutils_packages/civicutils/civicutils/data/example_cnv.txt")

# Query input genes in CIVIC
varMap = query_civic(list(snvData.keys()), identifier_type="entrez_symbol")
# Filter undesired evidences to avoid matching later on
varMap = filter_civic(varMap, evidence_status_in = ['ACCEPTED'], var_origin_not_in = ['GERMLINE'], output_empty=False)

# Match input SNV variants in CIVIC, pick highest tier available per input gene+variant
# Tier hierarchy: 1 > 1b > 2 > 3 > 4
(matchMap,matchedIds,varMap) = match_in_civic(snvData, dataType="CNV", identifier_type="entrez_symbol", select_tier="highest", varMap=varMap)

# Annotate matched CIVIC evidences with cancer specificity of the associated diseases
disease_name_not_in = []
disease_name_in = ['bladder']
alt_disease_names = ['solid tumor']
annotMap = annotate_ct(varMap, disease_name_not_in, disease_name_in, alt_disease_names)

# Filter CIVIC evidences to pick only those for the highest cancer specificity available
# ct hierarchy: ct > gt > nct
annotMap = filter_ct(annotMap,select_ct="highest")

# Get custom dictionary of support from data.yml (provided within the package)
# This defines how each combination of evidence direction + clinical significance in CIVIC is classified in terms of drug support (eg. sensitivity, resistance, unknown, etc.)
supportDict = get_dict_support()

# Process drug support of the matched variants using the annotated CIVIC evidences
annotMatch = process_therapy_support(matchMap,annotMap,supportDict)

# Write to output
# Do not report the CT classification of each disease, and write column with the overall drug support of the match for each available CT class
write_match(annotMatch, annotMap, rawData, extraHeader, dataType="CNV", outfile="/Users/hanns/Documents/Work/2023_05_CIViCutils_packages/try_example_cnv.tsv", hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False)





##### Read Expr

import sys
import os

sys.path.append("/Users/hanns/Documents/Work/2023_05_CIViCutils_packages/civicutils/civicutils")

# Import relevant functions
from read_and_write import readInExpr,get_dict_support,write_match
from query import query_civic
from filtering import filter_civic
from match import match_in_civic,annotate_ct,filter_ct,process_therapy_support

# Read in file of input SNV variants
(rawData,snvData,extraHeader) = readInExpr("/Users/hanns/Documents/Work/2023_05_CIViCutils_packages/civicutils/civicutils/data/example_expr.txt")

# Query input genes in CIVIC
varMap = query_civic(list(snvData.keys()), identifier_type="entrez_symbol")
# Filter undesired evidences to avoid matching later on
varMap = filter_civic(varMap, evidence_status_in = ['ACCEPTED'], var_origin_not_in = ['GERMLINE'], output_empty=False)

# Match input SNV variants in CIVIC, pick highest tier available per input gene+variant
# Tier hierarchy: 1 > 1b > 2 > 3 > 4
(matchMap,matchedIds,varMap) = match_in_civic(snvData, dataType="EXPR", identifier_type="entrez_symbol", select_tier="highest", varMap=varMap)

# Annotate matched CIVIC evidences with cancer specificity of the associated diseases
disease_name_not_in = []
disease_name_in = ['bladder']
alt_disease_names = ['solid tumor']
annotMap = annotate_ct(varMap, disease_name_not_in, disease_name_in, alt_disease_names)

# Filter CIVIC evidences to pick only those for the highest cancer specificity available
# ct hierarchy: ct > gt > nct
annotMap = filter_ct(annotMap,select_ct="highest")

# Get custom dictionary of support from data.yml (provided within the package)
# This defines how each combination of evidence direction + clinical significance in CIVIC is classified in terms of drug support (eg. sensitivity, resistance, unknown, etc.)
supportDict = get_dict_support()

# Process drug support of the matched variants using the annotated CIVIC evidences
annotMatch = process_therapy_support(matchMap,annotMap,supportDict)

# Write to output
# Do not report the CT classification of each disease, and write column with the overall drug support of the match for each available CT class
write_match(annot_match, annot_map, raw_data, extra_header, data_type="SNV", outfile="/cluster/work/nexus/antoine/Projects/2023_05_Molecular_Profil_CIViCutils/try_example_snv.tsv", has_support=True, has_ct=True, write_ct=False, write_support=True, write_complete=False)

match_map =annot_match
var_map =annot_map
raw_map =raw_data
header =extra_header
data_type="SNV"
outfile="/cluster/work/nexus/antoine/Projects/2023_05_Molecular_Profil_CIViCutils/try_example_snv.tsv"
has_support=True
has_ct=True
write_ct=False
write_support=True
write_complete=False

sorted_evidence_types = ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']
evidence_type = "PREDICTIVE"
special_cases = ["NON_SNV_MATCH_ONLY","NON_CNV_MATCH_ONLY","NON_EXPR_MATCH_ONLY"]
sorted_cts = ["ct","gt","nct"]
varmap_entries_variant = ['name','hgvs','types']
    
from utils import check_match_before_writing,check_keys,check_keys_not,check_data_type,check_dict_entry
check_match_before_writing(match_map, var_map, raw_map, has_support, has_ct, write_ct, write_support, write_complete)
check_data_type(data_type)
outfile = open(outfile,'w')

(out_header, clean_header, write_impact, write_exon) = write_header_line(data_type, header, write_support) 
outfile.write(out_header + "\n")

n_line="0"
line_list = raw_map[n_line]
extra_line = []

gene = line_list[0]
c_var = line_list[1]
p_var = line_list[2]
impact = line_list[3]
exon = line_list[4]
comb_id = c_var + "|" + p_var + "|" + impact + "|" + exon + "|" + str(n_line)

main_line = gene + "\t" + c_var + "\t" + p_var

tier = "tier_3"

gene_scores = []
gene_var_types = []
therapy_support = []
result_map = {}
write_line = False
all_variants = []
therapySupport = match_map[gene][comb_id][tier]["therapy_support"]
for tmp_var in match_map[gene][comb_id][tier]["matched"]:
        all_variants.append(tmp_var)

if write_support:
        for i in therapySupport:
                therapy_support.append(i.upper())

var_id = all_variants[0]
variant = var_map[gene][var_id]["name"]

gene_var_types.append(gene + ":" + variant + ":" + ",".join(var_map[gene][var_id]["types"]))
molecular_profile_ids = set(list(var_map[gene][var_id].keys())) ^ set(varmap_entries_variant)
molecular_profil_id = "1590"
gene_scores.append(gene + ':' + variant + ':' + molecular_profil_id + ':' + str(var_map[gene][var_id][molecular_profil_id]['civic_score']))

for evidence_type in sorted_evidence_types:
        if evidence_type in var_map[gene][var_id][molecular_profil_id]["evidence_items"].keys():
                if evidence_type not in result_map.keys():
                        result_map[evidence_type] = []
                        write_therapy = False
                if evidence_type == evidence_type:
                        write_therapy=True
                if has_ct:
                        check_keys(list(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys()),"var_map",sorted_cts,matches_all=True)
                        for ct in var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys():
                                if write_ct:
                                        results = write_evidences(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type][ct], write_therapy=write_therapy, write_ct=ct, write_complete=write_complete)
                                else:
                                        results = write_evidences(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type][ct], write_therapy=write_therapy, write_ct=None, write_complete=write_complete)
                                for x in results:
                                        result_map[evidence_type].append(gene + ":" + variant + ":" + molecular_profil_id + ':' + x)
                else:
                        check_keys_not(list(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type].keys()), "var_map", sorted_cts)
                        if write_ct:
                                raise ValueError("Option 'write_ct' cannot be selected when 'has_ct'=False!")
                        results = write_evidences(var_map[gene][var_id][molecular_profil_id]["evidence_items"][evidence_type], write_therapy=write_therapy, write_ct=None, write_complete=write_complete)
                        for x in results:
                                result_map[evidence_type].append(gene + ":" + variant + ":" + molecular_profil_id + ':' + x)






