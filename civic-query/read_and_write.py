import sys
import yaml
import json

# TODO: how to specify? or can they be imported as modules are instead?
# data.yml


# FIXME: should this go in init.py?
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


# TODO
# FIXME: allow several variants per row?

# Assumes header and that relevant info is contained in the following columns: Variant_dna, Variant_prot, Gene. Optional columns: Variant_impact, Variant_exon
# For further info, see docs of processSnvHeader()
def readInSnvs(infile):
    # dict lineNumber -> [dna,prot,gene,impact,exon]
    rawData = {}
    snvData = {}
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

        # Process rawData to have gene-centered dict
        # Returns dict of gene -> [var1,var2,..,varN], where a given var="lineNumber|dna|prot|impact|exon"
        if gene not in snvData.keys():
            snvData[gene] = {}
        # Collapse variant info separated with "|"
        # Keep track of what line each variant comes from
        variant = nLine + "|" + cVar + "|" + pVar + "|" + impact + "|" + exon
# FIXME
        if variant in snvData[gene].keys():
            print("Found duplicated variant '%s|%s' for gene '%s' in line '%s'!" %(cVar,pVar,gene,nLine))
            sys.exit(1)
        # This way, we ensure all variants will be unique in dict snvData
        snvData[gene][variant] = {}

    infile.close()
    return (rawData,snvData)


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
