import os
import sys
import yaml
import json

# TODO: how to specify? or can they be imported as modules are instead?
BinPath = os.path.split(os.path.realpath(__file__))[0]


# FIXME: should this go in init.py?
def get_dict_aminoacids():
    f = BinPath + "/data.yml"
    entry_name = "aminoacids"
    with open(f, 'r') as infile:
        try:
            data = yaml.safe_load(infile)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(1)

    # Sanity check that expected entry is contained in the yml file
    if entry_name not in data.keys():
        raise ValueError("Please provide a dictionary of one-letter to three-letter aminoacid codes via the '%s' entry in %s!" %(entry_name,f))

    # TODO add sanity checks for all aa being present?
    # TODO add sanity check for contained values eg. special cases like "*"
    dict_codes = data[entry_name]
    return dict_codes


# FIXME: should this go in init.py?
def get_dict_support():
    f = BinPath + "/data.yml"
    entry_name = "drug_support"
    with open(f, 'r') as infile:
        try:
            data = yaml.safe_load(infile)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(1)

    # Sanity check that expected entry is contained in the yml file
    if entry_name not in data.keys():
        raise ValueError("Please provide a custom dictionary of drug support for CIVIC evidences via the '%s' entry in %s!" %(entry_name,f))

    # TODO add sanity checks for all aa being present?
    # TODO add sanity check for contained values eg. special cases like "*"
    supportDict = data[entry_name]
    return supportDict


# Retrieve a given column name from list of header fields
# Throw a warning when required column is not found
def checkHeaderField(name,headerSplit,isRequired=True):
    pos = None
    if name in headerSplit:
        pos = headerSplit.index(name)
    else:
        if isRequired:
            raise ValueError("Required column '%s' could not be found in header '%s'" %(name," ".join(headerSplit)))
    return pos


# Retrieve the following column names from the given SNV header:
# - Gene: gene affected by the variant (one gene per row, so variants affecting multiple genes will be reported as separate lines)
# - Variant_dna: HGVS c. annotation for the variant (one per row)
# - Variant_prot: HGVS p. annotation (if available) for the variant (one per row)
# - Variant_impact: optional. Impact of the variant.
# - Variant_exon: optional. Exon or intron of the variant.
def processSnvHeader(header):
    headerSplit = header.strip().split('\t')
    genePos = checkHeaderField("Gene",headerSplit,isRequired=True)
    cPos = checkHeaderField("Variant_dna",headerSplit,isRequired=True)
    pPos = checkHeaderField("Variant_prot",headerSplit,isRequired=True)
    impactPos = checkHeaderField("Variant_impact",headerSplit,isRequired=False)
    exonPos = checkHeaderField("Variant_exon",headerSplit,isRequired=False)
    return (genePos,cPos,pPos,impactPos,exonPos)


# TODO
# Assumes header and that relevant info is contained in the following columns: Variant_dna, Variant_prot, Gene. Optional columns: Variant_impact, Variant_exon
# For further info, see docs of processSnvHeader()
def readInSnvs(infile):
    # dict lineNumber -> [gene,dna,prot,(impact,exon)]
    rawData = {}
    # dict gene -> variant (dna|prot|impact|exon|lineNumber) -> null
    snvData = {}
    inFile = open(infile,'r')
    header = inFile.readline().strip()
    (genePos,cPos,pPos,impactPos,exonPos) = processSnvHeader(header)
    for nLine,line in enumerate(inFile):
        lineSplit = line.strip().split("\t")
        cVar = lineSplit[cPos].strip()
        pVar = lineSplit[pPos].strip()
        gene = lineSplit[genePos].strip()
        impact = ""
        if impactPos:
            impact = lineSplit[impactPos].strip()
        exon = ""
        if exonPos:
            exon = lineSplit[exonPos].strip()
        rawData[str(nLine)] = [gene,cVar,pVar,impact,exon]

        # Process rawData to have gene-centered dict
        # Returns dict of gene -> [var1,var2,..,varN], where a given var="dna|prot|impact|exon|lineNumber"
        if gene not in snvData.keys():
            snvData[gene] = {}
        # Collapse variant info separated with "|"
        # Keep track of what line each variant comes from
        variant = cVar + "|" + pVar + "|" + impact + "|" + exon + "|" + str(nLine)
        # NOTE: Variants can never be duplicated because of the different row numbers assigned
        # if variant in snvData[gene].keys():
        #     print("Found duplicated variant '%s|%s' for gene '%s' in line '%s'!" %(cVar,pVar,gene,str(nLine)))
        #     sys.exit(1)
        snvData[gene][variant] = None
    inFile.close()
    return (rawData,snvData)


# Retrieve the following column names from the given CNV header:
# - Gene: gene affected by the variant (one gene per row, so variants affecting multiple genes will be reported in different lines)
# - Variant_cnv: type of CNV variant (should be one of the following terms: 'GAIN', 'DUPLICATION', 'DUP', 'AMPLIFICATION', 'AMP', 'DELETION', 'DEL', 'LOSS')
def processCnvHeader(header):
    headerSplit = header.strip().split('\t')
    genePos = checkHeaderField("Gene",headerSplit,isRequired=True)
    cnvPos = checkHeaderField("Variant_cnv",headerSplit,isRequired=True)
    return (genePos,cnvPos)


# Assumes header and that relevant info is contained in the following columns: Gene, Variant_cnv
# For further info, see docs of processCnvHeader()
# TODO
# FIXME: allow several variants per row?
def readInCnvs(infile):
    # dict lineNumber -> [gene,cnv]
    rawData = {}
    # dict gene -> variant (cnv|lineNumber) -> null
    cnvData = {}
    inFile = open(infile,'r')
    header = inFile.readline().strip()
    (genePos,cnvPos) = processCnvHeader(header)
    for nLine,line in enumerate(inFile):
        lineSplit = line.strip().split("\t")
        gene = lineSplit[genePos].strip()
        cnv = lineSplit[cnvPos].strip()
        rawData[str(nLine)] = [gene,cnv]
        # Process rawData to have gene-centered dict
        # Returns dict of gene -> [var1,var2,..,varN], where a given var="cnv|lineNumber"
        if gene not in cnvData.keys():
            cnvData[gene] = {}
        # Collapse variant info separated with "|"
        # Keep track of what line each variant comes from
        variant = cnv + "|" + str(nLine)
        # NOTE: Variants can never be duplicated because of the different row numbers assigned
        # if variant in cnvData[gene].keys():
        #     print("Found duplicated variant '%s' for gene '%s' in line '%s'!" %(cnv,gene,str(nLine)))
        #     sys.exit(1)
        cnvData[gene][variant] = None
    inFile.close()
    return (rawData,cnvData)


# Retrieve the following column names from the given EXPR header:
# - Gene: gene affected by the expression change (one gene per row)
# - logFC: log fold-change of the gene expression (will be translated into 'OVEREXPRESSION' or 'UNDEREXPRESSION')
def processExprHeader(header):
    headerSplit = header.strip().split('\t')
    genePos = checkHeaderField("Gene",headerSplit,isRequired=True)
    logfcPos = checkHeaderField("logFC",headerSplit,isRequired=True)
    return (genePos,logfcPos)


# Assumes header and that relevant info is contained in the following columns: Gene, logFC
# For further info, see docs of processExprHeader()
# TODO
# FIXME: allow several variants per row?
def readInExpr(infile):
    # dict lineNumber -> [gene,logFC]
    rawData = {}
    # dict gene -> expression (logFC|lineNumber) -> null
    exprData = {}
    inFile = open(infile,'r')
    header = inFile.readline().strip()
    (genePos,logfcPos) = processExprHeader(header)
    for nLine,line in enumerate(inFile):
        lineSplit = line.strip().split("\t")
        gene = lineSplit[genePos].strip()
        logFC = lineSplit[logfcPos].strip()
        # Sanity check on valid type and number
        check_logFC(logFC,gene)

        rawData[str(nLine)] = [gene,logFC]
        # Process rawData to have gene-centered dict
        # Returns dict of gene -> [var1,var2,..,varN], where a given var="logFC|lineNumber"
        if gene not in exprData.keys():
            exprData[gene] = {}
        # Collapse expression info separated with "|"
        # Keep track of what line each variant comes from
        expression = str(logFC) + "|" + str(nLine)
        # NOTE: Expression values can never be duplicated because of the different row numbers assigned
        # if expression in exprData[gene].keys():
        #     print("Found duplicated expression value '%s' for gene '%s' in line '%s'!" %(logFC,gene,str(nLine)))
        #     sys.exit(1)
        exprData[gene][expression] = None
    inFile.close()
    return (rawData,exprData)


# TODO: add more options from json.dump?
def write_to_json(myDict, outfile, indent=4):
    with open(outfile, 'w') as f:
        json.dump(myDict, f, ensure_ascii=False, indent=indent)
    return None


def write_match(matchMap, varMap, outfile, rawMap=None):

#     outFile = open(args.outfile,'w')
#     outHeader = header
#     outHeader += "\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\tCIViC_Drug_Support\tCIViC_Predictive\tCIViC_Diagnostic\tCIViC_Prognostic\tCIViC_Predisposing"
#     outFile.write(outHeader + "\n")



# TODO

    return None


def write_civic(varMap, outfile):

    return None
