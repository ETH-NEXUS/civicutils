import os
import sys
import yaml
import json

# TODO: how to specify? or can they be imported as modules are instead?
BinPath = os.path.split(os.path.realpath(__file__))[0]


# FIXME: should this go in init.py?
def get_dict_aminoacids():
    f = BinPath + "/data/data.yml"
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
    f = BinPath + "/data/data.yml"
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
def processSnvHeader(headerSplit):
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
    headerSplit = header.strip().split('\t')
    (genePos,cPos,pPos,impactPos,exonPos) = processSnvHeader(headerSplit)
    extraHeader = []
    if impactPos:
        extraHeader.append("Variant_impact")
    if exonPos:
        extraHeader.append("Variant_exon")
    extraPos = []
    for pos,x in enumerate(headerSplit):
        if (pos == genePos) or (pos == cPos) or (pos == pPos):
            continue
        if impactPos:
            if pos == impactPos:
                continue
        if exonPos:
            if pos == exonPos:
                continue
        extraHeader.append(x)
        extraPos.append(pos)
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
        for p in extraPos:
            rawData[str(nLine)].append(lineSplit[p].strip())

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
    return (rawData,snvData,extraHeader)


# Retrieve the following column names from the given CNV header:
# - Gene: gene affected by the variant (one gene per row, so variants affecting multiple genes will be reported in different lines)
# - Variant_cnv: type of CNV variant (should be one of the following terms: 'GAIN', 'DUPLICATION', 'DUP', 'AMPLIFICATION', 'AMP', 'DELETION', 'DEL', 'LOSS')
def processCnvHeader(headerSplit):
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
    headerSplit = header.strip().split('\t')
    (genePos,cnvPos) = processCnvHeader(headerSplit)
    extraHeader = []
    extraPos = []
    for pos,x in enumerate(headerSplit):
        if (pos == genePos) or (pos == cnvPos):
            continue
        extraHeader.append(x)
        extraPos.append(pos)

    for nLine,line in enumerate(inFile):
        lineSplit = line.strip().split("\t")
        gene = lineSplit[genePos].strip()
        cnv = lineSplit[cnvPos].strip()
        rawData[str(nLine)] = [gene,cnv]
        for p in extraPos:
            rawData[str(nLine)].append(lineSplit[p].strip())
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
    return (rawData,cnvData,extraHeader)


# Retrieve the following column names from the given EXPR header:
# - Gene: gene affected by the expression change (one gene per row)
# - logFC: log fold-change of the gene expression (will be translated into 'OVEREXPRESSION' or 'UNDEREXPRESSION')
def processExprHeader(headerSplit):
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
    headerSplit = header.strip().split('\t')
    (genePos,logfcPos) = processExprHeader(headerSplit)
    extraHeader = []
    extraPos = []
    for pos,x in enumerate(headerSplit):
        if (pos == genePos) or (pos == logfcPos):
            continue
        extraHeader.append(x)
        extraPos.append(pos)
    for nLine,line in enumerate(inFile):
        lineSplit = line.strip().split("\t")
        gene = lineSplit[genePos].strip()
        logFC = lineSplit[logfcPos].strip()

        rawData[str(nLine)] = [gene,logFC]
        for p in extraPos:
            rawData[str(nLine)].append(lineSplit[p].strip())
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
    return (rawData,exprData,extraHeader)


def write_to_json(inDict, outfile, indent=1):
    with open(outfile, 'w') as f:
        json.dump(inDict, f, ensure_ascii=False, indent=indent)
    return None


def write_to_yaml(inDict, outfile):
    with open(outfile, 'w') as f:
        # Preserve the original order of the entries in the manifest template
        yaml.dump(inDict, f, default_flow_style=False, sort_keys=False)
    return None


def write_header_line(dataType,header,writeSupport):
    sorted_evidence_types = ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']

    # Variables only relevant for dataType="SNV"
    writeImpact = False
    writeExon = False

    if dataType == "SNV":
        mainHeader = "Gene\tVariant_dna\tVariant_prot"
        if "Variant_impact" in header:
            mainHeader += "\tVariant_impact"
            writeImpact = True
        if "Variant_exon" in header:
            mainHeader += "\tVariant_exon"
            writeExon = True

    if dataType == "CNV":
        mainHeader = "Gene\tVariant_cnv"

    if dataType == "EXPR":
        mainHeader = "Gene\tlogFC"

    clean_header = []
    if header:
        for tmp in header:
            if (tmp != "Variant_impact" and tmp != "Variant_exon"):
                clean_header.append(tmp)
    if clean_header:
        mainHeader += "\t%s" %("\t".join(clean_header))

    if writeSupport:
        outHeader = "%s\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\tCIViC_Drug_Support\t%s" %(mainHeader,"\t".join(["CIViC_" + x for x in sorted_evidence_types]))
    else:
        outHeader = "%s\tCIViC_Tier\tCIViC_Score\tCIViC_VariantType\t%s" %(mainHeader,"\t".join(["CIViC_" + x for x in sorted_evidence_types]))

    return (outHeader,clean_header,writeImpact,writeExon)


def write_output_line(tier,mainLine,geneScores,geneVarTypes,drugSupport,resultMap,writeSupport):
    sorted_evidence_types = ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']

    # Remove the "tier" tag from the assined tier
    if tier.startswith("tier_"):
        tier = tier.replace("tier_", "")

    if geneScores:
        outLine = mainLine + "\t" + tier + "\t" + ";".join(geneScores)
    else:
        outLine = mainLine + "\t" + tier + "\t."

    if geneVarTypes:
        outLine += "\t" + ";".join(geneVarTypes)
    else:
        outLine += "\t."

    if writeSupport:
        if drugSupport:
            outLine += "\t" + ";".join(drugSupport)
        else:
            outLine += "\t."

    for evidence_type in sorted_evidence_types:
        # if evidence_type not in resultMap.keys():
        #     raise ValueError("Evidence type '%s' is not found in provided 'resultMap'!" %(evidence_type))
        if evidence_type in resultMap.keys():
            if resultMap[evidence_type]:
                outLine += "\t" + ";".join(resultMap[evidence_type])
            else:
                outLine += "\t."
        else:
            outLine += "\t."

    return outLine


# Write information about a single cancer type item into one or more structured strings
# Report drug names in the string for 'predictive' evidence (writeDrug=True)
# i.e. DISEASE[|DRUG1,DRUG2..](direction, significance(level(PMID,..,PMID),level(..)));
def write_evidences(item, writeDrug=False, writeCt=None, writeComplete=False):
    evidences = []
    # For each disease found in the provided item
    for disease in item.keys():
        # For each drug associated with the given cancer type
        for drug in item[disease].keys():
            # For each evidence associated with the given drug
            # Evidences are simplified by using the combined form 'direction:significance'
            for evidence in item[disease][drug].keys():
                # If drug=True, write drug information, i.e. DISEASE|DRUG(..)
                if writeDrug:
                    # Always one drug (single or combination with '+')
                    if writeCt:
                        outString = disease + '|' + writeCt.upper() + '|' + drug + '('
                    else:
                        outString = disease + '|' + drug + '('
                else:
                    if writeCt:
                        outString = disease + '|' + writeCt.upper() + '('
                    else:
                        outString = disease + '('
                # Split the evidence direction and clinical significance
                evidenceArr = evidence.strip().split(':')
                if (len(evidenceArr) != 2):
                    raise ValueError("Unexpected format of evidence '%s'! Please provide string as 'EVIDENCE_DIRECTION:CLINICAL_SIGNIFICANCE'." %(evidence))
                direction = evidenceArr[0]
                clin_signf = evidenceArr[1]
                outString += direction + ',' + clin_signf + '('
                # There may be several levels grouped per evidence
                levels = []
                for level in item[disease][drug][evidence].keys():
                    pmids = []
                    # There may be several publications (i.e. PMIDs) grouped per level
                    for z in item[disease][drug][evidence][level]:
                        if writeComplete:
                            pmids.append(z)
                        else:
                            zSplit = z.strip().split(":")
                            if (len(zSplit) != 5):
                                raise ValueError("Unexpected format of evidence item '%s'! Please provide string as 'EVIDENCE_ID:EVIDENCE_STATUS:SOURCE_STATUS:VARIANT_ORIGIN:RATING'." %(z))
                            pubId = zSplit[0].strip()
                            pmids.append(pubId)
                    levels.append(level + '(' + ','.join(pmids) + ')')
                outString += ','.join(levels) + '))'
                evidences.append(outString)
    return evidences


def write_match(matchMap, varMap, rawMap, header, dataType, outfile, hasSupport=True, hasCt=True, writeCt=False, writeSupport=True, writeComplete=False):
    # NOTE: uppercase is critical for matching!
    sorted_evidence_types = ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']
    evidenceType = "PREDICTIVE"
    special_cases = ["NON_SNV_MATCH_ONLY","NON_CNV_MATCH_ONLY","NON_EXPR_MATCH_ONLY"]
    sorted_cts = ["ct","gt","nct"]

    from utils import check_match_before_writing,check_keys,check_keys_not,check_data_type,check_dict_entry
    check_match_before_writing(matchMap,varMap,rawMap,hasSupport,hasCt,writeCt,writeSupport,writeComplete)
    check_data_type(dataType)
# TODO: sort line number numerically to ensure correct order of lines in output
    outfile = open(outfile,'w')

    # Retrieve the output header given the argument selection
    (outHeader,cleanHeader,writeImpact,writeExon) = write_header_line(dataType,header,writeSupport)
    outfile.write(outHeader + "\n")

    for nLine in rawMap.keys():
        lineArr = rawMap[nLine]
        extraLine = []
        if dataType == "SNV":
            if (len(lineArr) < 5):
                raise ValueError("Must provide at least 5 elements to describe a SNV variant (even if some can be empty): gene,dna,[prot],[impact],[exon],..'")
            gene = lineArr[0]
            cVar = lineArr[1]
            pVar = lineArr[2]
            impact = lineArr[3]
            exon = lineArr[4]
            combId = cVar + "|" + pVar + "|" + impact + "|" + exon + "|" + str(nLine)
            # Extract any additional fields that might be present for this line
            for pos in range(5,len(lineArr)):
                extraLine.append(lineArr[pos])
            # Build line string for writing to output
            mainLine = gene + "\t" + cVar + "\t" + pVar
            if writeImpact:
                mainLine += "\t" + impact
            if writeExon:
                mainLine += "\t" + exon

        if dataType == "CNV":
            if (len(lineArr) < 2):
                raise ValueError("Must provide at least 2 elements to describe a CNV variant: gene,cnv,..'")
            gene = lineArr[0]
            cnv = lineArr[1]
            combId = cnv + "|" + str(nLine)
            # Extract any additional fields that might be present for this line
            for pos in range(2,len(lineArr)):
                extraLine.append(lineArr[pos])
            # Build line string for writing to output
            mainLine = gene + "\t" + cnv

        if dataType == "EXPR":
            if (len(lineArr) < 2):
                raise ValueError("Must provide at least 2 elements to describe a EXPR variant: gene,logFC,..'")
            gene = lineArr[0]
            logFC = lineArr[1]
            combId = str(logFC) + "|" + str(nLine)
            # Extract any additional fields that might be present for this line
            for pos in range(2,len(lineArr)):
                extraLine.append(lineArr[pos])
            # Build line string for writing to output
            mainLine = gene + "\t" + str(logFC)

        # Sanity check that as many data fields were provided as in the header
        if len(extraLine) != len(cleanHeader):
            raise ValueError("Number of fields available does not match provided header!")
        # Add extra fields to the current line build
        for extra in extraLine:
            mainLine += "\t" + extra

        # Check if matchMap contains the provided input variants
        if gene not in matchMap.keys():
            raise ValueError("Provided gene '%s' is not contained in 'matchMap'." %(gene))
        if combId not in matchMap[gene].keys():
            raise ValueError("Provided variant '%s' is not contained in 'matchMap' of gene '%s'." %(combId,gene))
        for tier in matchMap[gene][combId].keys():
            geneScores = []
            geneVarTypes = []
            drugSupport = []
            resultMap = {}
            writeLine = False
            if tier != "tier_4":
                allVariants = []
                if hasSupport:
                    drug_support = matchMap[gene][combId][tier]["drug_support"]
                    for tmpVar in matchMap[gene][combId][tier]["matched"]:
                        allVariants.append(tmpVar)
                    if writeSupport:
                        for i in drug_support:
                            drugSupport.append(i.upper())
                else:
                    if writeSupport:
                        raise ValueError("Option 'writeSupport' cannot be selected when 'hasSupport'=False!")
                    for tmpVar in matchMap[gene][combId][tier]:
                        allVariants.append(tmpVar)

                for varId in allVariants:
                    # NOTE: check for special case when tier3 but no matching variant returned for the given data type
                    # This is a dummy tag and not an actual variant record from CIVICdb, so skip checking in varMap
                    if varId.upper() in special_cases:
                        # In this case, current line will be associated with tier3, but all columns will be empty with "."
                        for evidence_type in sorted_evidence_types:
                            if evidence_type not in resultMap.keys():
                                resultMap[evidence_type] = []
                        continue

                    variant = varMap[gene][varId]['name']
                    geneScores.append(gene + ':' + variant + ':' + str(varMap[gene][varId]['civic_score']))
                    geneVarTypes.append(gene + ':' + variant + ':' + ','.join(varMap[gene][varId]['types']))
                    for evidence_type in sorted_evidence_types:
                        if evidence_type in varMap[gene][varId]["evidence_items"].keys():
                            if evidence_type not in resultMap.keys():
                                resultMap[evidence_type] = []
                            writeDrug = False
                            if evidence_type == evidenceType:
                                writeDrug=True
                            if hasCt:
                                check_keys(list(varMap[gene][varId]["evidence_items"][evidence_type].keys()),"varMap",sorted_cts,matches_all=True)
                                for ct in varMap[gene][varId]["evidence_items"][evidence_type].keys():
                                    if writeCt:
                                        results = write_evidences(varMap[gene][varId]["evidence_items"][evidence_type][ct], writeDrug=writeDrug, writeCt=ct, writeComplete=writeComplete)
                                    else:
                                        results = write_evidences(varMap[gene][varId]["evidence_items"][evidence_type][ct], writeDrug=writeDrug, writeCt=None, writeComplete=writeComplete)
                                    for x in results:
                                        resultMap[evidence_type].append(gene + ":" + variant + ":" + x)
                            else:
                                check_keys_not(list(varMap[gene][varId]["evidence_items"][evidence_type].keys()),"varMap",sorted_cts)
                                if writeCt:
                                    raise ValueError("Option 'writeCt' cannot be selected when 'hasCt'=False!")
                                results = write_evidences(varMap[gene][varId]["evidence_items"][evidence_type], writeDrug=writeDrug, writeCt=None, writeComplete=writeComplete)
                                for x in results:
                                    resultMap[evidence_type].append(gene + ":" + variant + ":" + x)

                # Only write line current tier when there was at least one match for it
                if allVariants:
                    writeLine = True

            else:
                if hasSupport:
                    if matchMap[gene][combId][tier]["matched"]:
                        writeLine = True
                else:
                    if writeSupport:
                        raise ValueError("Option 'writeSupport' cannot be selected when 'hasSupport'=False!")
                    if matchMap[gene][combId][tier]:
                        writeLine = True

            if writeLine:
                outLine = write_output_line(tier,mainLine,geneScores,geneVarTypes,drugSupport,resultMap,writeSupport)
                outfile.write(outLine + "\n")

    return None
