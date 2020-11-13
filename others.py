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
                        print("Error! Could not find direction %s in support dictionary." %(direction))
                        sys.exit(1)
                    if clin_signf not in supportDict[direction].keys():
                        print("Error! Could not find clinical significance %s in support dictionary." %(clin_signf))
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
            print("Error! Unexpected ct case for gene %s." %(gene))
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
            print("Error! Unexpected support case for gene %s." %(gene))
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
            print("Error! Unexpected support case for gene %s." %(gene))
            sys.exit(1)

        ## Build support string for current drug(+gene in line)
        ## Format: DRUG:CT:SUPPORT
        drugSupport = drug + ':' + thisCT.upper() + ':' + tempSupport
        supportStrings.append(drugSupport)

    return(supportStrings)


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
            print("Error! Unexpected annotation format \'{}\'".format(annotation))
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
                print("Error! Unexpected annotation format \'{}\'".format(annotation))
                sys.exit(1)

            ## Special annotation case for drugs ie. 'DRUG1:annotation;..'
            if (indx == index_column+5) or (indx == index_column+6):
                # Retrieve drug name from field annotation
                drug = elements[0].upper()
                # Retrieve the second element (contains the relevant ClinTrial information)
                singleAnnot = elements[1]
                if drug not in drugDict.keys():
                    ## Drug should be present in drugDict already, as we should have parsed drug column previously
                    print("Error! Encountered unknown drug in {} from ClinTrial columns".format(annotation))
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

    return lineRoot, lineMap
