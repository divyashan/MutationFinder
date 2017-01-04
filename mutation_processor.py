import sys
import os
import os.path
import string
import operator
import json
import lxml.html
import simplejson
import requests
from tempfile import mkstemp
sys.path.insert(0, 'pdfminer-20140328/tools/')
from pdf2txt import convert
from mutation_finder import amino_acid_three_to_one_letter_map
from mutation_finder import extract_mutations_from_lines_to_file
from mutation_finder import mutation_finder_from_regex_filepath
import urllib3
import certifi

http = urllib3.PoolManager(
    cert_reqs='CERT_REQUIRED', # Force certificate check.
    ca_certs=certifi.where(),  # Path to the Certifi bundle.
)

PATH_TO_AUTH_PW = '../../Documents/MutationFinder/secrets.txt'

# Given a gene name, finds all relevant papers in the Google Drive and
# mines them for allele names in relation to the gene. This data is then sent
# to the Django database corresponding to links between papers and alleles.
def process_gene_links(gene_name):
    print "Converting PDF to txt files..."
    if not os.path.isdir(gene_name + "_text_files"):
        print 'processing pdf files'
        process_pdf_files(gene_name)
    text_file_path = gene_name + "_text_files"


    print "Converting text files to MutationFinder friendly input..."
    process_text_files(text_file_path)

    print "Beginning to look through each paper for mutations..."
    find_mutations("mutation_finder_input.txt")

    print "Sending data to Django database..."
    process_mutation_finder_output('output.txt', gene_name)


# gene_name: The gene to be evaluated. Currently, the path is hardcoded
#            because the script is normally run locally.
# Populates an output folder with text files corresponding to each pdf file
def process_pdf_files(gene_name):
    home = os.path.join("/Users", "divya")
    papers_by_gene = os.path.join(home, "Google Drive", "Curation Papers",
                                  "Curation Papers by Gene", gene_name)
    if not os.path.exists(gene_name+"_text_files"):
        os.makedirs(gene_name + '_text_files')
    text_file_folder = gene_name + '_text_files'

    for paper in os.listdir(papers_by_gene):
        if not os.path.isdir(os.path.join(papers_by_gene, paper)):
            try:
                pmid = int(paper.split("_")[0])
            except ValueError:
                continue

            output_path = os.path.join(text_file_folder, paper[:-4] + '.txt')
            file_path = os.path.join(papers_by_gene, paper)
            try:
                convert(['./pdf2txt.py', '-o', output_path, file_path])
            except:
                print "Text extraction not allowed :("


# textPDFPath: Filepath to the folder containing all the text files. Each text
#              file's name should be underscore delimited, with the first field
#              being the PubMed ID of the paper. If these conditions aren't
#              satisfied, the file will be skipped.
# Writes a file properly formatted for MutationFinder.
def process_text_files(text_file_path):
    textString = ""
    paper_text_map = dict()
    for pdf in os.listdir(text_file_path):
        with open(os.path.join(text_file_path, pdf), 'r') as textFile:
            try:
                pmid = int(pdf.split("_")[0])
            except ValueError:
                continue

            fileText = textFile.read().replace('\n', ' ').replace('\r', ' ').replace('\t', ' ')
            print str(pdf.split("_")[0])
            stringToAdd = ''.join([str(pdf.split("_")[0]),
                                  "\t", str(fileText), "\n"])
            stringToAdd = filter(lambda x: x in string.printable, stringToAdd)
            textString += stringToAdd
    text_file = open("mutation_finder_input.txt", "w")
    text_file.write(textString)
    text_file.close()


# input_file: Filepath to the textfile containing two tab-deliminted fields;
#             the first contains the UID (often a PubMed ID) and the second
#             contains the text to be processed in relation to that ID
# mutation_finder_home: Delimits where the mutation_finder files are located
# Finds the mutations within the file and outputs a file with two fields in
# each row where the first contains the UID of the text processed and the
# second contains a list of mutations found by MutationFinder
def find_mutations(input_filepath, mutation_finder_home='./'):
    try:
        input_file = open(input_filepath)
    except IOError:
        print 'Cannot open specified input file: ', input_filepath

    output_filepath = "output.txt"
    regular_expression_filepath = ''.join([mutation_finder_home, 'regex.txt'])
    mutation_extractor = mutation_finder_from_regex_filepath(
            regular_expression_filepath)
    extract_mutations_from_lines_to_file(input_file, output_filepath,
                                         mutation_extractor, store_spans=True)


# output_file: File in which the MutationFinder output is located
# gene_name: name of gene that these papers are related to
# Processes output of MutationFinder to output a dictionary mapping papers to
# alleles (identified by p names)
def process_mutation_finder_output(output_file, gene_name):
    lines = open(output_file, 'r').readlines()

    pw = open(PATH_TO_AUTH_PW).read()
    gene_url = "https://counsyl.com/api/internal/curation/alleles?gene_name=" + gene_name + "&format=json"
    pw = open(PATH_TO_AUTH_PW).read()
    r = requests.get(gene_url, auth=('api', pw), headers={'host': 'www.counsyl.com'})
    try:
        curated_alleles = json.loads(r.text)['results']
    except Exception:
        print "Ran into JSON parsing error with allele api request"
        return

    # Creating a dictionary of papers (by PubMed ID) to alleles they contain
    allele_paper_info = dict()
    for line in lines:
        split_array = line.split('\t')
        if len(split_array) < 2:
            # This paper has no alleles mentioned within
            continue
        else:
            # Getting rid of the newline at the end of the last element
            split_array[-1] = split_array[-1][:-1]
            alleles = split_array[1:]
            pmid = split_array[0]
            paper_alleles = []
            # Iterating through all alleles for a given paper
            for allele_str in alleles:
                allele_info = allele_str.split(':')
                allele_id = name_to_allele_id(allele_info[0], gene_name, curated_alleles)
                if allele_id is not -1:
                    try:
                        pmid_url = 'http://www.ncbi.nlm.nih.gov/pubmed/?term=' + pmid
                        parsed_html = lxml.html.parse(pmid_url)
                        title = parsed_html.find(".//title").text[:-15]
                    except IOError:
                        title = "Unknown"

                    allele_paper_temp = dict()
                    allele_paper_temp["p_name"] = allele_info[0]
                    allele_paper_temp["allele_id"] = allele_id
                    allele_paper_temp["context"] = allele_info[1]
                    allele_paper_temp["gene_name"] = gene_name
                    allele_paper_temp["pmid"] = pmid
                    allele_paper_temp["paper_title"] = title
                    allele_paper_temp["mention"] = allele_info[2]
                    allele_paper_temp["added_by"] = 'MutationFinder'
                    paper_alleles.append(allele_paper_temp)
            allele_paper_info[int(pmid)] = paper_alleles
    f = open('jsonmf.txt', 'w')
    f.write(json.dumps(allele_paper_info))
    f.close()
    headers = {'content-type': 'application/json'}
    r = requests.post("https://counsyl.com/curation/paper_links/add",
                      json=allele_paper_info)
    print r, '\n'

    return json.dumps(allele_paper_info)

# missed_alleles: A set of all alleles recognized by MF for a certain gene
#                 and not matched to any curated alleles
# curated_lleles: A set of all alleles currently curated for a certain gene
# returns: The number of connections to the Counsyl database MF has made that
#          curators haven't (yet).
# Measures the number of MutationFinder alleles that exist in the curated
# database.
def measure_MF_curator_overlap(missed_alleles, curated_alleles):
    overlap = 0
    for allele in missed_alleles:
        base_changes = find_poss_base_changes(allele)
        # Only takes care of situations where the possible base change is
        # narrowed down to one
        if len(base_changes) == 1:
            try:
                amino_acid_position = int(allele[1:-1])
            except ValueError:
                continue

            # This is assuming every mutation MutationFinder has found is BRCA1
            standard_name = '(BRCA1):c.' + str(amino_acid_position*3-1)+ base_changes[0][0] + '>' + base_changes[0][1]

            for curated_allele in curated_alleles:
                if standard_name in curated_allele:
                    overlap += 1
    print overlap


# gene_name: Name of the gene to filter curations by
# Returns a dictionary of papers keyed to alleles (of a certain gene)
# curators have identified in each of these papers.
def process_allele_curations_by_gene(gene_name):
    # Used to create a list of JSON objects where
    # a is a list of SeqAlleleCurations
    # that fit some criteria
    import operator
    import json
    import simplejson
    from counsyl.product.sequencing.utils.hgvs_names import get_aa_change
    gene_allele_curations = SeqAlleleCuration.objects.filter(
        allele__gene__name__contains=gene_name, needs_recuration=False,
        approved=True, autocurated=False)

    # Creating an array of JSON representations for each allele curation
    get_pmid = operator.attrgetter("pmid")
    print gene_allele_curations
    alleleReferenceAA = []
    for alleleCuration in gene_allele_curations:
        temp = dict()
        paperCurations = alleleCuration.get_paper_curations()
        pmids = map(get_pmid, paperCurations)

        # Uncomment this line and comment the next if you'd like to key alleles
        # by their standard names instead of p.names
        # temp["allele_name"] = alleleCuration.standard_name
        temp["allele_name"] = get_aa_change(alleleCuration.allele,
                                            alleleCuration.allele.canonical_transcript)
        temp["allele_id"] = alleleCuration.allele_id
        temp["pmids"] = pmids
        temp["curation_id"] = alleleCuration.id
        alleleReferenceAA.append(json.dumps(temp))

    # Creates a dictionary that maps alleles (represented by allele_name) to
    # PubMed IDs of papers that reference them
    allele_dict = dict()
    for json in alleleReferenceAA:
        curationData = simplejson.loads(json)
        allele_name = curationData["allele_name"]
        if allele_name is not None and allele_name != "":
            pmid_list = allele_dict.get(allele_name, [])
            pmid_list.extend(curationData["pmids"])
            allele_dict[allele_name] = pmid_list

    # Flips dictionary created above to be paper PubMed IDs whose values are
    # the allele names that the paper references
    paper_dict = dict()
    for allele in allele_dict.keys():
        # Getting rid of curations with no allele name
        if allele != '' or allele is not None:
            for paper_id in allele_dict[allele]:
                allele_list = paper_dict.get(paper_id, [])
                allele_list.extend(allele)
                paper_dict[paper_id] = allele_list
    return paper_dict


def get_gene_dicts(gene_name):
    p_name_dict = dict()
    c_name_dict = dict()
    r = requests.get("https://counsyl.com/api/internal/curation/alleles?gene_name=CFTR&format=json")

    print json.loads(r.text)
    poss_alleles = {}
    print "generating curated alleles dicts"
    for curated_allele in poss_alleles:
        comp_p_name = get_aa_change(curated_allele.allele,
                                    curated_allele.allele.
                                    canonical_transcript)
        standard_name = curated_allele.allele.standard_name
        print standard_name
        if comp_p_name is None:
            continue
        if '*' in comp_p_name:
            alternate_name = comp_p_name[:-1] + 'X'
            p_name_dict[alternate_name] = curated_allele.allele
        p_name_dict[comp_p_name] = curated_allele.allele
        c_name_dict[standard_name] = curated_allele.allele


def name_to_allele_id(name, gene_name, curated_alleles):
    # Currently only manages situations where there's only one possible
    # base pair change; also does this by testing all possible frame shifts
    # in amino acid position
    standard_names = []
    if 'c' in name:
        standard_names.append(name)
    else:
        base_changes = find_poss_base_changes(name)
        if len(base_changes) != 1:
            return -1
        try:
            amino_acid_position = int(name[1:-1])
        except ValueError:
            return -1

        beginningStr = '(' + gene_name + '):c.'
        endStr = base_changes[0][0] + '>' + base_changes[0][1]

        standard_names.extend([beginningStr + str(amino_acid_position*3-1) + endStr,
                               beginningStr + str(amino_acid_position*3) + endStr,
                               beginningStr + str(amino_acid_position*3-2) + endStr])

    for curated_allele in curated_alleles:
        for standard_name in standard_names:
            if curated_allele['latest_approved_curation']:
                curated_standard_name = curated_allele['latest_approved_curation']['standard_name']
                if curated_standard_name and standard_name in curated_standard_name:
                    return curated_allele['id']
    return -1



# mf_paper_dict: Mapping of papers (represented by PubMed ID) to alleles found
#              by MF(represented by p names)
# counsyl_paper_dict: Mapping of papers (represented by PubMed ID) to alleles
#                   found by curators
# Performs comparison between MF mapping of papers to alleles to
# mapping generated from curator data. Returns missed alleles in a dictionary
# keyed by PubMed ID as well as leftover alleles in extraClassifications
def compare_MF_alleles(counsyl_paper_dict, mf_paper_dict):
    print "# Papers (Counsyl Curators): ", len(counsyl_paper_dict.keys())
    print "# Papers (Mutation Finder): ", len(mf_paper_dict.keys())

    # Testing for sensitivity and specificity
    # and recording which alleles are missed by MutationFinder
    # missed_alleles: Describes alleles in every paper that are curated but NOT found by MF
    # captured_alleles: Describes alleles in common between MF and curators
    # mf_leftover_alleles: Describes alleles in every paper that are found by MF but NOT curated
    totalCounsylClassified = 0
    totalMFClassified = 0
    missed_alleles = dict()
    captured_alleles = dict()
    mf_leftover_alleles = dict()
    totalSensitivity = 0.
    totalSpecificity = 0.
    papers_in_common = 0
    all_curated_alleles = []
    allMFAlleles = []
    numLessThanOne = 0
    for paper in counsyl_paper_dict.keys():
        curated_alleles = sorted(set(counsyl_paper_dict[paper]))
        missed_alleles[paper] = []
        captured_alleles[paper] = []
        mf_leftover_alleles[paper] = []
        if paper in mf_paper_dict.keys():
            papers_in_common += 1
            mutation_finder_alleles = sorted(set(mf_paper_dict[paper]))
            leftoverMFAlleles = mutation_finder_alleles[:]
            numCounsyl = len(curated_alleles)
            numMF = len(mutation_finder_alleles)

            totalCounsylClassified += numCounsyl
            totalMFClassified += numMF

            all_curated_alleles.extend(curated_alleles)
            allMFAlleles.extend(mutation_finder_alleles)

            inCommon = 0
            for allele in curated_alleles:
                if allele[-1] == '*':
                    compAllele = allele[:-1] + 'X'
                else:
                    compAllele = allele
                if compAllele in mutation_finder_alleles:
                    inCommon += 1
                    captured_alleles[paper].append(allele)
                    leftoverMFAlleles.remove(compAllele)
                else:
                    missed_alleles[paper].append(allele)
            # Percentage of the Curated Alleles that were found
            sensitivity = float(inCommon)/float(len(curated_alleles))
            # Percentage of the MF Alleles that are incorrect
            specificity = float(len(mutation_finder_alleles)-inCommon)/float(len(mutation_finder_alleles))

            totalSensitivity += sensitivity
            totalSpecificity += numMF * specificity
            mf_leftover_alleles[paper].extend(leftoverMFAlleles)

        else:
            missed_alleles[paper].extend(curated_alleles)
    all_curated_alleles = sorted(set(all_curated_alleles))
    allMFAlleles = sorted(set(allMFAlleles))
    print '\n\n\n', numLessThanOne
    print "Total alleles classified by curators: ", len(all_curated_alleles)
    print "Total alleles classified by MF: ", len(allMFAlleles)
    print "Total papers by Counsyl", len(counsyl_paper_dict)
    print "Total papers by MF ", len(mf_paper_dict)
    print "Sensitivity: ", (totalSensitivity/papers_in_common)
    print "specificity: ", (totalSpecificity/totalMFClassified)
    print "Papers in Common: ", papers_in_common

    extraClassifications = []
    for paper in mf_leftover_alleles.keys():
        extraClassifications.extend(mf_leftover_alleles[paper])
    extraClassifications = sorted(set(extraClassifications))

    for paper in missed_alleles:
        if len(captured_alleles[paper]) > 0:
            print "Paper: ", paper, "\tCaptured: ", captured_alleles[paper],
            '\tMissed: ', missed_alleles[paper], '\n'

    return (missed_alleles, extraClassifications)

# mutation: a string representation of the mutation where
#           the first letter represents the beginning amino acid
#           and the last letter represents the ending amino acid
# Returns a set of all possible nucleotide changes that would result
# in the specified amino acid transition.
def find_poss_base_changes(mutation):
    amino_acid_three_to_one_letter_map = \
            dict([('ALA','A'),('GLY','G'),('LEU','L'),('MET','M'),\
             ('PHE','F'),('TRP','W'),('LYS','K'),('GLN','Q'),('GLU','E'),('SER','S'),\
             ('PRO','P'),('VAL','V'),('ILE','I'),('CYS','C'),('TYR','Y'),('HIS','H'),\
             ('ARG','R'),('ASN','N'),('ASP','D'),('THR','T'),('XAA','X'),('GLX','Z'),\
             ('ASX','B')])

    amino_acid_to_codons = {"PHE":["TTT", "TTC"], "LEU":["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "SER":["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
                                "TYR":["TAT", "TAC"], "PRO":["CCT", "CCC", "CCA", "CCG"], "HIS":["CAT", "CAC"],
                                "GLN":["CAA", "CAG"], "ARG":["CGT", "CGC", "CCA", "CCG", "AGA", "AGG"], "ILE":["ATT", "ATC", "ATA"], "MET":["ATG"],
                                "ASN":["AAT", "AAC"], "THR":["ACT", "ACC", "ACA", "ACG"], "LYS":["AAA", "AAG"],
                                "ASP":["GAT", "GAC"], "VAL":["GTT", "GTC", "GTA", "GTG"], "ALA":["GCT", "GCC", "GCA", "GCG"],
                                "GLU":["GAA", "GAG"], "GLY":["GGT", "GGC", "GGA", "GGG"], "GLX":["GAA", "GAG", "CAA", "CAG"], "TRP":["TGG"],
                                "ASX":["AAT", "AAC", "GAT", "GAC"], "CYS":["TGT", "TGC"]}
    one_letter_to_amino_acid_map = dict()
    for amino_acid in amino_acid_three_to_one_letter_map.keys():
        one_letter_to_amino_acid_map[amino_acid_three_to_one_letter_map[amino_acid]] = amino_acid

    startAA = one_letter_to_amino_acid_map[mutation[0]]
    endAA = one_letter_to_amino_acid_map[mutation[-1]]
    if startAA is 'XAA' or endAA is 'XAA':
        return []
    start_codons = amino_acid_to_codons[startAA]
    end_codons = amino_acid_to_codons[endAA]
    poss_transitions = []
    for start in start_codons:
        start_ascii = map(lambda x: ord(x), start)
        for end in end_codons:
            end_ascii = map(lambda x: ord(x), end)
            subtracted = [a_i-b_i for a_i,b_i in zip(start_ascii, end_ascii)]
            filtered = filter(lambda x: x is not 0, subtracted)

            if len(filtered) == 1:
                # This means only one base pair change is needed
                changed_index = map(map_nonzeroes, subtracted).index(1)
                poss_transitions.append((start[changed_index],
                                        end[changed_index]))
    return sorted(set(poss_transitions))


# Maps all nonzero values of x to 1 and 0 otherwise
def map_nonzeroes(x):
    return 1 if x is not 0 else 0


home = os.path.join("/Users", "divya")
papers_by_gene_folder = os.path.join(home, "Google Drive", "Curation Papers",
                                     "Curation Papers by Gene")
finished = ['ABCC8', 'ABCD1', 'ACADM']
problems = []
seenStart = False

process_mutation_finder_output('output.txt', 'ABCD1')
for gene in os.listdir(papers_by_gene_folder):
    print gene
    if gene not in finished:
        if os.path.isdir(os.path.join(papers_by_gene_folder, gene)):
            process_gene_links(gene)
