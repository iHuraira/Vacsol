from Bio.Blast import NCBIWWW as ncbi
from  Bio.Blast import  NCBIXML as xml_opener
import regex as re
import pandas as pd

#file_read = input("Enter File name: ")

def file_creation():
    with open("GCF_000063585.1_ASM6358v1_protein.faa") as fasta_file:
        pattern = "\w{2}_\d{1,}\.\d{1,}"
        regex_expression = re.findall(pattern, fasta_file.read())[:5]
        total_length = len(regex_expression)
        accession_col = pd.DataFrame({"Accession" : regex_expression})
        acid = accession_col["Accession"]
        return acid

acid = file_creation()

def blast_search(acid):
    xml_list = []
    for x in acid:
        blast_records = ncbi.qblast("blastp", "nr" , x, hitlist_size=100)
        xml_list.append(blast_records)
    print(len(xml_list))
    return xml_list

blast_record = blast_search(acid)

def read_file(xml_doc):
    score_list = []
    e_value_list = []
    accession_list = []
    scientific_name = []
    title_list = []
    percentage_identity_list = []
    subject_list = []

    for xml_element in xml_doc:
        for record in xml_opener.parse(xml_element):
            for alignment in record.alignments:
                pattern = "(?<=\[).+?(?=\])"
                regex_expression = re.search(pattern, alignment.title)
                scientific_name.append(regex_expression.group())
                title_list.append(alignment.title)
                for hsp in alignment.hsps:
                    score_list.append(hsp.score)
                    e_value_list.append(hsp.expect)
                    accession_list.append(alignment.hit_id)
                    identity = "{:.0%}".format(hsp.identities / hsp.align_length)
                    percentage_identity_list.append(identity)
                    subject_list.append(hsp.sbjct)

    results_parsed  =  pd.DataFrame({"ACCESSION" : accession_list,
                      "E_VALUE" : e_value_list,
                    "SCORE" : score_list,
                                     "SCIENTIFIC NAME" : scientific_name,
                                     "TITLE" : title_list,
                                    "PERCENTAGE_IDENTITY" : percentage_identity_list,
                                     "SUBJECT SEQUENCE" : subject_list})
    results_parsed.to_csv("results_500.csv")

read_file(blast_record)

