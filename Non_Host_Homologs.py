import pandas as pd
from Bio import Align

aligner = Align.PairwiseAligner()

def read_human_genome_file(file_name):
    with open(file_name) as fasta_file:
        Fasta_list = fasta_file.read().split(">")[1:]
        for i in Fasta_list:
            REFINED_LIST = []
            for i in Fasta_list:
                REFINED_LIST.append({"TITLE": "".join(i.split("\n")[:1]).split(" ")[0],
                                     "FASTA": "".join(i.split("\n")[1:])})
            return REFINED_LIST

fasta_human = read_human_genome_file("./Human_Reference_Genome/Human_Reference_Genome.faa")

fasta_bacteria = pd.read_csv("./Helicobacter/HP_LOCALIZATION.csv")
print(fasta_bacteria.columns)

PERCENT_LIST = []
SCORE_LIST = []
LOC_LIST = []
ID_LIST = []
QUERY_LIST = []
SUBJECT_LIST = []


for y_loc, y_id, y_proteome in zip(fasta_bacteria['LOC'],fasta_bacteria['ID'],fasta_bacteria['FASTA_PROTEOME']):
    for x in fasta_human:
        if len(x["FASTA"]) == len(y_proteome):
            #print(len(x["FASTA"]),len(y))
            aligner = Align.PairwiseAligner(match_score=1.0)
            score = aligner.score(x["FASTA"], y_proteome)
            percentage = score / len(x["FASTA"]) * 100
            if percentage < 30:
                #print(int(score), len(x["FASTA"]), len(y) ,f"{int(percentage)}%")
                PERCENT_LIST.append(f"{int(percentage)}%")
                SCORE_LIST.append(int(score))
                LOC_LIST.append(y_loc)
                ID_LIST.append(y_id)
                SUBJECT_LIST.append(y_proteome)
                QUERY_LIST.append(x["FASTA"])
                #print(int(percentage), int(score) ,y_loc, y_id, y_proteome, x["FASTA"])

FINAL_DATA = pd.DataFrame({
    "SCORE" : SCORE_LIST,
    "PERCENTAGE" : PERCENT_LIST,
    "ID" : ID_LIST,
    "QUERY" : QUERY_LIST,
    "SUBJECT" : SUBJECT_LIST
}).drop_duplicates()

FINAL_DATA.to_csv("./HELICOBACTER/HELICOBACTER_NON_HOST_HOMOLOGS.csv")