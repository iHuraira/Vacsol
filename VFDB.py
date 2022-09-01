import pandas as pd
from Bio import Align
import pip

aligner = Align.PairwiseAligner()

VFDB_DATA = pd.read_csv("./VFDB/VFDB_EXTRACTED.csv")
NON_HOMOLOGS = pd.read_csv("./Helicobacter/HELICOBACTER_NON_HOST_HOMOLOGS.csv")

LIST_SCORE = []
LIST_NON_HOMOLOG_SEQUENCE = []
LIST_VFDB_SEQUENCE = []
LIST_PERCENTAGE = []
LIST_IDS = []

for non_homolog, ids in zip(NON_HOMOLOGS["SUBJECT"],NON_HOMOLOGS["ID"]):
    for single_vfdb_seq in VFDB_DATA["FASTA"]:
        if len(non_homolog) == len(single_vfdb_seq):
            score = aligner.score(non_homolog, single_vfdb_seq)
            percentage = score / len(non_homolog) * 100
            if int(score) > 80:
                if percentage > 80:
                    LIST_SCORE.append(int(score))
                    LIST_PERCENTAGE.append(f"{int(percentage)}%")
                    LIST_IDS.append(ids)
                    LIST_VFDB_SEQUENCE.append(single_vfdb_seq)
                    LIST_NON_HOMOLOG_SEQUENCE.append(non_homolog)

DATA_FRAME = pd.DataFrame({
    "SCORE" : LIST_SCORE,
    "PERCENTAGE" : LIST_PERCENTAGE,
    "ID" : LIST_IDS,
    "QUERY" : LIST_NON_HOMOLOG_SEQUENCE,
    "SUBJECT" : LIST_VFDB_SEQUENCE
}).drop_duplicates()

DATA_FRAME.to_csv("./Helicobacter/VFDB.csv")