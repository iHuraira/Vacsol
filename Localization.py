from Bio import Entrez, SeqIO
import requests as r
from io import StringIO  # Python 3
import regex as re
import os
import pandas as pd
import numpy as np
from time import sleep
from Bio import Align

aligner = Align.PairwiseAligner()

string_fasta = """>WP_003355768.1 MULTISPECIES: tRNA (guanosine(46)-N7)-methyltransferase TrmB [Clostridium]
MRLRKKWWARPEIEASDKFAEEPKELRGKWNKEFNNNNDIHLELGCGRGGFISQLVEKNKDINYVGIDLKDEVIVYAIRK
VKEKEEEVKREFKNIKFVTMNIMGIAEVFDKNEISKIYINFCNPWPKERHNKRRLTHTKLLTEYKKFLKPNTEIWFKTDD
KELFEDSQEYFKESGFNIEYITYDLHNSDFKENIKTEYETKFETMGMKIMFLKARLL"""

#Fasta_list = []

#read_submitted_fasta_file() function will read a file provided to it and it will return a dictionary with fasta and protein ID.

def read_submitted_fasta_file(file_name):
    with open(file_name) as fasta_file:
        Fasta_list = fasta_file.read().split(">")[1:]
        for i in Fasta_list:
            REFINED_LIST = []
            for i in Fasta_list:
                REFINED_LIST.append({"TITLE": "".join(i.split("\n")[:1]).split(" ")[0],
                                     "FASTA": "".join(i.split("\n")[1:])})
            return REFINED_LIST

fasta = read_submitted_fasta_file("./Helicobacter/HP_PROTEOME.faa")
print(len(fasta))

#read_submitted_fasta function will read a fasta sequence provided to it and it will return a dictionary with fasta and protein ID.

def read_submitted_fasta(string_fasta):
    Fasta_list = string_fasta.split(">")[1:]
    for i in Fasta_list:
        Refined_Object = {"TITLE": "".join(i.split("\n")[:1]).split(" ")[0],
                          "FASTA": "".join(i.split("\n")[1:])}
        print(Refined_Object)


def get_single_strain(single):
    if single[0][0].isupper() == True:
        word = single[0] +" "+ single[1]
        files_list = os.listdir("./Bacteria")
        for i in files_list:
            x = i.replace("_", ".").split(".")[1:][:2]
            if x[0] + " " + x[1] == word:
                return i
    else:
        word = single[1] + " " + single[0]
        files_list = os.listdir("./Bacteria")
        for i in files_list:
            x = i.replace("_", ".").split(".")[1:][:2]
            if x[0] + " " + x[1] == word:
                return i



# get_names() function will fetch all of the names of the bacteria present in the directory "Bacteria"
def get_names():
    files_list = os.listdir("./Bacteria")
    cleaned_list = [x.replace("_", ".").split(".")[1:][:2] for x in files_list][:-1]
    refined_list = []
    for i in cleaned_list:
        refined_list.append({
            "First_Name" : i[0],
            "Second_Name" : i[1]
        })
    return refined_list

genus = get_names()

#short variable represents total number of bacteria genus.
shorted_list = []
for i in genus:
    shorted_list.append(i["First_Name"])
short = list(set(shorted_list))
print("Choose Bacteria Genus: ", short)
print(len(short))




#search for the bacteria first name.
final_list = []
search_genus = input("Enter Bacteria Genus: ")
second_name_filtered = []

for x in genus:
    if search_genus in x["First_Name"]:
        second_name_filtered.append(x["Second_Name"])
        final_list.append(x["First_Name"])
print(second_name_filtered)

#search for the second name.

search_second_name = input("Enter Bacteria Second Name: ")
for y in second_name_filtered:
    if y == search_second_name:
        final_list.append(y)

single = list(set(final_list))
first, second = single[0], single[1]

file_name = get_single_strain(single)

LOCATION_LIST = []
BIT_SCORE_LIST = []
PROTEIN_ID_LIST = []

with open("./Bacteria/"+file_name, "r") as fasta_file:
    for line in fasta_file.readlines():
        without_hastag = re.findall("^#", line)
        if not  without_hastag:
            list_sequence = line.split()
            #print(list_sequence[2:3][0])
            if list_sequence[2:3][0] == "inner":
                LOCATION_LIST.append(list_sequence[2:3][0] + " membrance")
            elif list_sequence[2:3][0] == "outer":
                LOCATION_LIST.append(list_sequence[2:3][0] + " membrance")
            else:
                LOCATION_LIST.append(list_sequence[2:3][0])

            BIT_SCORE_LIST.append(int(list_sequence[1]))
            PROTEIN_ID_LIST.append(list_sequence[0])



LOCATION_DATA = pd.DataFrame({
    "LOCATION" : LOCATION_LIST,
    "BIT SCORE" : BIT_SCORE_LIST,
    "PROTEIN ID" : [PROTEIN_ID.split("|")[1] for PROTEIN_ID in PROTEIN_ID_LIST],
    #"PROTEIN FASTA" : PROTEIN_FAST_LIST
}).reset_index(drop=True)

BIT_SCORE_MEAN = int(np.mean(LOCATION_DATA["BIT SCORE"]))

mean_filtered = LOCATION_DATA[LOCATION_DATA["BIT SCORE"] > 70]
remove_list = ["secreted","outer membrane","fimbrium" ]
NEW_DATA_FRAME = mean_filtered[mean_filtered["LOCATION"].isin(remove_list)].reset_index(drop=True)


PROTEIN_ID_LIST_SECOND = []
PROTEIN_LOC = []
PROTEIN_BIT_SCORE = []
PROTEIN_PROTEOME = []
PROTEIN_UNIPROT = []
PROTEIN_PERC = []

for single_id, loc, bit_score in zip(NEW_DATA_FRAME["PROTEIN ID"], NEW_DATA_FRAME['LOCATION'], NEW_DATA_FRAME['BIT SCORE']):
    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl + single_id + ".fasta"
    response = r.post(currentUrl)
    cData = ''.join(response.text)
    Seq = StringIO(cData)
    pSeq = list(SeqIO.parse(Seq, 'fasta'))
    for i in pSeq:
        for single_fasta in fasta:
            if len(single_fasta["FASTA"]) == len(i.seq):
                aligner = Align.PairwiseAligner(match_score=1.0)
                score = aligner.score(i.seq, single_fasta["FASTA"])
                percentage = score / len(single_fasta["FASTA"]) * 100
                if percentage > 80:
                    print(int(score), len(single_fasta["FASTA"]), len(i.seq), f"{int(percentage)}%")
                    PROTEIN_PROTEOME.append(single_fasta["FASTA"])
                    PROTEIN_UNIPROT.append(i.seq)
                    PROTEIN_ID_LIST_SECOND.append(single_id)
                    PROTEIN_LOC.append(loc)
                    PROTEIN_BIT_SCORE.append(bit_score)
                    PROTEIN_PERC.append(f"{int(percentage)}%")
                    sleep(3)

FASTA_INCLUDED = pd.DataFrame({
    "ID" : PROTEIN_ID_LIST_SECOND,
    "LOC" : PROTEIN_LOC,
    "BIT SCORE" : PROTEIN_BIT_SCORE,
    "PERCENTAGE" : PROTEIN_PERC,
    "FASTA_PROTEOME": PROTEIN_PROTEOME,
    "FASTA_UNIPROT" : PROTEIN_UNIPROT
})

print(FASTA_INCLUDED)

FASTA_INCLUDED.to_csv("./Helicobacter/HP_LOCALIZATION.csv")

