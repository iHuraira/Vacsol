import pandas as pd
from Bio import SeqIO
from tkinter import *
from tkinter.filedialog import askopenfile
from tkinter import messagebox
from tkinter import ttk
import os
import re
from time import sleep
import requests as r
from io import StringIO
from Bio import Align
from ttkwidgets.autocomplete import AutocompleteCombobox
from tkinter import filedialog
import threading
from PIL import ImageTk,Image
import sys

window = Tk()
window.title("VacSol")
window.geometry("995x700")
window.maxsize(height=700, width=995)

window.iconbitmap("./Images/Logo_Sep.ico")

aligner = Align.PairwiseAligner()

First_Color = Second_Color = "#EAEAEA"
Buttons_Colors = Label_Color = "#646C8F"
Button_Text = "#F5F5F5"
font_family = "Bahnschrift"

Submission_method_frame = Frame(window, width=1200, height=330, bg=First_Color)
Submission_method_frame.grid(row=0, column=0)

Job_Status = Label(Submission_method_frame, text="Job Status: Not Running",bg=First_Color, font=(font_family,(12)), fg=Buttons_Colors)
Job_Status.place(relx=0.65, rely=0.1)

Logo_Label_Image = ImageTk.PhotoImage(Image.open("./Images/Logo.png").resize((230,60)))
Logo_Label = Label(Submission_method_frame, image=Logo_Label_Image, bg=First_Color)
Logo_Label.place(relx=0.018, rely=0.07)

Description_Label = Label(Submission_method_frame, text="Vaccine Candidate Prediction Pipeline", bg=First_Color, font=(font_family,(16)), fg=Buttons_Colors)
Description_Label.place(relx=0.28, rely=0.16)

Submission_label = Label(Submission_method_frame, text="Choose a Submission Method", bg=First_Color, font=(font_family,(14)), fg=Buttons_Colors)
Submission_label.place(relx=0.018, rely=0.27)

def read_submitted_fasta_file(file_name):
    title_list = []
    fasta_list = []
    for seq_record in SeqIO.parse(file_name, "fasta"):
        title_list.append(seq_record.id)
        fasta_list.append(seq_record.seq)
    return title_list,fasta_list

def read_submitted_fasta_sequence(sequence):
    title_list = []
    fasta_list = []
    fasta_io = StringIO(sequence)
    records = SeqIO.parse(fasta_io, "fasta")
    for record in records:
        title_list.append(record.id)
        fasta_list.append(record.seq)
    return title_list,fasta_list

string_fasta = """>WP_003355768.1
MRLRKKWWARPEIEASDKFAEEPKELRGKWNKEFNNNNDIHLELGCGRGGFISQLVEKNKDINYVGIDLKDEVIVYAIRK
VKEKEEEVKREFKNIKFVTMNIMGIAEVFDKNEISKIYINFCNPWPKERHNKRRLTHTKLLTEYKKFLKPNTEIWFKTDD
KELFEDSQEYFKESGFNIEYITYDLHNSDFKENIKTEYETKFETMGMKIMFLKARLL"""

fasta_sequence = read_submitted_fasta_sequence(string_fasta)



Fasta_Entry_Var = StringVar()
Fasta_Entry = Text(Submission_method_frame, height=13, width=87,padx=8, pady=8)
Fasta_Entry.place(rely=0.4, relx=0.025)

readed_file_data = []

filename_var = StringVar()

def FASTA_BROWSE():

    try:
        filename = askopenfile(parent=window, initialdir="/Desktop", title="Choose Fasta File", filetypes=[("all files", "*.*")])
        filename_var.set(filename.name)
        Fasta_Entry.delete('1.0', END)
        FASTA_BROWSE_LABEL.config(text=filename.name.split("/")[-1])
    except AttributeError:
        messagebox.showwarning("Warning", "No File Selected")


Fasta_Browse = Button(Submission_method_frame, text="BROWSE", command=FASTA_BROWSE, height=2, width=20, bd=0, background=Buttons_Colors, font = (font_family,12), fg=Button_Text)
Fasta_Browse.place(relx=0.65, rely=0.4)

FASTA_BROWSE_LABEL = Label(Submission_method_frame, text="Select File", height=2, width=20, bd=0, background=First_Color, font = (font_family,12), fg=Buttons_Colors)
FASTA_BROWSE_LABEL.place(relx=0.65, rely=0.28)


def FILL_EXAMPLE():
    Fasta_Entry.insert(END, string_fasta)

def MUST_EVALUATE():
    if Must_Evalutate_var.get() == "True":
        messagebox.showinfo("Info", "Must Evaluate Selected")
        Localization.select()
        Non_Homologs_check_btn.select()
        VFDB_check_btn.select()
        Epitope_check_btn.select()

    else:
        Localization.deselect()
        Non_Homologs_check_btn.deselect()
        VFDB_check_btn.deselect()
        Epitope_check_btn.deselect()


SAVE_FILE_Var = StringVar()

def SAVE_FILE():
    folder_selected = filedialog.askdirectory()
    SAVE_FILE_Var.set(folder_selected)
    print(folder_selected)



def FASTA_SUBMIT_BTN():

    data_list = []

    input = Fasta_Entry.get("1.0",END).upper()
    if len(input) == 1 and filename_var.get() != "":

        title, fasta = read_submitted_fasta_file(filename_var.get())
        data_list.append({
                    "TITLE" : title,
                    "FASTA" : fasta
                })
        FASTA_BROWSE_LABEL.config(text="Select A File")

    elif len(input) != 0 and filename_var.get() != "":
        messagebox.showwarning("Warning", "Choose One Submission Method")
        filename_var.set("")
        FASTA_BROWSE_LABEL.config(text="Select A File")

    else:
        Amino_acid_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
                           "V"]
        fasta_io = StringIO(input)
        records = SeqIO.parse(fasta_io, "fasta")
        sure_var = ""
        for record in records:
            for single_amino_input in list(record.seq):
                if single_amino_input in Amino_acid_list:
                    sure_var = True
                else:
                    sure_var = False
                    break
            if sure_var == True:
                title, fasta = read_submitted_fasta_sequence(input)
                data_list.append({
                    "TITLE" : title,
                    "FASTA" : fasta
                })
                break

            else:
                messagebox.showerror("Error", "Mistake in Fasta Sequence")
                Fasta_Entry.delete('1.0', END)

    if len(data_list) == 0 and filename_var.get() == "":
        messagebox.showerror("Warning", "Choose a Submission Method")
    else:

        if Must_Evalutate_var.get() != "True":

            Job_Status.config(text="Job Status: Running")

            if SAVE_FILE_Var.get() != "":

                if Localization_var.get() == "True":
                    data = StringIO(LOC_TREE_var.get())
                    df = pd.read_csv(data, sep='\s+')[:-1]
                    df["LOCATION"] = df["BIT"]
                    df.drop(columns=["BIT", "ID"], inplace=True)

                    PROTEIN_ID_LIST_SECOND = []
                    PROTEIN_LOC = []
                    PROTEIN_BIT_SCORE = []
                    PROTEIN_UNIPROT = []

                    for single_id, loc, bit_score in zip(df["PROTEIN"], df['LOCATION'], df['SCORE']):
                        baseUrl = "http://www.uniprot.org/uniprot/"
                        currentUrl = baseUrl + single_id + ".fasta"
                        response = r.post(currentUrl)
                        cData = ''.join(response.text)
                        Seq = StringIO(cData)
                        pSeq = list(SeqIO.parse(Seq, 'fasta'))
                        for i in pSeq:
                            PROTEIN_UNIPROT.append(i.seq)
                            PROTEIN_ID_LIST_SECOND.append(single_id)
                            PROTEIN_LOC.append(loc)
                            PROTEIN_BIT_SCORE.append(bit_score)
                            print(single_id, " DONE")
                            sleep(5)

                    Localization_Data = pd.DataFrame({
                        "ID": PROTEIN_ID_LIST_SECOND,
                        "LOC": PROTEIN_LOC,
                        "SEQ": PROTEIN_UNIPROT,
                        "SCORE": PROTEIN_BIT_SCORE
                    })

                    PROTEIN_PROTEOME = []
                    PROTEIN_PERC = []
                    PROTEIN_EXTRACTED_LOCTREE = []
                    PROTEIN_EXTRACTED_LOCATION = []

                    for extracted_seq, extracted_loc in zip(Localization_Data["SEQ"], Localization_Data["LOC"]):
                        for single_Localization in data_list:
                            for single_Localization_fasta in single_Localization["FASTA"]:
                                aligner = Align.PairwiseAligner(match_score=1.0)
                                score = aligner.score(extracted_seq, single_Localization_fasta)
                                percentage = score / len(single_Localization_fasta) * 100
                                if Localization_identity_var.get() != 0:
                                    if percentage > Localization_identity_var.get():
                                        PROTEIN_PROTEOME.append(single_Localization_fasta)
                                        PROTEIN_PERC.append(f"{int(percentage)}%")
                                        PROTEIN_EXTRACTED_LOCTREE.append(extracted_seq)
                                        PROTEIN_EXTRACTED_LOCATION.append(extracted_loc)

                    FASTA_INCLUDED = pd.DataFrame({
                        "LOC": PROTEIN_EXTRACTED_LOCATION,
                        "PERCENTAGE": PROTEIN_PERC,
                        "FASTA_SUBMITTED": PROTEIN_PROTEOME,
                        "FASTA_UNIPROT": PROTEIN_EXTRACTED_LOCTREE
                    }).drop_duplicates()

                    if len(FASTA_INCLUDED) == 0:
                        messagebox.showerror("Error", "No Relevant Data found")
                    else:
                        FASTA_INCLUDED.to_csv(f"{SAVE_FILE_Var.get()}/FIRST_LOCALIZATION.csv")

                elif Non_Homologs_var.get() == "True":
                    print(True)

                    HUMAN_GENOME = pd.read_csv("./Data/Human_Reference_Genome.csv")["FASTA"]
                    PERCENT_LIST = []
                    SCORE_LIST = []
                    QUERY_LIST = []
                    SUBJECT_LIST = []

                    for single in data_list:
                        for single_fasta in single["FASTA"]:
                            for x in HUMAN_GENOME:
                                if len(x) >= len(single_fasta):
                                    aligner = Align.PairwiseAligner(match_score=1.0)
                                    score = aligner.score(x, single_fasta)
                                    percentage = score / len(x) * 100
                                    if Non_Homologs_identity_var.get() != 0:
                                        if percentage < Non_Homologs_identity_var.get():
                                            PERCENT_LIST.append(f"{int(percentage)}%")
                                            SCORE_LIST.append(int(score))
                                            SUBJECT_LIST.append(single_fasta)
                                            QUERY_LIST.append(x)

                        FINAL_DATA = pd.DataFrame({
                            "SCORE": SCORE_LIST,
                            "PERCENTAGE": PERCENT_LIST,
                            "HUMAN_PROTEINS": QUERY_LIST,
                            "SUBMITTED_SEQUENCE": SUBJECT_LIST
                        }).drop_duplicates()

                        if len(FINAL_DATA) == 0:
                            messagebox.showerror("Error", "No Relevant Data found")
                        else:
                            FINAL_DATA.to_csv(f"{SAVE_FILE_Var.get()}/NON_HOMOLOGS.csv")

                elif VFDB_var.get() == "True":

                    VFDB_DATA = pd.read_csv("./Data/VFDB.csv")

                    LIST_NON_HOMOLOG_SEQUENCE = []
                    LIST_VFDB_SEQUENCE = []
                    LIST_PERCENTAGE = []

                    aligner = Align.PairwiseAligner()

                    for single in data_list:
                        for single_fasta in single["FASTA"]:
                            for single_vfdb_seq, title in zip(VFDB_DATA["FASTA"], VFDB_DATA["TITLE"]):
                                score = aligner.score(single_fasta, single_vfdb_seq)
                                percentage = score / len(single_fasta) * 100
                                if VFDB_identity_var.get() != 0:
                                    if percentage > VFDB_identity_var.get():
                                        LIST_PERCENTAGE.append(f"{int(percentage)}%")
                                        LIST_VFDB_SEQUENCE.append(single_vfdb_seq)
                                        LIST_NON_HOMOLOG_SEQUENCE.append(single_fasta)

                    DATA_FRAME = pd.DataFrame({
                        "PERCENTAGE": LIST_PERCENTAGE,
                        "QUERY": LIST_NON_HOMOLOG_SEQUENCE,
                        "SUBJECT": LIST_VFDB_SEQUENCE
                    }).drop_duplicates()

                    if len(DATA_FRAME) == 0:
                        messagebox.showerror("Error", "No Relevant Data found")
                    else:
                        DATA_FRAME.to_csv(f"{SAVE_FILE_Var.get()}/FIRST_VFDB.csv")

                elif Epitope_var.get() == "True":

                    aligner = Align.PairwiseAligner()

                    B_CELL_DATA_LINEAR = pd.read_csv("./Data/Linear.csv")["FASTA"]
                    B_CELL_DATA_CHAINS = pd.read_csv("./Data/Chains.csv")["FASTA"]
                    T_CELL_DATA_CD4 = pd.read_csv("./Data/CD4_EPITOPES.csv")["FASTA"]
                    T_CELL_DATA_CD8 = pd.read_csv("./Data/CD8_EPITOPES.csv")["FASTA"]

                    MATCH_VALUE = 45

                    LINEAR_EPI = []
                    LINEAR_SEQ = []
                    LINEAR_PERC = []

                    CHAIN_EPI = []
                    CHAIN_SEQ = []
                    CHAIN_PERC = []

                    CD4_EPI = []
                    CD4_SEQ = []
                    CD4_PERC = []

                    CD8_EPI = []
                    CD8_SEQ = []
                    CD8_PERC = []

                    for single in data_list:
                        for single_fasta in single["FASTA"]:
                            FASTA_CHUNKS = [single_fasta[i: j] for i in range(len(single_fasta)) for j in
                                            range(i + 1, len(single_fasta) + 1) if len(single_fasta[i:j]) == Epitope_window_var.get()]
                            for y in FASTA_CHUNKS:
                                for B_CELL_EPI_LINEAR in B_CELL_DATA_LINEAR:
                                    if len(y) == len(B_CELL_EPI_LINEAR):
                                        scored_B_CELL_LINEAR = aligner.score(y, B_CELL_EPI_LINEAR)
                                        aligned_perc_B_CELL_LINEAR = int(scored_B_CELL_LINEAR) / len(y) * 100
                                        if aligned_perc_B_CELL_LINEAR >= Epitope_identity_B_cell_var.get():
                                            # print(f"{int(aligned_perc_B_CELL_LINEAR)}%", y, B_CELL_EPI_LINEAR, "LINEAR EPITOPE")
                                            LINEAR_EPI.append(B_CELL_EPI_LINEAR)
                                            LINEAR_SEQ.append(y)
                                            LINEAR_PERC.append(aligned_perc_B_CELL_LINEAR)

                                for B_CELL_EPI_CHAIN in B_CELL_DATA_CHAINS:
                                    if len(y) == len(B_CELL_EPI_CHAIN):
                                        scored_B_CELL_CHAIN = aligner.score(y, B_CELL_EPI_CHAIN)
                                        aligned_perc_B_CELL_CHAIN = int(scored_B_CELL_CHAIN) / len(y) * 100
                                        if aligned_perc_B_CELL_CHAIN >= Epitope_identity_B_cell_var.get():
                                            # print(f"{int(aligned_perc_B_CELL_CHAIN)}%", y, B_CELL_EPI_CHAIN, "CHAIN EPITOPE")
                                            CHAIN_EPI.append(B_CELL_EPI_CHAIN)
                                            CHAIN_SEQ.append(y)
                                            CHAIN_PERC.append(aligned_perc_B_CELL_CHAIN)

                                for T_CELL_CD4 in T_CELL_DATA_CD4:
                                    if len(y) == len(T_CELL_CD4):
                                        scored_CD4 = aligner.score(y, T_CELL_CD4)
                                        aligned_perc_CD4 = int(scored_CD4) / len(y) * 100
                                        if aligned_perc_CD4 >= Epitope_identity_T_cell_var.get():
                                            # print(f"{int(aligned_perc_CD4)}%", y, T_CELL_CD4, "T CELL CD4")
                                            CD4_EPI.append(T_CELL_CD4)
                                            CD4_SEQ.append(y)
                                            CD4_PERC.append(aligned_perc_CD4)

                                for T_CELL_CD8 in T_CELL_DATA_CD8:
                                    if len(y) == len(T_CELL_CD8):
                                        scored_CD8 = aligner.score(y, T_CELL_CD8)
                                        aligned_perc_CD8 = int(scored_CD8) / len(y) * 100
                                        if aligned_perc_CD8 >= Epitope_identity_T_cell_var.get():
                                            # print(f"{int(aligned_perc_CD8)}%", y, T_CELL_CD8, "T CELL CD8")
                                            CD4_EPI.append(T_CELL_CD8)
                                            CD4_SEQ.append(y)
                                            CD4_PERC.append(aligned_perc_CD8)

                    LINEAR_EXTRACTED = pd.DataFrame({
                        "LINEAR PERC": LINEAR_PERC,
                        "LINEAR EPI": LINEAR_EPI,
                        "LINEAR SEQ": LINEAR_SEQ
                    }).drop_duplicates()

                    CHAIN_EXTRACTED = pd.DataFrame({
                        "CHAIN PERC": CHAIN_PERC,
                        "CHAIN EPI": CHAIN_EPI,
                        "CHAIN SEQ": CHAIN_SEQ
                    }).drop_duplicates()

                    CD4_EXTRACTED = pd.DataFrame({
                        "CD4 PERC": CD4_PERC,
                        "CD4 EPI": CD4_EPI,
                        "CD4 SEQ": CD4_SEQ
                    }).drop_duplicates()

                    CD8_EXTRACTED = pd.DataFrame({
                        "CD8 PERC": CD8_PERC,
                        "CD8 EPI": CD8_EPI,
                        "CD8 SEQ": CD8_SEQ
                    }).drop_duplicates()

                    print(LINEAR_EXTRACTED)
                    print(CHAIN_EXTRACTED)
                    print(CD4_EXTRACTED)
                    print(CD8_EXTRACTED)

                    LINEAR_EXTRACTED.to_csv(f"{SAVE_FILE_Var.get()}/LINEAR_EXTRACTED.csv")
                    CHAIN_EXTRACTED.to_csv(f"{SAVE_FILE_Var.get()}/CHAIN_EXTRACTED.csv")
                    CD4_EXTRACTED.to_csv(f"{SAVE_FILE_Var.get()}/CD4_EXTRACTED.csv")
                    CD8_EXTRACTED.to_csv(f"{SAVE_FILE_Var.get()}/CD8_EXTRACTED.csv")

                else:
                    messagebox.showerror("Error", "Select a VacSol Method")
            else:
                messagebox.showwarning("Warning", "No Folder Selected To Place Data")

        else:

            if SAVE_FILE_Var.get() != "":

                Job_Status.config(text="Job Status: Running")

                try:
                    Localization_Result = []

                    data = StringIO(LOC_TREE_var.get())
                    df = pd.read_csv(data, sep='\s+')[:-1]
                    df["LOCATION"] = df["BIT"]
                    df.drop(columns=["BIT", "ID"], inplace=True)

                    PROTEIN_ID_LIST_SECOND = []
                    PROTEIN_LOC = []
                    PROTEIN_BIT_SCORE = []
                    PROTEIN_UNIPROT = []

                    for single_id, loc, bit_score in zip(df["PROTEIN"], df['LOCATION'], df['SCORE']):
                        baseUrl = "http://www.uniprot.org/uniprot/"
                        currentUrl = baseUrl + single_id + ".fasta"
                        response = r.post(currentUrl)
                        cData = ''.join(response.text)
                        Seq = StringIO(cData)
                        pSeq = list(SeqIO.parse(Seq, 'fasta'))
                        for i in pSeq:
                            PROTEIN_UNIPROT.append(i.seq)
                            PROTEIN_ID_LIST_SECOND.append(single_id)
                            PROTEIN_LOC.append(loc)
                            PROTEIN_BIT_SCORE.append(bit_score)
                            print(single_id, " DONE")
                            sleep(5)

                    Localization_Data = pd.DataFrame({
                        "ID" : PROTEIN_ID_LIST_SECOND,
                        "LOC" : PROTEIN_LOC,
                        "SEQ" : PROTEIN_UNIPROT,
                        "SCORE" : PROTEIN_BIT_SCORE
                    })

                    PROTEIN_PROTEOME = []
                    PROTEIN_PERC = []
                    PROTEIN_EXTRACTED_LOCTREE = []
                    PROTEIN_EXTRACTED_LOCATION = []

                    for extracted_seq,extracted_loc in zip(Localization_Data["SEQ"], Localization_Data["LOC"]):
                        for single_Localization in data_list:
                            for single_Localization_fasta in single_Localization["FASTA"]:
                                aligner = Align.PairwiseAligner(match_score=1.0)
                                score = aligner.score(extracted_seq, single_Localization_fasta)
                                percentage = score / len(single_Localization_fasta) * 100
                                if Localization_identity_var.get() != 0:
                                    if percentage > Localization_identity_var.get():
                                        PROTEIN_PROTEOME.append(single_Localization_fasta)
                                        PROTEIN_PERC.append(f"{int(percentage)}%")
                                        PROTEIN_EXTRACTED_LOCTREE.append(extracted_seq)
                                        PROTEIN_EXTRACTED_LOCATION.append(extracted_loc)

                    FASTA_INCLUDED = pd.DataFrame({
                        "LOC": PROTEIN_EXTRACTED_LOCATION,
                        "PERCENTAGE": PROTEIN_PERC,
                        "FASTA_SUBMITTED": PROTEIN_PROTEOME,
                        "FASTA_UNIPROT": PROTEIN_EXTRACTED_LOCTREE
                    }).drop_duplicates()

                    if len(FASTA_INCLUDED) == 0:
                        messagebox.showerror("Error", "No Relevant Data found Localization")
                        sys.exit()
                    else:
                        Localization_Result.append(FASTA_INCLUDED)
                        FASTA_INCLUDED.to_csv(f"{SAVE_FILE_Var.get()}/FIRST_LOCALIZATION.csv")

                    sleep(2)

                    Non_Homolog_Result = []

                    HUMAN_GENOME = pd.read_csv("./Data/Human_Reference_Genome.csv")["FASTA"]
                    PERCENT_LIST = []
                    SCORE_LIST = []
                    QUERY_LIST = []
                    SUBJECT_LIST = []

                    for single_Non_Homolog in Localization_Result:
                        for single_Non_Homolog_fasta in single_Non_Homolog["FASTA_SUBMITTED"].drop_duplicates():
                            for x in HUMAN_GENOME:
                                if len(x) >= len(single_Non_Homolog_fasta):
                                    aligner = Align.PairwiseAligner(match_score=1.0)
                                    score = aligner.score(x, single_Non_Homolog_fasta)
                                    percentage = score / len(x) * 100
                                    if Non_Homologs_identity_var.get() != 0:
                                        if percentage < Non_Homologs_identity_var.get():
                                            PERCENT_LIST.append(f"{int(percentage)}%")
                                            SCORE_LIST.append(int(score))
                                            SUBJECT_LIST.append(single_Non_Homolog_fasta)
                                            QUERY_LIST.append(x)

                        FINAL_DATA = pd.DataFrame({
                            "SCORE": SCORE_LIST,
                            "PERCENTAGE": PERCENT_LIST,
                            "HUMAN_PROTEINS": QUERY_LIST,
                            "SUBMITTED_SEQUENCE": SUBJECT_LIST
                        }).drop_duplicates()

                        if len(FINAL_DATA) == 0:
                            messagebox.showerror("Error", "No Relevant Data found Non Homologs")
                            sys.exit()
                        else:
                            Non_Homolog_Result.append(FINAL_DATA)
                            FINAL_DATA.to_csv(f"{SAVE_FILE_Var.get()}/NON_HOMOLOGS.csv")

                    sleep(2)

                    VFDB_DATA = pd.read_csv("./Data/VFDB.csv")

                    LIST_NON_HOMOLOG_SEQUENCE = []
                    LIST_VFDB_SEQUENCE = []
                    LIST_PERCENTAGE = []

                    print(VFDB_identity_var.get())

                    aligner = Align.PairwiseAligner()

                    for single_VFDB in Non_Homolog_Result:
                        for single_VFDB_fasta in single_VFDB["SUBMITTED_SEQUENCE"].drop_duplicates():
                            for single_vfdb_seq, title in zip(VFDB_DATA["FASTA"], VFDB_DATA["TITLE"]):
                                score = aligner.score(single_VFDB_fasta, single_vfdb_seq)
                                percentage = score / len(single_VFDB_fasta) * 100
                                if VFDB_identity_var.get() != 0:
                                    if percentage > VFDB_identity_var.get():
                                        LIST_PERCENTAGE.append(f"{int(percentage)}%")
                                        LIST_VFDB_SEQUENCE.append(single_vfdb_seq)
                                        LIST_NON_HOMOLOG_SEQUENCE.append(single_VFDB_fasta)

                    DATA_FRAME = pd.DataFrame({
                        "PERCENTAGE": LIST_PERCENTAGE,
                        "QUERY": LIST_NON_HOMOLOG_SEQUENCE,
                        "SUBJECT": LIST_VFDB_SEQUENCE
                    }).drop_duplicates()

                    if len(DATA_FRAME) == 0:
                        messagebox.showerror("Error", "No Relevant Data found")
                        sys.exit()
                    else:
                        DATA_FRAME.to_csv(f"{SAVE_FILE_Var.get()}/VFDB.csv")

                    sleep(2)

                    B_CELL_DATA_LINEAR = pd.read_csv("./Data/Linear.csv")["FASTA"]
                    B_CELL_DATA_CHAINS = pd.read_csv("./Data/Chains.csv")["FASTA"]
                    T_CELL_DATA_CD4 = pd.read_csv("./Data/CD4_EPITOPES.csv")["FASTA"]
                    T_CELL_DATA_CD8 = pd.read_csv("./Data/CD8_EPITOPES.csv")["FASTA"]


                    LINEAR_EPI = []
                    LINEAR_SEQ = []
                    LINEAR_PERC = []

                    CHAIN_EPI = []
                    CHAIN_SEQ = []
                    CHAIN_PERC = []

                    CD4_EPI = []
                    CD4_SEQ = []
                    CD4_PERC = []

                    CD8_EPI = []
                    CD8_SEQ = []
                    CD8_PERC = []

                    for single in data_list:
                        for single_fasta in single["FASTA"]:
                            FASTA_CHUNKS = [single_fasta[i: j] for i in range(len(single_fasta)) for j in
                                            range(i + 1, len(single_fasta) + 1) if
                                            len(single_fasta[i:j]) == Epitope_window_var.get()]
                            for y in FASTA_CHUNKS:
                                for B_CELL_EPI_LINEAR in B_CELL_DATA_LINEAR:
                                    if len(y) == len(B_CELL_EPI_LINEAR):
                                        scored_B_CELL_LINEAR = aligner.score(y, B_CELL_EPI_LINEAR)
                                        aligned_perc_B_CELL_LINEAR = int(scored_B_CELL_LINEAR) / len(y) * 100
                                        if aligned_perc_B_CELL_LINEAR >= Epitope_identity_B_cell_var.get():
                                            # print(f"{int(aligned_perc_B_CELL_LINEAR)}%", y, B_CELL_EPI_LINEAR, "LINEAR EPITOPE")
                                            LINEAR_EPI.append(B_CELL_EPI_LINEAR)
                                            LINEAR_SEQ.append(y)
                                            LINEAR_PERC.append(aligned_perc_B_CELL_LINEAR)

                                for B_CELL_EPI_CHAIN in B_CELL_DATA_CHAINS:
                                    if len(y) == len(B_CELL_EPI_CHAIN):
                                        scored_B_CELL_CHAIN = aligner.score(y, B_CELL_EPI_CHAIN)
                                        aligned_perc_B_CELL_CHAIN = int(scored_B_CELL_CHAIN) / len(y) * 100
                                        if aligned_perc_B_CELL_CHAIN >= Epitope_identity_B_cell_var.get():
                                            # print(f"{int(aligned_perc_B_CELL_CHAIN)}%", y, B_CELL_EPI_CHAIN, "CHAIN EPITOPE")
                                            CHAIN_EPI.append(B_CELL_EPI_CHAIN)
                                            CHAIN_SEQ.append(y)
                                            CHAIN_PERC.append(aligned_perc_B_CELL_CHAIN)

                                for T_CELL_CD4 in T_CELL_DATA_CD4:
                                    if len(y) == len(T_CELL_CD4):
                                        scored_CD4 = aligner.score(y, T_CELL_CD4)
                                        aligned_perc_CD4 = int(scored_CD4) / len(y) * 100
                                        if aligned_perc_CD4 >= Epitope_identity_T_cell_var.get():
                                            # print(f"{int(aligned_perc_CD4)}%", y, T_CELL_CD4, "T CELL CD4")
                                            CD4_EPI.append(T_CELL_CD4)
                                            CD4_SEQ.append(y)
                                            CD4_PERC.append(aligned_perc_CD4)

                                for T_CELL_CD8 in T_CELL_DATA_CD8:
                                    if len(y) == len(T_CELL_CD8):
                                        scored_CD8 = aligner.score(y, T_CELL_CD8)
                                        aligned_perc_CD8 = int(scored_CD8) / len(y) * 100
                                        if aligned_perc_CD8 >= Epitope_identity_T_cell_var.get():
                                            # print(f"{int(aligned_perc_CD8)}%", y, T_CELL_CD8, "T CELL CD8")
                                            CD4_EPI.append(T_CELL_CD8)
                                            CD4_SEQ.append(y)
                                            CD4_PERC.append(aligned_perc_CD8)

                    LINEAR_EXTRACTED = pd.DataFrame({
                        "LINEAR PERC": LINEAR_PERC,
                        "LINEAR EPI": LINEAR_EPI,
                        "LINEAR SEQ": LINEAR_SEQ
                    })

                    CHAIN_EXTRACTED = pd.DataFrame({
                        "CHAIN PERC": CHAIN_PERC,
                        "CHAIN EPI": CHAIN_EPI,
                        "CHAIN SEQ": CHAIN_SEQ
                    })

                    CD4_EXTRACTED = pd.DataFrame({
                        "CD4 PERC": CD4_PERC,
                        "CD4 EPI": CD4_EPI,
                        "CD4 SEQ": CD4_SEQ
                    })

                    CD8_EXTRACTED = pd.DataFrame({
                        "CD8 PERC": CD8_PERC,
                        "CD8 EPI": CD8_EPI,
                        "CD8 SEQ": CD8_SEQ
                    })

                    sleep(2)
                    messagebox.showinfo("Completed", "Epitope Mapping Completed")

                    LINEAR_EXTRACTED.to_csv(f"{SAVE_FILE_Var.get()}/LINEAR.csv")
                    CHAIN_EXTRACTED.to_csv(f"{SAVE_FILE_Var.get()}/CHAIN.csv")
                    CD4_EXTRACTED.to_csv(f"{SAVE_FILE_Var.get()}/CD4.csv")
                    CD8_EXTRACTED.to_csv(f"{SAVE_FILE_Var.get()}/CD8.csv")

                    Job_Status.config(text="Job Status: Not Running", fg=First_Color)

                except:
                    messagebox.showerror("Error", "Make Sure You Are Connected To The Internet\n and Choose The Parameters Correctly")

            else:
                messagebox.showwarning("Warning", "No Folder Selected To Place Data")


Save_file_btn = Button(Submission_method_frame, text="SAVE FILE", command=SAVE_FILE, height=2, width=20, bd=0, background=Buttons_Colors, font = (font_family,12), fg=Button_Text)
Save_file_btn.place(relx=0.65, rely=0.625)

Fasta_Entry_Btn = Button(Submission_method_frame, text="SUBMIT", command=FASTA_SUBMIT_BTN, height=2, width=20, bd=0, background=Buttons_Colors, font = (font_family,12), fg=Button_Text)
Fasta_Entry_Btn.place(relx=0.65, rely=0.85)

Must_Evalutate_Frame = Frame(window, width=1200, height=50, bg=Second_Color, bd=0)
Must_Evalutate_Frame.grid(row=1, column=0)

Must_Evalutate_var = StringVar()
Must_Evalutate = Checkbutton(Must_Evalutate_Frame, variable=Must_Evalutate_var,text="Must Evaluate" ,onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, command=MUST_EVALUATE, font=(font_family, 11))
Must_Evalutate.deselect()
Must_Evalutate.place(relx=0.015, rely=0.4)

Must_Evalutate_Details = Label(Must_Evalutate_Frame, text="If Must Evaluate is Checked Do not Uncheck Any Method",bg=Second_Color, fg="Red", font=(font_family, 10))
Must_Evalutate_Details.place(relx=0.12, rely=0.47)

Example_Btn = Button(Submission_method_frame, text="Press Here For Example Sequence",bg=Second_Color, fg=Buttons_Colors, font=(font_family, 10), bd=0,background=First_Color, command=FILL_EXAMPLE)
Example_Btn.place(relx=0.45, rely=0.3)

Localization_Details = Label(Must_Evalutate_Frame, text="Choose From Greater Than Given Values Reliability Score and Percentage",bg=Second_Color, fg="Red", font=(font_family, 10))
Localization_Details.place(relx=0.41, rely=0.47)

Localization_Frame = Frame(window, width=1200, height=80, bg=Second_Color, bd=0)
Localization_Frame.grid(row=2, column=0)

Localization_var = StringVar()
Localization = Checkbutton(Localization_Frame, variable=Localization_var,text="Localization", onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11))
Localization.deselect()
Localization.place(relx=0.015, rely=0.38)


Localization_file_name = StringVar()

def get_single_strain(single):
    if single[0][0].isupper() == True:
        word = single[0] +" "+ single[1]
        files_list = os.listdir("./Data/Bacteria")
        for i in files_list:
            x = i.replace("_", ".").split(".")[1:][:2]
            if x[0] + " " + x[1] == word:
                return i
    else:
        word = single[1] + " " + single[0]
        files_list = os.listdir("./Data/Bacteria")
        for i in files_list:
            x = i.replace("_", ".").split(".")[1:][:2]
            if x[0] + " " + x[1] == word:
                return i

# get_names() function will fetch all of the names of the bacteria present in the directory "Bacteria"
def get_names():
    files_list = os.listdir("./Data/Bacteria")
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

LOC_TREE_var = StringVar()

def Localization_select(event):

    final_list = []
    genus_selected = Localization_combo_box_1_var.get()
    second_name_filtered = []

    for x in genus:
        if genus_selected in x["First_Name"]:
            second_name_filtered.append(x["Second_Name"])
            final_list.append(x["First_Name"])

    Localization_combo_box_2.config(values=second_name_filtered)

    for y in second_name_filtered:
        if y == Localization_combo_box_2.get():
            final_list.append(y)

    single = list(set(final_list))
    if len(single) == 2:
        file_name = get_single_strain(single)

        LOCATION_LIST = []
        BIT_SCORE_LIST = []
        PROTEIN_ID_LIST = []

        with open("./Data/Bacteria/" + file_name, "r") as fasta_file:
            for line in fasta_file.readlines():
                without_hastag = re.findall("^#", line)
                if not without_hastag:
                    list_sequence = line.split()
                    # print(list_sequence[2:3][0])
                    if list_sequence[2:3][0] == "inner":
                        LOCATION_LIST.append(list_sequence[2:3][0] + " membrance")
                    elif list_sequence[2:3][0] == "outer":
                        LOCATION_LIST.append(list_sequence[2:3][0] + " membrance")
                    else:
                        LOCATION_LIST.append(list_sequence[2:3][0])

                    BIT_SCORE_LIST.append(int(list_sequence[1]))
                    PROTEIN_ID_LIST.append(list_sequence[0])

        LOCATION_DATA = pd.DataFrame({
            "LOCATION": LOCATION_LIST,
            "BIT SCORE": BIT_SCORE_LIST,
            "PROTEIN ID": [PROTEIN_ID.split("|")[1] for PROTEIN_ID in PROTEIN_ID_LIST],
            # "PROTEIN FASTA" : PROTEIN_FAST_LIST
        }).reset_index(drop=True)

        if Reliability_var.get() != 0:

            #BIT_SCORE_MEAN = int(np.mean(LOCATION_DATA["BIT SCORE"]))

            mean_filtered = LOCATION_DATA[LOCATION_DATA["BIT SCORE"] > Reliability_var.get()]
            remove_list = ["secreted", "outer membrane", "fimbrium"]
            NEW_DATA_FRAME = mean_filtered[mean_filtered["LOCATION"].isin(remove_list)].reset_index(drop=True)

            if len(NEW_DATA_FRAME) == 0:
                messagebox.showwarning("Warning", "No Data Available")
            else:
                Relaibility_based_proteins.config(text=str((len(NEW_DATA_FRAME["PROTEIN ID"]))))
                LOC_TREE_var.set(NEW_DATA_FRAME)



        Localization_file_name.set(file_name)
    elif len(single) == 1:
        Localization_combo_box_2.current(0)
    else:
        messagebox.showwarning("Warning","Select Strain")

Localization_combo_box_1_var = StringVar()
Localization_combo_box_1 = AutocompleteCombobox(Localization_Frame, completevalues=short, textvariable=Localization_combo_box_1_var)
Localization_combo_box_1.current(0)
Localization_combo_box_1.bind("<<ComboboxSelected>>", Localization_select)
Localization_combo_box_1.place(relx=0.12, rely=0.45)

Localization_combo_box_Label = Label(Localization_Frame, text="Select Bacteria Genus", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Localization_combo_box_Label.place(relx=0.12, rely=0.1)


Localization_combo_box_2_var = StringVar()
Localization_combo_box_2 = ttk.Combobox(Localization_Frame, values=[""], textvariable=Localization_combo_box_2_var)
Localization_combo_box_2.bind("<<ComboboxSelected>>", Localization_select)
Localization_combo_box_2.place(relx=0.265, rely=0.45)

Localization_combo_box_Label_Second = Label(Localization_Frame, text="Select Bacteria Strain", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Localization_combo_box_Label_Second.place(relx=0.265, rely=0.1)

Reliability_var = IntVar()
Reliability = ttk.Combobox(Localization_Frame, values=[reliable for reliable in range(10,100, 5)], textvariable=Reliability_var)
Reliability.bind("<<ComboboxSelected>>", Localization_select)
Reliability.place(relx=0.41, rely=0.45)

Reliability_Label = Label(Localization_Frame, text="Select Reliability Score", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Reliability_Label.place(relx=0.41, rely=0.1)

Localization_identity_var = IntVar()
Localization_identity = ttk.Combobox(Localization_Frame, values=[Localization_identity_values for Localization_identity_values in range(10,110, 10)], textvariable=Localization_identity_var)
Localization_identity.bind("<<ComboboxSelected>>", Localization_select)
Localization_identity.place(relx=0.555, rely=0.45)

Localization_identity_Label = Label(Localization_Frame, text="Select Identity Percentage", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Localization_identity_Label.place(relx=0.555, rely=0.1)

Relaibility_based_proteins = Label(Localization_Frame, text="0", bg=Second_Color, font=((8)), fg=Label_Color)
Relaibility_based_proteins.place(relx=0.7, rely=0.1)

Relaibility_based_proteins_Label = Label(Localization_Frame, text="Reliable Proteins Found", bg=Second_Color, fg=Label_Color)
Relaibility_based_proteins_Label.place(relx=0.7, rely=0.45)

Non_Homologs_Frame = Frame(window, width=1200, height=80, bg=Second_Color, bd=0)
Non_Homologs_Frame.grid(row=3, column=0)

Non_Homologs_Details = Label(Non_Homologs_Frame, text="Choose Non Homologs Percentage Less Than Given Values",bg=Second_Color, fg="Red", font=(font_family, 10))
Non_Homologs_Details.place(relx=0.35, rely=0.35)

Non_Homologs_var = StringVar()
Non_Homologs_check_btn = Checkbutton(Non_Homologs_Frame, variable=Non_Homologs_var,text="Non Host Homologs" ,onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11))
Non_Homologs_check_btn.deselect()
Non_Homologs_check_btn.place(relx=0.015, rely=0.4)

Non_Homologs_identity_var = IntVar()
Non_Homologs_identity = ttk.Combobox(Non_Homologs_Frame, values=[Non_Homologs_identity_values for Non_Homologs_identity_values in range(1,21)], textvariable=Non_Homologs_identity_var)
Non_Homologs_identity.place(relx=0.18, rely=0.4)

Non_Homologs_identity_Label = Label(Non_Homologs_Frame, text="Select Identity Percentage", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Non_Homologs_identity_Label.place(relx=0.18, rely=0.015)

VFDB_Frame = Frame(window, width=1200, height=80, bg=Second_Color, bd=0)
VFDB_Frame.grid(row=4, column=0)

VFDB_Details = Label(VFDB_Frame, text="Choose Virulence Factors Percentage Greater Than Given Values",bg=Second_Color, fg="Red", font=(font_family, 10))
VFDB_Details.place(relx=0.35, rely=0.35)

VFDB_var = StringVar()
VFDB_check_btn = Checkbutton(VFDB_Frame, variable=VFDB_var,text="Virulence Factors",onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11))
VFDB_check_btn.deselect()
VFDB_check_btn.place(relx=0.015, rely=0.4)

VFDB_identity_var = IntVar()
VFDB_identity = ttk.Combobox(VFDB_Frame, values=[VFDB_identity_values for VFDB_identity_values in range(10,110, 10)], textvariable=VFDB_identity_var)
VFDB_identity.place(relx=0.18, rely=0.4)

VFDB_identity_Label = Label(VFDB_Frame, text="Select Identity Percentage", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
VFDB_identity_Label.place(relx=0.18, rely=0.015)

Epitope_Frame = Frame(window, width=1200, height=80, bg=Second_Color, bd=0)
Epitope_Frame.grid(row=5, column=0)

Epitope_Details = Label(Epitope_Frame, text="Choose Epitope\nPercentage Greater Than Given Values",bg=Second_Color, fg="Red", font=(font_family, 10))
Epitope_Details.place(relx=0.615, rely=0.2)

Epitope_var = StringVar()
Epitope_check_btn = Checkbutton(Epitope_Frame, variable=Epitope_var,text="Epitope Mapping" , onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11))
Epitope_check_btn.deselect()
Epitope_check_btn.place(relx=0.015, rely=0.4)

Epitope_identity_Label_B_cell = Label(Epitope_Frame, text="Select B CELL Identity", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_identity_Label_B_cell.place(relx=0.325, rely=0.015)

Epitope_identity_Label_T_cell = Label(Epitope_Frame, text="Select T CELL Identity", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_identity_Label_T_cell.place(relx=0.4725, rely=0.015)

Epitope_window_Label = Label(Epitope_Frame, text="Select Window Length", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_window_Label.place(relx=0.18, rely=0.015)

Epitope_window_var = IntVar()
Epitope_window = ttk.Combobox(Epitope_Frame, values=[Epitope_window_values for Epitope_window_values in range(12,21)], textvariable=Epitope_window_var)
Epitope_window.place(relx=0.18, rely=0.4)

Epitope_identity_B_cell_var = IntVar()
Epitope_identity_B_cell = ttk.Combobox(Epitope_Frame, values=[Epitope_identity_B_cell_values for Epitope_identity_B_cell_values in range(10,110, 10)], textvariable=Epitope_identity_B_cell_var)
Epitope_identity_B_cell.place(relx=0.325, rely=0.4)

Epitope_identity_T_cell_var = IntVar()
Epitope_identity_T_cell = ttk.Combobox(Epitope_Frame, values=[Epitope_identity_T_cell_values for Epitope_identity_T_cell_values in range(10,110, 10)], textvariable=Epitope_identity_T_cell_var)
Epitope_identity_T_cell.place(relx=0.4725, rely=0.4)


window.mainloop()

