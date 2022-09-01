import pandas as pd
from Bio import Align

aligner = Align.PairwiseAligner()

VFDB_DATA = pd.read_csv("./Helicobacter/VFDB.csv")['QUERY'].drop_duplicates()
print(VFDB_DATA)

B_CELL_DATA = pd.read_csv("./B-CELL/B_CELL.csv")["Sequence"]
#print(B_CELL_DATA)

for x in VFDB_DATA:
    n = 16
    FASTA_CHUNKS = [x[i:i + n] for i in range(0, len(x), n)]
    for y in FASTA_CHUNKS:
        for i in B_CELL_DATA:
            if len(y) == len(i):
                score = aligner.score(y, i)
                percentage = score / len(y)
                if percentage > 0.4:
                    print(score, len(y),len(i), percentage, y, i)
