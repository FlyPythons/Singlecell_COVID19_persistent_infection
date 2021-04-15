#!/usr/bin/env python
import sys
import pandas as pd
import gzip
import numpy

meta_file = sys.argv[1]     # sample information file
barcode_file = sys.argv[2]  # barcodes.tsv.gz download from GSE158055
matrix_file = sys.argv[3]   # matrix.mtx download from GSE158055

meta = pd.read_csv(meta_file, sep="\t")
sample_dict = {}
samples = {"S-HC008":"HD01","S-HC009":"HD02", "S-HC010":"HD03", "S-HC011":"HD04", "S-HC012": "HD05",
    "S-HC018-2": "HD06", "S-HC019-2": "HD07", "S-HC020-2":"HD08", "S-HC001": "HD09", "S-HC002": "HD10",
    "S-HC003" :"HD11", "S-HC004": "HD12", "S-HC005": "HD13", "S-HC006": "HD14", "S-HC007": "HD15",
    "S-HC021" :"HD16", "S-HC022": "HD17", "S-HC023": "HD18", "S-HC024": "HD19", "S-HC025": "HD20",
    "S-M043-1": "SD01", "S-M048": "SD02", "S-S091": "SD03", "S-M019": "SD04", "S-M064": "SD05",
    "S-M065": "SD06", "S-M018": "SD07", "S-M035-1": "SD08", "S-M023": "SD09", "S-M061-2": "SD10",
    "S-M044-2": "SD11", "S-M017": "SD12", "S-M024": "SD13", "S-M034": "SD14", "S-M058-2": "SD15",
    "S-M041-2": "SD16",
    "S-S073-2": "LD01", "S-S075-2": "LD02"
}

for i in meta[["Sample name", "characteristics: Meta sample"]].iterrows():
    if i[1][0] in samples:
        sample_dict[i[1][1]] = i[1][0]

line_dict = {}
groups_fh = open("groups.tsv", "w")
bc_fh = open("barcodes.tsv", "w")
j = 1
max_line = 0
for n, line in enumerate(gzip.open(barcode_file, "rb").readlines()):
    line = line.decode().strip()
    parts = line.split("_")
    #print(parts)
    _sample = "_".join(parts[:-1])
    #print(_sample)
    if _sample in sample_dict:
        max_line = n+1
        line_dict[str(n+1)] = str(j)
        j += 1
        groups_fh.write("%s\n" % (samples[sample_dict[_sample]][:2]))
        bc_fh.write("%s\n" % (samples[sample_dict[_sample]]+"_"+parts[-1]))
bc_fh.close()
groups_fh.close()

mtx_fh = open("matrix.mtx", "w")
datas = []
counts = [27943, len(line_dict), 0]
print(max_line)
for line in open(matrix_file):
    line = line.strip()
    parts = line.split()
    if parts[0][0] == "%" or parts[0]=="1462702":
        continue
    if int(parts[0]) > max_line:
        break
    if parts[0] in line_dict:
        bc, gene, count = parts
        datas.append(" ".join([gene, line_dict[bc], count]))
        counts[2] += 1

mtx_fh.write('''\
%%%%MatrixMarket matrix coordinate real general
%%
%s %s %s
''' %(counts[0], counts[1], counts[2]))
mtx_fh.write("\n".join(datas)+"\n")
