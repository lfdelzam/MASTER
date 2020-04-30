#!/usr/bin/env python3

import os
import argparse

usage = 'python extract_genes_with_pfam_best_hit.py [options]'
description = 'This program extracts genes with pfam hit (best hit). Creates genes.fna, proteins.faa and annotation.gff files'

parser = argparse.ArgumentParser(description=description, usage=usage)
parser.add_argument('-i', dest='inf', help='input file, --tblout from hmmsearch', required=True)
parser.add_argument('-a', dest='a', help='input file, <Genome_name>.gff', required=True)
parser.add_argument('-p', dest='p', help='input file, <Genome_name>.faa', required=True)
parser.add_argument('-n', dest='n', help='input file, <Genome_name>.fna', required=True)
parser.add_argument('-f', dest='f', help='Genome contig sequences file, <Genome_name>.fa', required=True)
parser.add_argument('-o', dest='o', help='dataset name', required=True)
parser.add_argument('-g', dest='g', help='genome name', required=True)

args = parser.parse_args()

best_hit = {}
scores = {}

# functions

def extracting(filein, fileout):

    with open(filein, "r") as finp, open(fileout, "w") as foutp:
        non_first_line = False
        for line in finp:
            line = line.rstrip()
            if line.startswith(">"):
                copy = False
                id = line.split()[0][1:]
                remove = line.split()[0]
                rest = line.replace(remove, "")
                if id in best_hit:
                    v = best_hit[id]
                    if non_first_line:
                        print("", file=foutp)
                    print(">{} Pfam_hit_{}_{}_Score_{} Prodigal_description{}".format(id, v[1], v[2], v[3], rest), file=foutp)
                    copy = True
            else:
                if copy:
                    print(line, end="", file=foutp)
                    non_first_line = True
        print("", file=foutp)


def extracting_gff(filein, fileout, genome):
    selected_contigs = []
    with open(filein, "r") as finp, open(genome, "r") as fing, open(fileout, "w") as foutp:
        print("##gff-version 3", file=foutp)
        for line in finp:
            line = line.rstrip()
            if not line.startswith("#"):
                    target = line.split()[0]
                    idline = line.split()[8]
                    pre1 = idline.split(";")[0]
                    idgff = pre1.replace("ID=", "")
                    target_name = target+"_"+idgff.split("_")[1]
                    if target_name in best_hit:
                        if best_hit[target_name][0] == idgff:
                                selected_contigs.append(target)
                                print(line, file=foutp)
#                            line = line.split()
#                            text = "\t".join(line[1:])
#                            print(target_name, text, sep="\t", file = foutp)

        print("##FASTA", file=foutp)
        Fst_time = True
        for line in fing:
            line = line.rstrip()
            if line.startswith(">"):
                copy = False
                id = line.split()[0][1:]
                if id in selected_contigs:
                    copy = True
                    if not Fst_time:
                        print("", file=foutp)
#                    print(line, file = foutp)
                    print(">"+id.rstrip(), file=foutp)
                    Fst_time = False
            else:
                if copy:
                    print(line, end="", file=foutp)
        print("", file=foutp)

# MAIN

with open(args.inf, "r") as fin:
    for line in fin:
        line = line.rstrip()
        if not line.startswith("#"):
            line = line.split()
            target_name = line[0]
            query_name = line[2]
            query_accession = line[3]
            score = line[5]
            id = line[25].split(";")[0].replace("ID=", "")
            if target_name in scores:
                if float(score) > float(scores[target_name]):
                    scores[target_name] = score
                    best_hit[target_name] = (id, query_name, query_accession, score)
            else:
                scores[target_name] = score
                best_hit[target_name] = (id, query_name, query_accession, score)

# Creating output directory

outdir = os.path.join("Gene_calling", args.o, args.g, "best_hit_pfam")
outdirgff = os.path.join("GFF_files", args.o)
if not os.path.exists(outdir):
    os.mkdir(outdir)

if not os.path.exists(outdirgff):
    os.mkdir(outdirgff)

# creating outputfiles in the New directory

file_p = os.path.join(outdir, args.g+".faa")
file_n = os.path.join(outdir, args.g+".fna")
file_g = os.path.join(outdirgff, args.g+".gff")

extracting(args.p, file_p)
extracting(args.n, file_n)
extracting_gff(args.a, file_g, args.f)
