#!/usr/bin/env python3

import argparse
import tempfile
import re
import os
import numpy as np
from pathlib import Path
from biotite.sequence.phylo import neighbor_joining
from collections import defaultdict

from magnumopus.needleman_wunsch import needleman_wunsch
from magnumopus.ispcr import ispcr
from magnumopus.mapping import minimap
import magnumopus.phylo_inference as ph

#usage:
#magop.py [-a ASSEMBLY [ASSEMBLY ...]] -p PRIMERS [-r READS [READS ...]] [-s REF_SEQS]


def parse_arguments():
    map = argparse.ArgumentParser(description="Extracts 16S sequences from different data types(reads and assemblies) and returns a newick phylogenetic tree")

    map.add_argument("-s","--REF_SEQS", type=str, help="path to a fasta file containing reference sequences", required=False)
    map.add_argument("-p","--PRIMERS", type=str, help="a path to a fasta file containing the primers", required=True)
    map.add_argument("-a","--ASSEMBLY",nargs="+", type=str, help="a space separated list of paths to fasta files containing assemblies or path to a directory containing your assemblies", required=False)
    map.add_argument("-r","--READS", nargs="+", type=str, help="a space separated list of paths to fastq files containing illumina reads", required=False)
    return map.parse_args()
    


def main():
    m = parse_arguments()
    amps = []
    seq_names = []

    if m.ASSEMBLY:
        assemblies = ph.parse_assembly(m.ASSEMBLY)
        assembly_amps, assembly_names = ph.get_assem_amps(assemblies, m.PRIMERS)
        amps.extend(assembly_amps)
        seq_names.extend(assembly_names)

    if m.READS==True and m.REF_SEQS==False:
        raise ValueError(f"Missing reference sequences")

    if m.READS:
        paired_reads = ph.parse_reads(m.READS) #list of tuples of strs
        read_amps, read_names = ph.get_read_amps(paired_reads, m.PRIMERS, m.REF_SEQS)
        amps.extend(read_amps)
        seq_names.extend(read_names)

    if m.REF_SEQS:
        ref_amps, ref_names = ph.get_ref_amps(m.REF_SEQS, m.PRIMERS)
        amps.extend(ref_amps)
        seq_names.extend(ref_names)

    # make tree 
    dist_matrix = ph.get_pairwise_distance(amps)
    my_tree = ph.make_tree(dist_matrix, seq_names)
    print(my_tree)


if __name__ == "__main__":
    main()