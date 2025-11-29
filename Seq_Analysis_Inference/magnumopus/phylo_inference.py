import tempfile
import re
import os
import numpy as np
from pathlib import Path
from biotite.sequence.phylo import neighbor_joining
from collections import defaultdict

from .needleman_wunsch import needleman_wunsch
from .ispcr import ispcr
from .mapping import minimap


                            #######################################################
                            ##### phylogenetic inference and main functions #######
                            #######################################################

def parse_assembly(assembly_input: list[str]) -> list[str]:
    """
    returns list of paths to assembly files
    """
    assemblies = []
    for item in assembly_input:
        if " " in item or "*" not in item: #file(s)
            assemblies.extend(item.split())
        else:
            p = Path(item)
            if p.is_dir(): #directory
                assemblies.extend(str(f) for f in p.iterdir() if f.is_file())
            else: #glob
                assemblies.extend(str(f) for f in Path().glob(item))
        
    return assemblies


def get_assem_amps(assemblies: list[str], primers:str) -> tuple[list]:
    """
    finds one amplicon per assembly and returns amplicons and assembly names 
    """
    assembly_amps = []
    assembly_names = []
    for assem in assemblies:
        amp = ispcr(primers, assem, 2000)
        assembly_amps.append(amp)

        with open(assem) as f:
            first_line = f.readline().strip()
            if first_line.startswith(">"):
                a_name = first_line[1:].split()[0]
                assembly_names.append(a_name)

    return assembly_amps, assembly_names


def parse_reads(reads_input: list) -> list[tuple[str]]:
    """
    ensures each read has a mate pair and places them together in a list
    """
    pattern = re.compile(r"(.+?)_([12])\.fastq$")
    pairs = defaultdict(dict)

    for f in reads_input:
        match = pattern.match(f)
        sample_id, read_num = match.groups()
        pairs[sample_id][read_num] = f

    paired_list = [] #list of tuples of strs
    for sample, reads in pairs.items():
        if "1" not in reads or "2" not in reads:
            raise ValueError(f"Missing mate pair for sample {sample}: {reads}")
        paired_list.append((reads["1"], reads["2"]))

    return paired_list


def get_read_amps(paired_reads: list[tuple[str]], primers:str, ref_seqs:str) -> tuple[list[str]]:
    """
    Maps reads to reference to find consensus to produce amplicons
    Finds one amplicon per consensus and returns amplicons and read names
    """
    read_amps = []
    read_names = []
    con = []
    for preads in paired_reads:
        consensus = minimap(ref_seqs, preads[0], preads[1])
        con.append(consensus)
        
        #find read name
        pattern = re.compile(r'^(?P<sample>.+?)_[12]\.fastq$')   
        name = Path(preads[0]).name                     
        m = pattern.match(name)
        if m:
            r_name = m.group("sample")
            read_names.append(r_name)

        # write consensus to temp file
        with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as temp:
            temp.write(f">{r_name}_consensus\n{consensus}")
            con_path = temp.name
        
        #find amplicon
        amp = ispcr(primers, con_path, 2000)
        read_amps.append(amp)

    return read_amps, read_names


def get_ref_amps(references:str, primers:str)-> tuple[list[str]]:
    """
    parses references from given fasta files and returns list of amplicons and names for each reference
    """
    pattern = re.compile(r"^>([^:]+)")
    ref_names = []
    ref_amps = []

    current_name = None
    seq_lines = []

    with open(references) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_name is not None:
                # If we already collected a sequence, process it
                    # Write to temp FASTA
                    tmp = tempfile.NamedTemporaryFile(
                        mode="w", delete=False, suffix=".fasta"
                    )
                    tmp.write(f">{current_name}\n{''.join(seq_lines)}\n")
                    tmp_ref = tmp.name
                    tmp.close()

                    # Run and store ispcr
                    amp = ispcr(primers, tmp_ref, 2000)
                    ref_names.append(current_name)
                    ref_amps.append(amp)

                    # Delete temporary file
                    os.remove(tmp_ref)

                # Start a new sequence
                m = pattern.match(line)
                current_name = m.group(1) if m else None
                seq_lines = []

            else:
                seq_lines.append(line)

        # Process the last sequence
        if current_name is not None:

            tmp = tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".fasta"
            )
            tmp.write(f">{current_name}\n{''.join(seq_lines)}\n")
            tmp_ref = tmp.name
            tmp.close()

            amp = ispcr(primers, tmp_ref, 2000)
            ref_names.append(current_name)
            ref_amps.append(amp)

            os.remove(tmp_ref)

    return ref_amps, ref_names


def rev_comp(amp1: str) -> str:
    """
    finds reverse complement
    """
    seq = amp1.upper()          
    #make complementary strand
    complement = ""
    for nuc in seq:
        if nuc == 'A':
            complement += 'T'
        elif nuc == 'T':
            complement += 'A'
        elif nuc == 'C':
            complement += 'G'
        elif nuc == 'G':
            complement += 'C'
        else:
            complement += 'N'  # for unknown bases
                
    #reverse complement
    rc_amp1 = complement[::-1] 
    
    return rc_amp1


def get_pairwise_distance(amps: list[str]) -> list[list[int]]:
    """
    aligns amplicons and finds best alignment to find pairwise distances
    """
    n = len(amps)
    # initialize n x n matrix of zeros
    dist_matrix = [[0] * n for _ in range(n)]

    # compute pairwise distances
    for i in range(n):
        for j in range(i + 1, n):  # only compute upper triangle
            amp1 = ''.join(line.strip() for line in amps[i].splitlines() if not line.startswith('>'))
            amp2 = ''.join(line.strip() for line in amps[j].splitlines() if not line.startswith('>'))

            alignments1, score1, pairwise_distance1 = needleman_wunsch(amp1, amp2, 1, -1, -1)

            rc_amp2 = rev_comp(amp2)
            alignments2, score2, pairwise_distance2 = needleman_wunsch(amp1, rc_amp2, 1, -1, -1)

            if score1 > score2:
                dist_matrix[i][j] = pairwise_distance1
                dist_matrix[j][i] = pairwise_distance1 # symmetric matrix
            else:
                dist_matrix[i][j] = pairwise_distance2 
                dist_matrix[j][i] = pairwise_distance2 


    return dist_matrix


def make_tree(pd_matrix:list[list[int]], seq_names: list[str]) -> str:
    """
    creates newick tree from pairwise-distance matrix
    """
    dists = np.array(pd_matrix)
    tree = neighbor_joining(dists)
    my_tree = tree.to_newick(labels=seq_names,round_distance=2)
    return my_tree