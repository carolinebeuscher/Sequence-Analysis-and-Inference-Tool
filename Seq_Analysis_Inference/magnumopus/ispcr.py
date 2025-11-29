import subprocess
import tempfile
import os


                            ###########################
                            ##### ispcr functions #####
                            ###########################

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:
    """
    Performs in silico PCR (isPCR) and returns an amplicon (amplified sequence)
    """
    hits = step_one(primer_file, assembly_file)
    pairs = step_two(hits, max_amplicon_size)
    amps = step_three(pairs, assembly_file)
    return amps 


def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    """
    Identifies locations where primers would anneal to the target sequence
    """
    blast_out = run_blastn(primer_file, assembly_file)
    filtered_hits = filter_blastn(blast_out)
    sorted_hits = sort_blastn(filtered_hits)
    return sorted_hits 


def step_two(sorted_hits: list[str], max_amplicon_size: int) -> list[tuple[list[str]]]:
    """
    Identifies pairs of primer annealing sites which would yield an amplicon 
    """
    hit_pairs = primer_pairs(sorted_hits, max_amplicon_size)
    return hit_pairs 


def step_three(hit_pairs: list[tuple[list[str]]], assembly_file: str) -> str:
    """
    Extracts amplified sequences
    """
    bed_string = amplicons(hit_pairs)

    #writes bed_string to temp file to run seqtk in subprocess
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_bed:
        temp_bed.write(bed_string)
        temp_bed_sp = temp_bed.name

    bed_sp = ["seqtk", "subseq", assembly_file, temp_bed_sp]
    amp_seqs = subprocess.run(bed_sp, capture_output=True, text=True, check=True)

    os.unlink(temp_bed_sp)
    return amp_seqs.stdout

## step one ##

def run_blastn(primer_file, assembly_file):
    """
    runs BLASTN to search for sequence matches between the provided primer file and assembly file
    """
    blastn_sp = ["blastn", "-query", primer_file, "-subject", assembly_file, "-task", "blastn-short", "-word_size", "6", "-penalty", "-2", "-outfmt", "6 std qlen"]
    blast = subprocess.run(blastn_sp, capture_output=True, text=True)
    return blast.stdout


def filter_blastn(blast_out):
    """
    sort BLAST output to only keep full length hits with a percent identity of at least 80%
    """
    filter_sp = "awk '$3>=80 && $4==$13'"
    sorted_hits = subprocess.run(filter_sp, input=blast_out, capture_output=True, text=True, shell = True)
    return sorted_hits.stdout


def sort_blastn(filtered_hits) -> list[list[str]]:
    """
    split BLAST lines into lists of columns and sort by position
    """ 
    sorted_hits = [line.split('\t') for line in filtered_hits.strip().split('\n') if line]
    sorted_hits.sort(key=lambda x: (x[1], int(x[9]))) 
    return sorted_hits

## step two ##

def primer_pairs(sorted_hits, max_amplicon_size) -> list[tuple[list[str]]]:
    """
    identifies BLAST hits that are less than the given amplicon size apart and pointing towards one another
    """
    paired_hits = []
    for counter, first in enumerate(sorted_hits):
        first_or = f_or_r(first)
        for second in sorted_hits[counter+1:]:
            if first[1] != second[1]: #if contigs arent the same move on to next
                break 

            second_or = f_or_r(second)
            if first_or + second_or != "FR":
                continue

            dist = abs(int(first[9]) - int(second[8]))
            if dist > max_amplicon_size:
                break #break inner loop bc no use cont to check 

            paired_hits.append((first, second))

    return paired_hits
  

def f_or_r(hit: list[str]) -> str:
    """
    identifies the directionality of a blast hit
    """
    start, end = int(hit[8]), int(hit[9])
    if start < end:
        return "F"
    else:
        return "R"

## step three ##

def amplicons(hit_pairs):
    """
    amplifies the sequences produced by the PCR reaction and returns first amplicon
    """
    #can edit if you want to return list of amplicons instead of just one amplicon
    #amplicons = [] 
    bed = []
    for hit1, hit2 in hit_pairs:
        contig = hit1[1]
        start = int(hit1[9])
        end = int(hit2[9])-1
        bed.append(f"{contig}\t{start}\t{end}\n")
    #bed_string = "".join(bed)
    bed_string = bed[0]

    return bed_string