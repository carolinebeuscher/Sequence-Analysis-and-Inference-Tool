# MagnumOpus - DNA Sequence Analysis Tool

MagnumOpus is a bioinformatics tool that constructs a newick phylogenetic tree from provided reads and assemblies.

![alt text](/Seq_Analysis_Inference/ex_tree.svg)

## Features

### 1. In silico PCR 
- Analyzes PCR primers against DNA assemblies
- Handles multiple assembly files simultaneously
- Identifies potential amplicons within specified size limits

### 2. Sequence Alignment
- Implements Needleman-Wunsch algorithm for global sequence alignment
- Supports both forward and reverse complement sequence alignment
- Customizable scoring parameters (match, mismatch, gap penalties)
- Returns aligned sequences, alignment score, and pairwise distance

### 3. Read Mapping
- Maps sequencing reads to reference sequences using minimap2
- Generates SAM format output
- Supports paired-end read mapping
- Includes consensus sequence generation

### 4. Phylogenetic Inference
- Generates Newick phylogenetic tree from pairwise distance matrix
- Use https://itol.embl.de/upload.cgi to visualize your tree if desired (see example above)

### 5. File Format Support
- FASTA file parsing and manipulation
- SAM/BAM file handling
- FASTQ read file processing
- Support for multiple input formats

## Dependencies

### Python Packages
- numpy (>=1.21.0)
- pandas (>=1.3.0)
- biopython (>=1.79)
- setuptools (>=45.0.0)
- wheel (>=0.37.0)

### External Tools
- minimap2 (for read mapping)
- samtools (for SAM/BAM file manipulation)
- blastn (for primer analysis)


## Usage

The main script magnop.py can be used with the following arguments:
   ```bash
   python magop.py -p primers.fasta [-a assemblies/*.fasta] [-r reads/*.fastq] [-s reference.fasta]
   ```

### Arguments: 
- -p, --primers: Path to primers file (required)
- -a, --assemblies: Path to assembly files (optional)
- -r, --reads: Path to read files (optional)
- -s, --reference: Path to reference sequence (optional)

## Contributing
This is a course project repository.


