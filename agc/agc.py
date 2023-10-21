#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw
import numpy as np
np.int = int
__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, 'rt') as fasta_file:
        current_sequence = None
        sequence_lines = []

        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                
                if current_sequence is not None and len("".join(sequence_lines)) >= minseqlen:
                    yield "".join(sequence_lines)
                current_sequence = line
                sequence_lines = []
            else:
                
                sequence_lines.append(line)

        
        if current_sequence is not None and len("".join(sequence_lines)) >= minseqlen:
            yield "".join(sequence_lines)


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    sequence_counts = {}  
    
    sequences = read_fasta(amplicon_file, minseqlen)

    
    for sequence in sequences:
        if sequence in sequence_counts:
            sequence_counts[sequence] += 1
        else:
            sequence_counts[sequence] = 1

    
    sorted_sequences = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)

    
    filtered_sequences = filter(lambda x: x[1] >= mincount, sorted_sequences)

    
    for sequence, count in filtered_sequences:
        yield [sequence, count]

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    if len(alignment_list) != 2:
        raise ValueError("L'alignement doit contenir exactement deux séquences.")

    
    sequence1, sequence2 = alignment_list

    
    identical_count = sum(1 for nt1, nt2 in zip(sequence1, sequence2) if nt1 == nt2)

    
    alignment_length = len(sequence1)
    identity_percentage = (identical_count / alignment_length) * 100

    return identity_percentage

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int=0, kmer_size: int=0) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """

    
    unique_sequences = dereplication_fulllength(amplicon_file, minseqlen, mincount)

    
    OTU_list = []

    for seq1 in unique_sequences:
        added_to_otu = False

        for otu in OTU_list:
            seq2, count = otu
            
            align = nw.global_align(seq1[0], seq2[0], gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))

            identity = get_identity(align)
            if identity > 97:
                added_to_otu = True
                break

        if not added_to_otu:
            
            OTU_list.append(seq1)

    return OTU_list


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, 'w') as file:
        for i, (sequence, occurrence) in enumerate(OTU_list, start=1):
            file.write(f">OTU_{i} occurrence:{occurrence}\n")
            formatted_sequence = textwrap.fill(sequence, width=80)
            file.write(formatted_sequence + "\n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    # Lecture du fichier d'amplicons
    print("Lecture du fichier d'amplicons...")
    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    output_file = args.output_file

    
    print("Dé-duplication des séquences...")
    unique_sequences = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))

    if not unique_sequences:
        print("Aucune séquence trouvée après dé-duplication.")
        sys.exit(1)

    
    print("Regroupement glouton des séquences en OTUs...")
    OTU_list = abundance_greedy_clustering(amplicon_file, minseqlen, mincount)

    if not OTU_list:
        print("Aucune OTU trouvée.")
        sys.exit(1)

    
    print(f"Écriture des OTUs dans {output_file}...")
    write_OTU(OTU_list, output_file)

    print("Le processus est terminé.")


if __name__ == '__main__':
    main()

