import mappy as mp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import argparse
from collections import defaultdict
import gzip
from functools import partial

def get_options():
    description = "Cuts out loci based on alignment of reference sequences"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python loci_cutter.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--infiles',
                    required=True,
                    help='List of file paths to cut, one per line')
    IO.add_argument('--query',
                    help='Fasta file of sequences to align and cut. Each cut will be conducted with a pair of sequences paired using the same key.')
    IO.add_argument('--cutoff',
                    type=float,
                    default=0.7,
                    help='Cutoff of alignment length to confirm match. '
                         'Default = 0.7')
    IO.add_argument('--outpref',
                    default="result",
                    help='Output prefix ')
    return parser.parse_args()

def get_best_map(index, pair_id, seq_index, sequence, cutoff):
    a = mp.Aligner(index, preset="asm10")

    best_map = (None, None, None, 0, 0, 0)

    for hit in a.map(sequence):
        if not hit.is_primary:
            continue
        query_hit = hit.blen

        # set cutoff for minimum alignment length
        if query_hit < cutoff * len(sequence):
            continue

        if query_hit > best_map[-1]:
            best_map = (pair_id, seq_index, hit.ctg, hit.r_st, hit.r_en, query_hit)

    return best_map

def cut_loci(infiles, seq_pair_dict, cutoff):
    file_list = []

    cut_records = []

    with open(infiles, "r") as f:
        for line in f:
            file_list.append(line.strip())

    for file in file_list:
        # check if FASTA is gzipped
        gzipped = False
        with open(file, 'rb') as test_f:
            gzipped = True if test_f.read(2) == b'\x1f\x8b' else False

        _open = partial(gzip.open, mode='rt') if gzipped == True else open
        
        for pair_id, seq_pair in seq_pair_dict.items():
        
            best_map_pair = [None, None]
            # find best match for each sequence
            for seq_index, seq in enumerate(seq_pair):
                best_map_pair[seq_index] = get_best_map(file, pair_id, seq_index, seq, cutoff)

            seq1_valid = best_map_pair[0][0] != None
            seq2_valid = best_map_pair[1][0] != None

            # make sure at least one sequences found a hit
            if seq1_valid or seq2_valid:
                with _open(file) as handle:
                    fasta_sequences = SeqIO.parse(handle, 'fasta')
                    for fasta in fasta_sequences:
                        contig_id, sequence = fasta.id, str(fasta.seq)
                        seq1_match = contig_id == best_map_pair[0][2]
                        seq2_match = contig_id == best_map_pair[1][2]
                        locus_1 = None
                        locus_2 = None
                        detail = ""

                        # both sequences match contig
                        if seq1_match and seq2_match:
                            # work out which way round to cut
                            locus_1 = min(best_map_pair[0][3], best_map_pair[1][3])
                            locus_2 = max(best_map_pair[0][4], best_map_pair[1][4])
                            detail = "complete"
                        # first sequence matches
                        elif seq1_match and not seq2_match:
                            locus_1 = best_map_pair[0][3]
                            locus_2 = best_map_pair[0][4]
                            detail = "1_only"
                        # second sequence matches
                        elif not seq1_match and seq2_match:
                            locus_1 = best_map_pair[1][3]
                            locus_2 = best_map_pair[1][4]
                            detail = "2_only"
                        
                        # ensure match found
                        if locus_1 != None and locus_2 != None:
                            cut = sequence[locus_1:locus_2]

                            cut_records.append(SeqRecord(Seq(cut), id=contig_id + "_" + pair_id + "_" + detail,
                                                        description=fasta.description))                        
    return cut_records

def main():
    options = get_options()
    infiles = options.infiles
    query = options.query
    cutoff = options.cutoff
    outpref = options.outpref
    
    # check if FASTA is gzipped
    gzipped = False
    with open(query, 'rb') as test_f:
        gzipped = True if test_f.read(2) == b'\x1f\x8b' else False

    _open = partial(gzip.open, mode='rt') if gzipped == True else open
    
    # get pairs of query sequences
    seq_pair_dict = defaultdict(list)

    with _open(query) as handle:
        fasta_sequences = SeqIO.parse(handle, 'fasta')
        for fasta in fasta_sequences:
            id, sequence = fasta.id, str(fasta.seq)
            seq_pair_dict[id].append(sequence)
    
    # check that each sequence pair only has two sequences
    for pair_id, seq_pair in seq_pair_dict.items():
        assert len(seq_pair) == 2
    
    # generate list of cut loci
    cut_records = cut_loci(infiles, seq_pair_dict, cutoff)

    # write cut loci
    SeqIO.write(cut_records, outpref + "_cut.fa", "fasta")

    return 0

if __name__ == "__main__":
    main()