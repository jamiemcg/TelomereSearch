#!/usr/bin/env python

# Find telomeric repeats in contigs. Scans the start/end of contig sequences to
# detect telomeric repeats TTAGGG/CCCTAA allowing for slight variation
# 
# Jamie McGowan, 2021 <jamie.mcgowan@earlham.ac.uk>
# 
# Usage: python TelomereSearch.py -i assembly.fasta

import argparse
import re
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type = str, help = "Input FASTA file", required = True)
parser.add_argument("-l", "--length", type = int, help = "Scan the first/last N bp (default = 200 bp)", required = False)
parser.add_argument("-t", "--threshold", type = float, help = "A telomere is reported if >= X%% of bp scanned is composed of the telomeric repeats (default = 0.4). Cannot be used with '-m'", required = False)
parser.add_argument("-m", "--min_copies", type = int, help = "A telomere is reported if at least N copies of the telomeric repeat sequence is identified. Cannot be used with '-t'", required = False)
parser.add_argument("-f", "--forward", type = str, help = "Regular expression to search start of contig (default = C{2,4}T{1,2}A{1,3}")
parser.add_argument("-r", "--reverse", type = str, help = "Regular expression to search end of contig (default = T{1,3}A{1,2}G{2,4}")

args = parser.parse_args()

# Examine the first/last N bp of the contig (N's excluded below)
if args.length:
    window_size = int(args.length)
else:
    window_size = 300

if args.threshold and args.min_copies:
    print("Error: You can only specify one of '-t THRESHOLD' or '-m MIN_COPIES'")
    sys.exit(-1)

# Min proportion of nucleotides classified as hits for telomere to be called
if args.threshold:
    threshold = float(args.threshold)
else:
    threshold = 0.4

if args.min_copies:
    min_copies = int(args.min_copies)

fasta_file = args.input

if not args.forward and not args.reverse:
    # Default regular expressions for TTAGGG/CCCTAA, allowing some variation
    # Modify these for other organisms
    telomere_F = "C{2,4}T{1,2}A{1,3}"
    telomere_R = "T{1,3}A{1,2}G{2,4}"
elif args.forward and args.reverse:
    telomere_F = args.forward
    telomere_R = args.reverse
else:
    print("Error!")
    print("If specifying a custom regular expression, you mush supply BOTH forward and reverse sequences using -f and -r parameters")
    sys.exit(-1)

# Keep track of how many contigs have hits
N_sequences = 0
N_both = []
N_start_only = []
N_end_only = []

print("\t".join(["Contig", "Length", "First " + str(window_size) + " bp", "Start Repeats", "Start Repeats Count", "Start Repeats Length", "Last " + str(window_size) + " bp", "End Repeats", "End Repeats Cound", "End Repeats Length", "Start Telomere", "End Telomere"]))

for record in SeqIO.parse(fasta_file, "fasta"):
    N_sequences += 1
    sequence_id = str(record.description)

    sequence = str(record.seq)

    # Convert to upper case
    sequence = sequence.upper()

    # Strip N's from start and end of sequence
    sequence = sequence.strip("N")

    start_telomere_found = False
    end_telomere_found  = False

    contig_start = sequence[:window_size]
    contig_end = sequence[(window_size * -1):]

    # Look at start of contig for telomeric repeats
    hits_start = re.findall(telomere_F, contig_start)
    hits_start_length = sum(len(hit) for hit in hits_start)

    if not args.min_copies:
        if (hits_start_length / (window_size * 1.0)) >= threshold:
            start_telomere_found = True
    else:
        if len(hits_start) >= min_copies:
            start_telomere_found = True

    hits_end = re.findall(telomere_R, contig_end)
    hits_end_length = sum(len(hit) for hit in hits_end)

    if not args.min_copies:
        if (hits_end_length / (window_size * 1.0)) >= threshold:
            end_telomere_found = True
    else:
        if len(hits_end) >= min_copies:
            end_telomere_found = True

    if start_telomere_found and end_telomere_found:
        N_both.append(sequence_id)
    elif start_telomere_found:
        N_start_only.append(sequence_id)
    elif end_telomere_found:
        N_end_only.append(sequence_id)

    print("\t".join(map(str, [sequence_id, len(sequence), contig_start, hits_start, len(hits_start), hits_start_length, contig_end, hits_end, len(hits_end), hits_end_length, start_telomere_found, end_telomere_found])))


# Print summary

print()
print("Input filename:", fasta_file)
print()
print("Number of sequence searched:", N_sequences)
print()
print("Number of sequences with telomeric repeats at both ends (passing threshold):", len(N_both))
for i in N_both:
    print(i)

print()
print("Number of sequences with telomeric repeats at start only (passing threshold):", len(N_start_only))
for i in N_start_only:
    print(i)

print()
print("Number of sequences with telomeric repeats at end only (passing threshold):", len(N_end_only))
for i in N_end_only:
    print(i)

print()
print("Parameters Used:")
print("Target sequences:", telomere_F, "and", telomere_R)
print("Window length:", window_size)
print("Threshold:", threshold)
