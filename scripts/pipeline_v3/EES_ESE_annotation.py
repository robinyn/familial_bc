#!/usr/bin/env python3

# Title: Variant consequence/ESE/ESS/dRSCU annotation script for VCF files
#
# Description: This script is a modified version of the custom_annotation.py and annotation_script.sh scripts by
#   Arthur Boffelli Castro. The two scripts were combined and modularized into functions so that they could be
#   accessed from a wrapper script.
#
#   The script annotates variant consequences using Variant Effect Predictor (VEP) from Ensembl.
#   Exonic splicing enhancer/silencer motifs overlapping a given variant, as well as the
#   change in relative synonymous codon usage (dRSCU) values between the reference and alternate codons for
#   synonymous variants are annotated using custom Python functions.
#
# Date: 19/Mar/2024
# Author: Euisuk (Robin) Han
# Github: https://www.github.com/robinyn/familial_bc
#
# Original scripts by Arthur can be found in https://www.github.com/aboffelli/variant_annotation

import re
import os
import subprocess
import shlex
import sys

def initialize(resources_directory):
    # RSCU file
    rscu_file = '{resources}/RSCU_table.tsv'.format(resources=resources_directory)

    # ESE/ESS list files
    ese_file = '{resources}/RESCUE-ESE_hexamers_200703.txt'.format(resources=resources_directory)
    ess_file = '{resources}/ESS_hexamers_200824.txt'.format(resources=resources_directory)

    # Load the ESS/ESE hexamers and RSCU values from the tables
    ese_set = set()
    ess_set = set()
    rscu_table = dict()

    with open(ese_file, "r") as ese, open(ess_file, "r") as ess, open(rscu_file, "r") as rscu:
        for line in ese:
            ese_set.add(line.strip())

        for line in ess:
            ess_set.add(line.strip())

        for line in rscu:
            line = line.strip().split("\t")
            rscu_table[line[0]] = float(line[1])

    return ese_set, ess_set, rscu_table

def reverse_complement(sequence):
    complement_dictionary = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    rev_comp = ''.join(complement_dictionary.get(base, base) for base in reversed(sequence))
    return rev_comp

def fix_headers(line):

    new_csq_header = line

    if line.startswith('##INFO=<ID=CSQ,Number=.,Type=String,'
                        'Description="Consequence annotations from '
                        'Ensembl VEP. Format:'):
        new_csq_header = ('##INFO=<ID=REDCSQ,Number=.,Type=String,'
                        'Description="Consequence annotations from '
                        'Ensembl VEP. Format:AF|EUR_AF|SweGen_AF|'
                        'gnomADe_AF|gnomADe_NFE_AF|PhyloP|GERP">\n'
                        '##INFO=<ID=CSQ,Number=.,Type=String,'
                        'Description="Consequence annotations from '
                        'Ensembl VEP. Format:Allele|Gene|SYMBOL|Feature|'
                        'BIOTYPE|STRAND|EXON|INTRON|Consequence|Codons">\n')

    elif line.startswith('##INFO=<ID=Encode,Number='):
        new_csq_header = ('##INFO=<ID=ENCODE,Number=.,Type=String,'
        'Description="Consequence annotations from '
        'Ensembl VEP. Format:ProteinName_CellLine:Strand:'
        'Log2FoldChange:NegLog10Value">\n')

    elif line.startswith('##INFO=<ID=TargetScan,Number='):
        new_csq_header = ('##INFO=<ID=TARGETSCAN,Number=.,Type=String,'
        'Description="Consequence annotations from '
        'Ensembl VEP. Format:SYMBOL:miRNA:STRAND:'
        'StartCoordinate:EndCoordinate">\n')

    elif line.startswith('##INFO=<ID=ClinVar,Number='):
        new_csq_header = ('##INFO=<ID=CLINVAR,Number=.,Type=String,'
        'Description="Consequence annotations from '
        'Ensembl VEP. Format:VARID|ALLELEID|CLNDISDB|'
        'CLNDN|CLNHGVS|CLNREVSTAT|CLNSIG|CLNSIGCONF|'
        'CLNSIGINCL|CLNVC|CLNVCSO|CLNVI|DBVARID|GENEINFO|MC|ORIGIN">\n')

    elif line.startswith('##INFO=<ID=ClinVar_') or \
    line.startswith('##INFO=<ID=SweGen') or \
    line.startswith('##INFO=<ID=GERP') or \
    line.startswith('##INFO=<ID=PhyloP'):
        new_csq_header = ""

    line = new_csq_header

    return line

def round_values(csq):

    csq = csq.split("|")

    if csq[18]:
        print(csq[18])
        encode = []
        for RBP in csq[18].split("&"):
            RBP = RBP.split(":")

            RBP[3] = str(round(float(RBP[3]), 4))
            RBP[4] = str(round(float(RBP[4]), 4))

            RBP.pop(1)

            encode.append(":".join(RBP))

        encode = "&".join(encode)
        csq[18] = encode

    if csq[16]:
        csq[16] = str(round(float(csq[16]), 4))

    if csq[17]:
        csq[17] = str(round(float(csq[17]), 4))

    csq = "|".join(csq)

    return csq

def find_flanking_seq(line):

    flanking_sequence = re.search(r"FSEQ=(\S[A-Z]+\[.+\/.+\][A-Z]+)",line[7])

    if flanking_sequence:
        flanking_sequence = flanking_sequence.group(1)
    else:
        ref_allele = line[3]
        alt_allele = line[4]

        left_seq = re.search(r"LSEQ=(\w+)", line[7])
        right_seq = re.search(r"RSEQ=(\w+)", line[7])

        if left_seq and right_seq:
            flanking_sequence = "{}[{}/{}]{}".format(left_seq.group(1), ref_allele, alt_allele, right_seq.group(1))
        else:
            print("Flanking sequence annotations cannot be found!")
            raise(Exception)

    return flanking_sequence

def calculate_dRSCU(transcript_csq, rscu_table):

    variant_type = transcript_csq[7]
    codon = transcript_csq[8]
    dRSCU=""

    if "synonymous_variant" in variant_type:
        ref_codon = codon[:3].upper()
        alt_codon = codon[-3:].upper()
        dRSCU = str(round(rscu_table[alt_codon] - rscu_table[ref_codon], 2))

    return [dRSCU]

def annotate_ESE_ESS(transcript_csq, flanking_sequence, ese_set, ess_set):

    ese_ref = []
    ese_alt = []
    ess_ref = []
    ess_alt = []

    flanking_sequence = re.search(r'(\w{5})\[(\S+)\/(\S+)\](\w{5})', flanking_sequence)

    ref_allele = flanking_sequence.group(2)
    alt_allele = flanking_sequence.group(3)
    left_sequence = flanking_sequence.group(1)
    right_sequence = flanking_sequence.group(4)

    exon = transcript_csq[6]

    if exon and (ref_allele!="-" and alt_allele!="-" and len(ref_allele + alt_allele) == 2):
        ref_seq = left_sequence + ref_allele + right_sequence
        alt_seq = left_sequence + alt_allele + right_sequence

        if transcript_csq[5] == '-1':
            ref_seq = reverse_complement(ref_seq)
            alt_seq = reverse_complement(alt_seq)

        for i in range(len(ref_seq) - 5):
            hexamer = ref_seq[i:i+6]
            if hexamer in ese_set:
                ese_ref.append(hexamer)
            if hexamer in ess_set:
                ess_ref.append(hexamer)

        for i in range(len(alt_seq) - 5):
            hexamer = alt_seq[i:i+6]
            if hexamer in ese_set:
                ese_alt.append(hexamer)
            if hexamer in ess_set:
                ess_alt.append(hexamer)

    ese_ref = ":".join(ese_ref)
    ese_alt = ":".join(ese_alt)
    ess_ref = ":".join(ess_ref)
    ess_alt = ":".join(ess_alt)

    return [ese_ref, ese_alt, ess_ref, ess_alt]

def ESE_ESS_annotation(file, output_directory, resources_directory):
    ese_set, ess_set, rscu_table = initialize(resources_directory)

    # Create a new directory to store the output files if it does not exist.
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    newfile = "{}/custom_{}".format(output_directory, file.split("/")[-1])
    with open(file, "r") as input_vcf, open(newfile, "w") as output_vcf:
        for line in input_vcf:
            if line.startswith("#"):
                line = fix_headers(line)
            else:
                line = line.split("\t")

                flanking_sequence = find_flanking_seq(line)
                csq = re.search(r'CSQ=(\S+)', line[7]).group(1).split(",")

                csq[0] = round_values(csq[0])

                fixed_csq = csq[0].split("|")[11:18]
                fixed_csq = "|".join(fixed_csq)

                encode = csq[0].split("|")[18]
                targetscan = csq[0].split("|")[19]

                clinvar = csq[0].split("|")[20:]
                clinvar = "|".join(clinvar)

                var_csq = []

                for transcript_csq in csq:
                    transcript_csq = transcript_csq.split("|")

                    dRSCU = calculate_dRSCU(transcript_csq, rscu_table)
                    ESE_EES = annotate_ESE_ESS(transcript_csq, flanking_sequence, ese_set, ess_set)
                    transcript_csq = transcript_csq[0:11] + dRSCU + ESE_EES
                    var_csq.append("|".join(transcript_csq))

                var_csq = ",".join(var_csq)

                csq = "REDCSQ={};CSQ={};ENCODE={};TARGETSCAN={};CLINVAR={}".format(fixed_csq, var_csq, encode, targetscan, clinvar)

                line = "\t".join(line)
                line = re.sub(r'CSQ=\S+', csq, line)

            output_vcf.write(line)

if __name__ == "__main__":
    file=sys.argv[1].strip().removesuffix("/")
    output_directory=sys.argv[2].strip().removesuffix("/")
    resources_directory=sys.argv[3].strip().removesuffix("/")

    ESE_ESS_annotation(file, output_directory, resources_directory)