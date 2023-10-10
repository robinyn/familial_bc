#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Custom VCF annotation for BRIDGES files.

Description: The script annotates the presence of Exonic Splicing
    Enhancer/Silencer (ESE/ESS) hexamers overlapped by the variant, the
    difference in the values of the Relative Synonymous Codon Usage (RSCU)
    between the reference codon and altered codon, and changes the format of
    the lines to remove redundant information of CSQ.
    The header is also updated.

Created on: 2021-12-10
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation

Updated on: 2023-10-10
Author: Euisuk Robin Han

Update notes:
	- Added command-line arguments to allow the script to be included
	  into the automated pipeline script.
	- Added processing of miRNA target site annotation
"""

import re
import glob
import time
import sys

start_time = time.time()


def reverse_complement(sequence):
    """
    Function that creates the reverse complement of a DNA sequence.

    :param sequence: sequence of DNA
    :return: reverse complement
    """
    rev_comp = ''.join(
        complement_dictionary.get(base, base) for base in reversed(sequence))
    return rev_comp

data_directory=sys.argv[1].strip().removesuffix("/")
output_directory=sys.argv[2].strip().removesuffix("/")
resources_directory=sys.argv[3].strip().removesuffix("/")

ese_file = '{resources}/RESCUE-ESE_hexamers_200703.txt'.format(resources=resources_directory)
ess_file = '{resources}/ESS_hexamers_200824.txt'.format(resources=resources_directory)

list_of_files = glob.glob("{data}/**/*.vcf".format(data=data_directory), recursive=True)
for file in list_of_files.copy():
    if '/targetscan_encode_vep' not in file:
        list_of_files.remove(file)

# Complement dictionary
complement_dictionary = {'A': 'T', 'C': 'G',
                         'T': 'A', 'G': 'C'}
# Load RSCU dictionary
rscu_dictionary = {"TTT": 0.92, "TTC": 1.08, "TTA": 0.46, "TTG": 0.77,
                   "CTT": 0.80, "CTC": 1.16, "CTA": 0.43,
                   "CTG": 2.38, "ATT": 1.09, "ATC": 1.39, "ATA": 0.52,
                   "ATG": 1.00, "GTT": 0.73, "GTC": 0.95,
                   "GTA": 0.47, "GTG": 1.85, "TCT": 1.12, "TCC": 1.30,
                   "TCA": 0.92, "TCG": 0.33, "CCT": 1.15,
                   "CCC": 1.30, "CCA": 1.09, "CCG": 0.46, "ACT": 1.00,
                   "ACC": 1.41, "ACA": 1.15, "ACG": 0.45,
                   "GCT": 1.05, "GCC": 1.61, "GCA": 0.91, "GCG": 0.44,
                   "TAT": 0.89, "TAC": 1.11, "TAA": 0.84,
                   "TAG": 0.66, "CAT": 0.84, "CAC": 1.16, "CAA": 0.53,
                   "CAG": 1.47, "AAT": 0.95, "AAC": 1.05,
                   "AAA": 0.88, "AAG": 1.12, "GAT": 0.93, "GAC": 1.07,
                   "GAA": 0.85, "GAG": 1.15, "TGT": 0.92,
                   "TGC": 1.08, "TGA": 1.50, "TGG": 1.00, "CGT": 0.48,
                   "CGC": 1.11, "CGA": 0.65, "CGG": 1.23,
                   "AGT": 0.90, "AGC": 1.44, "AGA": 1.26, "AGG": 1.27,
                   "GGT": 0.64, "GGC": 1.35, "GGA": 1.00,
                   "GGG": 1.01}

# Load the ESS and ESE hexamers from the tables
ese_set = set()
ess_set = set()

with open(ese_file) as ese, open(ess_file) as ess:
    for line in ese:
        ese_set.add(line.strip())

    for line in ess:
        ess_set.add(line.strip())

file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)
    new_file = re.sub(r'5_TargetScanAnnotation(\S*/)(targetscan)',
                      r'6_CustomAnnotation\1custom_\2', file)

    with open(file, 'r') as vcf, open(new_file, 'w') as out_vcf:
        for line in vcf:
            # Update the header lines.
            if line.startswith("#"):
                line = line.strip()

                # Update the CSQ fields to the new format.
                if line.startswith('##INFO=<ID=CSQ,Number=.,Type=String,'
                                   'Description="Consequence annotations from '
                                   'Ensembl VEP. Format:'):
                    new_csq_header = ('Existing_variation|AF|EUR_AF|SweGen_AF|'
                                     'gnomAD_AF|gnomAD_NFE_AF|PhyloP|GERP,'
                                     'Gene|SYMBOL|Feature|STRAND|EXON|INTRON|'
                                     'Consequence|Codons|Encode|delta_RSCU|'
                                     'ESEs_REF|ESEs_ALT|ESSs_REF|ESSs_ALT|TargetScan')

                    # Change the line using regex substitution.
                    line = re.sub(r'(.*Format:).*(">)', r'\g<1> ' +
                                  new_csq_header + r'\g<2>', line)
                    print(line, file=out_vcf)

                # Update the Encode info line.
                elif line.startswith('##INFO=<ID=Encode,Number=.,Type=String,'
                                     'Description="Consequence annotations from'
                                     ' Ensembl VEP. Format: Feature|Encode">'):
                    line = line.rstrip('Feature|Encode">')
                    line += ('ProteinName_CellLine:Strand:Log2FoldChange:'
                             'NegLog10Value">')
                    print(line, file=out_vcf)

                else:  # All other header lines.
                    print(line, file=out_vcf)

            # Variants lines.
            else:
                line = line.strip()
                split_line = line.split('\t')
                ref_base = split_line[3]
                alt_base = split_line[4]
                line_info = split_line[7]

                # Retrieve the flanking sequence.
                left_seq = re.search(r"LSEQ=(\w+)", line_info).group(1)
                right_seq = re.search(r"RSEQ=(\w+)", line_info).group(1)

                # Retrieve the transcripts between CSQ and Encode.
                transcripts = re.search(r'CSQ=(\S+);Encode=', line_info).group(1).split(',')
                # Retrieve the Encode info.
                encode_line = re.search(r'Encode=(\S+);TargetScan=', line_info).group(1).split(',')
                # Retrieve the TargetScan annotations
                targetscan = re.search(r'TargetScan=(\S+)',line_info).group(1).split(',')

                # Start fixed_csq variable that will be composed by the
                # information that does not change for the transcripts
                # (redundant CSQ).
                fixed_csq = ''
                for index, transcript in enumerate(transcripts.copy()):
                    csq = transcript.split('|')
                    strand = csq[4]
                    exon = csq[5]
                    consequence = csq[7]
                    codon = csq[8].upper()
                    delta_rscu = ''
                    ref_motifs_ese = []
                    ref_motifs_ess = []
                    alt_motifs_ese = []
                    alt_motifs_ess = []
                    encode_info = encode_line[index].split('|')[1]

                    # Round encode numbers
                    if encode_info:
                        # If encode has more than one protein, split it.
                        if '&' in encode_info:
                            encode_info = encode_info.split('&')
                        # If only one protein, transform it into a list item.
                        else:
                            encode_info = [encode_info]
                        # Loop through the encode_info list.
                        for number, protein in enumerate(encode_info):
                            protein = protein.split(':')
                            # Round the numbers
                            protein[3] = str(round(float(protein[3]), 2))
                            protein[4] = str(round(float(protein[4]), 2))
                            # Remove the number 1000.
                            protein.pop(1)
                            # Join the list in the right index.
                            encode_info[number] = ':'.join(protein)
                        # Join all the proteins again
                        encode_info = '&'.join(encode_info)

                    # Round PhyloP and GERP to four numbers
                    phylop = csq[14]
                    gerp = csq[15]
                    # Change the csq items to join them in the end.
                    if phylop:
                        csq[14] = str(round(float(phylop), 4))
                    if gerp:
                        csq[15] = str(round(float(gerp), 4))
                    if not fixed_csq:
                        fixed_csq = [csq[3]]
                        for inf in csq[9:]:
                            fixed_csq.append(inf)

                    # Calculate the deltaRSCU for synonymous variants:
                    if 'synonymous_variant' in consequence:
                        ref_codon = codon[:3]
                        alt_codon = codon[-3:]
                        delta_rscu = round(rscu_dictionary[alt_codon]
                                           - rscu_dictionary[ref_codon], 2)

                    # Check for ESE and ESS in exonic SNV variants
                    if exon and (len(split_line[3]) + len(split_line[4]) == 2):
                        # Retrieve the surrounding sequence of the variant. As
                        # we want to check hexamers, we only need 5 more bases
                        # on each side of the variant to have hexamers that
                        # include the variant.
                        ref_seq = left_seq[-5:] + ref_base + right_seq[:5]
                        alt_seq = left_seq[-5:] + alt_base + right_seq[:5]

                        # Check if the variant is located in the reverse strand
                        # and get the reverse complement.
                        if strand == '-1':
                            ref_seq = reverse_complement(ref_seq)
                            alt_seq = reverse_complement(alt_seq)

                        # Check all possible hexamers in the reference sequence
                        # and if they are present in the the ESE and ESS sets.
                        for i in range(len(ref_seq) - 5):
                            hexamer = ref_seq[i:i + 6]
                            if hexamer in ese_set:
                                ref_motifs_ese.append(hexamer)
                            if hexamer in ess_set:
                                ref_motifs_ess.append(hexamer)

                        # Check all possible hexamers in the altered sequence
                        # and if they are present in the the ESE and ESS sets.
                        for i in range(len(alt_seq) - 5):
                            hexamer = alt_seq[i:i + 6]
                            if hexamer in ese_set:
                                alt_motifs_ese.append(hexamer)
                            if hexamer in ess_set:
                                alt_motifs_ess.append(hexamer)

                    # Parse TargetScan annotation
                    if not targetscan[index].split('|')[1]=="":
                        targetscan_info = targetscan[index].split('|')[1].split(':')
                        targetscan_info.pop(0)
                    else:
                        targetscan_info = ""

                    # Remove fixed info from the original csq.
                    csq.pop(3)
                    csq = csq[0:8]

                    # Append the new info to csq and join it together again.
                    csq.append(encode_info)
                    csq.append(str(delta_rscu))
                    csq.append(';'.join(ref_motifs_ese))
                    csq.append(';'.join(alt_motifs_ese))
                    csq.append(';'.join(ref_motifs_ess))
                    csq.append(';'.join(alt_motifs_ess))
                    csq.append(':'.join(targetscan_info))
                    csq = '|'.join(csq)
                    # Put the new csq in the right transcript
                    transcripts[index] = csq

                # Join the fixed csq and transcripts together again
                fixed_csq = '|'.join(fixed_csq)
                transcripts = ','.join(transcripts)

                # Replace the old CSQ by the new CSQ in info column of the line.
                split_line[7] = re.sub(r'(CSQ=)\S.+', r'\g<1>' + fixed_csq +
                                       ',' + transcripts, line_info)

                # Join the line
                line = '\t'.join(split_line)
                # Print the line in the file
                print(line, file=out_vcf)
    file_count += 1

# Print the run time.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
