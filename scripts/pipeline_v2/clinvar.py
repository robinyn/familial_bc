#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Custom ClinVar Annotation.

Description: Script that annotates the clinical information of the variant
    based on the file variant_summary.txt obtained from ClinVar at
    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz.

Created on: 2021-12-10
Author: Arthur Boffelli Castro

GitHub: https://github.com/aboffelli/variant_annotation

Updated on: 2023-10-10
Author: Euisuk Robin Han

Update notes:
	- Added command-line arguments to allow the script to be included
	  into the automated pipeline script.
	- Modified to work with both SWEA and BRIDGES datasets
"""
import os
import time
import sys
import glob
import re

output_directory=sys.argv[1]
resources_directory=sys.argv[2]
data_type=sys.argv[3]

start_time = time.time()

# Put all files in a list, removing anything that is not a vcf file.
files_directory = r"{root}/6_CustomAnnotation/".format(root=output_directory)
list_of_files = glob.glob("{}/**/*.vcf".format(files_directory), recursive=True)

for file in list_of_files.copy():
    if '/custom' not in file:
        list_of_files.remove(file)

# Initiate a dictionary to store the ClinVar file info.
clinvar_dict = {}
with open(r'{resources}/variant_summary_GRCh37.txt'.format(resources=resources_directory),
          'r') as variant:
    print("Loading ClinVar info...")
    for line in variant:
        # Use the chromosome numbers as keys for a nested dictionary.

        if data_type=="swea":
            chrom = line.split('\t')[18]
        else:
            chrom = 'chr' + line.split('\t')[18]

        position = line.split('\t')[31]
        if chrom not in clinvar_dict:
            clinvar_dict[chrom] = {}
        # Add the line as a set to each position.
        if position not in clinvar_dict[chrom]:
            clinvar_dict[chrom][position] = {line}
        else:
            clinvar_dict[chrom][position].add(line)
    print("Done.")

# New line to add in the vcf header.
new_info = ('##INFO=<ID=ClinVar,Number=.,Type=String,Description='
            '"ClinVar annotation from python script clinvar.py. Format: '
            'AlleleID|Type|Name|GeneID|GeneSymbol|HGNC_ID|'
            'ClinicalSignificance|ClinSigSimple|LastEvaluated|'
            'RS#(dbSNP)|nsv/esv(dbVar)|RCVaccession|PhenotypeIDS|'
            'PhenotypeList|Origin|OriginSimple|Assembly|ChromosomeAccession|'
            'Chromosome|Start|Stop|ReferenceAllele|AlternateAllele|'
            'Cytogenetic|ReviewStatus|NumberSubmitters|Guidelines|'
            'TestedInGTR|OtherIDs|SubmitterCategories|VariationID|'
            'PositionVCF|ReferenceAlleleVCF|AlternateAlleleVCF">\n')

directory = '{root}/7_ClinVarAnnotation/'.format(root=output_directory)
if not os.path.exists(directory):
    os.makedirs(directory)

# File count that will be printed in the screen.
file_count = 1
for file in list_of_files:
    print(f"{file_count}/{len(list_of_files)}", end='\r', flush=True)

    new_file = re.sub(r'6_CustomAnnotation(\S*/)(custom)',
                      r'7_ClinVarAnnotation\1clinvar_\2', file)

    with open(file, 'r') as vcf, open(new_file, 'w') as outvcf:
        for vcf_line in vcf:
            # Non header lines
            if not vcf_line.startswith("#"):
                vcf_line = vcf_line.split('\t')
                # Get chromosome number, position and ref and alt bases.
                chrom = vcf_line[0]
                pos = vcf_line[1]
                ref = vcf_line[3]
                alt = vcf_line[4]

                # The search in the ClinVar file is based on the position,
                # reference, and altered bases.
                search = f'{pos}\t{ref}\t{alt}'
                # Search only in the respective chromosome set.
                if pos in clinvar_dict[chrom]:
                    for clinvar_line in clinvar_dict[chrom][pos]:

                        # If the position and bases are the same as the end of
                        # the Clinvar line.
                        if search == '\t'.join(clinvar_line.strip().split(
                                '\t')[-3:]):
                            # Change the '|' in the clinvarline to ',', and the
                            # spaces to '_'. Join the ClinVar line with '|'.
                            clinvar_line = clinvar_line.replace('|', ',')
                            clinvar_line = clinvar_line.replace(' ', '_')
                            clinvar_line = '|'.join(clinvar_line.split('\t'))

                            # Add the clinvar line to the info section of the
                            # vcf line
                            vcf_line[7] = vcf_line[7] + ';ClinVar=' + \
                                          clinvar_line.strip()

                            break
                # Write the line to the file.
                outvcf.write('\t'.join(vcf_line))

            # Header lines
            else:
                if vcf_line.startswith('#CHROM'):
                    # Add the new line in the header, just before the last
                    # header line.
                    outvcf.write(new_info)
                # Write the last line of the header.
                outvcf.write(vcf_line)

        # Raise the file count
        file_count += 1

# Print the run time in the screen.
print('\nRun time: {:.2f} seconds'.format(time.time() - start_time))
