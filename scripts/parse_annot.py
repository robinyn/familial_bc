# Name: parse_annot.py
# Version 1.0
# Author: Euisuk Robin Han
# Date: 26 Sept 2023
# Description: A script to parse annotations from SWEA/BRIDGES annotation pipeline V2

import argparse
import os
import re

def init_args():
    parser = argparse.ArgumentParser(prog = "SWEA/BRIDGES annotation parser", \
                                     description="A script to parse annotations from SWEA/BRIDGES annotation pipeline V2")
    parser.add_argument('-i', '--input', default=".", nargs="?", \
                        help="Directory of input files (Default = current directory)")
    parser.add_argument('-o', '--output', default="output.tsv", nargs="?", \
                        help="Directory and name of the output file. (Default = ./output.tsv)")

    args = parser.parse_args()
    input_dir = args.input
    output_file = args.output

    return input_dir, output_file

def parse_vcf(input_dir, output_file):
    header=("#sample_name\tchromosome\tposition\tvariant_id\tref_allele\talt_allele\tquality\tfilter\tAC\tAF"
    "\tAN\tBaseQRankSum\tClippingRankSum\tDP\tExcessHet\tFS\tMLEAC\tMLEAF"
    "\tMQ\tMQRankSum\tNDA\tQD\tReadPosRankSum\tSOR\tFSEQ\tCSQ\tClinVar\tformat\tadd_info\n")

    annot_dict={}

    for item in header.removeprefix("#").removesuffix("\n").split("\t"):
        annot_dict[item]="NA"

    try:
        files_list = os.listdir(input_dir)
        for file in files_list:
            if '.vcf' not in file:
                files_list.remove(file)

        print("{num_file} VCF files to parse".format(num_file=len(files_list)))

        with open(output_file, "w") as output_tsv:
            output_tsv.write(header)
            for index, file in enumerate(files_list):
                print("\rProcessing {} ({}/{})\n".format(file, index, len(files_list)))
                with open("{dir}/{name}".format(dir=input_dir, name=file), "r") as input_vcf:
                    for line in input_vcf:
                        # Disregard metadata/header lines
                        if line.startswith("##"):
                            continue
                        elif line.startswith("#CHROM"):
                            annot_dict["sample_name"] = line.split("\t")[-1].strip()
                            continue

                        line = line.strip().split("\t")

                        annot_dict["chromosome"] = line[0]
                        annot_dict["position"] = line[1]
                        annot_dict["variant_id"] = line[2]
                        annot_dict["ref_allele"] = line[3]
                        annot_dict["alt_allele"] = line[4]
                        annot_dict["quality"] = line[5]
                        annot_dict["filter"] = line[6]
                        info = line[7]
                        annot_dict["format"] = line[8]
                        annot_dict["add_info"] = line[9]

                        parsed_info = re.findall(r"(\w+)=([\w.,+\-\[\]/\&\|:]+)(?:;|$)", info)

                        if parsed_info[-1][0]=="ClinVar":
                            parsed_info.pop(-1)

                        clinvar = re.findall(r"(ClinVar)=(\S+)", info)

                        if clinvar:
                            parsed_info.append(clinvar[0])

                        for item in parsed_info:
                            annot_dict[item[0]]=item[1]

                        output_tsv.write("\t".join(annot_dict.values()) + "\n")

                        for item in header.removesuffix("\n").split("\t")[1:]:
                            annot_dict[item]="NA"

    except Exception as e:
        print("ERROR: {error}".format(error=e))

input_dir, output_file = init_args()
parse_vcf(input_dir, output_file)