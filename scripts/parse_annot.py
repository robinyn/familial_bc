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
    parser.add_argument('-t', '--type', default="swea", nargs="?", \
                        help="Type of data to parse (swea/bridges)")

    args = parser.parse_args()
    input_dir = args.input
    output_file = args.output
    data_type = args.type

    return input_dir, output_file, data_type

def parse_vcf(input_dir, output_file, data_type):

    if data_type=="swea":
        header=("#sample_name\tchromosome\tposition\tvariant_id\tref_allele\talt_allele\tquality\tfilter\tCSQ\tClinVar\tformat\tadd_info\tAC\tAF"
        "\tAN\tBaseQRankSum\tClippingRankSum\tDP\tExcessHet\tFS\tMLEAC\tMLEAF"
        "\tMQ\tMQRankSum\tNDA\tQD\tReadPosRankSum\tSOR\tFSEQ\n")
    elif data_type=="bridges":
        header=("#sample_name\tchromosome\tposition\tvariant_id\tref_allele\talt_allele\tquality\tfilter\tCSQ\tClinVar\tformat\tadd_info\tTYPE\tDP"
        "\tVD\tAF\tBIAS\tREFBIAS\tVARBIAS\tPMEAN\tPSTD\tQUAL"
        "\tQSTD\tSBF\tODDRATIO\tMQ\tSN\tHIAF\tADJAF\tSHIFT3\tMSI\tMSILEN\tNM\tHICNT\tHICOV\tLSEQ\tRSEQ\tGDAMP\tTLAMP\tNCAMP\tAMPFLAG\n")
    else:
        print("ERROR: Invalid data type")
        exit

    annot_dict={}

    for item in header.removeprefix("#").removesuffix("\n").split("\t"):
        annot_dict[item]="NA"

    try:
        # files_list = os.listdir(input_dir)
        # for file in files_list:
        #     if '.vcf' not in file:
        #         files_list.remove(file)

        files_list=[]

        input_dir=input_dir.removesuffix("/")

        for directory_item in os.walk(input_dir):
            for file in directory_item[2]:
                if file.endswith('.vcf'):
                    files_list.append(directory_item[0]+"/"+file)


        print("{num_file} VCF files to parse".format(num_file=len(files_list)))

        with open(output_file, "w") as output_tsv:
            output_tsv.write(header)
            for index, file in enumerate(files_list):
                print("Processing {} ({}/{})".format(file, index+1, len(files_list)))
                with open(file, "r") as input_vcf:
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

                        if parsed_info[-1][0]=="CSQ":
                            parsed_info.pop(-1)

                        clinvar = re.findall(r"(ClinVar)=(\S+)", info)
                        csq = re.findall(r"(CSQ)=(\S+?)(?:;ClinVar|$)", info)

                        if clinvar:
                            parsed_info.append(clinvar[0])

                        if csq:
                            parsed_info.append(csq[0])

                        for item in parsed_info:
                            annot_dict[item[0]]=item[1]

                        if data_type=="bridges":
                            del annot_dict["SAMPLE"]

                        output_tsv.write("\t".join(annot_dict.values()) + "\n")

                        for item in header.removesuffix("\n").split("\t")[1:]:
                            annot_dict[item]="NA"
            print()
    except Exception as e:
        print("ERROR: {error}".format(error=e))

input_dir, output_file, data_type = init_args()
parse_vcf(input_dir, output_file, data_type)
