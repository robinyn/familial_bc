# Name: parse_variant_type.py
# Version 1.0
# Author: Euisuk Robin Han
# Date: 28 Sept 2023
# Description: A script to re-format parsed annotations from SWEA/BRIDGES annotation pipeline V2

import argparse

def init_args():
    parser = argparse.ArgumentParser(prog = "SWEA/BRIDGES annotation parser", \
                                     description="A script to re-format parsed annotations from SWEA/BRIDGES annotation pipeline V2")
    parser.add_argument('-i', '--input', default=".", nargs="?", \
                        help="Directory and name of the input file (Default = ./input.tsv)")
    parser.add_argument('-o', '--output', default="./output/", nargs="?", \
                        help="Directory for the output files. (Default = ./output/)")

    args = parser.parse_args()
    input_file = args.input
    output_dir = args.output

    if output_dir.strip().endswith("/"):
        output_dir = output_dir.removesuffix("/")

    return input_file, output_dir

def parse_variant(input_file, output_dir):
    variant_dict=dict()
    filter_type=dict()

    try:
        with open(input_file, "r") as input_table:
            for index, line in enumerate(input_table):
                print("Processing line {}".format(index+1), end="\r", flush=True)
                line = line.strip().split("\t")
                genes_at_position = set()
                variant_types_at_position = set()
                miRNA_target_at_position = set()

                # Disregard header line
                if line[0].startswith("#"):
                    continue
                # If variant did not pass filtering, add filter type to dict,
                # increase count and skip to next iteration

                if line[7] != "PASS":
                    if line[7] not in filter_type:
                        filter_type[line[7]]=1
                    else:
                        filter_type[line[7]]+=1
                    continue

                # Parse variant info/consequence
                sample_name = line[0]
                chromosome = line[1]
                position = line[2]
                ref_allele = line[4]
                alt_allele = line[5]

                fixed_csq = line[25].split(",")[0].split("|")
                variable_csq = line[25].split(",")[1:]

                clinvar = line[26]

                # Create variant ID (chr-pos-ref-alt)
                variant = "{}-{}-{}-{}".format(chromosome, position, ref_allele, alt_allele)

                # Parse fixed csq
                known = fixed_csq[0]
                global_allele_freq = fixed_csq[1]
                eur_allele_freq = fixed_csq[2]
                swe_allele_freq = fixed_csq[3]

                # Loop through all transcripts overlapping the position of the variant and
                # identify all genes/variant types
                for transcript in variable_csq:
                    transcript = transcript.split("|")
                    genes_at_position.add(transcript[1])

                    for type in transcript[6].split("&"):
                        variant_types_at_position.add(type)

                    if transcript[14] != "":
                        for i, targetscan in enumerate(transcript[14].split("&")):
                            if i == 0:
                                miRNA_target_at_position.add(":".join(targetscan.split(":")[0:3]))
                            else:
                                miRNA_target_at_position.add(":".join(targetscan.split(":")[1:4]))

                if variant in variant_dict:
                    variant_dict[variant][-1]+="|{}".format(sample_name)
                else:
                    variant_dict[variant]=[known, "|".join(genes_at_position), "|".join(variant_types_at_position),
                                           global_allele_freq, eur_allele_freq, swe_allele_freq, "|".join(miRNA_target_at_position), clinvar, sample_name]
        print()
        per_variant_summary_file = "{}/per_variant_summary.tsv".format(output_dir)
        filter_summary_file ="{}/filter_summary.tsv".format(output_dir)

        with open(per_variant_summary_file, "w") as outfile:
            outfile.write("variant\tknown_variation\tgene\ttype\tAF\tEUR_AF\tSwe_AF\tmiRNA\tClinVar\tsamples\n")
            for var in variant_dict.keys():
                outfile.write("{}\t{}\n".format(var, "\t".join(variant_dict[var])))

        with open(filter_summary_file, "w") as outfile:
            outfile.write("filter_type\tcount\n")
            for filter in filter_type.keys():
                outfile.write("{}\t{}\n".format(filter, filter_type[filter]))

    except Exception as e:
        print("ERROR: {}".format(e))

input_file, output_dir = init_args()
parse_variant(input_file,output_dir)