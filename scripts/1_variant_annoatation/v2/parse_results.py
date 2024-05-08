# Version 1.0
# Author: Euisuk Robin Han
# Date: 27 Oct 2023
# Description: A script to re-format parsed annotations from SWEA/BRIDGES annotation pipeline V2
#
# CHANGE LOG
# Integrated additional filtering for BRIDGES data (13 Dec 23)
#

import argparse

def init_args():
    parser = argparse.ArgumentParser(prog = "SWEA/BRIDGES annotation parser", \
                                     description="A script to re-format parsed annotations from SWEA/BRIDGES annotation pipeline V2")
    parser.add_argument('-i', '--input', default=".", nargs="?", \
                        help="Directory and name of the input file (Default = ./input.tsv)")
    parser.add_argument('-o', '--output', default="./output/", nargs="?", \
                        help="Directory for the output files. (Default = ./output/)")
    parser.add_argument('-t', '--type', default="swea", nargs="?", \
                        help="Type of data to parse (swea/bridges)")

    args = parser.parse_args()
    input_file = args.input
    output_dir = args.output
    data_type=args.type

    if output_dir.strip().endswith("/"):
        output_dir = output_dir.removesuffix("/")

    return input_file, output_dir, data_type


def parse_variant(input_file, output_dir, data_type):
    variant_dict=dict()
    filter_type=dict()
    sample_dict=dict()
    transcript_set = set()


    try:
        with open(input_file, "r") as input_table:
            for index, line in enumerate(input_table):
                print("Processing line {}".format(index+1), end="\r", flush=True)
                line = line.strip().split("\t")

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

                # Additional filtering for BRIDGES data
                if data_type=="bridges":
                    DP = float(line[13])
                    AF = float(line[15])
                    QUAL = float(line[21])
                    MQ = float(line[25])
                    NM = float(line[32])

                    add_filter = []

                    if QUAL < 30:
                        add_filter.append("q30")
                    if AF < 0.2:
                        add_filter.append("f0.2")
                    if MQ < 60:
                        add_filter.append("Q60")
                    if NM > 2:
                        add_filter.append("NM2.0")
                    if AF*DP < 7.5:
                        add_filter.append("fd7.5")

                    add_filter = ";".join(add_filter)

                    if add_filter:
                        if add_filter not in filter_type:
                            filter_type[add_filter]=1
                        else:
                            filter_type[add_filter]+=1
                        continue

                # Parse variant info/consequence
                sample_name = line[0]
                chromosome = line[1]
                position = line[2]
                ref_allele = line[4]
                alt_allele = line[5].split(",")

                genotype = line[11].split(":")[0]

                for i, allele in enumerate(alt_allele, start=1):
                    genotype = genotype.replace(str(i), allele)

                genotype = genotype.replace("0", ref_allele).replace("/", "")

                fixed_csq = line[8].split(",")[0].split("|")
                variable_csq = line[8].split(",")[1:]

                clinvar = line[9]

                # Create variant ID (chr-pos-ref-alt)
                variant = "{}-{}-{}-{}".format(chromosome, position, ref_allele, ",".join(alt_allele))

                # Parse fixed csq
                known = fixed_csq[0]
                global_allele_freq = fixed_csq[1]
                eur_allele_freq = fixed_csq[2]
                swe_allele_freq = fixed_csq[3]
                gnome_AD_AF = fixed_csq[4]
                gnome_AD_NFE_AF = fixed_csq[5]
                phyloP = fixed_csq[6]
                gerp = fixed_csq[7]

                gene_symbol = ""
                strand = None
                variant_type =""
                exon = ""
                intron = ""
                codon = ""
                encode = ""
                rscu = ""
                ref_ese = ""
                alt_ese = ""
                ref_ess = ""
                alt_ess = ""

                variant_types_at_pos = set()
                genes_at_pos = set()

                # Loop through all transcripts overlapping the position of the variant and
                # identify all genes/variant types
                for transcript in variable_csq:
                    transcript = transcript.split("|")
                    gene_symbol = transcript[1]
                    transcript_id = transcript[2]
                    strand = transcript[3]
                    exon = transcript[4]
                    intron = transcript[5]
                    variant_type = transcript[6]
                    codon = transcript[7]
                    encode = transcript[8]
                    rscu = transcript[9]
                    ref_ese = transcript[10]
                    alt_ese = transcript[11]
                    ref_ess = transcript[12]
                    alt_ess = transcript[13]
                    miRNA_target = transcript[14]

                    for type in transcript[6].split("&"):
                        variant_types_at_pos.add(type)

                    transcript_set.add((variant, gene_symbol, transcript_id, strand, exon, intron, variant_type, codon, encode, rscu, ref_ese, alt_ese, ref_ess, alt_ess, miRNA_target))

                    genes_at_pos.add(gene_symbol)

                if sample_name in sample_dict:
                    sample_dict[sample_name].append([variant, genotype])
                else:
                    sample_dict[sample_name]=[[variant, genotype]]

                if variant in variant_dict:
                    variant_dict[variant][-1]+="|{}".format(sample_name)
                else:
                    variant_dict[variant]=[known, "|".join(genes_at_pos), "|".join(variant_types_at_pos), global_allele_freq, eur_allele_freq, swe_allele_freq, gnome_AD_AF, gnome_AD_NFE_AF, phyloP, gerp, clinvar, sample_name]
        print()
        per_variant_summary_file = "{}/per_variant_summary.tsv".format(output_dir)
        filter_summary_file ="{}/filter_summary.tsv".format(output_dir)
        per_transcript_summary_file = "{}/per_transcript_summary.tsv".format(output_dir)
        per_sample_summary_file = "{}/per_sample_summary.tsv".format(output_dir)

        with open(per_variant_summary_file, "w") as outfile:
            outfile.write("variant\tknown_variation\tgene\ttype\tAF\tEUR_AF\tSwe_AF\tgnome_AD_AF\tgnome_AD_NFE_AF\tphyloP\tgerp\tClinVar\tsamples\n")
            for var in variant_dict.keys():
                outfile.write("{}\t{}\n".format(var, "\t".join(variant_dict[var])))

        with open(per_transcript_summary_file, "w") as outfile:
            outfile.write("variant\tgene\ttranscript_id\tstrand\texon\tintron\tvariant_type\tcodon\tencode\trscu\tref_ese\talt_ese\tref_ess\talt_ess\tmiRNA_target\n")
            for item in transcript_set:
                outfile.write("{}\n".format("\t".join(item)))

        with open(per_sample_summary_file, "w") as outfile:
            outfile.write("sampleID\tvariant\tgenotype\n")
            for sample in sample_dict.keys():
                for item in sample_dict[sample]:
                    outfile.write("{}\t{}\t{}\n".format(sample, item[0], item[1]))

        with open(filter_summary_file, "w") as outfile:
            outfile.write("filter_type\tcount\n")
            for filter in filter_type.keys():
                outfile.write("{}\t{}\n".format(filter, filter_type[filter]))

    except Exception as e:
        print("ERROR: {}".format(e))

input_file, output_dir, data_type = init_args()
parse_variant(input_file,output_dir, data_type)
