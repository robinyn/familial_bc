import argparse
import time
import os
import subprocess
import glob
import re
import shlex
import annotation_functions as custom

def init_args():
    parser = argparse.ArgumentParser(prog = "Variant Annotation Pipeline V3")
    subparser = parser.add_subparsers(required=True, dest="command")

    filter_parser = subparser.add_parser("filter", help="Filter variants using GATK VariantFiltration command")
    filter_parser.add_argument('--snp-filter', help=("Directory and name of the text file containing filter names"
                                                     "and values to be used to filter SNPs. Follow GATK arguments file"
                                                     "format."))
    filter_parser.add_argument('--indel-filter', help=("Directory and name of the text file containing filter names"
                                                     "and values to be used to filter indels. Follow GATK arguments file"
                                                     "format."))
    filter_parser.add_argument('-i', '--input-directory', required=True,
                                 help=("Directory of the VCF files to be filtered. If a file is specified, only the file "
                                       "will be filtered"))
    filter_parser.add_argument('-o', '--output-directory', default="output",
                                 help="Directory for the filtered output files.")

    annotate_parser = subparser.add_parser("annotate", help="Annotate VCF files")
    annotate_parser.add_argument('-d', '--data-directory', required=True,
                                 help=("Directory of the raw VCF files. All VCF files present in the "
                                 "provided directory will be annotated. If a file is specified, only the file will be annotated."))
    annotate_parser.add_argument('-r', '--resources-directory', required=True,
                                 help="Directory of the resource files.")
    annotate_parser.add_argument('-o', '--output-directory', default="output",
                                 help="Directory for the output files.")

    setup_parser = subparser.add_parser("setup", help=("Download and install necessary dependencies, "
                                        "download publically available resources for annotation."))
    setup_parser.add_argument('-d', '--directory', required=True,
                              help=("Directory to download the resource files to. The downloaded dependencies will "
                                    "be installed in a folder named 'bin' inside this directory."))

    args = parser.parse_args()

    if args.command == "annotate":
        if not os.path.isdir(args.data_directory):
            print("ERROR: Provided directory for the raw data does not exist.")
            exit()
        if not os.path.isdir(args.resources_directory):
            print("ERROR: Provided directory for the resources files does not exist.")
            exit()
        if os.path.isdir(args.output_directory):
            response = ask_yn("The provided output directory already exists. Do you want to continue? (y/n): ")

            if response in ["n","N"]:
                exit()

    elif args.command == "setup":
        if os.path.isdir(args.directory):
            response = ask_yn("The provided directory already exists. Do you want to continue? (y/n): ")

            if response in ["n","N"]:
                exit()

    return args

def ask_yn(message):
    response = input(message)

    while response not in ["N","n","Y","y"]:
        response = input("Invalid response. Please try again. (y/n): ")

    return response

def download_dependencies(directory):

    if not os.path.isdir(directory):
        os.mkdir(directory)

    directory = "{}/bin".format(directory)

    if not os.path.isdir(directory):
        os.mkdir(directory)

    os.chdir(directory)

    command = "conda create -f"

def filter_variants(input_directory, output_directory):



    return 0

def run_pipeline(input_directory, output_directory, resources_directory):

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    base_output_directory = output_directory

    output_directory = "{}/1_VEPAnnotations".format(base_output_directory)

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    custom.vep_annotation(input_directory, output_directory, resources_directory)

    input_directory = output_directory
    output_directory = "{}/2_CustomAnnotations".format(base_output_directory)

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    custom.ESE_ESS_annotation(input_directory, output_directory, resources_directory)

    return 0


args = init_args()
if args.command == "annotate":
    run_pipeline(args.data_directory, args.output_directory, args.resources_directory)

elif args.command == "setup":
    download_dependencies(directory=args.directory)