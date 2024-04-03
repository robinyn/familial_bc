import argparse
import os
import glob
import sys
import subprocess
import shlex
import multiprocessing
import EES_ESE_annotation as custom

def init_args():
    parser = argparse.ArgumentParser(prog = "Variant Annotation Pipeline V3")
    subparser = parser.add_subparsers(required=True, dest="command")

    annotate_parser = subparser.add_parser("annotate", help="Annotate VCF files")
    annotate_parser.add_argument('-d', '--data-directory', required=True,
                                 help=("Directory of the raw VpipCF files. All VCF files present in the "
                                 "provided directory will be annotated. If a file is specified, only the file will be annotated."))
    annotate_parser.add_argument('-r', '--resources-directory', required=True,
                                 help="Directory of the resource files.")
    annotate_parser.add_argument('-o', '--output-directory', default="output",
                                 help="Directory for the output files.")
    annotate_parser.add_argument('-c', '--cores', default=1, type=int,
                                 help="Number of cores to use for multiprocessing.")
    annotate_parser.add_argument('--snp-filter', default="resources/snp_filters.txt",
                                 help=("Directory and name of the text file containing filter names"
                                       "and values to be used to filter SNPs. Follow GATK arguments file"
                                       "format."))
    annotate_parser.add_argument('--indel-filter', default="resources/indel_filters.txt",
                                 help=("Directory and name of the text file containing filter names"
                                       "and values to be used to filter indels. Follow GATK arguments file"
                                       "format."))
    annotate_parser.add_argument('--vep', action="store_true", default=False)
    annotate_parser.add_argument('--custom', action="store_true", default=False)
    annotate_parser.add_argument('--fs', action="store_true", default=False)
    annotate_parser.add_argument('--filt', action="store_true", default=False)

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

def unpack(args):
    return custom.ESE_ESS_annotation(*args)

def poolProcesses(cores, func, args):

    p = multiprocessing.Pool(cores)
    p.map(func, args)
    p.close()

    return 0

def run_pipeline(input_directory, output_directory, resources_directory, \
                 cores, filt, fs, vep, custom, snp_filter, indel_filter):

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    if not (filt or fs or vep or custom):
        filt = fs = vep = custom = True

    script_path = os.path.dirname(os.path.realpath(sys.argv[0]))

    base_output_directory = output_directory.rstrip("/")

    io_dictionary = {"none":input_directory,
                    "filt":"{}/1_Filtering/Merged".format(base_output_directory),
                    "fs":"{}/2_FlankingSequence".format(base_output_directory),
                    "vep":"{}/3_VEPAnnotations".format(base_output_directory),
                    "custom":"{}/4_CustomAnnotations".format(base_output_directory)}

    prev_operation = "none"

    if filt:

        print("======================== 1. Filtering variants ========================")

        input_directory = io_dictionary[prev_operation]
        output_directory = "{}/1_Filtering".format(base_output_directory)


        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)
            os.mkdir("{}/SNVs".format(output_directory))
            os.mkdir("{}/Indels".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data=input_directory), recursive=True)

        print("{} VCF files to process.".format(len(list_of_files)))

        print("Separating SNVs and indels.")

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/filtering.sh separateVariants {} {} {} {}".format(script_path, file, output_directory, snp_filter, indel_filter)))

        poolProcesses(cores, subprocess.run, command)

        print("Filtering SNVs.")

        os.mkdir("{}/FilteredSNVs".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data="{}/SNVs".format(output_directory)), recursive=True)

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/filtering.sh filterSNVs {} {} {} {}".format(script_path, file, output_directory, snp_filter, indel_filter)))

        poolProcesses(cores, subprocess.run, command)

        print("Filtering indels.")

        os.mkdir("{}/FilteredIndels".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data="{}/Indels".format(output_directory)), recursive=True)

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/filtering.sh filterIndels {} {} {} {}".format(script_path, file, output_directory, snp_filter, indel_filter)))

        poolProcesses(cores, subprocess.run, command)

        print("Merging filtered variants.")

        os.mkdir("{}/Merged".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data="{}/FilteredIndels".format(output_directory)), recursive=True)

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/filtering.sh mergeVariants {} {} {} {}".format(script_path, file, output_directory, snp_filter, indel_filter)))

        poolProcesses(cores, subprocess.run, command)

        prev_operation = "filt"

        print("Filtering complete.")

    if fs:
        print("================== 2. Retrieving flanking sequences ===================")

        input_directory = io_dictionary[prev_operation]
        output_directory = io_dictionary["fs"]

        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

        os.mkdir("{}/temp".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data=input_directory), recursive=True)

        print("{} VCF files to process.".format(len(list_of_files)))

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/flanking_sequences.sh {} {}".format(script_path, file, output_directory)))

        poolProcesses(cores, subprocess.run, command)

        print("Flanking sequence retrieval complete.")

        prev_operation = "fs"

    if vep:

        print("========================== 3. VEP annotations =========================")

        input_directory = io_dictionary[prev_operation]
        output_directory = io_dictionary["vep"]

        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

        # List all VCF files to annotate in the provided directory
        list_of_files = glob.glob("{data}/**/*.vcf".format(data=input_directory), recursive=True)

        print("{} VCF files to process.".format(len(list_of_files)))

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/vep_annotation.sh -i {} -o {} -r {}".format(script_path, file, output_directory, resources_directory)))

        poolProcesses(cores, subprocess.run, command)

        prev_operation = "vep"

    if custom:

        print("======================== 4. Custom annotations ========================")

        input_directory = io_dictionary[prev_operation]
        output_directory = io_dictionary["custom"]

        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

        # List all VCF files to annotate in the provided directory
        list_of_files = glob.glob("{data}/**/*.vcf".format(data=input_directory), recursive=True)

        print("{} VCF files to process.".format(len(list_of_files)))

        command = []

        for file in list_of_files:
            command.append((file, output_directory, resources_directory))

        poolProcesses(cores, unpack, command)

        prev_operation = "custom"

    return 0

if __name__ == "__main__":
    args = init_args()
    if args.command == "annotate":
        run_pipeline(args.data_directory, args.output_directory, \
                     args.resources_directory, args.cores, \
                     args.filt, args.fs, args.vep, args.custom, \
                     args.snp_filter, args.indel_filter)

    elif args.command == "setup":
        download_dependencies(directory=args.directory)