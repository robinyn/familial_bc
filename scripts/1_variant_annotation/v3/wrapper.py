import argparse
import os
import shutil
import glob
import sys
import subprocess
import shlex
import multiprocessing
import tqdm
import EES_ESE_annotation as custom

#################### FUNCTION DEFINITIONS ####################

def init_args():
    '''
    Initialize and parse command-line arguments using the Argparse module
    and verify their validity.

    Arguments:
        None.
    Returns:
        args (object): Parsed arguments of various types.
    '''

    # Initialize an argparse.ArgumentParser object and parse command-line arguments.
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
    annotate_parser.add_argument('-q', '--quiet', action="store_true", default=False)
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

    # Verify parsed arguments.
    if args.command == "annotate":

        # Check if input data directory exists
        if not os.path.isdir(args.data_directory):
            print("ERROR: Provided directory for the raw data does not exist.")
            exit()

        # Check if resources directory exists
        if not os.path.isdir(args.resources_directory):
            print("ERROR: Provided directory for the resources files does not exist.")
            exit()

        # Check if the SNP filters file exists
        if not os.path.isfile(args.snp_filter):
            print("ERROR: SNP filter filter not found. Please check that the resources/snp_filters.txt "
                  "file exists or provide the directory and file name of the filters using the option '--snp-filter'.")
            exit()

        # Check if the indel filters file exists
        if not os.path.isfile(args.indel_filter):
            print("ERROR: Indel filter filter not found. Please check that the resources/indel_filters.txt "
                  "file exists or provide the directory and file name of the filters using the option '--indel-filter'.")
            exit()

        # If specified output directory already exists, verify overwrite.
        # If the directory doesn't exist, create directory.
        if os.path.isdir(args.output_directory):
            cont = ask_yn("The provided output directory already exists."
                          "The pipeline will remove any existing files.\n"
                          "Do you want to continue? (y/n): ")
            if not cont:
                exit()

            list_of_dir = []

            for file in glob.glob("{}/**".format(args.output_directory), recursive=True):
                if os.path.isdir(file):
                    list_of_dir.append(file)
                else:
                    os.remove(file)

            list_of_dir.pop(0)

            for directory in reversed(list_of_dir):
                os.rmdir(directory)

        else:
            os.mkdir(args.output_directory)

    elif args.command == "setup":

        # If specified resources directory already exists, verify overwrite.
        # If the directory doesn't exist, create directory.
        if os.path.isdir(args.directory):
            cont = ask_yn("The provided resources directory already exists."
                          "The pipeline will remove any existing files.\n"
                          "Do you want to continue? (y/n): ")
            if not cont:
                exit()

            shutil.rmtree(args.directory)
            os.mkdir(args.directory)

        else:
            os.mkdir(args.directory)

    return args

def ask_yn(message):
    '''
    Promt a Y/N question and continue asking until the user response
    is either Y/y/N/n.

    Arguments:
        message (str): the promt to be printed.
    Returns:
        response (bool): the user response converted into a boolean.
                         True for yes, False for no.
    '''

    # Dictionary to convert user response to boolean values.
    answer_dict = {"N":False, "n":False, "Y":True, "y":True}

    # Promt the question and record the response.
    response = input(message)

    # Continue asking until the user response is either Y/y/N/n.
    while response not in ["N","n","Y","y"]:
        response = input("Invalid response. Please try again. (y/n): ")

    # Convert response to boolean.
    response = answer_dict[response.strip()]

    return response


def download_resources(directory):
    '''
    Wrapper function to download all publically available resource files using the setup.sh script.
    The script will download the most up-to-date files, whenever possible.
    SweGen allele frequency data is not publically available and will not be downloaded.

    Arguments:
        directory (str): the directory to download the resource files to,
    Returns:
        None.
    '''

    # Retrieve current script path.
    script_path = os.path.dirname(os.path.realpath(sys.argv[0]))

    # Formulate shell command to run setup.sh script.
    command = shlex.split("bash ${script}/setup.sh".format(script=script_path))

    # Invoke a subprocess to run the setup script.
    subprocess.run(command)

def unpack(args):
    '''
    Unpacks arguments to pass to the ESE_ESS_annotation function. Required for
    passing multiple ESE_ESS_annotation processes to the multiprocessing pool.

    Arguments:
        args (str:list): list of output/resource directories and file name of
                         the VCF to annotate.
    Returns:
        custom.ESE_ESS_annotation(*args)
    '''

    return custom.ESE_ESS_annotation(*args)

def poolProcesses(cores, func, args, quiet):
    '''
    Pool processes for parallele processing using the Multiprocessing module and
    initiate processing. The function will also display a progress bar using the
    tqdm module.

    Arguments:
        cores (int): number of cores to use for parallel processing.
        func (func): name of the function to call.
        args (str:list): list of arguments to pass to the function.
        quiet (bool): suppresses the progress bar if set to true.
    Returns:
        None.
    '''

    # Create a pool with the specified number of workers.
    p = multiprocessing.Pool(cores)

    # Map processes to the pool for parallel processing. Display progress bar unless the
    # --quiet option was specified.
    if not quiet:
        list(tqdm.tqdm(p.imap(func, args), total=len(args), ncols=term_width-10, ascii=" =", \
                       bar_format='{percentage:3.0f}%|{bar}|  [ {n_fmt}/{total_fmt} {elapsed} ]'))
    else:
        list(p.imap(func, args))

    # Close the worker pool.
    p.close()

def run_pipeline(input_directory, output_directory, resources_directory, \
                 cores, filt, fs, vep, custom, snp_filter, indel_filter, quiet):

    '''
    Wrapper function to run the pipeline using the different filtering/annotation
    scripts. The function will run the entire pipeline, unless options to only run
    specific parts are provided.

    Arguments:
        input_directory (str): directory of the input VCF files to be annotated.
        output_directory (str): directory for the annotated output VCF files.
        resources_directory (str): directory where the resource files are located.
        cores (int): number of cores to use.
        filt (bool): option to run the filtering step.
        fs (bool): option to run the flanking sequence retrieval step.
        vep (bool): option to run the VEP annotation step.
        custom (bool): option to run the custom annotation step.
        snp_filter (str): directory and name of the SNP filtering options file.
        indel_filter (str): directory and name of the indel filtering options file.
        quiet (bool): suppresses the progress bar if set to true.
    Returns:
        None.
    '''

    # If no particular parts of the pipeline were specified, run the whole pipeline.
    if not (filt or fs or vep or custom):
        filt = fs = vep = custom = True

    # Retrieve current script path.
    script_path = os.path.dirname(os.path.realpath(sys.argv[0]))

    # Base output directory.
    base_output_directory = output_directory.rstrip("/")

    # Subdirectories for outputs from different parts of the pipeline.
    io_dictionary = {"none":input_directory,
                    "filt":"{}/1_Filtering/Merged".format(base_output_directory),
                    "fs":"{}/2_FlankingSequence".format(base_output_directory),
                    "vep":"{}/3_VEPAnnotations".format(base_output_directory),
                    "custom":"{}/4_CustomAnnotations".format(base_output_directory)}

    # Set the previous operation to "none"
    prev_operation = "none"

    #============================ 1. FILTERING ============================
    if filt:

        # Print line to show the beginning of the filtering step.
        print("{s:=^{n}}".format(n=term_width, s=" 1. Filtering variants "))

        # Set input/output directory
        input_directory = io_dictionary[prev_operation]
        output_directory = "{}/1_Filtering".format(base_output_directory)

        os.mkdir(output_directory)
        os.mkdir("{}/SNVs".format(output_directory))
        os.mkdir("{}/Indels".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data=input_directory), recursive=True)

        print("{} VCF files to process.".format(len(list_of_files)))

        print("Separating SNVs and indels.")

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/filtering.sh separateVariants {} {} {} {}".format(script_path, \
                                                                                                  file, output_directory, snp_filter, indel_filter)))

        poolProcesses(cores, subprocess.run, command, quiet)

        print("Filtering SNVs.")

        os.mkdir("{}/FilteredSNVs".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data="{}/SNVs".format(output_directory)), recursive=True)

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/filtering.sh filterSNVs {} {} {} {}".format(script_path, \
                                                                                            file, output_directory, snp_filter, indel_filter)))

        poolProcesses(cores, subprocess.run, command, quiet)

        print("Filtering indels.")

        os.mkdir("{}/FilteredIndels".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data="{}/Indels".format(output_directory)), recursive=True)

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/filtering.sh filterIndels {} {} {} {}".format(script_path, \
                                                                                              file, output_directory, snp_filter, indel_filter)))

        poolProcesses(cores, subprocess.run, command, quiet)

        print("Merging filtered variants.")

        os.mkdir("{}/Merged".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data="{}/FilteredIndels".format(output_directory)), recursive=True)

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/filtering.sh mergeVariants {} {} {} {}".format(script_path, \
                                                                                               file, output_directory, snp_filter, indel_filter)))

        poolProcesses(cores, subprocess.run, command, quiet)

        prev_operation = "filt"

        print("Filtering complete.")

    if fs:
        print("{s:=^{n}}".format(n=term_width, s=" 2. Retrieving flanking sequences "))

        input_directory = io_dictionary[prev_operation]
        output_directory = io_dictionary["fs"]

        os.mkdir(output_directory)

        os.mkdir("{}/temp".format(output_directory))

        list_of_files = glob.glob("{data}/**/*.vcf".format(data=input_directory), recursive=True)

        print("{} VCF files to process.".format(len(list_of_files)))

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/flanking_sequences.sh {} {}".format(script_path, file, output_directory)))

        poolProcesses(cores, subprocess.run, command, quiet)

        print("Flanking sequence retrieval complete.")

        prev_operation = "fs"

    if vep:
        print("{s:=^{n}}".format(n=term_width, s=" 3. VEP annotations "))

        input_directory = io_dictionary[prev_operation]
        output_directory = io_dictionary["vep"]

        os.mkdir(output_directory)

        # List all VCF files to annotate in the provided directory
        list_of_files = glob.glob("{data}/**/*.vcf".format(data=input_directory), recursive=True)

        print("{} VCF files to process.".format(len(list_of_files)))

        command = []

        for file in list_of_files:
            command.append(shlex.split("bash {}/vep_annotation.sh -i {} -o {} -r {}".format(script_path, \
                                                                                            file, output_directory, resources_directory)))

        poolProcesses(cores, subprocess.run, command, quiet)

        prev_operation = "vep"

    if custom:
        print("{s:=^{n}}".format(n=term_width, s=" 4. Custom annotations "))

        input_directory = io_dictionary[prev_operation]
        output_directory = io_dictionary["custom"]

        os.mkdir(output_directory)

        # List all VCF files to annotate in the provided directory
        list_of_files = glob.glob("{data}/**/*.vcf".format(data=input_directory), recursive=True)

        print("{} VCF files to process.".format(len(list_of_files)))

        command = []

        for file in list_of_files:
            command.append((file, output_directory, resources_directory))

        poolProcesses(cores, unpack, command, quiet)

        prev_operation = "custom"


######################### MAIN LOGIC #########################

if __name__ == "__main__":

    term_width = os.get_terminal_size().columns - 20

    print("Variant Annotation Pipeline V3\n")
    print("{s:=^{n}}\n".format(n=term_width, s=" Variant Annotation Pipeline V3 "))
    print("Initializing...")

    args = init_args()


    if args.command == "annotate":

        print("Data directory: {}".format(args.data_directory))
        print("Resources directory: {}".format(args.resources_directory))
        print("Output directory: {}".format(args.output_directory))
        print("Number of cores: {}".format(args.cores))
        print("SNP filter file: {}".format(args.snp_filter))
        print("INDEL filter file: {}\n".format(args.indel_filter))

        run_pipeline(args.data_directory, args.output_directory, \
                     args.resources_directory, args.cores, \
                     args.filt, args.fs, args.vep, args.custom, \
                     args.snp_filter, args.indel_filter, args.quiet)

    elif args.command == "setup":
        download_resources(directory=args.directory)