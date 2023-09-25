#! /bin/sh

function progressBar {
    let _current_progress=(${1}*100/${2}*100)/100
    let _completed=(${_current_progress}*5)/10
    let _incomplete=(50-$_completed)

    _fill=$(printf "%${_completed}s")
    _empty=$(printf "%${_incomplete}s")

    printf "\rProgress: [${_fill// /#}${_empty// /-}] ${_current_progress}%% (${1}/${2})"

    if [ ${_current_progress} -eq 100 ]
    then
        printf "\n"
    fi
}

while getopts d:r:s:f: flag
do
    case "${flag}" in
        d) data_directory=${OPTARG};;
        r) root_directory=${OPTARG};;
        s) resources_directory=${OPTARG};;
        f) start_flag=${OPTARG};;
    esac
done

printf "SWEA/BRIDGES Variant Annotation Pipeline V2\n"
printf "Initializing arguments\n"

if [[ -z ${data_directory} || ! -d ${data_directory} ]]
then
    printf "\nERROR: Invalid directory for raw data\n"
    exit 0
fi

if [[ -z ${root_directory} || ! -d ${root_directory} ]] 
then
    printf "\nERROR: Invalid root directory\n"
    exit 0
fi

if [[ -z ${start_flag} ]]
then
    start_flag="1"
fi

data_directory=$(readlink -f $(echo $data_directory | sed 's/\/$//'))
root_directory=$(readlink -f $(echo $root_directory | sed 's/\/$//'))
resources_directory=$(readlink -f $(echo $resources_directory | sed 's/\/$//'))

export -f progressBar

printf "\nRaw data directory: ${data_directory}\n"
printf "Project root directory: ${root_directory}\n"
printf "Project resources directory: ${resources_directory}\n\n"

padding="------------------------------------------------------"

case "${start_flag}" in
    1)
        text="Filtering raw data"
        printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

        bash ${root_directory}/scripts/hard_filtering_script.sh ${data_directory} ${root_directory}
    ;;
    2)
        text="Retrieving flanking sequences"
        printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

        bash ${root_directory}/scripts/flanking_sequences.sh ${root_directory} ${resources_directory}
    ;;
    3)
        text="Creating base annotations"
        printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

        bash ${root_directory}/scripts/annotation_script.sh ${root_directory} ${resources_directory}
    ;;
    4)
        text="Adding ENCODE RBS annotations"
        printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

        bash ${root_directory}/scripts/rna_binding_annotation.sh ${root_directory} ${resources_directory}
    ;;
    5)
        text="Adding gene name annotations"
        printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

        bash ${root_directory}/scripts/gene_name_annotation.sh ${root_directory}
    ;;
    6)
        text="Adding miRNA target site annotations"
        printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

        bash ${root_directory}/scripts/miRNA_target_site_annotation.sh ${root_directory} ${resources_directory}
    ;;
    7)
        text="Adding ESE/ESS, RSCU annotations"
        printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

        python3 ${root_directory}/scripts/custom_annotation.py ${data_directory} ${root_directory} ${resources_directory}
    ;;
    8)
        text="Adding ClinVar annotations"
        printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

        python3 ${root_directory}/scripts/clinvar.py ${root_directory} ${resources_directory}
    ;;
esac
printf "Pipeline complete.\n"
