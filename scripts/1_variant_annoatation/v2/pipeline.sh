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

while getopts d:r:s:f:t:o: flag
do
    case "${flag}" in
        d) data_directory=${OPTARG};;
        s) resources_directory=${OPTARG};;
        f) start_flag=${OPTARG};;
		t) data_type=${OPTARG};;
		o) output_directory=${OPTARG};;
    esac
done

printf "SWEA/BRIDGES Variant Annotation Pipeline V2\n"
printf "Initializing arguments\n"

if [[ -z ${data_directory} || ! -d ${data_directory} ]]
then
    printf "\nERROR: Invalid directory for raw data\n"
    exit 0
fi

if [[ -z ${output_directory} ]]
then
	printf "\nNo output directory provided.\nPipeline will output to root directory\n"
	output_directory="."
fi

if [[ -z ${output_directory} || ! -d ${output_directory} ]]
then
    printf "\nERROR: Invalid output directory\n"
    exit 0
fi

if [[ -z ${data_type} ]]
then
	printf "\nData type not provided.\nPipeline will assume SWEA dataset"
	data_type="swea"
fi

if [[ -z ${start_flag} && ${data_type} == "swea" ]]
then
    start_flag="1"
fi

if [[ -z ${start_flag} && ${data_type} == "bridges" ]]
then
    start_flag="3"
fi

if [[ ${data_type} == "bridges" && ${start_flag} -lt 3 ]]
then
	printf "\nERROR: BRIDGES data is pre-filtered and already contain flanking sequence annotations.\nPipeline will run from step 3 (base annotations)\n"
    start_flag="3"
fi

data_directory=$(readlink -f $(echo $data_directory | sed 's/\/$//'))
output_directory=$(readlink -f $(echo $output_directory | sed 's/\/$//'))
resources_directory=$(readlink -f $(echo $resources_directory | sed 's/\/$//'))
scripts_directory=$(dirname $(readlink -f $0))

export -f progressBar

printf "\nAnnotating ${data_type} data\n"
printf "Raw data directory: ${data_directory}\n"
printf "Output directory: ${output_directory}\n"
printf "Project resources directory: ${resources_directory}\n\n"

padding="------------------------------------------------------"

while [[ ${start_flag} -lt 9 ]]
do
	case "${start_flag}" in
	    1)
		text="Filtering raw data"
		printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

		bash ${scripts_directory}/hard_filtering_script.sh ${data_directory} ${output_directory}
	    ;;
	    2)
		text="Retrieving flanking sequences"
		printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

		bash ${scripts_directory}/flanking_sequences.sh ${output_directory} ${resources_directory}
	    ;;
	    3)
		text="Creating base annotations"
		printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

		bash ${scripts_directory}/annotation_script.sh ${output_directory} ${resources_directory} ${data_type} ${data_directory}
	    ;;
	    4)
		text="Adding ENCODE RBS annotations"
		printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

		bash ${scripts_directory}/rna_binding_annotation.sh ${output_directory} ${resources_directory}
	    ;;
	    5)
		text="Adding miRNA target site annotations"
		printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

		bash ${scripts_directory}/miRNA_target_site_annotation.sh ${output_directory} ${resources_directory}
	    ;;
	    6)
		text="Adding ESE/ESS, RSCU annotations"
		printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

		rsync -a -f"+ */" -f"- *" ${output_directory}/5_TargetScanAnnotation/ ${output_directory}/6_CustomAnnotation/

		python3 ${scripts_directory}/custom_annotation_${data_type}.py ${output_directory}/5_TargetScanAnnotation/ ${output_directory} ${resources_directory}
	    ;;
	    7)
		text="Adding ClinVar annotations"
		printf "%s %s %s\n" "${padding:${#text}}" "$text" "${padding:${#text}}"

		rsync -a -f"+ */" -f"- *" ${output_directory}/6_CustomAnnotation/ ${output_directory}/7_ClinVarAnnotation/

		python3 ${scripts_directory}/clinvar.py ${output_directory} ${resources_directory} ${data_type}
	    ;;
	esac
	start_flag=$(($start_flag+1))
done

printf "Pipeline complete.\n"
