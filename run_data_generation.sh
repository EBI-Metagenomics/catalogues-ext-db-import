#!/usr/bin/env bash

set -e

VERSION="1.0.0"

# The script runs RNACentral and UniProt file generation, validation and some output reorganisation

function Usage {
    echo "Usage: $0 [-f ftp-folder-name] [-v catalogue-version] [-r /path/to/results/folder]"
    echo "Options:"
    echo "-f          FTP name of the catalogue, for example, human-oral or non-model-fish-gut"
    echo "-v          Catalogue version, for example, v1.0"
    echo "-r          Full path to nextflow pipeline results folder"
    echo "-j          Full path to the previous catalogue JSON, if this is an update (skip if this is not an update)"
    echo "--version   Show script version"
    exit 1
}

print_version() {
    echo "Post-processing data generation script version $VERSION"
}

case "$1" in
    --version)
        print_version
        exit 0
        ;;
esac


GET_REPS() {
    cut -f14 "${RESULTS_PATH}"/genomes-all_metadata.tsv | grep -v "Species" | sort -u
}


function GenerateRNACentralJSON {
    echo "Generating RNAcentral JSON. This could take a while..."
    mkdir -p "${RESULTS_PATH}"/additional_data/rnacentral
    mkdir -p "${RESULTS_PATH}"/additional_data/rnacentral/GFFs

    echo "Copying GFFs"
    for R in $(GET_REPS)
    do
        mkdir -p "${RESULTS_PATH}/all_genomes/${R::-2}/${R}/genomes1/"
        python3 /nfs/production/rdf/metagenomics/pipelines/prod/genomes-pipeline/helpers/move_files.py --command move --dir-from "${RESULTS_PATH}/all_genomes/${R::-2}/${R}/" --dir-to "${RESULTS_PATH}/all_genomes/${R::-2}/${R}/genomes1/"
        python3 /nfs/production/rdf/metagenomics/pipelines/prod/genomes-pipeline/helpers/move_files.py --command copy --dir-from "${RESULTS_PATH}/all_genomes/${R::-2}/${R}/genomes1/" --dir-to "${RESULTS_PATH}/additional_data/rnacentral/GFFs/" --filename "${R}.gff"
    done

    echo "Running JSON generation"
    mitload miniconda && conda activate pybase
    rna_cmd="python3 /nfs/production/rdf/metagenomics/pipelines/prod/catalogues-ext-db-import/helpers/database_import_scripts/rnacentral/generate_rnacentral_json.py \
    -r /nfs/production/rdf/metagenomics/pipelines/prod/catalogues-ext-db-import/helpers/database_import_scripts/rnacentral/rfam_model_lengths_15.0.txt \
    -m '${RESULTS_PATH}/genomes-all_metadata.tsv' -o '${RESULTS_PATH}/additional_data/rnacentral/${CATALOGUE_FOLDER}-rnacentral.json' \
    -d '${RESULTS_PATH}/additional_data/ncrna_deoverlapped_species_reps/' -g '${RESULTS_PATH}/additional_data/rnacentral/GFFs/'\
     -f '${RESULTS_PATH}/additional_data/mgyg_genomes/'"

    if [[ -n $PREV_JSON_PATH ]]; then
        rna_cmd="$rna_cmd --previous-json $PREV_JSON_PATH"
    fi

    eval $rna_cmd


    echo "Removing GFFs"
    rm -r "${RESULTS_PATH}/additional_data/rnacentral/GFFs/"

    if [[ ! -f "${RESULTS_PATH}/additional_data/rnacentral/${CATALOGUE_FOLDER}-rnacentral.json" ]]
    then
        echo "Did not generate the RNAcentral JSON successfully"
        exit 1
    fi

}


function CheckRNACentralErrors {
    # Read the file line by line, check specific lines, and extract second column
    value_3=$(awk -F '\t' 'NR==3 {print $2}' "${RESULTS_PATH}/additional_data/rnacentral/${CATALOGUE_FOLDER}-rnacentral.json.report")
    value_5=$(awk -F '\t' 'NR==5 {print $2}' "${RESULTS_PATH}/additional_data/rnacentral/${CATALOGUE_FOLDER}-rnacentral.json.report")
    value_6=$(awk -F '\t' 'NR==6 {print $2}' "${RESULTS_PATH}/additional_data/rnacentral/${CATALOGUE_FOLDER}-rnacentral.json.report")

    # Check if any of the values are greater than 0
    if (( value_3 > 0 || value_5 > 0 || value_6 > 0 )); then
        echo "Warning: the RNAcentral script had to skip some records. Make sure this is what is expected. If not,
        abort execution of the script and investigate."
        cat "${RESULTS_PATH}/additional_data/rnacentral/${CATALOGUE_FOLDER}-rnacentral.json.report"

        while true; do
            read -r -p "Do you want to abort? (yes/no): " answer
            case "$answer" in
                no)
                    break
                    ;;
                yes)
                    echo "Aborted."
                    exit 1
                    ;;
                *)
                    echo "Invalid response. Please enter 'yes' or 'no'."
                    ;;
            esac
        done
    fi
}


function RunRNACentralValidator {
    echo "Running RNAcentral validator"
    mitload miniconda && conda activate pybase
    cd /nfs/production/rdf/metagenomics/pipelines/prod/rnacentral-data-schema/
    python3 /nfs/production/rdf/metagenomics/pipelines/prod/rnacentral-data-schema/validate.py \
    "${RESULTS_PATH}/additional_data/rnacentral/${CATALOGUE_FOLDER}-rnacentral.json" > \
    "${RESULTS_PATH}/additional_data/rnacentral/validator_output.txt"

    cd "${RESULTS_PATH}"

    if [ -s "${RESULTS_PATH}/additional_data/rnacentral/validator_output.txt" ]; then
    echo "RNAcentral validator found issues. Aborting."
    cat "${RESULTS_PATH}/additional_data/rnacentral/validator_output.txt"
    exit 1
    fi
}


function GenerateUniprotFiles {
    echo "Converting taxonomy for Uniprot"
    mkdir -p ${RESULTS_PATH}/additional_data/uniprot ${RESULTS_PATH}/additional_data/uniprot/uniprot-files
    if [[ -f ${RESULTS_PATH}/additional_data/gtdbtk_results.tar.gz && ! -d ${RESULTS_PATH}/additional_data/gtdbtk_results ]]
    then
        tar -xf ${RESULTS_PATH}/additional_data/gtdbtk_results.tar.gz -C ${RESULTS_PATH}/additional_data/
    fi

    mitload miniconda && conda activate pybase
    python3 /nfs/production/rdf/metagenomics/pipelines/prod/catalogues-ext-db-import/helpers/database_import_scripts/uniprot/preprocess_taxonomy_for_uniprot.py \
    -g ${RESULTS_PATH}/additional_data/gtdbtk_results/ -v "2" \
    -m ${RESULTS_PATH}/genomes-all_metadata.tsv --species-level-taxonomy -t 4 \
    -o ${RESULTS_PATH}/additional_data/uniprot/preprocessed_taxonomy.tsv

    echo "Generating UniProt files"
    ACCS=$(ls ${RESULTS_PATH}/additional_data/prokka_gbk_species_reps/*.gbk | rev | cut -d '/' -f1 | rev | sed "s/\.gbk//")

    for F in $ACCS; do python3 /nfs/production/rdf/metagenomics/pipelines/prod/catalogues-ext-db-import/helpers/database_import_scripts/uniprot/convert_gbk.py \
    -g ${RESULTS_PATH}/additional_data/prokka_gbk_species_reps/${F}.gbk \
    -o ${RESULTS_PATH}/additional_data/uniprot/uniprot-files/${F}_uniprot.gbk \
    -t ${RESULTS_PATH}/additional_data/uniprot/preprocessed_taxonomy.tsv; done

    echo "Generating UniProt metadata"
    python3 /nfs/production/rdf/metagenomics/pipelines/prod/catalogues-ext-db-import/helpers/database_import_scripts/uniprot/generate_uniprot_metadata.py \
    -m ${RESULTS_PATH}/genomes-all_metadata.tsv -o ${RESULTS_PATH}/additional_data/uniprot/${CATALOGUE_FOLDER}_${CATALOGUE_VERSION}_uniprot_metadata.tsv \
    -p ${RESULTS_PATH}/additional_data/uniprot/preprocessed_taxonomy.tsv

    echo "Adding version"
    echo "$VERSION" > "${RESULTS_PATH}/additional_data/uniprot/VERSION.txt"

    echo "UniProt cleanup"
    # gzip the gtdb directory
    cd ${RESULTS_PATH}/additional_data/ || {
    echo "Error: Cannot change to directory ${RESULTS_PATH}/additional_data/. Exiting." >&2
    exit 1
    }
    if [[ -d gtdbtk_results ]]; then
      if [[ ! -f gtdbtk_results.tar.gz ]]; then
          if tar -czvf gtdbtk_results.tar.gz gtdbtk_results; then
              rm -r gtdbtk_results
          else
              echo "GTDB-Tk compression failed. Directory not removed." >&2
          fi
      else
          echo "GTDB-Tk archive exists. Deleting directory without backing up."
          rm -r gtdbtk_results
      fi
    fi
}


while getopts 'f:v:r:j:' flag; do
    case "${flag}" in
        f) export CATALOGUE_FOLDER=$OPTARG ;;
        v) export CATALOGUE_VERSION=$OPTARG ;;
        r) export RESULTS_PATH=$OPTARG ;;
        j) PREV_JSON_PATH=$OPTARG ;;
        *) Usage exit 1 ;;
    esac
done

if [[ -z $CATALOGUE_FOLDER ]] || [[ -z $RESULTS_PATH ]] || [[ -z $CATALOGUE_VERSION ]]; then
  echo 'Not all of the arguments are provided'
  Usage
fi

if [[ -n "$PREV_JSON_PATH" ]]; then
    export PREV_JSON_PATH
fi

cd "${RESULTS_PATH}"
GenerateRNACentralJSON
CheckRNACentralErrors
RunRNACentralValidator
GenerateUniprotFiles
