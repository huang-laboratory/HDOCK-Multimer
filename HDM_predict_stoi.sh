# Script for running HDM assembly strategy for predicting stoichiometry from given ranges

trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR
set -e

# ------------------------------------------------------------------------------
# Parse arguments from command line
NUM_CPUS=20
NUM_STEPS=0
RANKING="itscore"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -stois|--stois)
            FILE_CANDIDATE_STOIS=$(readlink -f "$2")
            shift 2;;
        -ncpu|--ncpu)
            NUM_CPUS=$2
            shift 2;;
        -steps|--steps)
            NUM_STEPS=$2
            shift 2;;
        -ranking|--ranking)
            RANKING=$2
            shift 2;;
        -h|--help)
            echo "Usage: $0 -stois <candidate_stois.txt> [-ncpu N] [-steps N] [-ranking itscore|confidence]"
            echo ""
            echo "Required arguments:"
            echo "    -stois    : Path to TXT file containing candidate stoichiometries and corresponding subcomponent directories"
            echo ""
            echo "Optional arguments:"
            echo "    -ncpu     : Maximum number of CPUs used in parallel (default: 20)"
            echo "    -steps    : The number of energy minimization steps during structural relaxation (default: 0)"
            echo "    -ranking   : The score used for ranking the predicted models (default: itscore, options: itscore|confidence)"
            exit 0;;
        *)
            echo "[ERROR] Wrong command argument: $1"
            echo "Type '$0 -h' for help."
            exit 1;;
    esac
done


# ------------------------------------------------------------------------------
# Validate required arguments
if [[ -z "$FILE_CANDIDATE_STOIS" ]]; then
    echo "[ERROR] Missing required arguments"
    echo "Usage: $0 -stois <candidate_stois.txt> [-ncpu N] [-steps N] [-ranking itscore|confidence]"
    exit 1
fi

if [[ ! "$RANKING" =~ ^(itscore|confidence)$ ]]; then
    echo "[ERROR] Invalid ranking method: $RANKING"
    echo "Options for ranking: itscore, confidence"
    exit 1
fi

declare -A complex_by_copy subunit_by_copy stoi_by_copy dir_by_copy

while IFS= read -r line || [[ -n $line ]]; do
    line=$(echo $line | tr -d "r")
    complex=$(echo $line | awk '{print $1}')
    subunit=$(echo $line | awk '{print $2}')
    copy_number=$(echo $line | awk '{print $3}')
    file_stoi=$(echo $line | awk '{print $4}')
    subcomponent_dir=$(echo $line | awk '{print $5}')

    complex_by_copy["$copy_number"]=$complex
    subunit_by_copy["$copy_number"]=$subunit
    stoi_by_copy["$copy_number"]=$file_stoi
    dir_by_copy["$copy_number"]=$subcomponent_dir
done < "$FILE_CANDIDATE_STOIS"

echo ""
echo "[INFO] Candidate stoichiometry file:          $FILE_CANDIDATE_STOIS"
echo "[INFO] Ranking method:                        $RANKING"
echo "[INFO] Number of candidate stoichiometries:   ${#complex_by_copy[@]}"
echo -e "\n[INFO] Start running"

START_TIME=$(date +%s)


# ------------------------------------------------------------------------------
# Start modeling
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/" && pwd)"
WORK_DIR=$(dirname "$FILE_CANDIDATE_STOIS")

declare -A confidence_by_copy itscore_by_copy
for k in "${!stoi_by_copy[@]}"; do
    echo ""
    echo "================================================================================"
    echo "[INFO] Processing copy number: $k"
    echo "--------------------------------------------------------------------------------"
    echo "[INFO] Complex:           ${complex_by_copy[$k]}"
    echo "[INFO] Subunit:           ${subunit_by_copy[$k]}"
    echo "[INFO] Stoichiometry:     ${stoi_by_copy[$k]}"
    echo "[INFO] Subcomponent dir:  ${dir_by_copy[$k]}"
    echo "================================================================================"

    ASSEMBLY_PATH="$WORK_DIR/$k"
    mkdir "$ASSEMBLY_PATH"

    cd "$ASSEMBLY_PATH"
    file_stoi="$WORK_DIR/${stoi_by_copy[$k]}"
    subcomponent_dir="$WORK_DIR/${dir_by_copy[$k]}"
    python ${BASE_DIR}/scripts/assemble/dataprocess.py \
        --pdbdir $subcomponent_dir \
        --stoi $file_stoi
    python ${BASE_DIR}/scripts/assemble/predict_stoi.py \
        --pdbdir $subcomponent_dir \
        --nmax 10 \
        --nround 500 \
        --rmsd 5.0 \
        --population 100 \
        --workers $NUM_CPUS \
        --md $NUM_STEPS \
        --stoi $file_stoi

    log_file="$ASSEMBLY_PATH/aa2.log"
    top1_confidence=""
    top1_itscore=""
    while IFS= read -r line; do
        if [[ -z "$top1_confidence" && "$line" =~ ^Confidence:[[:space:]]+initial[[:space:]]+model[[:space:]]+1[[:space:]]+confidence_score:[[:space:]]*([+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?) ]]; then
            top1_confidence="${BASH_REMATCH[1]}"
            continue
        fi
        if [[ -z "$top1_itscore" && "$line" =~ ^Normalized_ITscore:[[:space:]]+initial[[:space:]]+model[[:space:]]+1[[:space:]]+normalized[[:space:]]+it_score:[[:space:]]*([+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?) ]]; then
            top1_itscore="${BASH_REMATCH[1]}"
            continue
        fi
    done < "$log_file"
    if [[ -n "$top1_confidence" && -n "$top1_itscore" ]]; then
        confidence_by_copy[$k]=$top1_confidence
        itscore_by_copy[$k]=$top1_itscore
        echo "[INFO] Top 1 model confidence score: $top1_confidence"
        echo "[INFO] Top 1 model normalized ITscore: $top1_itscore"
    else
        echo "[WARNING] Failed to assemble a full complex for copy number $k"
    fi
done


# ------------------------------------------------------------------------------
# Select the best stoichiometry based on the specified ranking method 
if (( ${#confidence_by_copy[@]} == 0 )) || (( ${#itscore_by_copy[@]} == 0 )); then
    echo "All candidate stoichiometries fail to assemble a full complex"
    echo "HDM prediction failed"
    exit 1
fi

best_copy_number=""
if [[ $RANKING == "itscore" ]]; then 
    best_copy_number=$(
        for k in "${!itscore_by_copy[@]}"; do
            printf "%s\t%s\n" "$k" "${itscore_by_copy[$k]}"
        done | sort -k2,2g | head -n 1 | cut -f1
    )
else
    best_copy_number=$(
        for k in "${!confidence_by_copy[@]}"; do
            printf "%s\t%s\n" "$k" "${confidence_by_copy[$k]}"
        done | sort -k2,2gr | head -n 1 | cut -f1
    )
fi 

END_TIME=$(date +%s)
RUN_TIME=$(($END_TIME - $START_TIME))

echo -e "\n[INFO] HDM program finish, total cost $RUN_TIME s"

complex=${complex_by_copy[$best_copy_number]}
subunit=${subunit_by_copy[$best_copy_number]}
echo ""
echo "Final predicted stoichiometry by HDM:"
echo "${complex}, subunit ${subunit}, copy number: $best_copy_number"
