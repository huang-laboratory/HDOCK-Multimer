# Script for running HDM assembly strategy with crosslinking restraints

trap 'echo "[ERROR] Script failed at line $LINENO"; exit 1' ERR
set -e


# ------------------------------------------------------------------------------
# Parse arguments from command line
NUM_CPUS=20
NUM_STEPS=0
NMAX=10

while [[ $# -gt 0 ]]; do
    case "$1" in
        -stoi|--stoi)
            FILE_STOI=$(readlink -f "$2")
            shift 2;;
        -sub_dir|--sub_dir)
            SUBCOMPONENT_DIR=$(readlink -f "$2")
            shift 2;;
        -crosslink|--crosslink)
            FILE_CROSSLINK=$(readlink -f "$2")
            shift 2;;
        -nmax|--nmax)
            NMAX=$2
            shift 2;;
        -ncpu|--ncpu)
            NUM_CPUS=$2
            shift 2;;
        -steps|--steps)
            NUM_STEPS=$2
            shift 2;;
        -h|--help)
            echo "Usage: $0 -stoi <stoi.json> -sub_dir <subcomponent_dir/> -crosslink <crosslink.txt> [-nmax N] [-ncpu N]"
            echo ""
            echo "Required arguments:"
            echo "    -stoi     : Path to stoichiometry JSON file"
            echo "    -sub_dir  : Directory containing subcomponent structure PDB files"
            echo "    -crosslink: Path to crosslinking restraints TXT file"
            echo ""
            echo "Optional arguments:"
            echo "    -nmax     : Maximum number of output models (default: 10)"
            echo "    -ncpu     : Maximum number of CPUs used in parallel (default: 20)"
            echo "    -steps    : The number of energy minimization steps during structural relaxation (default: 0)"
            exit 0;;
        *)
            echo "[ERROR] Wrong command argument: $1"
            echo "Type '$0 -h' for help."
            exit 1;;
    esac
done


# ------------------------------------------------------------------------------
# Validate required arguments
if [[ -z "$FILE_STOI" || -z "$SUBCOMPONENT_DIR" || -z "$FILE_CROSSLINK" ]]; then
    echo "[ERROR] Missing required arguments"
    echo "Usage: $0 -stoi <stoi.json> -sub_dir <subcomponent_dir/> -crosslink <crosslink.txt> [-nmax N] [-ncpu N] [-steps N]"
    exit 1
fi

echo ""
echo "[INFO] Stoichiometry file:            $FILE_STOI"
echo "[INFO] Subcomponent directory:        $SUBCOMPONENT_DIR"
echo "[INFO] Crosslinking restraints file:  $FILE_CROSSLINK"
echo -e "\n[INFO] Start running"

START_TIME=$(date +%s)


# ------------------------------------------------------------------------------
# Start modeling
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/" && pwd)"
WORK_DIR=$(dirname "$FILE_STOI")

ASSEMBLY_PATH="$WORK_DIR/assemble_crosslinking"
mkdir "$ASSEMBLY_PATH"

cd "$ASSEMBLY_PATH"
python ${BASE_DIR}/crosslinks/assemble/dataprocess.py \
    --pdbdir $SUBCOMPONENT_DIR \
    --stoi $FILE_STOI
python ${BASE_DIR}/crosslinks/assemble/xl_assemble.py \
    --pdbdir $SUBCOMPONENT_DIR \
    --output assemble.pdb \
    --nmax $NMAX \
    --nround 500 \
    --rmsd 5.0 \
    --population 100 \
    --workers $NUM_CPUS \
    --generations 50 \
    --md $NUM_STEPS \
    --stoi $FILE_STOI \
    --xlinks $FILE_CROSSLINK \
    --xl_threshold1 0.7 \
    --xl_threshold2 0.8

RESULTS_PATH="$WORK_DIR/results/"
mkdir "$RESULTS_PATH"

RESULTS_FILE="$RESULTS_PATH/results.log"
while IFS= read -r line; do
    [[ "$line" == Iteration* ]] || continue

    if [[ "$line" =~ ^Iteration[[:space:]]+([0-9]+)[[:space:]]+model[[:space:]]+([0-9]+)[[:space:]]+XL-satisfaction:[[:space:]]*([^,[:space:]]+),[[:space:]]+XL-scaled[[:space:]]+it_score:[[:space:]]*([^[:space:]]+)[[:space:]]*$ ]]; then
        iteration="${BASH_REMATCH[1]}"
        model_id="${BASH_REMATCH[2]}"
        xl_satisfaction="${BASH_REMATCH[3]}"
        scaled_it_score="${BASH_REMATCH[4]}"
        xl_fmt=$(awk -v x="$xl_satisfaction" 'BEGIN{printf "%.4f", x+0}')
        it_fmt=$(awk -v x="$scaled_it_score"   'BEGIN{printf "%.4f", x+0}')
        echo "Model ${model_id} Crosslink_Satisfaction: ${xl_fmt} Scaled_IT-score: ${it_fmt}" >> "$RESULTS_FILE"
        cp "$ASSEMBLY_PATH/assemble${iteration}_${model_id}.pdb" "$RESULTS_PATH/model_${model_id}.pdb"
    else
        echo "[WARN] Unparsed line: $line" >&2
    fi
done < "$ASSEMBLY_PATH/aa2.log"

END_TIME=$(date +%s)
RUN_TIME=$(($END_TIME - $START_TIME))

echo -e "\n[INFO] HDM program finish, total cost $RUN_TIME s"
