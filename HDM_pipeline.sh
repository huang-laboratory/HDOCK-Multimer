# Script for running the entire HDM program

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
        -mono_dir|--mono_dir)
            MONO_DIR=$(readlink -f "$2")
            shift 2;;
        -sub_dir|--sub_dir)
            SUBCOMPONENT_DIR=$(readlink -f "$2")
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
            echo "Usage: $0 -stoi <stoi.json> -mono_dir <mono_dir/> -sub_dir <subcomponent_dir/> [-nmax N] [-ncpu N]"
            echo ""
            echo "Required arguments:"
            echo "    -stoi     : Path to stoichiometry JSON file"
            echo "    -mono_dir : Directory containing monomer structure PDB files"
            echo "    -sub_dir  : Directory containing subcomponent structure PDB files"
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
if [[ -z "$FILE_STOI" || -z "$MONO_DIR" || -z "$SUBCOMPONENT_DIR" ]]; then
    echo "[ERROR] Missing required arguments"
    echo "Usage: $0 -stoi <stoi.json> -mono_dir <mono_dir/> -sub_dir <subcomponent_dir/> [-nmax N] [-ncpu N] [-steps N]"
    exit 1
fi

echo ""
echo "[INFO] Stoichiometry file:        $FILE_STOI"
echo "[INFO] Monomer dir:               $MONO_DIR"
echo "[INFO] Subcomponent file dir:     $SUBCOMPONENT_DIR"
echo -e "\n[INFO] Start running"

START_TIME=$(date +%s)


# ------------------------------------------------------------------------------
# Determinate modeling strategies
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/" && pwd)"
chmod +x ${BASE_DIR}/tools/transchain ${BASE_DIR}/tools/ppscore
JQ=${BASE_DIR}/tools/jq
WORK_DIR=$(dirname "$FILE_STOI")

if [[ ! -f "$WORK_DIR/res_indices_D.txt" ]]; then
    python ${BASE_DIR}/scripts/filter_monofiles.py $FILE_STOI $MONO_DIR
fi

entry_count=$($JQ length "$FILE_STOI")

declare -A stoi_chain_num

hetero_flag=0
if (( entry_count > 1 )); then
    hetero_flag=1
    chain_num=0
    for ((i=0; i<entry_count; i++)); do
        key=$($JQ -r ".[$i].asym_ids[0]" "$FILE_STOI")
        val=$($JQ ".[$i].asym_ids | length" "$FILE_STOI")
        stoi_chain_num[$key]=$val
        ((chain_num += val))
    done
elif (( entry_count == 1 )); then
    key=$($JQ -r ".[0].asym_ids[0]" "$FILE_STOI")
    val=$($JQ ".[0].asym_ids | length" "$FILE_STOI")
    stoi_chain_num[$key]=$val
    chain_num=$val
else
    echo "[ERROR] Wrong stoi file, please check $FILE_STOI" >&2
    exit 1
fi

if [[ ! -f "$WORK_DIR/build_types.txt" ]]; then
    python ${BASE_DIR}/scripts/check_build_type.py $FILE_STOI
fi 
if [[ ! -f "$WORK_DIR/build_types.txt" ]]; then
    echo "[ERROR] Determination of modeling strategies failed" >&2
    exit 1
fi

modeling_strategies=()
while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%"${line##*[![:space:]]}"}"
    line="${line#"${line%%[![:space:]]*}"}"
    modeling_strategies+=("$line")
done < "$WORK_DIR/build_types.txt"


# ------------------------------------------------------------------------------
# Start modeling
for method in "${modeling_strategies[@]}"; do
    if [[ $method == "asym" ]]; then
        echo -e "\n************************************************************"
        echo -e "[INFO] Start multi-body docking"

        METHOD_PATH="$WORK_DIR/$method/"
        cd $METHOD_PATH
        if (( hetero_flag )); then
            chains=""
            copies=""
            for key in "${!stoi_chain_num[@]}"; do
                chains+=" ${key}.pdb"
                copies+=" ${stoi_chain_num[$key]}"
            done
            if [[ -d "pairwise_out/" ]]; then
                python ${BASE_DIR}/scripts/multidock/multi_dock.py $chains -s $copies -i -p 2 -m 1 -r $chain_num -j $NUM_CPUS -d -b pairwise_out/
            else
                python ${BASE_DIR}/scripts/multidock/multi_dock.py $chains -s $copies -i -p 2 -m 1 -r $chain_num -j $NUM_CPUS -d
            fi
        else
            for file in A*.pdb; do
                filename="${file##*/}"
                if [[ "$file" == *_* ]]; then
                    prefix="${filename%%_*}"
                    csym="${prefix#A}"
                else
                    csym=1
                fi
            done
            copy=$(($chain_num / $csym))
            if [[ -d "pairwise_out" ]]; then
                python ${BASE_DIR}/scripts/multidock/multi_dock.py $filename -s $copy -i -p 2 -m 1 -r $chain_num -j $NUM_CPUS -d -b pairwise_out/
            else
                python ${BASE_DIR}/scripts/multidock/multi_dock.py $filename -s $copy -i -p 2 -m 1 -r $chain_num -j $NUM_CPUS -d
            fi
        fi
        python ${BASE_DIR}/scripts/multidock/refine_multi_ids.py $FILE_STOI

    elif [[ $method == "assemble" ]]; then
        echo -e "\n************************************************************"
        echo -e "[INFO] Start assembly"

        METHOD_PATH="$WORK_DIR/$method/"
        cd $METHOD_PATH
        if [[ -f "interfaces.csv" && -f "res_indices_A.txt" ]]; then
            python ${BASE_DIR}/scripts/assemble/assemble.py --pdbdir $SUBCOMPONENT_DIR --output assemble.pdb --nmax $NMAX --nround 500 --rmsd 5.0 --population 100 --workers $NUM_CPUS --generations 50 --md $NUM_STEPS --stoi $FILE_STOI
        else
            python ${BASE_DIR}/scripts/assemble/dataprocess.py --pdbdir $SUBCOMPONENT_DIR --stoi $FILE_STOI
            python ${BASE_DIR}/scripts/assemble/assemble.py --pdbdir $SUBCOMPONENT_DIR --output assemble.pdb --nmax $NMAX --nround 500 --rmsd 5.0 --population 100 --workers $NUM_CPUS --generations 50 --md $NUM_STEPS --stoi $FILE_STOI
        fi

        if [[ -d "dock/" ]]; then
            cd dock/
            if [[ -f "warns.log" ]]; then
                echo "[ERROR] Split subunits are not connected successfully"
            else
                python ${BASE_DIR}/scripts/multidock/multi_dock.py *.pdb -i -p 2 -m 1 -r $chain_num -j $NUM_CPUS -d
            fi
        fi
    
    elif [[ $method == "cubicT" ]]; then
        echo -e "\n************************************************************"
        echo -e "[INFO] Start Tetrahedral symmetric docking"

        METHOD_PATH="$WORK_DIR/$method/"
        cd $METHOD_PATH
        python ${BASE_DIR}/scripts/symmdock/Tsym/cubicTsym.py --n2 10 --n3 20 --nmax $NMAX --workers $NUM_CPUS --output Tsym.pdb --md 0

    elif [[ $method == "cubicO" ]]; then
        echo -e "\n************************************************************"
        echo -e "[INFO] Start Octahedral symmetric docking"

        METHOD_PATH="$WORK_DIR/$method/"
        cd $METHOD_PATH
        python ${BASE_DIR}/scripts/symmdock/Osym/cubicOsym.py --n2 10 --n3 10 --n4 10 --nmax $NMAX --workers $NUM_CPUS --output Osym.pdb --md $NUM_STEPS

    elif [[ $method == c* ]]; then
        csym=${method:1}
        echo -e "\n************************************************************"
        echo -e "[INFO] Start C$csym symmetric docking"

        METHOD_PATH="$WORK_DIR/$method/"
        cd $METHOD_PATH
        python ${BASE_DIR}/scripts/symmdock/csym/csymdock.py $chain_num $csym $NMAX

    elif [[ $method == d* ]]; then
        dsym=${method:1}
        echo -e "\n************************************************************"
        echo -e "[INFO] Start D$dsym symmetric docking"

        METHOD_PATH="$WORK_DIR/$method/"
        cd $METHOD_PATH
        python ${BASE_DIR}/scripts/symmdock/dsym/dsymdock.py $dsym $NMAX

    else
        echo "[ERROR] Unsupported modeling strategy !!!" >&2
        exit 1
    fi
done 

# ------------------------------------------------------------------------------
# Models clustering
python ${BASE_DIR}/scripts/clustermodels.py $FILE_STOI $NMAX


END_TIME=$(date +%s)
RUN_TIME=$(($END_TIME - $START_TIME))

echo -e "\n[INFO] HDM program finish, total cost $RUN_TIME s"

