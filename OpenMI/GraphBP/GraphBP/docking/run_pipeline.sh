#!/bin/bash
set -e  # Exit on error
set -u  # Error on undefined variables


# === Source configuration file ===
source "$(dirname "${BASH_SOURCE[0]}")/dock_with_vina.sh"
echo "$(dirname "${BASH_SOURCE[0]}")/dock_with_vina.sh"


# === Set script directory ===
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "Script directory: $SCRIPT_DIR"

BASE_DIR="$SCRIPT_DIR/$RECEPTOR_PREFIX/experiment_epoch_${EPOCH}_mols_${MOLS}_bs_${BS}_aurora_${AURORA}"
LIGAND_DIR="$BASE_DIR/ligands"
PREPARED_LIG_DIR="$BASE_DIR/prepared_ligands"
PREPARED_REC_DIR="$SCRIPT_DIR/$RECEPTOR_PREFIX"
OUTPUT_DIR="$BASE_DIR/vina_outputs"
RESULTS_FILE="$BASE_DIR/vina_results.csv"

mkdir -p "$PREPARED_LIG_DIR" "$OUTPUT_DIR"

# === Prepare receptor ===
echo "Preparing receptor..."
mk_prepare_receptor.py -i "$RECEPTOR" -o "$PREPARED_REC_DIR/$RECEPTOR_PREFIX" -p -v -g --box_size $SIZE --box_center $CENTER --allow_bad_res


echo "Running autogrid..."
cd "$PREPARED_REC_DIR" || exit 1
autogrid4 -p "$RECEPTOR_PREFIX.gpf" -l "$RECEPTOR_PREFIX.glg"
cd -  # go back to the previous directory

# === Start results CSV ===
echo "ligand,affinity_kcal/mol" > "$RESULTS_FILE"

# === Process each ligand ===
echo "Processing ligands..."
shopt -s nullglob
sdf_files=("$LIGAND_DIR"/*.sdf)
if [ ${#sdf_files[@]} -eq 0 ]; then
    echo "No .sdf files found in $LIGAND_DIR"
    exit 1
fi
for sdf in "${sdf_files[@]}"; do
    base=$(basename "$sdf" .sdf)
    echo "Docking ligand: $base"

    scrubbed="$PREPARED_LIG_DIR/${base}_scrubbed.sdf"
    pdbqt="$PREPARED_LIG_DIR/${base}_scrubbed.pdbqt"
    out="$OUTPUT_DIR/${base}_out.pdbqt"

    # Scrub ligand
    scrub.py "$sdf" -o "$scrubbed"

    # Prepare ligand
    if ! mk_prepare_ligand.py -i "$scrubbed" -o "$pdbqt"; then
        echo "Failed to prepare $base"
        continue
    fi

    # Run Vina
    if ! vina --ligand "$pdbqt" --maps "$PREPARED_REC_DIR/$RECEPTOR_PREFIX" --scoring ad4 --exhaustiveness "$EXHAUST" --out "$out"; then
        echo "Vina failed for $base"
        continue
    fi

    # if ! vina --receptor "$PREPARED_REC_DIR/$RECEPTOR_PREFIX".pdbqt --ligand "$pdbqt" \
    #    --config "$PREPARED_REC_DIR/$RECEPTOR_PREFIX".box.txt \
    #    --exhaustiveness=8 --out "$out"; then
    #     echo "Vina failed for $base"
    #     continue
    # fi

    # Extract best binding score
    score=$(grep "^REMARK VINA RESULT:" "$out" | awk '{print $4}')
    echo "$base,$score" >> "$RESULTS_FILE"
    echo "$base: $score kcal/mol"
done

echo "All docking complete! Results saved to $RESULTS_FILE"
