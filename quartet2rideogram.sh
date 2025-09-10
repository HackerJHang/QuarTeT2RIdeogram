#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.

# --- Help Message ---
function show_help() {
    echo "
Usage: $(basename "$0") -f <fasta> -g <gff> -c <candidates_dir> -o <prefix> [--telotype <type>]

A pipeline to generate chromosome ideograms from FASTA, GFF, and Quartet TR results.

Required arguments:
  -f, --fasta         Path to the chromosome-level genome assembly FASTA file.
  -g, --gff           Path to the genome annotation GFF3 file.
  -c, --candidates    Path to the Quartet 'Candidates' directory containing *.candidate files.
  -o, --output_prefix The prefix for all output files (e.g., Hap1_ideogram).

Optional arguments:
  --telotype          Type of telomere repeat to search for ('plant' or 'animal'). Default: plant.
  -h, --help          Display this help message and exit.
"
}

# --- Argument Parsing ---
if [ $# -eq 0 ]; then
    show_help
    exit 1
fi

TELO_TYPE="plant" # Default value

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--fasta)
            FASTA_FILE="$2"
            shift 2
            ;;
        -g|--gff)
            GFF_FILE="$2"
            shift 2
            ;;
        -c|--candidates)
            CANDIDATES_DIR="$2"
            shift 2
            ;;
        -o|--output_prefix)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        --telotype)
            TELO_TYPE="$2"
            if [[ "$TELO_TYPE" != "plant" && "$TELO_TYPE" != "animal" ]]; then
                echo "Error: --telotype must be 'plant' or 'animal'." >&2
                exit 1
            fi
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# --- Validate Inputs ---
if [ -z "$FASTA_FILE" ] || [ -z "$GFF_FILE" ] || [ -z "$CANDIDATES_DIR" ] || [ -z "$OUTPUT_PREFIX" ]; then
    echo "Error: All required arguments (-f, -g, -c, -o) must be provided."
    show_help
    exit 1
fi

if [ ! -f "$FASTA_FILE" ]; then echo "Error: Fasta file not found: $FASTA_FILE"; exit 1; fi
if [ ! -f "$GFF_FILE" ]; then echo "Error: GFF file not found: $GFF_FILE"; exit 1; fi
if [ ! -d "$CANDIDATES_DIR" ]; then echo "Error: Candidates directory not found: $CANDIDATES_DIR"; exit 1; fi


# --- Setup Environment ---
SCRIPT_DIR=$(dirname "$(realpath "$0")")
INTERMEDIATE_DIR="_intermediate_files"
mkdir -p "$INTERMEDIATE_DIR"
echo "Intermediate files will be stored in: $INTERMEDIATE_DIR"

# Define file paths
BASE_KARYOTYPE="${INTERMEDIATE_DIR}/${OUTPUT_PREFIX}_base_karyotype.txt"
CHROME_MAPPING="${INTERMEDIATE_DIR}/${OUTPUT_PREFIX}_chr_mapping.txt"
CENTROMERE_SUMMARY="${INTERMEDIATE_DIR}/${OUTPUT_PREFIX}_centromere_summary.txt"
FINAL_KARYOTYPE="${INTERMEDIATE_DIR}/${OUTPUT_PREFIX}_final_karyotype.txt"
GAP_LABELS="${INTERMEDIATE_DIR}/${OUTPUT_PREFIX}_gap_labels.txt"
TELOMERE_LABELS="${INTERMEDIATE_DIR}/${OUTPUT_PREFIX}_telomere_labels.txt"
COMBINED_LABELS="${INTERMEDIATE_DIR}/${OUTPUT_PREFIX}_combined_labels.txt"


# --- Pipeline Execution ---
echo -e "\n--- [Step 1/7] Generating base karyotype from FASTA ---"
python3 "${SCRIPT_DIR}/scripts/01_generate_karyotype.py" "$FASTA_FILE" "$BASE_KARYOTYPE"

echo -e "\n--- [Step 2/7] Creating chromosome name mapping (e.g. scaffold_1 -> CHR1) ---"
# Uses natural sort (-V) to handle chr1, chr2, chr10 correctly.
cut -f1 "$BASE_KARYOTYPE" | tail -n +2 | sort -V | awk -v OFS='\t' '{print $0, "CHR" NR}' > "$CHROME_MAPPING"
echo "Mapping file created at: $CHROME_MAPPING"

echo -e "\n--- [Step 3/7] Finding centromeres using GFF and TR data ---"
python3 "${SCRIPT_DIR}/scripts/02_find_centromeres.py" "$BASE_KARYOTYPE" "$GFF_FILE" "$CANDIDATES_DIR" "$CENTROMERE_SUMMARY"

echo -e "\n--- [Step 4/7] Merging centromere data into final karyotype ---"
python3 "${SCRIPT_DIR}/scripts/05_merge_karyotype.py" "$BASE_KARYOTYPE" "$CENTROMERE_SUMMARY" "$FINAL_KARYOTYPE"

echo -e "\n--- [Step 5/7] Identifying gaps and telomeres directly from FASTA ---"
python3 "${SCRIPT_DIR}/scripts/03_parse_gaps.py" "$FASTA_FILE" "$GAP_LABELS"
python3 "${SCRIPT_DIR}/scripts/04_parse_telomeres.py" "$FASTA_FILE" "$TELOMERE_LABELS" --telotype "$TELO_TYPE"

# Combine label files for RIdeogram
cat "$GAP_LABELS" <(tail -n +2 "$TELOMERE_LABELS" 2>/dev/null) > "$COMBINED_LABELS"
echo "Combined label file created at: $COMBINED_LABELS"

echo -e "\n--- [Step 6/7] Calculating gene density for plotting ---"
# This step is handled inside the R script which calls GFFex

echo -e "\n--- [Step 7/7] Plotting the final ideogram ---"
# The R script will use the mapping file to rename chromosomes and will use default colors.
Rscript "${SCRIPT_DIR}/scripts/06_plot_ideogram.R" \
  --karyotype "$FINAL_KARYOTYPE" \
  --gff "$GFF_FILE" \
  --labels "$COMBINED_LABELS" \
  --mapping "$CHROME_MAPPING" \
  --output "${OUTPUT_PREFIX}_ideogram.svg"

echo -e "\n--- Pipeline Finished Successfully! ---"
echo "Output SVG: ${OUTPUT_PREFIX}_ideogram.svg"
echo "Output PNG: ${OUTPUT_PREFIX}_ideogram.png"
echo "-----------------------------------------"

