Quartet2RIdeogram

<img width="1144" height="1058" alt="image" src="https://github.com/user-attachments/assets/422849c6-fc2f-4f89-95fe-eb35cdfb063a" />



A bioinformatic pipeline to automate the creation of publication-quality chromosome ideograms from Quartet genome assembly results.

This pipeline takes a chromosome-level FASTA file, a GFF annotation file, and the Tandem Repeat (TR) candidates directory from Quartet. It intelligently identifies centromeres based on gene-poor regions and high TR coverage, parses gaps and telomeres directly from the FASTA sequence, and generates a beautiful, publication-ready chromosome map using the R package RIdeogram.

Features
	•	Automated Workflow: A single command runs the entire pipeline from start to finish.
	•	Intelligent Centromere Prediction: Combines gene density information (from GFF) and Tandem Repeat coverage (from Quartet) to accurately locate centromeric regions.
	•	Direct FASTA Analysis: No need for Quartet’s .stat or .telo.info files. Gaps (Ns) and telomeres are identified directly from your genome assembly.
	•	Customizable Telomere Types: Supports both plant (TTTAGGG) and animal (TTAGGG) telomere repeat patterns.
	•	Standardized Naming: Automatically renames chromosomes from any input format (e.g., scaffold_1, group1) to a standard CHR1, CHR2, etc., format for clear visualization.
	•	Publication-Ready Output: Generates SVG (vector), PNG, and PDF formats.

Installation

This pipeline is designed to run on a Linux-based system with Conda installed.

Clone the repository
**
cd Quartet2RIdeogram**

Create a Conda environment and install dependencies:

The provided **install.sh** script will create a dedicated Conda environment named q2r and install all necessary Python and R packages.

# Give the script execution permission
**chmod +x install.sh**

# Run the installer
**./install.sh**

Usage

The entire pipeline is run using the quartet2rideogram.sh script.

Command

bash quartet2rideogram.sh -f <genome.fasta> -g <annotation.gff> -c <candidates_dir> -o <output_prefix> [--telotype <plant|animal>]

Parameters
	•	-f, --fasta (Required): Path to the chromosome-level genome FASTA file.
	•	-g, --gff (Required): Path to the corresponding GFF3 annotation file.
	•	-c, --candidates (Required): Path to the Quartet Candidates directory, which contains the tandem repeat .candidate files.
	•	-o, --output (Required): Prefix for all output files (e.g., Hap1_map).
	•	--telotype (Optional): The type of telomere repeat sequence to search for. Can be plant (default) or animal.

Example

# Activate the conda environment first
**conda activate q2r**

# Run the pipeline
**bash quartet2rideogram.sh \
  -f /path/to/your.fasta \
  -g /path/to/your.gff \
  -c /path/to/your/quartet_results/hap1/Candidates/ \
  -o _Chromosome_Map \
  --telotype plant**

This will generate Hap1_Chromosome_Map_ideogram.svg, .png, and .pdf in the current directory, along with an _intermediate_files folder containing all the intermediate data.

Pipeline Steps

The main script automates the following steps:
	1.	01_generate_karyotype.py: Creates a basic karyotype file (Chr, Start, End) from the FASTA headers.
	2.	02_find_centromeres.py: The core logic. It first finds gene-poor regions from the GFF, then cross-references them with high Tandem Repeat coverage regions from the Quartet .candidate files to define centromeres.
	3.	03_parse_gaps.py: Scans the FASTA file for stretches of Ns and defines them as gaps.
	4.	04_parse_telomeres.py: Scans the ends of each chromosome in the FASTA file for telomeric repeat sequences.
	5.	05_merge_karyotype.py: Combines the basic karyotype with the identified centromere locations.
	6.	06_plot_ideogram.R: Takes all the generated data files and uses RIdeogram to plot the final chromosome map.

License

This project is licensed under the GNU GPLv3 License - see the LICENSE file for details.
