Quartet2RIdeogram
A pipeline to generate publication-quality chromosome ideograms from Quartet assembly finishing results using the RIdeogram R package.

This tool automates the process of parsing Tandem Repeat (TR) data from Quartet, identifying gaps and telomeres from the genome FASTA itself, and integrating them with genome annotation data (GFF) to produce a comprehensive and visually appealing chromosomal map.

(This is a placeholder image. The actual output will be a detailed chromosome map.)

Features
Automated Workflow: A single command runs the entire pipeline from raw inputs to final plot.

Intelligent Centromere Detection: Combines low gene density regions (from GFF) with high tandem repeat (TR) coverage regions (from Quartet) to accurately identify centromeres.

Direct FASTA analysis: Gaps (stretches of 'N's) and telomeric repeats are identified directly from the genome assembly, removing dependencies on other tools.

Comprehensive Visualization: Plots key chromosomal features including:

Chromosome length

Centromere positions

Telomere locations

Assembly gaps

Gene density overlaid as a heatmap

Customizable & Extendable: The modular script-based design makes it easy to modify or extend the pipeline.

Installation
1. Clone the repository
git clone [https://github.com/your-username/quartet2rideogram.git](https://github.com/your-username/quartet2rideogram.git)
cd quartet2rideogram

2. Create a Conda Environment (Recommended)
It is highly recommended to use Conda to manage dependencies.

# Create and activate a new conda environment
conda create -n q2r python=3.9 r-base=4.2 -y
conda activate q2r

3. Run the installation script
The install.sh script will install all required Python and R packages.

bash install.sh

This will install:

Python packages: pandas

R packages: RIdeogram, stringr, getopt

Usage
The main pipeline is executed via the quartet2rideogram.sh script.

Command
bash quartet2rideogram.sh \\
    -f <genome.fasta> \\
    -g <annotation.gff> \\
    -c <candidates_dir> \\
    -o <output_prefix> \\
    --telotype <plant|animal>

Arguments
-f, --fasta: Path to the chromosome-level genome assembly FASTA file. (Required)

-g, --gff: Path to the genome annotation GFF3 file. (Required)

-c, --candidates: Path to the Quartet 'Candidates' directory containing *.candidate files for TR analysis. (Required)

-o, --output_prefix: The prefix for all output files (e.g., Hap1_ideogram). (Required)

--telotype: Type of telomere repeat to search for ('plant' or 'animal'). Default is plant. (Optional)

-h, --help: Display the help message.

Example
bash quartet2rideogram.sh \\
  -f /data/my_genome/hap1.fasta \\
  -g /data/my_genome/hap1.gff3 \\
  -c /data/quartet_results/hap1/Candidates/ \\
  -o Hap1_ChrMap \\
  --telotype plant

Output Files
After a successful run, the following files will be generated in your working directory:

Hap1_ChrMap_ideogram.svg: The final chromosome ideogram in SVG format (editable).

Hap1_ChrMap_ideogram.png: The final chromosome ideogram in PNG format (publication-ready).

_intermediate_files/: A directory containing all intermediate files.

Workflow Overview
Generate Base Karyotype: Chromosome names and lengths are extracted from the input FASTA file.

Identify Centromeres: Intersects low gene-density regions (from GFF) and high-TR regions (from Quartet *.candidate files) to define centromeres.

Find Gaps & Telomeres: The pipeline directly scans the FASTA file to locate assembly gaps (Ns) and terminal telomeric repeat sequences.

Prepare RIdeogram Inputs: Formats all extracted data into files compatible with RIdeogram.

Plot Ideogram: An R script generates the final SVG and PNG plots.

License
This project is licensed under the MIT License.