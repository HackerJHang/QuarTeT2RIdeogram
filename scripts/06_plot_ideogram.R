#!/usr/bin/env Rscript

# Load necessary libraries
# These packages should be installed via the install.sh script.
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(RIdeogram))
suppressPackageStartupMessages(library(stringr))
# 显式加载 rsvg 包，以便更好地控制PDF转换
suppressPackageStartupMessages(library(rsvg))

# --- Argument Specifications ---
# Defines the command-line arguments the script accepts.
spec <- matrix(c(
  'karyotype', 'k', 1, "character", "Path to the final karyotype file with original names.",
  'gff',       'g', 1, "character", "Path to the GFF annotation file.",
  'labels',    'l', 1, "character", "Path to the combined labels file (gaps, telomeres).",
  'mapping',   'm', 1, "character", "Path to the chromosome name mapping file (original -> new).",
  'output',    'o', 1, "character", "Path for the output SVG file.",
  'help',      'h', 0, "logical",   "Display this help message."
), byrow=TRUE, ncol=5)

# --- Parse Arguments ---
# This block reads arguments from the command line.
tryCatch({
  opt <- getopt(spec)
}, error = function(e) {
  message("Error: Missing required arguments. Use -h or --help for usage information.")
  message(getopt(spec, usage=TRUE))
  quit(status=1)
})

# Display help message and exit if -h or --help is used.
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  quit(status=0)
}

# --- Load Data ---
print("--- Step 1: Loading all data files into R ---")
karyotype_data <- read.table(opt$karyotype, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
label_data     <- read.table(opt$labels, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
name_mapping   <- read.table(opt$mapping, sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("original", "new"))

# --- Calculate Gene Density ---
# GFFex needs original chromosome names to match the GFF file.
print("--- Step 2: Extracting gene density from GFF ---")
gene_density <- GFFex(
  input     = opt$gff,
  karyotype = opt$karyotype, # Use original karyotype path for GFFex
  feature   = "gene",
  window    = 100000
)

# --- Rename Chromosomes for Plotting ---
print("--- Step 3: Standardizing chromosome names for plotting (e.g., scaffold_1 -> CHR1) ---")

# A helper function to apply the renaming based on the mapping file.
rename_chr <- function(df, mapping) {
  df$Chr <- mapping$new[match(df$Chr, mapping$original)]
  df <- df[!is.na(df$Chr), ] # Remove rows that couldn't be mapped
  return(df)
}

# Apply renaming to all data frames that will be used for plotting.
karyotype_data <- rename_chr(karyotype_data, name_mapping)
label_data     <- rename_chr(label_data, name_mapping)
gene_density   <- rename_chr(gene_density, name_mapping)

# The karyotype data must be sorted numerically by the new CHR names for correct plotting order.
karyotype_data <- karyotype_data[order(as.numeric(gsub("CHR", "", karyotype_data$Chr))), ]

print("Chromosome renaming complete.")

# --- Plot Final Ideogram ---
print("--- Step 4: Plotting the ideogram (using default colors for gene density)... ---")
ideogram(
  karyotype  = karyotype_data,
  overlaid   = gene_density,
  label      = label_data,
  label_type = "marker",
  # By not specifying a 'colorset1' parameter, RIdeogram uses its default heatmap colors.
  output     = opt$output,
  
  # --- 更新参数：将图例固定在左上角 ---
  Lx = 0,
  Ly = 0
)

# --- Convert to PNG ---
print("--- Step 5: Converting SVG to PNG format ---")
convertSVG(
  svg    = opt$output,
  device = "png",
  dpi    = 300
)

# --- 修正步骤：转换成 PDF 并指定正确的文件名 ---
print("--- Step 6: Converting SVG to PDF format ---")
# 明确定义输出的PDF文件名
pdf_output_path <- sub("\\.svg$", ".pdf", opt$output)
# 直接调用rsvg_pdf来确保文件名被正确设置，避免生成默认的Rplots.pdf
rsvg_pdf(svg = opt$output, file = pdf_output_path)


print(paste("--- Plotting finished successfully! ---"))

