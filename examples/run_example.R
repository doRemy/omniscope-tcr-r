# =============================================================================
# Example: Processing Omniscope TCR sequencing data
# Author: Rémy Pétremand
#
# This script demonstrates the full processing pipeline from raw Omniscope
# TCR output files to cell-wise and contig-wise tables.
#
# Replace the file paths below with your own Omniscope output files.
# =============================================================================

# Source the pipeline functions (run from the repo root)
source("R/omniscope_tcr_r.R")

# -----------------------------------------------------------------------------
# 1. Set file paths
#    Inputs are CSVs exported from an Omniscope TCR run.
#    Each argument also accepts a directory path — CSV files will be
#    auto-discovered recursively.
# -----------------------------------------------------------------------------

file.name.contigs <- "path/to/contigs.csv"   # cell-level contig table
file.name.counts  <- "path/to/counts.csv"    # contig-level UMI counts
file.name.pairs   <- "path/to/pairs.csv"     # alpha-beta chain pairing

# -----------------------------------------------------------------------------
# 2. Load and QC-filter TCR data
#
#    Returns a list with three elements:
#      $cellwise   — one row per cell × chain
#      $contigwise — one row per unique contig_id (with counts and pairing)
#      $summary    — QC counts (raw, filtered, merged)
# -----------------------------------------------------------------------------

TCR.table.list <- load_omniscope_tcr(
  input.contigs           = file.name.contigs,
  input.counts            = file.name.counts,
  input.pairs             = file.name.pairs,
  filter_nonproductive    = TRUE,
  filter_stop_codon       = TRUE,
  filter_out_of_frame     = TRUE,
  require_complete_vdj    = FALSE,
  remove_multi_chain      = FALSE,
  resolve_ambiguous_genes = "none",
  verbose                 = TRUE
)

# Retrieve the output tables
df.cellwise   <- TCR.table.list$cellwise    # one row per cell × chain
df.contigwise <- TCR.table.list$contigwise  # one row per unique contig
df.qc         <- TCR.table.list$summary     # QC summary

# Quick look
head(df.cellwise)
head(df.contigwise)
print(df.qc)
