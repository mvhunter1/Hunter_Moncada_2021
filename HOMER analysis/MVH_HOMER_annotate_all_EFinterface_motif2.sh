#!/bin/bash

# *******************************************************************************

# Experiment Parameters:
INDIR="/Users/hunterm/Documents/bash_scripts/homer/HOMER_singlecell/EF_interface/"
INFILE="EF_interface_markers_all_homer.txt"

MOTIF2="/Users/hunterm/Documents/bash_scripts/homer/HOMER_singlecell/EF_interface/homer_results_20200721_112932_EFinterface_all/homerResults/motif5.motif"


OUTDIR="/Users/hunterm/Documents/bash_scripts/homer/HOMER_singlecell/EF_interface/" # base directory for outputs
TSTAMP=$(date +"homer_results_%Y%m%d_%H%M%S_EFinterface_all_motif5/")
OUTDIR="$OUTDIR$TSTAMP"

mkdir "$OUTDIR"

# *******************************************************************************

# Tool locations:
HOMER="/Applications/homer/bin/"


# *******************************************************************************

echo ""
echo "Searching zebrafish genome for motif5 target genes..."

$HOMER/findMotifs.pl "$INDIR$INFILE" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-find "$MOTIF2" \
	&> "$OUTDIR/out_m5all.txt"

echo "Done!"

