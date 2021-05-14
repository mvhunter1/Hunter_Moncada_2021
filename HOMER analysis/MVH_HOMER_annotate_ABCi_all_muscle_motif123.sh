#!/bin/bash

# *******************************************************************************

# Experiment Parameters:
INDIR="/Users/hunterm/Documents/bash_scripts/homer/HOMER_Visium/ABC/ABC_interface/"
INFILE="interface_markers_all_ABC.txt"
MOTIF1="/Users/hunterm/Documents/bash_scripts/homer/HOMER_Visium/ABC/ABC_interface/compare_w_muscle/homer_results_20200414_114354_ABCi_all_muscle/homerResults/motif1.motif"
MOTIF2="/Users/hunterm/Documents/bash_scripts/homer/HOMER_Visium/ABC/ABC_interface/compare_w_muscle/homer_results_20200414_114354_ABCi_all_muscle/homerResults/motif2.motif"
MOTIF3="/Users/hunterm/Documents/bash_scripts/homer/HOMER_Visium/ABC/ABC_interface/compare_w_muscle/homer_results_20200414_114354_ABCi_all_muscle/homerResults/motif3.motif"


OUTDIR="/Users/hunterm/Desktop/" # base directory for outputs
TSTAMP=$(date +"homer_results_%Y%m%d_%H%M%S_ABCi_all_muscle_motif123/")
OUTDIR="$OUTDIR$TSTAMP"

mkdir "$OUTDIR"

# *******************************************************************************

# Tool locations:
HOMER="/Applications/homer/bin/"


# *******************************************************************************


$HOMER/findMotifs.pl "$INDIR$INFILE" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-find "$MOTIF1" \
	&> "$OUTDIR/out_m1all.log"

$HOMER/findMotifs.pl "$INDIR$INFILE" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-find "$MOTIF2" \
	&> "$OUTDIR/out_m2all.log"

$HOMER/findMotifs.pl "$INDIR$INFILE" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-find "$MOTIF3" \
	&> "$OUTDIR/out_m3all.log"

