#!/bin/bash

# *******************************************************************************

# Experiment Parameters:
INDIR="/Users/hunterm/Documents/bash_scripts/homer/HOMER_Visium/ABC/ABC_interface/"
INFILE_D="interface_markers_down_ABC.txt"
INFILE_U="interface_markers_up_ABC.txt"
MOTIF1u="/Users/hunterm/Documents/bash_scripts/homer/HOMER_Visium/ABC/ABC_interface/homer_results_20200409_080953_ABCi_up/homerResults/motif1.motif"
MOTIF1d="/Users/hunterm/Documents/bash_scripts/homer/HOMER_Visium/ABC/ABC_interface/homer_results_20200409_091526_ABCi_down/homerResults/motif1.motif"

OUTDIR="/Users/hunterm/Desktop/" # base directory for outputs
TSTAMP=$(date +"homer_results_%Y%m%d_%H%M%S_ABCi_motif1/")
OUTDIR="$OUTDIR$TSTAMP"

mkdir "$OUTDIR"

# *******************************************************************************

# Tool locations:
HOMER="/Applications/homer/bin/"


# *******************************************************************************



$HOMER/findMotifs.pl "$INDIR$INFILE_D" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-find "$MOTIF1d" \
	&> "$OUTDIR/out_m1d_down.log"

$HOMER/findMotifs.pl "$INDIR$INFILE_U" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-find "$MOTIF1d" \
	&> "$OUTDIR/out_m1d_up.log"


$HOMER/findMotifs.pl "$INDIR$INFILE_D" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-find "$MOTIF1u" \
	&> "$OUTDIR/out_m1u_down.log"

$HOMER/findMotifs.pl "$INDIR$INFILE_U" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-find "$MOTIF1u" \
	&> "$OUTDIR/out_m1u_up.log"