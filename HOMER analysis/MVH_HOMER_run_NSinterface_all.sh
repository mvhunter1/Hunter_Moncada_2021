#!/bin/bash

# *******************************************************************************
# 
# (c) 2019, Nathaniel R. Campbell, Xavier and White Labs, MSKCC
# nac2026@med.cornell.edu | nathaniel.r.campbell@gmail.com
# 
# run HOMER findMotifs.pl on list of differentially expressed genes (upregulated only) from RNA-seq experiment
# 
# 
# *******************************************************************************

# Experiment Parameters:
INDIR="/Users/hunterm/Documents/bash_scripts/homer/HOMER_nucseq/NS_all/" # folder with input files
INFILE="NS_interface_markers_homer.txt" # file with DE genes of interest, ENSDARG ID's in first column

OUTDIR="/Users/hunterm/Documents/bash_scripts/homer/HOMER_nucseq/NS_all/" # base directory for outputs
TSTAMP=$(date +"homer_results_%Y%m%d_%H%M%S_interface_all/") # time stamp for outputs
OUTDIR="$OUTDIR$TSTAMP" # directory for outputs

mkdir "$OUTDIR" # make output directory

# *******************************************************************************

# Tool locations:
HOMER="/Applications/homer/bin/"


# *******************************************************************************


echo ""
echo "Running HOMER on $INFILE"
echo "Saving to $OUTDIR"
echo ""
echo ""

# Run HOMER findMotifs.pl and save stdout to out.log
$HOMER/findMotifs.pl "$INDIR$INFILE" zebrafish "$OUTDIR" \
	-start -500 \
	-end 500 \
	-len 8,10,16 \
	-p 1 \
	&> "$OUTDIR/out.log"

# add file info to out.log
echo "" >> "$OUTDIR/out.log"
echo "" >> "$OUTDIR/out.log"
echo "Ran HOMER on $INFILE" >> "$OUTDIR/out.log"
echo "Saving to $OUTDIR" >> "$OUTDIR/out.log"
echo "" >> "$OUTDIR/out.log"

echo "Done!"
