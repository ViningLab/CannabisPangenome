#! /bin/env bash

#Input/Output files:
#  --sr SR               Structural annotation mappings (syri.out) identified
#                        by SyRI (default: None)
#  --genomes GENOMES     File containing path to genomes (default: None)
#  -o O                  Output file name. Acceptable format: pdf, png, svg
#                        (default: plotsr.pdf)


#COMP="EH23a.chrX_AH3Ma.chrXrevcomp"
#COMP="AH3Ma.chrXrevcomp_AH3Mb.chrY"
#COMP="AH3Mb.chrY_BCMb.chrY"

# genomes.txt *MUST* be in the correct order and unused genomes are commented out.

#COMP="EH23a.chrXrevcomp_AH3Ma.chrX"
#COMP="AH3Ma.chrX_AH3Mb.chrYrevcomp"
COMP="AH3Mb.chrYrevcomp_BCMb.chrYrevcomp"

CMD="plotsr \
    --sr "$COMP"syri.out \
    --genomes genomes.txt \
    -o "$COMP"_plot.png"

#echo $CMD
#eval $CMD

#  --nosyn               Do not plot syntenic regions (default: False)
#  --noinv               Do not plot inversions (default: False)
#  --notr                Do not plot translocations regions (default: False)
#  --nodup               Do not plot duplications regions (default: False)
#  -s S                  minimum size of a SR to be plotted (default: 10000)


CMD="plotsr \
    -H 3 \
    -W 6 \
    --sr EH23a.chrXrevcomp_AH3Ma.chrXsyri.out \
    --sr AH3Ma.chrX_AH3Mb.chrYrevcompsyri.out \
    --sr AH3Mb.chrYrevcomp_BCMb.chrYrevcompsyri.out \
    --genomes genomes.txt \
    --noinv \
    --notr \
    --nodup \
    -o AH3MXY_plot.png"

#    -v \
#
echo $CMD
#
eval $CMD




CMD="plotsr \
    --sr CM011610.1_CM022973.1syri.out \
    --sr CM022973.1_CM010797.2syri.out \
    --genomes genomes2.txt \
    -o output_plot2.png"

#echo $CMD
#eval $CMD


CMD="plotsr \
    --sr NC_044378.1_CM028020.1syri.out \
    --sr CM028020.1_CM011610.1_complementsyri.out \
    --sr CM011610.1_complement_CM022973.1syri.out \
    --sr CM022973.1_CM010797.2syri.out \
    --genomes genomes3.txt \
    -o output_plot3.png"

#echo $CMD
#eval $CMD



CMD="plotsr \
    --sr GCF_000146045_GCA_000977955syri.out \
    --genomes genomes_yeast.txt \
    -o output_plot_yeast.png"

#echo $CMD
#eval $CMD


