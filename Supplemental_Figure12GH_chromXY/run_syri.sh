#!/usr/bin/env bash

#conda activate syri_env
# conda activate syri


#gA="EH23a.chrX"
#gB="AH3Ma.chrXrevcomp"
#gC="AH3Mb.chrY"
#gD="BCMb.chrY"

gA="EH23a.chrXrevcomp"
gB="AH3Ma.chrX"
gC="AH3Mb.chrYrevcomp"
gD="BCMb.chrYrevcomp"



#G1=$gA
#G2=$gB

#G1=$gB
#G2=$gC

#
G1=$gC
#
G2=$gD

fsuffix=".fasta"

# syri -c A_B.bam -r A.fa -q B.fa -F B --prefix A_B &

# -c INFILE [-r REF] [-q QRY] 
# [-F {T,S,B,P}] 
# [--prefix PREFIX]
# [--nc NCORES]


# CMD="syri -c $gA\_$gB.bam -r $gA$fsuffix -q "$gB""$fsuffix" -F B --prefix "$gA"_"$gB
#CMD="syri -c "$gB"_"$gC".bam -r $gB$fsuffix -q $gC$fsuffix -F B --prefix "$gB"_"$gC

CMD="syri --nosr -c "$G1"_"$G2".bam -r $G1$fsuffix -q "$G2$fsuffix" -F B --prefix "$G1"_"$G2

#CMD=$CMD" --tdgaplen 500000 -b 60"
#CMD=$CMD" --tdgaplen 10000 -b 10"




#
echo $CMD
#eval $CMD






