#!/bin/bash
year=$1
Date=$2

outdir=/eos/user/y/youzhiyong/WCDA/WCDARec/$year/$Date
indir=/eos/lhaaso/raw/wcda/$year/$Date

eos mkdir -p $outdir

#ifile=ES.77772.WCDA_EVENT.P110MC16_Z.es-1.20201230235937.246.reduced.root
ifile=$3

cd /workfs/ybj/youzhiyong/WCDA/WCDARec
./main $outdir/${ifile}.rec.root $indir/$ifile
