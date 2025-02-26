#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified March 11, 2024

Description:  Plots Illumina read positions and barcode hamming distance.

Usage:  plotreadposition.sh in=<file.fq> out=<file.tsv> expected=<barcodes.txt>

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx300m"
z2="-Xms300m"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
}
calcXmx "$@"

function plot() {
	local CMD="java $EA $EOOM $z $z2 -cp $CP hiseq.PlotReadPosition $@"
	echo $CMD >&2
	eval $CMD
}

plot "$@"
