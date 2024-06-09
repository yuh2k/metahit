#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified November 14, 2023
Description:  Shreds sequences into shorter, possibly overlapping sequences.

Usage: shred.sh in=<file> out=<file> length=<int>

in=<file>       Input sequences.
out=<file>      Destination of output shreds.
length=500      Desired length of shreds if a uniform length is desired.
minlen=-1       Shortest allowed shred.  The last shred of each input sequence
                may be shorter than desired length if this is not set.
maxlen=-1       Longest shred length.  If minlength and maxlength are both
                set, shreds will use a random flat length distribution.
overlap=0       Amount of overlap between successive shreds.
reads=-1        If nonnegative, stop after this many input sequences.
equal=f         Shred each sequence into subsequences of equal size of at most
                'length', instead of a fixed size.
qfake=30        Quality score, if using fastq output.

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

z="-Xmx4000m"
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

stats() {
	local CMD="java $EA $EOOM $z -cp $CP jgi.Shred $@"
#	echo $CMD >&2
	eval $CMD
}

stats "$@"
