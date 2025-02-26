#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 15, 2025

Description:  Grades metagenome bins for completeness and contamination.
The contigs should be labeled with their taxID; the header should contain
'tid_X' somewhere where X is a number.


Usage:  gradebins.sh ref=assembly bin*.fa
or
gradebins.sh ref=assembly.fa in=bin_directory

File parameters:
ref=<file>      The original assembly that was binned.
in=<directory>  Location of bin fastas.
hist=<file>     Histogram output.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

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

z="-Xmx2g"
z2="-Xms2g"
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

gradeBins() {
	local CMD="java $EA $EOOM $z -cp $CP bin.GradeBins $@"
	echo $CMD >&2
	eval $CMD
}

gradeBins "$@"
