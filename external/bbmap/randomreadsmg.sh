#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 9, 2025

Description:  Generates synthetic reads from a set of fasta assemblies.
Each assembly is assigned a random coverage level.

Usage:  randomreadsmg.sh *.fa out=reads.fq.gz

File parameters:
in=<file,file>  Assembly input.  Can be a single file, a directory of files,
                or comma-delimited list.  Unrecognized arguments with no '='
                sign will also be treated as input files.
out=<file>      Synthetic read output destination.
out2=<file>     Read 2 output if twin files are desired for paired reads.

Processing parameters:
paired=true     Generate paired reads.
mindepth=1      Minimum assembly average depth.
maxdepth=256    Maximum assembly average depth.
depth=          Sets minimum and maximum to the same level.
depthvariance=0.3  Contigs within an assembly will have coverage that varies
                   by plus or minus this fraction of average.
length=150      Read length.
avginsert=300   Average insert size; only affects paired reads.
threads=        Set the number of threads; default is logical core count.
seed=-1         If positive, use the specified RNG seed.  This will cause
                deterministic output if threads=1.

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

z="-Xmx4g"
z2="-Xms4g"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	setEnvironment
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

generate() {
	local CMD="java $EA $EOOM $z -cp $CP bin.RandomReadsMG $@"
	echo $CMD >&2
	eval $CMD
}

generate "$@"
