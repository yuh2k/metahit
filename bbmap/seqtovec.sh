#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified Sept 27, 2023

Description:  Generates one-hot 4-bit vectors from sequence.

Usage:  sectovec.sh in=<sequence data> out=<text vectors>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file>       Sequence data.
out=<file>      Vectors in tsv form, with the last column as the result.
overwrite=f     (ow) Set to false to force the program to abort rather than
                overwrite an existing file.

Processing parameters:
width=55        Maximum vector width, in bases; the actual vector size will be
                4+4*width+1 (where the +1 is the desired output).  For
                longer sequences, only the first 'width' bases will be used;
                shorter sequences will be zero-padded.
result=-1       Set to the desired output for all vectors.
parse=f         Set to true to parse the result from sequence headers,
                from a tab-delimited 'result=X' term.
rcomp=f         If true, also output vectors for the reverse-complement.

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

z="-Xmx1g"
z2="-Xms1g"
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

seqtovec() {
	local CMD="java $EA $EOOM $z -cp $CP ml.SequenceToVector $@"
	echo $CMD >&2
	eval $CMD
}

seqtovec "$@"
