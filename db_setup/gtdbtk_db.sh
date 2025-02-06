#!/usr/bin/env bash
#
# Function: Download and configure the GTDB-Tk database
# Usage: bash download_gtdbtk_db.sh [DB_DIR] [RELEASE_VERSION]
#       If DB_DIR is not specified, it will be downloaded to the current directory's "database" folder by default
#       If RELEASE_VERSION is not specified, version 207 will be downloaded by default
#

DB_DIR=${1:-database}
RELEASE_VERSION=${2:-207}

# Create the database directory
mkdir -p "${DB_DIR}"

# Download the compressed package and save it to the database directory
TAR_FILE="${DB_DIR}/gtdbtk_data.tar.gz"

echo "Starting to download the GTDB-Tk database (release ${RELEASE_VERSION})..."
wget "https://data.ace.uq.edu.au/public/gtdb/data/releases/release${RELEASE_VERSION}/${RELEASE_VERSION}.0/auxillary_files/gtdbtk_data.tar.gz" \
     -O "${TAR_FILE}"

# Unzip the database to the specified directory
tar -xzvf "${TAR_FILE}" -C "${DB_DIR}"

echo "The GTDB-Tk database has been downloaded and unzipped to: ${DB_DIR}"
echo "You can specify the database directory when using GTDB-Tk in the following way (or set an environment variable):"
echo "  gtdbtk classify_wf --data_path ${DB_DIR} --genome_dir your_genomes"
echo "Or set it in the environment variable:"
echo "  export GTDBTK_DATA_PATH=${DB_DIR}"

