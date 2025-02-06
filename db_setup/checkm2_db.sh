#!/usr/bin/env bash
#
# Function: Download and configure the CheckM2 database
# Usage: bash download_checkm2_db.sh [DB_DIR]
#        If DB_DIR is not specified, it will be downloaded to the "database" folder in the current directory by default
#

DB_DIR=${1:-database}

mkdir -p "${DB_DIR}"

# Use the built-in setup_database command of CheckM2 to automatically download and configure the database
checkm2 setup_database "${DB_DIR}"

echo "The CheckM2 database has been downloaded to: ${DB_DIR}"
echo "You can specify the database directory when using CheckM2 in the following way (can be written in a script or environment variable):"
echo "  checkm2 --db-path ${DB_DIR} [other commands]"