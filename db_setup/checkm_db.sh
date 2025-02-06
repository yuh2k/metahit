#!/usr/bin/env bash
#
# Function: Download and configure the CheckM database
# Usage: bash download_checkm_db.sh [DB_DIR]
#        If DB_DIR is not specified, it will be downloaded to the "database" folder in the current directory by default
#

# If the script has the first parameter, use it as the database directory, otherwise default to "database"
DB_DIR=${1:-database}

# Create the database directory if it does not exist
mkdir -p "${DB_DIR}"

# Set the CheckM database environment variable
export CHECKM_DATA_PATH="${DB_DIR}"

# Correct usage (no need to pass the URL unless you want a custom one)
checkm data download "${DB_DIR}"

echo "The CheckM database has been downloaded to: ${DB_DIR}"
echo "Please ensure that CHECKM_DATA_PATH is set before using CheckM, for example:"
echo "  export CHECKM_DATA_PATH=${DB_DIR}"