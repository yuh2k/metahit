#!/bin/bash

# Define default download directory
DB_DIR="$../databases"
mkdir -p "$DB_DIR"

# Function to download a file with progress
download_file() {
    local url=$1
    local dest=$2
    echo "Downloading: $url"
    wget -c "$url" -O "$dest"
}

# Function to extract tar.gz files
extract_tar_gz() {
    local file=$1
    local dest_dir=$2
    echo "Extracting: $file"
    tar -xvzf "$file" -C "$dest_dir"
}

# 1. Download CheckM Database
CHECKM_URL="https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
CHECKM_DB="$DB_DIR/checkm_data.tar.gz"

download_file "$CHECKM_URL" "$CHECKM_DB"
extract_tar_gz "$CHECKM_DB" "$DB_DIR"
checkm data setRoot "$DB_DIR"

# 2. Download CheckM2 Database
dat
checkm2 database --download --path "$DB_DIR"


# 3. Download GTDB-Tk Database
GTDBTK_URL="https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
GTDBTK_DB="$DB_DIR/gtdbtk_data.tar.gz"

download_file "$GTDBTK_URL" "$GTDBTK_DB"
extract_tar_gz "$GTDBTK_DB" "$DB_DIR"

conda activate gtdbtk-2.4.0
conda env config vars set GTDBTK_DATA_PATH="$DB_DIR"

# Cleanup downloaded files
rm -f "$CHECKM_DB" "$CHECKM2_DB" "$GTDBTK_DB"

echo "All databases downloaded and extracted to: $DB_DIR"