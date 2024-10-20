#!/bin/bash

# setup.sh
# This script sets up the necessary dependencies for the MetaHit pipeline.
# It downloads BWA (from GitHub, builds it), Samtools, FastQC, and BBTools
# into the "external" directory within the repository.

# Exit immediately if a command exits with a non-zero status
set -e

# Function to print informational messages
function echo_info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}

# Function to print error messages
function echo_error() {
    echo -e "\033[1;31m[ERROR]\033[0m $1" >&2
}

# Define the external and bin directory paths
EXTERNAL_DIR="$(pwd)/external"
BIN_DIR="${EXTERNAL_DIR}/bin"

# Create the external and bin directories if they don't exist
if [ ! -d "$EXTERNAL_DIR" ]; then
    echo_info "Creating 'external' directory."
    mkdir -p "$EXTERNAL_DIR"
else
    echo_info "'external' directory already exists."
fi

if [ ! -d "$BIN_DIR" ]; then
    echo_info "Creating 'bin' directory inside 'external'."
    mkdir -p "$BIN_DIR"
else
    echo_info "'bin' directory already exists inside 'external'."
fi

# Update package lists
echo_info "Updating package lists..."
sudo apt update

# Install system dependencies
echo_info "Installing system dependencies: build-essential, wget, unzip, openjdk-11-jdk, perl, git..."
sudo apt install -y build-essential wget unzip openjdk-11-jdk perl git

# Function to download and build BWA from GitHub
function install_bwa() {
    echo_info "Cloning BWA repository from GitHub..."
    cd "$EXTERNAL_DIR"
    
    # Remove existing 'bwa' directory if it exists to prevent conflicts
    if [ -d "bwa" ]; then
        echo_info "Removing existing 'bwa' directory to ensure a clean build."
        rm -rf bwa
    fi

    git clone https://github.com/lh3/bwa.git
    cd bwa
    echo_info "Building BWA..."
    make

    echo_info "Copying BWA binary to 'external/bin'."
    cp bwa "$BIN_DIR/"

    cd ..
    # Clean up
    rm -rf bwa
    echo_info "BWA built and installed successfully."
}

# Function to download, extract, and install Samtools
function install_samtools() {
    SAMTOOLS_VERSION="1.20"
    SAMTOOLS_TARBALL="samtools-${SAMTOOLS_VERSION}.tar.bz2"
    SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

    if [ ! -f "${EXTERNAL_DIR}/${SAMTOOLS_TARBALL}" ]; then
        echo_info "Downloading Samtools ${SAMTOOLS_VERSION}..."
        wget -O "${EXTERNAL_DIR}/${SAMTOOLS_TARBALL}" "${SAMTOOLS_URL}"
    else
        echo_info "Samtools tarball already exists, skipping download."
    fi

    if [ ! -f "${EXTERNAL_DIR}/samtools/bin/samtools" ]; then
        echo_info "Extracting Samtools..."
        tar -xjf "${EXTERNAL_DIR}/${SAMTOOLS_TARBALL}" -C "$EXTERNAL_DIR"
        cd "${EXTERNAL_DIR}/samtools-${SAMTOOLS_VERSION}"
        echo_info "Configuring Samtools..."
        ./configure --prefix="${EXTERNAL_DIR}/samtools"
        echo_info "Building Samtools..."
        make
        echo_info "Installing Samtools..."
        make install
        cd ../..

        echo_info "Creating a symbolic link to samtools binary in 'external/bin'..."
        ln -sf "${EXTERNAL_DIR}/samtools/bin/samtools" "${BIN_DIR}/samtools"

        echo_info "Cleaning up Samtools source directory..."
        rm -rf "${EXTERNAL_DIR}/samtools-${SAMTOOLS_VERSION}" "${EXTERNAL_DIR}/${SAMTOOLS_TARBALL}"

        echo_info "Samtools installed successfully."
    else
        echo_info "Samtools binary already exists in 'external/bin', skipping compilation."
    fi
}

# Function to download and install FastQC
function install_fastqc() {
    FASTQC_VERSION="0.11.9"
    FASTQC_ZIP="fastqc_v${FASTQC_VERSION}.zip"
    FASTQC_URL="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip"

    if [ ! -f "${EXTERNAL_DIR}/${FASTQC_ZIP}" ]; then
        echo_info "Downloading FastQC ${FASTQC_VERSION}..."
        wget -O "${EXTERNAL_DIR}/${FASTQC_ZIP}" "${FASTQC_URL}"
    else
        echo_info "FastQC zip file already exists, skipping download."
    fi

    if [ ! -f "${BIN_DIR}/fastqc" ]; then
        echo_info "Extracting FastQC..."
        unzip "${EXTERNAL_DIR}/${FASTQC_ZIP}" -d "${EXTERNAL_DIR}"
        mv "${EXTERNAL_DIR}/FastQC" "${EXTERNAL_DIR}/fastqc"
        chmod +x "${EXTERNAL_DIR}/fastqc/fastqc"

        echo_info "Creating a symbolic link to FastQC binary in 'external/bin'..."
        ln -sf "${EXTERNAL_DIR}/fastqc/fastqc" "${BIN_DIR}/fastqc"

        # Clean up
        rm "${EXTERNAL_DIR}/${FASTQC_ZIP}"
        echo_info "FastQC installed successfully."
    else
        echo_info "FastQC binary already exists in 'external/bin', skipping extraction."
    fi
}

# Function to download and install BBTools
function install_bbtools() {
    BBTOOLS_VERSION="39.10"  # Latest version as per user request
    BBTOOLS_TARBALL="Bbmap_${BBTOOLS_VERSION}.tar.gz"
    BBTOOLS_URL="https://sourceforge.net/projects/bbmap/files/latest/download"
    if [ ! -f "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" ]; then
        echo_info "Downloading BBTools ${BBTOOLS_VERSION}..."
        wget -O "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" "${BBTOOLS_URL}"
    else
        echo_info "BBTools tarball already exists, skipping download."
    fi

    if [ ! -f "${BIN_DIR}/bbmap.sh" ]; then
        echo_info "Extracting BBTools..."
        tar -xzf "${EXTERNAL_DIR}/${BBTOOLS_TARBALL}" -C "$EXTERNAL_DIR" || { echo_error "Failed to extract BBTools tarball."; exit 1; }

        echo_info "Creating symbolic links for BBTools binaries in 'external/bin'..."
        cd "${EXTERNAL_DIR}/bbmap"
        for bin in *.sh; do
            ln -sf "${EXTERNAL_DIR}/bbmap/${bin}" "${BIN_DIR}/${bin}"
        done

        cd ../..

        echo_info "BBTools installed successfully."
    else
        echo_info "BBTools binaries already exist in 'external/bin', skipping extraction."
    fi
}


# Install all dependencies
install_bwa
install_samtools
install_fastqc
install_bbtools  # BBTools installation

# Verify installations
echo_info "Verifying installations..."

if [ ! -f "${BIN_DIR}/bwa" ]; then
    echo_error "BWA binary not found in external/bin."
    exit 1
else
    echo_info "BWA installed successfully."
fi

if [ ! -f "${BIN_DIR}/samtools" ]; then
    echo_error "Samtools binary not found in external/bin."
    exit 1
else
    echo_info "Samtools installed successfully."
fi

if [ ! -f "${BIN_DIR}/fastqc" ]; then
    echo_error "FastQC binary not found in external/bin."
    exit 1
else
    echo_info "FastQC installed successfully."
fi

# Fix missing shared library for FastQC (libnsl.so.1)
echo_info "Checking for libnsl.so.1 library..."

if ldconfig -p | grep -q "libnsl.so.1"; then
    echo_info "libnsl.so.1 is already installed."
else
    echo_info "Installing libnsl-dev to provide libnsl.so.1..."
    sudo apt install -y libnsl-dev
fi

rm external/Bbmap_39.10.tar.gz

# Ensure all external binaries have execute permissions
echo_info "Ensuring all external binaries have execute permissions."
chmod +x "${BIN_DIR}/bwa"
chmod +x "${BIN_DIR}/samtools"
chmod +x "${BIN_DIR}/fastqc"
chmod +x "${BIN_DIR}"/*

# Optionally, add external/bin to PATH
echo_info "Dependencies have been successfully installed and configured."
echo_info "You can add the 'external/bin' directory to your PATH for easier access to these tools."
echo_info "Run the following command in your terminal, or add it to your ~/.bashrc or ~/.bash_profile:"
echo_info "export PATH=\"${BIN_DIR}:\$PATH\""
echo_info "Then, reload your shell configuration with:"
echo_info "source ~/.bashrc"

echo_info "Setup completed. Please ensure all scripts have been correctly updated and rerun your pipeline."
