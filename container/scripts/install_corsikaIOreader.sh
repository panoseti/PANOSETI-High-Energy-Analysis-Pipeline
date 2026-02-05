#!/bin/bash

# Install corsikaIOreader for PANOSETI Docker container
# This tool reads CORSIKA binary output files and converts to ROOT format

CLONE_DIR="corsikaIOreader"
INSTALL_BIN="/usr/local/bin"  # System-wide installation in container

# Check if ROOT is available
if ! command -v root-config &> /dev/null; then
    echo "Error: ROOT not found!"
    echo "ROOT must be installed before corsikaIOreader"
    exit 1
fi

echo "Found ROOT installation:"
root-config --version
echo ""

# Ensure ROOT is activated
source $(dirname $(dirname $(which root-config)))/bin/thisroot.sh

# Check if repository exists
if [ ! -d "${CLONE_DIR}" ]; then
    echo "Error: corsikaIOreader repository not found"
    echo "Run ./download_corsikaIOreader.sh first"
    exit 1
fi

echo ""
echo "Building corsikaIOreader..."
echo ""

cd "${CLONE_DIR}"

# Build corsikaIOreader with ROOT
make corsikaIOreader

if [ $? -eq 0 ]; then
    echo ""
    echo "corsikaIOreader built successfully!"
    echo ""

    # Install to system bin directory
    echo "Installing to ${INSTALL_BIN}..."
    mkdir -p "${INSTALL_BIN}"
    cp corsikaIOreader "${INSTALL_BIN}/"

    if [ $? -eq 0 ]; then
        echo ""
        echo "Installation complete!"
        echo ""
        echo "Binary installed to: ${INSTALL_BIN}/corsikaIOreader"
        echo "This location is in PATH and accessible from anywhere in the container"
        echo ""
        echo "Verify installation:"
        echo "  corsikaIOreader --help"
        echo ""
    else
        echo "Error: Failed to install corsikaIOreader"
        exit 1
    fi
    
    # Verify the binary is in PATH
    if command -v corsikaIOreader &> /dev/null; then
        echo "Installation verified successfully!"
    else
        echo "Warning: corsikaIOreader installed but not in PATH"
    fi
else
    echo "Error: corsikaIOreader build failed!"
    echo "Check that ROOT is properly installed and activated"
    exit 1
fi