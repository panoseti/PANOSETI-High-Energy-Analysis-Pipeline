#!/bin/bash

# Download corsikaIOreader from GitHub for Docker environment
# This tool reads CORSIKA binary output files and converts them to ROOT format

REPO_URL="https://github.com/nkorzoun/corsikaIOreader.git"
CLONE_DIR="corsikaIOreader"

if [ -d "${CLONE_DIR}" ]; then
    echo "corsikaIOreader repository already cloned at ${CLONE_DIR}"
    echo "Skipping download..."
    exit 0
else
    echo "Cloning corsikaIOreader repository from GitHub..."
    git clone "${REPO_URL}" "${CLONE_DIR}"

    if [ $? -eq 0 ]; then
        echo "Repository cloned: ${CLONE_DIR}"
        exit 0
    else
        echo "Error: Failed to clone corsikaIOreader repository"
        exit 1
    fi
fi