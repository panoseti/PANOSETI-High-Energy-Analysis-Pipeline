#!/bin/bash

# Download CORSIKA 7.7550 for Docker container
# This script attempts multiple download methods:
# 1. Uses CORSIKA_USER and CORSIKA_PASSWORD if provided (recommended)
# 2. Falls back to direct download if credentials are embedded in image
# 3. Skips if CORSIKA tarball is already present

USER="${CORSIKA_USER:-}"
PASSWORD="${CORSIKA_PASSWORD:-}"
CORSIKA_TAR="corsika-77550.tar.gz"

# Check if CORSIKA is already downloaded
if [ -f "${CORSIKA_TAR}" ] && [ -d "corsika-77550" ]; then
    echo "CORSIKA 7.7550 already downloaded and extracted"
    exit 0
fi

if [ -f "${CORSIKA_TAR}" ]; then
    echo "CORSIKA tarball found: ${CORSIKA_TAR}"
    echo "Proceeding with pre-downloaded tarball..."
else
    echo "Downloading CORSIKA 7.7550..."

    # Try to download with credentials if provided
    if [ -n "$USER" ] && [ -n "$PASSWORD" ]; then
        echo "Using provided credentials to download from official CORSIKA server"
        wget --user="$USER" --password="$PASSWORD" \
            https://web.iap.kit.edu/corsika/download/corsika-v770/${CORSIKA_TAR}

        if [ $? -eq 0 ]; then
            echo "Download completed successfully"
        else
            echo "Warning: Failed to download with credentials"
            echo "CORSIKA must be downloaded manually or with valid credentials"
            echo "Get credentials from https://www.iap.kit.edu/corsika/"
            echo ""
            echo "Alternative: Pre-download the tarball and place it in the build context"
            echo "  cd container"
            echo "  wget --user=YOUR_USER --password=YOUR_PASS \\"
            echo "    https://web.iap.kit.edu/corsika/download/corsika-v770/corsika-77550.tar.gz"
            echo "  docker build -t panoseti:latest ."
            exit 1
        fi
    else
        echo "No CORSIKA credentials provided (CORSIKA_USER and CORSIKA_PASSWORD)"
        echo ""
        echo "To continue, you must:"
        echo "  1. Register at https://www.iap.kit.edu/corsika/"
        echo "  2. Either:"
        echo "     a) Set environment variables and rebuild:"
        echo "        export CORSIKA_USER=your_username"
        echo "        export CORSIKA_PASSWORD=your_password"
        echo "        make build"
        echo ""
        echo "     b) OR pre-download the tarball:"
        echo "        cd container"
        echo "        wget --user=YOUR_USER --password=YOUR_PASS \\"
            echo "          https://web.iap.kit.edu/corsika/download/corsika-v770/corsika-77550.tar.gz"
        echo "        make build"
        exit 1
    fi
fi

# Extract CORSIKA
if [ -f "${CORSIKA_TAR}" ] && [ ! -d "corsika-77550" ]; then
    echo "Extracting CORSIKA..."
    tar -xzf "${CORSIKA_TAR}"
    if [ $? -eq 0 ]; then
        echo "CORSIKA extracted successfully"
    else
        echo "Error: Failed to extract CORSIKA"
        exit 1
    fi
fi
