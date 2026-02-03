#!/bin/bash

# Install CORSIKA 7.755 with PANOSETI requirements for Docker environment
# Compiling with: QGSJET-II-04, URQMD 1.3cr, IACT, CHERENKOV, ATMEXT, VOLUMEDET
# Base: Ubuntu 22.04 (from rootproject/root:6.28.10-ubuntu22.04)

if [ ! -d "corsika-77550" ]; then
    echo "Error: CORSIKA source directory not found"
    echo "Please run download_corsika.sh first"
    exit 1
fi

cd corsika-77550

echo "Cleaning previous installation..."
make distclean 2>/dev/null || true

echo "Extracting bernlohr package (required for IACT)..."
if [ -d "bernlohr" ]; then
    cd bernlohr
    if [ ! -d "bernlohr-1.67" ]; then
        tar -xzf bernlohr-1.67.tar.gz
    fi
    cd ..
fi

echo ""
echo "Installing CORSIKA 7.755 with PANOSETI configuration:"
echo "  - High energy model: QGSJET-II-04"
echo "  - Low energy model: URQMD 1.3cr"
echo "  - Detector geometry: Volume detector (VOLUMEDET)"
echo "  - Options:"
echo "    * IACT (Cherenkov telescope detector)"
echo "    * CHERENKOV (Cherenkov radiation)"
echo "    * ATMEXT (external atmosphere)"
echo ""

# Run coconut configuration with PANOSETI options
# This is non-interactive and follows standard PANOSETI configuration
# Options: QGSJET-II-04, URQMD, VOLUMEDET, IACT, CHERENKOV, ATMEXT
./coconut --no-cache << COCONUT_INPUTS
2
5
4
2
1b
z
y
f
COCONUT_INPUTS

if [ $? -eq 0 ]; then
    echo ""
    echo "CORSIKA 7.755 with PANOSETI options compiled successfully!"
    echo ""
    echo "Executable: /opt/panoseti/corsika-77550/run/corsika77550Linux_QGSII_urqmd"
    echo ""
    ls -lh run/corsika77550Linux_QGSII_urqmd
    echo ""
    echo "Configuration:"
    echo "  - QGSJET-II-04 (high energy hadronic model)"
    echo "  - URQMD 1.3cr (low energy hadronic model)"
    echo "  - VOLUMEDET (volume detector geometry)"
    echo "  - IACT (Cherenkov telescope detector - uses bernlohr)"
    echo "  - CHERENKOV (Cherenkov radiation generation)"
    echo "  - ATMEXT (external atmosphere functions)"

    # Create symlink for easier access
    if [ -f "run/corsika77550Linux_QGSII_urqmd" ]; then
        ln -sf /opt/panoseti/corsika-77550/run/corsika77550Linux_QGSII_urqmd /usr/local/bin/corsika
        echo "CORSIKA executable linked to /usr/local/bin/corsika"
    fi
else
    echo ""
    echo "Error: CORSIKA installation failed!"
    exit 1
fi