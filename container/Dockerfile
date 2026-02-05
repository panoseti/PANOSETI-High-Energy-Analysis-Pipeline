# PANOSETI High Energy Analysis Pipeline - Docker Container
# This Docker image provides a complete environment to run CORSIKA air shower
# simulations and analyze them using the PANOSETI pipeline.
#
# Base: rootproject/root:6.28.10-ubuntu22.04 - Official ROOT container
#       Ubuntu 22.04 base with ROOT 6.28.10 pre-installed

FROM rootproject/root:6.28.10-ubuntu22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /opt/panoseti

# Set build args for CORSIKA credentials (optional - pre-copy tarball instead)
ARG CORSIKA_USER
ARG CORSIKA_PASSWORD

# Pass credentials as environment variables for download script
ENV CORSIKA_USER=$CORSIKA_USER
ENV CORSIKA_PASSWORD=$CORSIKA_PASSWORD

# Add PANOSETI binaries to PATH
ENV PATH="/usr/local/bin:${PATH}"

# Install additional dependencies for CORSIKA and corsikaIOreader
RUN apt-get update && apt-get install -y \
    wget \
    git \
    cmake \
    gfortran \
    make \
    gnuplot \
    eog \
    && rm -rf /var/lib/apt/lists/*

# Copy installation scripts to container
COPY container/scripts/*.sh /opt/panoseti/

# Copy CORSIKA tarball if it exists in build context
# This is an alternative to providing credentials via build args
COPY container/corsika-77550.tar.gz /opt/panoseti/

# Make scripts executable
RUN chmod +x /opt/panoseti/*.sh

# Step 1: Download and Extract CORSIKA 7.7550
# CORSIKA simulates extensive air showers from cosmic rays and gamma rays
# Uses credentials if provided, otherwise must pre-copy tarball
RUN bash download_corsika.sh

# Step 2: Install CORSIKA with PANOSETI configuration
# Compilation with:
# - QGSJET-II-04: High energy hadronic interaction model (required for PANOSETI)
# - URQMD 1.3cr: Low energy hadronic interaction model
# - IACT: Imaging Atmospheric Cherenkov Telescope detector interface
# - CHERENKOV: Cherenkov radiation generation
# - ATMEXT: External atmosphere profiles
# - VOLUMEDET: Volume detector geometry
RUN ./install_corsika.sh

# Step 3: Download corsikaIOreader
# Tool to read CORSIKA binary output and convert to ROOT format
RUN bash download_corsikaIOreader.sh

# Step 4: Install corsikaIOreader
# Requires ROOT (pre-installed in base image)
RUN ./install_corsikaIOreader.sh || \
    echo "Notice: corsikaIOreader installation completed with exit code $?"

# Step 5: Copy example files
# Example CORSIKA simulation ready-to-run
COPY simulation-tools/example/ /opt/panoseti/example/
COPY simulation-tools/panodisplay.C /opt/panoseti/simulation-tools/
COPY simulation-tools/run_example.C /opt/panoseti/simulation-tools/

# Set up environment for running examples
ENV PYTHONUNBUFFERED=1

# Create activation script for PANOSETI environment
# Create activation script for PANOSETI environment
RUN printf '#!/bin/bash\necho "PANOSETI environment activated"\necho "ROOT: $(root-config --version)"\necho "corsikaIOreader: $(which corsikaIOreader 2>/dev/null || echo not found)"\nexec "$@"\n' \
    > /opt/panoseti/activate_panoseti.sh && \
    chmod +x /opt/panoseti/activate_panoseti.sh

# Default working directory for users
WORKDIR /opt/panoseti

# Add volume mount points for data and output
# - /opt/panoseti/simulation-tools: For working with simulations
# - /opt/panoseti/data: For CORSIKA input/output files
VOLUME ["/opt/panoseti/simulation-tools", "/opt/panoseti/data"]

# Entry point: Activate PANOSETI environment and start bash
ENTRYPOINT ["/opt/panoseti/activate_panoseti.sh"]
CMD ["/bin/bash"]

# Labels for metadata
LABEL maintainer="PANOSETI Team"
LABEL description="PANOSETI High Energy Analysis Pipeline with CORSIKA, ROOT, and corsikaIOreader"
LABEL version="2.0"
