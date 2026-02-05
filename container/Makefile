# Makefile for building PANOSETI Docker container

.PHONY: build test run shell clean example example-sim example-vis ci-example-test help

# Default target
all: build

# Variables
IMAGE_NAME ?= panoseti
IMAGE_TAG ?= latest
CONTAINER_DIR ?= $(shell pwd)

# Optional build args for CORSIKA download
BUILD_ARGS :=
ifneq ($(strip $(CORSIKA_USER)),)
BUILD_ARGS += --build-arg CORSIKA_USER=$(CORSIKA_USER)
endif
ifneq ($(strip $(CORSIKA_PASSWORD)),)
BUILD_ARGS += --build-arg CORSIKA_PASSWORD=$(CORSIKA_PASSWORD)
endif

help:
	@echo "PANOSETI Docker Container Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  make build         - Build the Docker image"
	@echo "  make test          - Test the Docker container (runs all checks)"
	@echo "  make run           - Run the container with shell"
	@echo "  make shell         - Run the container interactively (shortcut for run)"
	@echo "  make clean         - Remove the Docker image"
	@echo "  make help          - Show this help message"
	@echo ""
	@echo "Simulation & Visualization:"
	@echo "  make example       - Run example CORSIKA simulation (requires X11 for visualization)"
	@echo "  make example-sim    - Run example CORSIKA simulation only (batch mode)"
	@echo "  make example-vis    - Visualize existing example.root (requires X11)"
	@echo "  make ci-example-test - Run PANOSETI example in CI (writes test_output/)"
	@echo ""
	@echo "Variables:"
	@echo "  IMAGE_NAME         - Docker image name (default: panoseti)"
	@echo "  IMAGE_TAG          - Docker image tag (default: latest)"
	@echo ""
	@echo "Examples:"
	@echo "  make build"
	@echo "  make test"
	@echo "  make IMAGE_NAME=mypanoseti IMAGE_TAG=v1.0 build"
	@echo "  make run"
	@echo "  make example       # Runs simulation and visualization (requires xhost +local:docker)"
	@echo ""

test:
	@echo "Running container tests..."
	@docker run --rm $(IMAGE_NAME):$(IMAGE_TAG) bash -c "root-config --version && corsikaIOreader --help"
	@docker run --rm $(IMAGE_NAME):$(IMAGE_TAG) bash -c "cd /opt/panoseti/corsika-77550/run && timeout 30 ./corsika77550Linux_QGSII_urqmd < all-inputs 2>&1 | tail -20"

build:
	@echo "============================================="
	@echo "Building PANOSETI Docker Container"
	@echo "============================================="
	@echo ""
	@if ! command -v docker > /dev/null 2>&1; then \
		echo "Error: Docker is not installed"; \
		echo "Please install Docker from https://docs.docker.com/get-docker/"; \
		exit 1; \
	fi
	@if [ ! -f $(CONTAINER_DIR)/corsika-77550.tar.gz ]; then \
		if [ -z "$(CORSIKA_USER)" ] || [ -z "$(CORSIKA_PASSWORD)" ]; then \
			echo "Error: corsika-77550.tar.gz not found and CORSIKA_USER/CORSIKA_PASSWORD not set"; \
			echo "Either set credentials or pre-download the tarball into $(CONTAINER_DIR)"; \
			exit 1; \
		fi; \
		echo "Downloading CORSIKA tarball..."; \
		wget --user="$(CORSIKA_USER)" --password="$(CORSIKA_PASSWORD)" \
			https://web.iap.kit.edu/corsika/download/corsika-v7750/corsika-77550.tar.gz -O $(CONTAINER_DIR)/corsika-77550.tar.gz; \
	fi
	@echo "Docker found: $$(docker --version)"
	@echo ""
	@echo "Building image: $(IMAGE_NAME):$(IMAGE_TAG)"
	@echo "Base image: rootproject/root:6.28.10-ubuntu22.04"
	@echo ""
	@echo "This may take 30-60 minutes (CORSIKA compilation)..."
	@echo ""
	@docker build $(BUILD_ARGS) -t $(IMAGE_NAME):$(IMAGE_TAG) $(CONTAINER_DIR)
	@if [ $$? -eq 0 ]; then \
		echo ""; \
		echo "============================================="; \
		echo "Build completed successfully!"; \
		echo "============================================="; \
		echo ""; \
		docker images $(IMAGE_NAME):$(IMAGE_TAG); \
		echo ""; \
		echo "Quick start commands:"; \
		echo ""; \
		echo "  # Run container with shell:"; \
		echo "  make run"; \
		echo ""; \
		echo "  # Run example visualization (with X11):"; \
		echo "  xhost +local:docker"; \
		echo "  docker run -it --rm -e DISPLAY=$$${DISPLAY} -v /tmp/.X11-unix:/tmp/.X11-unix $(IMAGE_NAME):$(IMAGE_TAG)"; \
		echo ""; \
		echo "  # Run with volume mounts for data:"; \
		echo "  docker run -it --rm -v $$(pwd)/data:/opt/panoseti/data $(IMAGE_NAME):$(IMAGE_TAG)"; \
	else \
		echo ""; \
		echo "Error: Docker build failed!"; \
		exit 1; \
	fi

run shell:
	@echo "Starting PANOSETI container..."
	docker run -it --rm $(IMAGE_NAME):$(IMAGE_TAG)

clean:
	@echo "Removing Docker image: $(IMAGE_NAME):$(IMAGE_TAG)"
	-docker rmi $(IMAGE_NAME):$(IMAGE_TAG)

# Run example CORSIKA simulation with visualization
# Requires X11 forwarding: xhost +local:docker
example:
	@echo "============================================="
	@echo "Running Example CORSIKA Simulation"
	@echo "============================================="
	@echo ""
	@echo "Running CORSIKA simulation..."
	@docker run --rm \
		-v $(CONTAINER_DIR)/example:/opt/panoseti/example \
		-w /opt/panoseti/corsika-77550/run \
		-e DISPLAY=$${DISPLAY} \
		-v /tmp/.X11-unix:/tmp/.X11-unix \
		$(IMAGE_NAME):$(IMAGE_TAG) \
		bash -c "\
		./corsika77550Linux_QGSII_urqmd < /opt/panoseti/example/example.inp 2>&1 | tee /opt/panoseti/example/corsika-example.log && \
		echo \'\' && \
		echo \'=============================================\' && \
		echo \'CORSIKA simulation complete!\' && \
		echo \'=============================================\' && \
		echo \'\' && \
		if [ -n \$${DISPLAY} ]; then \
			echo \'Starting visualization...\'; \
			mkdir -p /opt/panoseti/example/data && \
			cp /opt/panoseti/corsika-77550/run/atmabs.dat /opt/panoseti/example/data/ && \
			cp /opt/panoseti/corsika-77550/run/mirreff.dat /opt/panoseti/example/data/ && \
			cp /opt/panoseti/corsika-77550/run/quanteff.dat /opt/panoseti/example/data/ && \
			cd /opt/panoseti/example && \
			corsikaIOreader -cors DAT000001 -histo example-output.root -abs CORSIKA && \
			cp example-output.root /opt/panoseti/simulation-tools/ && \
			cd /opt/panoseti/simulation-tools && \
			root -l run_example.C\$$ && \
			mv example-output.root /opt/panoseti/example/example-output.root; \
		else \
			echo \'Simulation complete. To visualize, run with X11:\'; \
			echo \'  make example-vis\'; \
		fi"

# Run example CORSIKA simulation only (batch mode, no visualization)
example-sim:
	@echo "============================================="
	@echo "Running CORSIKA Simulation (Batch Mode)"
	@echo "============================================="
	@docker run --rm \
		-v $(CONTAINER_DIR)/example:/opt/panoseti/example \
		-w /opt/panoseti/corsika-77550/run \
		$(IMAGE_NAME):$(IMAGE_TAG) \
		bash -c "\
		./corsika77550Linux_QGSII_urqmd < /opt/panoseti/example/example.inp 2>&1 | tee /opt/panoseti/example/corsika-example.log && \
		echo \'\' && \
		echo \'=============================================\' && \
		echo \'Simulation complete! Output saved to:\' && \
		echo \'  corsika-example.log (simulation log)\' && \
		echo \'  DAT*, CER* (CORSIKA output files)\' && \
		echo \'=============================================\' && \
		echo \'\' && \
		echo \'To convert to ROOT format:\'; \
		echo \'  cd /opt/panoseti/example\'; \
		echo \'  mkdir -p data && cp /opt/panoseti/corsika-77550/run/atmabs.dat mirreff.dat quanteff.dat data/\'; \
		echo \'  corsikaIOreader -cors DAT000001 -histo example-output.root -abs CORSIKA\'"

# Visualize existing example.root (requires X11)
example-vis:
	@echo "============================================="
	@echo "Visualizing Example Simulation"
	@echo "============================================="
	@if [ -z "$${DISPLAY}" ]; then \
		echo "Error: X11 display not set"; \
		echo ""; \
		echo "Run with X11 forwarding enabled:"; \
		echo "  xhost +local:docker"; \
		echo "  export DISPLAY=$${DISPLAY}"; \
		echo "  make example-vis"; \
		exit 1; \
	fi
	@echo "Enabling X11 connections..."
	@xhost +local:docker
	@echo "Running visualization..."
	@docker run -it --rm \
		-v $(CONTAINER_DIR)/example:/opt/panoseti/example \
		-w /opt/panoseti/example \
		-e DISPLAY=${DISPLAY} \
		-v /tmp/.X11-unix:/tmp/.X11-unix \
		$(IMAGE_NAME):$(IMAGE_TAG) \
		bash -c "\
		if [ -f example.root ]; then \
			cp example.root /opt/panoseti/simulation-tools/ && \
			cd /opt/panoseti/simulation-tools && \
			echo \'Visualizing example.root...\'; \
			root -l run_example.C\$$ && \
			mv example.root /opt/panoseti/example/example.root.vis; \
		elif [ -f DAT000001 ]; then \
			echo \'Converting CORSIKA output first...\'; \
			mkdir -p data && \
			cp /opt/panoseti/corsika-77550/run/atmabs.dat mirreff.dat quanteff.dat data/ && \
			corsikaIOreader -cors DAT000001 -histo example-output.root -abs CORSIKA && \
			cp example-output.root /opt/panoseti/simulation-tools/ && \
			cd /opt/panoseti/simulation-tools && \
			echo \'Visualizing example-output.root...\'; \
			root -l run_example.C\$$ && \
			mv example-output.root /opt/panoseti/example/example-output.root.vis; \
		elif [ -f example.root.vis ]; then \
			cp example.root.vis /opt/panoseti/simulation-tools/example.root && \
			cd /opt/panoseti/simulation-tools && \
			echo \'Visualizing example.root...\'; \
			root -l run_example.C\$$ && \
			mv example.root /opt/panoseti/example/example.root.vis; \
		elif [ -f example-output.root.vis ]; then \
			cp example-output.root.vis /opt/panoseti/simulation-tools/example-output.root && \
			cd /opt/panoseti/simulation-tools && \
			echo \'Visualizing example-output.root...\'; \
			root -l run_example.C\$$ && \
			mv example-output.root /opt/panoseti/example/example-output.root.vis; \
		else \
			echo \'Error: No example data found\'; \
			echo \'Run make example-sim first to generate simulation data\'; \
			exit 1; \
		fi"

# CI-only: run PANOSETI example and collect logs/output
ci-example-test:
	@mkdir -p $(CONTAINER_DIR)/test_output
	@docker run --rm \
		-v $(CONTAINER_DIR)/test_output:/output \
		-w /opt/panoseti/corsika-77550/run \
		$(IMAGE_NAME):$(IMAGE_TAG) \
		bash -c "./corsika77550Linux_QGSII_urqmd < /opt/panoseti/example/example.inp 2>&1 > /output/corsika.log || true"
	@ls -lh $(CONTAINER_DIR)/test_output/DAT* $(CONTAINER_DIR)/test_output/CER* 2>/dev/null || echo "No CORSIKA output files created"
	@tail -50 $(CONTAINER_DIR)/test_output/corsika.log
