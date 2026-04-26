# Makefile for BCAR project

# Check if we're in the conda environment
ifndef CONDA_PREFIX
$(error Please run 'conda activate bcar-env' before using make)
endif

# Compiler settings
CFLAGS = -O3 -std=gnu99 -Wextra -I$(CONDA_PREFIX)/include
LDFLAGS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib

# Source directory
SRC_DIR = src

# Output executables
BCAR       = bcar
#BCAR_DEBUG = bcar_debug

# Default target
#all: $(BCAR) $(BCAR_DEBUG)
all: $(BCAR)

# Combined pipeline: sort + merge + consensus
$(BCAR): $(SRC_DIR)/bcar_pipeline.c $(SRC_DIR)/sort_module.c $(SRC_DIR)/sort_module.h $(SRC_DIR)/seq_module.c $(SRC_DIR)/seq_module.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRC_DIR)/bcar_pipeline.c $(SRC_DIR)/sort_module.c $(SRC_DIR)/seq_module.c -o $@ -lz -lm -pthread

# Clustering debug tool (no consensus, evaluates against ground truth barcodes in headers)
#$(BCAR_DEBUG): $(SRC_DIR)/bcar_debug_cluster.c $(SRC_DIR)/sort_module.c $(SRC_DIR)/sort_module.h
#	$(CC) $(CFLAGS) $(LDFLAGS) $(SRC_DIR)/bcar_debug_cluster.c $(SRC_DIR)/sort_module.c -o $@ -lz -lm

# Clean target
clean:
	rm -f $(BCAR) $(BCAR_DEBUG)

# Phony targets
.PHONY: all clean
