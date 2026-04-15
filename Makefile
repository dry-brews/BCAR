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

# Output executable
BCAR = bcar

# Default target
all: $(BCAR)

# Combined pipeline: sort + merge
$(BCAR): $(SRC_DIR)/bcar_pipeline.c $(SRC_DIR)/sort_module.c $(SRC_DIR)/sort_module.h $(SRC_DIR)/seq_module.c $(SRC_DIR)/seq_module.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRC_DIR)/bcar_pipeline.c $(SRC_DIR)/sort_module.c $(SRC_DIR)/seq_module.c -o $@ -lz -lm -pthread

# Clean target
clean:
	rm -f $(BCAR)

# Phony targets
.PHONY: all clean
