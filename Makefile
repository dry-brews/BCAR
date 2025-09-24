# Makefile for BCAR2 project

# Check if we're in the conda environment
ifndef CONDA_PREFIX
$(error Please run 'conda activate bcar2-env' before using make)
endif

# Compiler settings
CXX = $(CONDA_PREFIX)/bin/g++
CC = $(CONDA_PREFIX)/bin/gcc
CXXFLAGS = -O3 -w -I$(CONDA_PREFIX)/include
CFLAGS = -O3 --std=gnu99 -w -I$(CONDA_PREFIX)/include
LDFLAGS = -L$(CONDA_PREFIX)/lib

# Source directory
SRC_DIR = src

# Output executables
FASTQ_SORTER = fastq_sorter
SEQ_MERGER = seq_merge

# Default target
all: $(FASTQ_SORTER) $(SEQ_MERGER) $(SEQ_MERGER_DEBUG)

# Rule for fastq_sorter
$(FASTQ_SORTER): $(SRC_DIR)/fastq_sorter_v02.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@ -lz

# Rule for seq_merge (bc_merger)
$(SEQ_MERGER): $(SRC_DIR)/bc_merger_v06.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@ -lm -pthread

# Clean target
clean:
	rm -f $(FASTQ_SORTER) $(SEQ_MERGER)

# Phony targets
.PHONY: all clean
