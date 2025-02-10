# Makefile for BCAR2 project

# Compiler settings
CXX = g++
CC = gcc
CXXFLAGS = -O3 --std=c++17 -w
CFLAGS = -O3 --std=gnu99 -w

# Source directory
SRC_DIR = src

# Output executables
FASTQ_SORTER = fastq_sorter
BC_MERGER = bc_merge

# Default target
all: $(FASTQ_SORTER) $(BC_MERGER)

# Rule for fastq_sorter
$(FASTQ_SORTER): $(SRC_DIR)/fastq_sorter.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ -lz

# Rule for seq_merge (bc_merger)
$(BC_MERGER): $(SRC_DIR)/bc_merger.c
	$(CC) $(CFLAGS) $< -o $@ -lm -pthread

# Clean target
clean:
	rm -f $(FASTQ_SORTER) $(BC_MERGER)

# Phony targets
.PHONY: all clean
