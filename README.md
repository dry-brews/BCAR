# BCAR: Barcode Collapse by Aligning Reads

BCAR is a C/C++ tool designed to merge raw sequencing reads based on shared barcode identity.
BCAR is designed to work on very large sequencing datasets that may not fit in memory.
It will correct sequencing errors and produce a single consensus read for each barcode in your library.
The Q-scores attached to the consensus read are based on the observed frequency of the consensus base among the raw reads at each position.
BCAR2 works in two phases to sort the fastq(s) and then merge the sorted reads by barcode.
The first script (fastq_sorter) uses a memory-efficient disk sorting algorithm to sort reads on the basis of their barcodes.
The second script (bc_merger) works with a pair of sorted fastq files to generate consensus reads for each barcode. 

## Features
- Uses a fast implementation of the Needleman-Wunsch alignment algorithm to handle indels (insertions and deletions) between reads.
- fastq_sorter accepts any number of input .fastq or .fastq.gz files (or combinations thereof). The files will be combined into a **single** barcode map.
- Your barcodes and their counts will appear in the headers of your consensus reads
- Bayesian estimation of post-merged quality scores
- Handling of single-end or paired-end reads (paired-end inputs will yield paired-end outputs)

## Requirements
- The barcode must occur at a **fixed location** in the **forward** read.
- If the barcode is on the reverse read, simply swap the order you pass the files to the scripts. 
- If the barcode position is not fixed but is adjacent to a constant sequence, you can preprocess your data using tools like [CutAdapt](https://cutadapt.readthedocs.io/) to trim your reads, ensuring that the barcodes are at a fixed location.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/dry-brews/BCAR.git

```

### 2. Build and Install
BCAR allows for simple compilation using a Makefile. Installation has been tested on Linux (Ubuntu) and MacOS. For unsupported systems, you may need to compile independently using gcc/g++ or clang.
```bash
cd BCAR
chmod +x install.sh
./install.sh
```

To test your installation:
```bash
./fastq_sorter
./bc_merger
```

### 3. Usage
```bash
Usage: ./fastq_sorter [options]
Options:
  --in1 input1.fastq,input2.fastq,...
  --in2 pairs1.fastq,pairs2.fastq,... (if using paired-end reads)
  --bc-start int   Barcode start position (Zero-indexed, default: 0)
  --bc-len int     Barcode length (default: 18)
  --out1 file      Output file for sorted read 1
  --out2 file      Output file for sorted read 2
  --temp dir       Temporary directory for storing chunk files (recommend .)

Usage: ./seq_merge [options]
Options:
  --in1 file       Input FASTQ file 1
  --in2 file       Input FASTQ file 2 (if using paired-end reads)
  --bc-start int   Barcode start position (Zero-indexed, default: 0)
  --bc-len int     Barcode length (default: 18)
  --out1 file      Output file for consensus read 1
  --out2 file      Output file for consensus read 2
  --threads int    Number of threads (default: 1)
```

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on the GitHub repository.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

For further information, contact [Bryan Andrews](mailto:andrewsb@uchicago.edu).


