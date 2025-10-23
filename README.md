# BCAR: Barcode Collapse by Aligning Reads

BCAR is a C/C++ tool designed to merge raw sequencing reads based on shared barcode identity.
BCAR is designed to work on very large sequencing datasets that may not fit in memory.
It will correct sequencing errors and produce a single consensus read for each barcode in your library.
The Q-scores of the consensus are calculated from the input scores of the input sequences.
BCAR2 works in two phases to sort the fastq(s) and then merge the sorted reads by barcode.
The first script (bcar_sort) uses a memory-efficient disk sorting algorithm to sort reads on the basis of their barcodes.
The second script (bcar_merge) works with a pair of sorted fastq files to generate consensus reads for each barcode. 

## Features
- Uses a fast banded implementation of the Needleman-Wunsch algorithm generate an optimal global alignment between each read associated with a barcode and the consensus read for that barcode.
- Your barcodes and their counts will appear in the headers of your consensus reads
- Bayesian estimation of post-merged quality scores
- Handles single-end or paired-end reads (paired-end inputs will yield paired-end outputs)

## Requirements
- The barcode must occur at a **fixed location** in the read (the forward read for paired-end inputs). You are recommended to put the barcode at the **front** of the read, or else indels occuring before the barcode may cause the barcode to be mis-read. If you have a 3` barcode, you will probably want to reverse-complement your reads with a tool like [seqtk](https://github.com/lh3/seqtk), running BCAR on the reverse-complemented reads. 
- If the barcode position is not fixed but is adjacent to a constant sequence, you can preprocess your data using tools like [CutAdapt](https://cutadapt.readthedocs.io/) to trim your reads, ensuring that the barcodes are at a fixed location.
- If you have very long reads (>100,000bp), pay special attention to the --max-len flag in bcar_merge. If you have even a single read that is longer than --max-len, it may cause unexpected behavior. You can adjust --max-len upwards, but you will eventually run into memory constraints with very long reads.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/dry-brews/BCAR.git
```

### 2. Build and Install
Use the provided installation script to compile from source. Installation may take several minutes. Installation has been tested on Linux (Ubuntu) and MacOS. 
```bash
cd BCAR
chmod +x install.sh
./install.sh
```

To test your installation:
```bash
chmod +x test_installation.sh
./test_installation.sh
```

### 3. Usage
```bash
Usage: ./bcar_sort [options]
Options:
  --in [file(s)]        input1.fastq,input2.fastq,...
  --in-pairs [file(s)]  pairs1.fastq,pairs2.fastq,... (if using paired-end reads)
  --bc-start [int]      Barcode start position (Zero-indexed, default: 0)
  --bc-len [int]        Barcode length (default: 18)
  --out [file]          Output file for sorted read 1
  --out-pairs [file]    Output file for sorted read 2
  --temp [dir]          Temporary directory for storing chunk files (recommend .)

Usage: ./bcar_merge [options]
Options:
  --in [file]           Input FASTQ file 1
  --in-pairs [file]     Input FASTQ file 2 (if using paired-end reads)
  --bc-start [int]      Barcode start position (Zero-indexed, default: 0)
  --bc-len [int]        Barcode length (default: 18)
  --max-len [int]       Maximum read length (default: 131072)
  --gap [float]         Gap penalty used during alignment (default:-3.0)
  --out [file]          Output file for consensus read 1
  --out-pairs [file]    Output file for consensus read 2
  --threads [int]       Number of threads (default: 1)
  --no-alignment        If added, skip alignment and make barcode with unaligned consensus 
```

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on the GitHub repository.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

For further information, contact [Bryan Andrews](mailto:andrewsb@uchicago.edu).


