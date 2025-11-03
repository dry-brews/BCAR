# BCAR: Barcode Collapse by Aligning Reads

BCAR is a C/C++ tool designed to merge raw sequencing reads based on shared barcode identity.
BCAR is designed to work on very large sequencing datasets that may not fit in memory.
It will correct sequencing errors and produce a single consensus read for each barcode in your library.
The Q-scores of the consensus are calculated from the quality scores of the input sequences.
BCAR2 works in two phases to sort the fastq(s) and then merge the sorted reads by barcode.
The first script (bcar_sort) uses a memory-efficient disk sorting algorithm to sort reads on the basis of their barcodes.
The second script (bcar_merge) works with a sorted fastq file to generate consensus reads for each barcode. 

## Features
- Uses a fast banded implementation of the Needleman-Wunsch algorithm generate an optimal global alignment between each read associated with a barcode and the consensus read for that barcode.
- Your barcodes and their counts will appear in the headers of your consensus reads
- Bayesian estimation of post-merged quality scores

## Requirements
- The barcode must occur at a **fixed location** in the read (the forward read for paired-end inputs). You are recommended to put the barcode at the **front** of the read, or else indels occuring before the barcode may cause the barcode to be mis-read. If you have a 3` barcode, you will probably want to reverse-complement your reads with a tool like [seqtk](https://github.com/lh3/seqtk), running BCAR on the reverse-complemented reads. 
- If the barcode position is not fixed but is adjacent to a constant sequence, you can preprocess your data using tools like [CutAdapt](https://cutadapt.readthedocs.io/) to trim your reads, ensuring that the barcodes are at a fixed location.

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
#### Example
```bash
~/path/to/BCAR/bcar_sort --in my_trimmed_reads.fastq --bc-start 0 --bc-len 25 --temp . --out my_sorted_reads.fastq

~/path/to/BCAR/bcar_merge --in my_sorted_reads.fastq --bc-start 0 --bc-len 25 --threads 12 --out my_merged_reads.fastq
```

```bash
Usage: ./bcar_sort [options]
Options:
  --in [file(s)]        input1.fastq,input2.fastq,...
  --bc-start [int]      Barcode start position (Zero-indexed, default: 0)
  --bc-len [int]        Barcode length (default: 18)
  --out [file]          Output file for sorted read 1
  --temp [dir]          Temporary directory for storing chunk files (recommend .)

Usage: ./bcar_merge [options]
Options:
  --in [file]           Input FASTQ file 1
  --out [file]          Output file for consensus reads
  --bc-start [int]      Barcode start position (Zero-indexed, default: 0)
  --bc-len [int]        Barcode length (default: 18)
  --max-len [int]       Maximum read length (default: 131072)
  --gap [float]         Gap penalty used during alignment (default:-1.0)
  --threads [int]       Number of threads (default: 1)
  --no-alignment        If added, skip alignment and make barcode with unaligned consensus 
```

### 4. Special Use Cases
#### Very long reads (>100kb)
If you have very long reads (>100kb), pay special attention to the --max-len flag in bcar_merge. If you have even a single read that is longer than --max-len, it may cause unexpected behavior.
You should either pre-trim your reads to a maximum length or determine the longest read length and set --max-len slightly higher than that value.
For reads of this length, you are likely to run into memory constraints.
Memory is shared across all active threads, so if you set fewer threads, more memory will be available to each.

#### Paired end data
Paired-end reads are not supported by the default bcar_merge build, because we generally see better behavior when running all pre-processing steps (e.g., FLASH) prior to running BCAR. 
However, if your experimental design necessitates running BCAR on paired data, you can install a version that does handle paired-end reads by running the following commands while in the BCAR directory:
```bash
conda activate bcar-env
${CONDA_PREFIX}/bin/gcc -O3 -w -I${CONDA_PREFIX}/include -L${CONDA_PREFIX}/lib -lm -pthread -o bcar_merge_paired ./src/bc_merger_v11.c
```
Use the flags --in-pairs and --out-pairs for both the sorting and merging scripts.

#### Very many reads (>1 billion)
In our experience, BCAR is sufficiently fast to handle hundreds of millions of reads within a few hours.
If you have more reads than that and many CPUs, you can split your reads into smaller files to be processed in parallel.
Suppose you sorted your dataset and determined it contains 10.5 million barcodes, and you have 10 CPUs to work with.
Run the following:
```bash
python ./src/split_sorted_reads_into_chunks.py --in1 my_sorted_reads.fastq --bc-start 0 --bc-len 18 --size 1050000
```
This will split your data into 10 files with 1,050,000 barcodes, each called my_sorted_reads_0001.fastq, etc., that you can process separately.
Afterwards, you can concatenate the BCAR output files.

#### Many missense mutations, few indels
BCAR defaults to a gap penalty of -1.0, which is appropriate when indel errors and missense errors are approximately equally abundant.
On simulated reads, BCAR is quite accurate with this gap score across a broad range of indel and missense errors.
However, if you have very frequent missense errors (>1%) and quite infrequent indel errors (<0.1%), then using a more severe gap penalty may improve accuracy.
A gap penalty of -3.0 is probably appropriate for many datasets. 

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on the GitHub repository.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

For further information, contact [Bryan Andrews](mailto:andrewsb@uchicago.edu).


