# BCAR: Barcode Collapse by Aligning Reads

BCAR is a Python/Cython tool designed to merge raw sequencing reads based on shared barcode identity. It will correct sequencing errors (mismatch by default and also indels with the --align flag) and produce a single consensus read for each barcode in your library. The Q-scores attached to the consensus read are based on the observed frequency of the consensus base among the raw reads at each position (post-alignment, if aligned).

## Features
- Uses a fast implementation of the Needleman-Wunsch alignment algorithm to handle indels (insertions and deletions) between reads.
- Accepts any number of input .fastq or .fastq.gz files (or combinations thereof). The files will be combined into a **single** barcode map.
- Your barcodes and their counts will appear in the headers of your consensus reads
- Empirical Bayes estimate of post-merged quality scores
- (New!) Works in single-read or paired end mode, depending on whether the --rev flag is included

## Limitations
- Does not consider the quality scores of the input reads

## Requirements
- The barcode must occur at a **fixed location** in the **forward** read.
- If the barcode is on the reverse read, simply pass the reverse read to the --fwd option and vice versa. 
- If the barcode position is not fixed but is adjacent to a constant sequence, you can preprocess your data using tools like [CutAdapt](https://cutadapt.readthedocs.io/) to trim your reads, ensuring that the barcodes are at a fixed location.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/dry-brews/BCAR.git
cd BCAR
```

### 2. Build and Install
BCAR relies on Cython for the slow steps of alignment and merging. Follow these steps to build and install:

Ensure that you have `Cython` installed in your Python environment:
```bash
conda install cython
```
or
```bash
pip install cython
```
Then run
```bash
python setup.py build
```
```bash
python setup.py install
```

To test your installation:
```bash
cd test
bcar --fwd Test_barcodes.fastq --rev Test_barcodes.fastq --align
```

## Troubleshooting
Ensure that you are using python 3.8 or later. Running setup.py may fail if you don't have the correct permissions. In this case, run "sudo python setup.py build". If you use a python -> python3 alias, you may need to run "sudo python3 setup.py build" instead.

## Usage
Once installed, BCAR can be invoked as a command-line script. The typical usage is:

```bash
bcar [-h] --fwd FWD_FASTQS [FWD_FASTQS ...] --rev REV_FASTQS [REV_FASTQS ...] [--out1 OUTPUT_FASTQ_FWD] [--out2 OUTPUT_FASTQ_REV] [--BC-start BC_START] [--BC-len BC_LEN] [--min-qscore MIN_QSCORE] [--min-count MIN_COUNT] [--align]
```

### Example Workflow
1. **Pre-trim Reads (if necessary)**:
   Use CutAdapt or another preprocessing tool to ensure that barcodes occur at a fixed location.
   ```bash
   cutadapt -g <constant_sequence> -o trimmed_reads.fastq input_reads.fastq
   ```

2. **Run BCAR**:
   ```bash
   bcar --fwd forward_reads_1.fastq.gz forward_reads_2.fastq.gz \
        --rev reverse_reads_1.fastq.gz reverse_reads_2.fastq.gz \
        --out1 consensus_forward.fastq --out2 consensus_reverse.fastq \
        --BC-start 10 --BC-len 16 --min-qscore 20 --min-count 5
   ```

3. **Merge consensus reads**:
   Use FLASH or something similar to merge your reads
   ```bash
   flash consensus_forward.fastq consensus_reverse.fastq -o my_merged_files
   ```

## How BCAR Works
BCAR reads both barcodes and raw reads as memory-efficient python integers, grouping all reads from each barcode into a separate list. It then (optionally) aligns the reads from each barcode and calls the most frequent base at each position, using an empirical Bayesian strategy to assign a Q-score to the basecall. It returns a .fastq file containing one read for each barcode, with the barcode and the count in the header of each read.

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on the GitHub repository.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

For further information, contact [Bryan Andrews](mailto:andrewsb@uchicago.edu).


