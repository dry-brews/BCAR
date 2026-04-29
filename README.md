# BCAR: Barcode Collapse by Aligning Reads

BCAR is a C/C++ tool designed to merge raw sequencing reads based on shared barcode identity.
BCAR is designed to work on very large sequencing datasets that may not fit in memory.
It will correct sequencing errors and produce a single consensus read for each barcode in your library.
The Q-scores of the consensus are calculated from the quality scores of the input sequences.

On simulated reads, BCAR outperforms comparable tools like PacRAT and fgbio when indel errors are reasonably common.
When sequencing errors (both indels and substitutions) are very rare, BCAR has similar performance compared to other tools.

If you find this code helpful, please cite our pre-print: https://doi.org/10.64898/2026.03.27.714882

## Features
- Automatically clusters barcodes to correct for read errors in the barcodes themselves
- Uses a fast banded implementation of the Needleman-Wunsch algorithm to globally align each read associated with a barcode and the consensus read for that barcode.
- Bayesian estimation of post-merged quality scores.
- Barcodes and diagnostic information appear in the headers of the consensus reads for easy post-filtering.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/dry-brews/BCAR.git
```

### 2. Build and Install
Use the provided installation script to compile from source.
Installation may take several minutes. Installation has been tested on Linux (Ubuntu) and MacOS.
Windows users are recommended to use Windows Subsystem for Linux (WSL).

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

## Usage
### Example
```bash
~/path/to/BCAR/bcar_sort --in my_trimmed_reads.fastq --bc-start 0 --bc-len 25 --tmp-dir . --out my_sorted_reads.fastq
```

```
Usage: ./bcar [options]

Required:
  --in <file[,file,...]>           FASTQ input(s), comma-separated, gzip OK

Barcode extraction (choose mode):
  --bc-start <int>                 Fixed position start (0-based, default: 0)
  --bc-len <int>                   Fixed barcode length (default: 15)
  --context <pattern>              Context pattern with Ns for barcode positions
  --max-context-mismatches <int>   Allowed mismatches in flanking seqs (default: 1)

Clustering:
  --max-mismatches <int>           Hamming distance radius (default: 1)
  --max-bc-indels <int>            Max indels in clustering/extraction (0-2, default: 0)
  --max-n <int>                    Max N bases in barcode (default: 2)

Consensus:
  --threads <int>                  Number of worker threads (default: 1)
  --gap <float>                    Gap penalty for alignment (default: -1.0)
  --max-len <int>                  Maximum read length (default: 100,000)
  --no-alignment                   Skip alignment

Output:
  --out <file>                     Output file (default: bcar_out.fastq)
  --fail <file>                    Reads where barcode not found (if context used)

Resources:
  --chunk-mem <bytes>              Sort buffer size (default: 4GB)
  --max-open-files <int>           Max files during merge (default: 240)
  --tmp-dir <path>                 Temp directory (default: .)

```

### Specifying a barcode
BCAR can detect barcodes through two methods: position-delimited or context.
When the context flag is not present, BCAR will default to position-delimited identification, and the barcode is the portion of your read staring at position --bc-start and extending for --bc-len bases.
If a context flag is used, it overrides --bc-start and --bc-len.
If your workflow allows you to use either position-delimited or context detection, we recommend to use context.

#### Context flags
Context can be provided in a grep-like format like this:
```
--context GAGTATGNNNNNTACGATC
--context ^NNNNNTACGATC
--context GAGTATGNNNNN$
```
The first call finds a 5-base barcode between two constant sequences.
The second and third calls find 5-base barcodes anchored at the start or end of the read, respectively, flanked by a single constant sequence.

### Inputs and outputs
In many cases, raw reads can be put into bcar in fastq format without any pre-processing.
However, if you have paired reads, you will likely want to merge your reads before passing them to BCAR, using a tool like [FLASH](https://github.com/ebiggers/flash).
Additionally, if your barcoded reads have sequence that is expected to vary within a single barcode group, such as UMIs, phasing regions, indices, or any other variable sequence, you want to first trim these sequences off using a tool like [CutAdapt](https://cutadapt.readthedocs.io/en/stable/)
Your input reads are likely to look something like this:
```
@AV240202:241115-32-PE300-RR-SW:2413661235:1:21304:2797:3178 2:N:0:GAGATTCC+GGCTCTGA
CAGGACCAACAATATGGATT...
+
NNMNNNNFNNNNNLNNNNNN...
@AV240202:241115-32-PE300-RR-SW:2413661235:1:11501:1469:2622 2:N:0:TCCGGAGA+ATAGAGGC
CGCGACGCACCGAAAAGGTT...
+
=(F?7LHH7JK.CB/GJNDM...
@AV240202:241115-32-PE300-RR-SW:2413661235:1:12104:1202:2303 2:N:0:ATTACTCG+TATAGCCT
GCTCGTGGTTCGTTAGGGTT...
+
LNMBNHMKM<N7MJMMLK>N...

```
The headers contain information about the sequencing run and the individual read.
BCAR will output another fastq that looks like this:
```
@bc=ATACGAGATCGGGTCGTG;bcid=464;count=62;min_q=40;mean_q=40;minor_frac_mean=0.004;minor_frac_max=0.049
ATACGAGATCGGGTCGTGTTATCAGTC...
+
IIIIIIIIIIIIIIIIIIIIIIIIIII...
@bc=AGGGCGTACTGACTGTCT;bcid=90;count=1;min_q=7;mean_q=42;minor_frac_mean=0.000;minor_frac_max=0.000
AGGGCGTACTGACTGTCTTTATCAGTC...
+
@8F(=*ILND(JLFFDLNNKNNNNNNN...
@bc=TTGGGCCCAGTGGCAGGA;bcid=2;count=95;min_q=40;mean_q=40;minor_frac_mean=0.009;minor_frac_max=0.326
TTGGGCCCAGTGGCAGGATTATCAGTC...
+
IIIIIIIIIIIIIIIIIIIIIIIIIII...
```
The headers here contain information about the barcode.
bcid=464 is a unique number assigned to each barcode group.
If barcode clustering is used, there may be multiple raw barcode sequences assigned to this number.
bc=TTGGGCCCAGTGGCAGGA is the barcode inferred from the conensus sequence using the same logic that was used to detect barcodes in the raw reads.
count=62 means that 62 raw reads were used to generate this consensus sequence.
min_q and mean_q refer to the quality scores assigned to each position the consensus read.
minor_frac_mean and minor_frac_max refer to the evidence for the second-most-common base at each position in the consensus read.
Take caution that while the headers always start with "@", quality score lines can also begin with "@".

Notice there are a few possible outcomes.
When you have many reads that all agree with each other (e.g., bcid=464), you will have consistently maxed-out quality scores of I=40 and minor_frac will be low.
When you have only a single read mapped to a given barcode (bcid=90), BCAR will return exactly that read with its original quality scores, and minor_frac_mean and minor_frac_max will always be zero. Note that quality scores >40 are possible in these reads, even though BCAR caps consensus quality scores at 40.

Sometimes, you may see consensus read with a high count and good quality scores, but minor_frac_max is high.
This may indicate that the barcode is associated with multiple variants in your library.
However, minor_frac_max can be surprisingly high even on simulated reads where no barcode misassociation is expected.
You should take high minor_frac_max values as an indication to look closer at your the possibility of barcode misassociation rather than a guarantee of such.

You will likely want to filter your consensus reads before using them to characterize your library.
How you do so is up to you, and the best way to do so depends on your experimental and sequencing strategies.
We have had success by filtering on a min_q >30, count >1.

### Threads
You should assign as many threads as you have available on your machine.
Reading, writing, and sorting read are each handled by a single thread, but the aligning and merging of each barcode group is a multithreaded process with speed approximately proportional to the number of threads available.

### 4. Special Use Cases
#### Very long reads (>100kb)
If you have very long reads (>100kb), pay special attention to the --max-len flag in bcar_merge. If you have even a single read that is longer than --max-len, it may cause unexpected behavior.
You should either pre-trim your reads to a maximum length or determine the longest read length and set --max-len slightly higher than that value.
For reads of this length, you are likely to run into memory constraints.
Memory is shared across all active threads, so if you set fewer threads, more memory will be available to each.

#### Many missense mutations, few indels
BCAR defaults to a gap penalty of -1.0, which is appropriate when indel errors and missense errors are approximately equally abundant.
On simulated reads, BCAR is quite accurate with this gap score across a broad range of indel and missense errors.
However, if you have very frequent missense errors (>1%) and quite infrequent indel errors (<0.1%), then using a more severe gap penalty may improve accuracy.
A gap penalty of -3.0 is probably appropriate for many such datasets.

If barcodes are expected to be extremely rare, BCAR can be run in --no-alignment mode, which is faster and may be more accurate when no indels are present.

## Contributing
Contributions are welcome! Feel free to open issues or submit pull requests on the GitHub repository.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

For further information, contact [Bryan Andrews](mailto:andrewsb@uchicago.edu).


