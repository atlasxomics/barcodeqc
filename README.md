# barcodeqc
`barcodeqc` is a lightweight command-line tool for the rapid evaluation of AtlasXomics epigenomic DBiT-seq experiments. The tool takes fastq files from read2 an input and analyzes the distribution of barcode sequences to screen for the top preliminary failure modes (low linker conservation, mismatched barcodes, missing/overabundent barcodes). 

The tool requires the reads with the barcoding schema described in [Zhang et al. 2023](https://www.nature.com/articles/s41586-023-05795-1) for Illumina short-read sequencing.

- read1: genomic sequence
- read2: linker1 | barcodeA | linker2 | barcodeB | genomic sequence

<div>
    <img src="./static/barcode_scheme.png" alt="scheme" width="400"/>
</div>

## Installation
We have tested this tool in Mac OS and Linux enviroments. We do not support running on Windows platforms.  This tool requires [python3.10](https://www.python.org/downloads/) or greater. We rely on GitHub for distributing this repository. Please ensure [git](https://github.com/git-guides/install-git) is installed locally before continuing. 

- [seqtk](https://github.com/lh3/seqtk)
    - Must be available on PATH for this tool to run.  We have tested with seqtk v1.4-r122.

## Steps 
By default, `barcodeqc` works with a random subsample of the original data to enable local processing on most machines.
1. **subsample**: reads from the read2 file are randomly downsampled to a value specified by the parameter SAMPLE_READS (default )

## Outputs
The tool outputs a results directory in the directory where the tool was run, named with the 'SAMPLE_NAME' supplied to the tool.  