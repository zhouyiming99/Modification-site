# Nanopore RNA Processing and Modification Analysis Pipeline

This repository contains a comprehensive suite of shell scripts for processing Oxford Nanopore Technologies (ONT) long-read sequencing data, specifically focused on RNA basecalling, alignment, and the detection of RNA modifications (m6A, Pseudouridine) and PolyA tail lengths.

## Scripts Overview

### 1. Core Processing
*   **`basecalling.sh`**: Performs basecalling using Guppy, merges FASTQ files, and filters out failed FAST5 reads.
*   **`minimap2.sh`**: Conducts splice-aware alignment using `minimap2` and generates sorted, indexed BAM files via `samtools`.

### 2. RNA Modification Detection
*   **`nanopsu.sh`**: A complete pipeline for **Pseudouridine ($\Psi$)** detection. It includes basecalling, alignment, intron gap removal, feature extraction, and machine learning-based prediction.
*   **`m6anet.sh`**: Detects **m6A** modifications using a combination of `nanopolish` (for signal-to-sequence alignment) and `m6anet` inference.
*   **`modkit.sh`**: A batch processing script for **PolyA tail length** extraction and general **base modification** calling using ONT `modkit`. It automatically calculates statistics for multiple BAM files in a directory.

## Prerequisites

The following tools must be installed and available in your `$PATH`:
*   [Guppy / Dorado](https://community.nanoporetech.com/downloads) (Basecalling)
*   [Minimap2](https://github.com/lh3/minimap2) (Alignment)
*   [Samtools](http://www.htslib.org/) (BAM manipulation)
*   [Nanopolish](https://github.com/jts/nanopolish) (Signal alignment)
*   [m6anet](https://github.com/GoekeLab/m6anet) (m6A prediction)
*   [NanoPsu](https://github.com/your-source/nanopsu) (psU prediction)
*   [Modkit](https://github.com/nanoporetech/modkit) (Modification pileup)

## Usage

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/your-username/your-repo-name.git
    cd your-repo-name
    ```

2.  **Configuration**:
    Open the scripts and update the directory paths (e.g., `REF`, `INPUT_DIR`, `guppy`) to match your local environment and reference genome/transcriptome.

3.  **Permissions**:
    ```bash
    chmod +x *.sh
    ```

4.  **Execution**:
    Run the specific script based on your analysis needs:
    ```bash
    # Basic Processing
    ./basecalling.sh
    ./minimap2.sh

    # Modification Analysis
    ./nanopsu.sh
    ./m6anet.sh
    ./modkit.sh
    ```

## Outputs
-   **`prediction.csv`** (from `nanopsu.sh`): psU probabilities per site.
-   **`run/`** (from `m6anet.sh`): m6A probability scores.
-   **`*_modifications.bed`** (from `modkit.sh`): Modification frequencies in BED format.
-   **`*_polya_lengths.txt`** (from `modkit.sh`): Estimated PolyA tail length per read.

## Notes
- Scripts are configured by default for **Direct RNA (SQK-RNA002)** or **RNA R9.4.1** chemistry. Adjust the `--config` or `--kit` parameters if using R10.4.1 or cDNA kits.
- Ensure your BAM files contain modification tags (`MM`/`ML`) for `modkit.sh` to work correctly.
