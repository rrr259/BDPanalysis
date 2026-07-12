# BDPanalysis

The project combines a Python data-processing pipeline with an R  downstream analysis to create a catalogue of bidirectional promoters. The data for this research was taken from the Gene Expression Omnibus (GEO) database via series GSE87821, provided by King and Klose 2017.

## Repository contents

- `DATA_PROCESSING.py` - interactive Python pipeline for RNA-seq data preparation, alignment, promoter-region processing, read counting, and optional IGV bedGraph generation.
- `Analysis.Rmd` - R Markdown analysis for normalization, ranking bidirectional promoter regions, DESeq2 differential expression analysis, visualization, catalogue generation, and functional enrichment analysis.

## Workflow overview

The analysis is split into two stages.

### 1. Data preparation in Python

`DATA_PROCESSING.py` is designed to:

1. Download RNA-seq FASTQ files from SRA.
2. Run FastQC quality control.
3. Optionally trim adapters with Cutadapt.
4. Optionally summarize QC results with MultiQC.
5. Download reference genome FASTA and GTF files.
6. Build a STAR genome index.
7. Align paired-end reads with STAR.
8. Prepare promoter-region BED files for bidirectional promoter analysis.
9. Convert BED regions to GTF format.
10. Count reads over custom promoter regions using featureCounts.
11. Optionally create positive- and negative-strand bedGraph files for IGV visualization.

### 2. Downstream analysis in R

`Analysis.Rmd` uses the `counts.txt` output from featureCounts to:

1. Load and prepare count data.
2. Normalize counts using RPM.
3. Select divergent promoter-region pairs ending in `_2` and `_3`.
4. Rank bidirectional promoter activity across sample groups.
5. Run DESeq2 differential expression analysis for OCT4 and BRG1 conditions.
6. Generate MA plots, PCA plots, and Venn diagrams.
7. Build a catalogue of bidirectional promoters.
8. Annotate promoter regions with Ensembl gene information.
9. Run GO-based functional enrichment analysis using ranked gene sets.

## Requirements

Python packages:

- `pandas`

R packages:

- `biomaRt`
- `dplyr`
- `data.table`
- `tidyverse`
- `DESeq2`
- `pheatmap`
- `viridis`
- `GO.db`
- `VennDiagram`
- `grid`
- `EnhancedVolcano`
- `ggrepel`

External command-line tools:

- NCBI EDirect tools: `esearch`, `efetch`
- SRA Toolkit: `fastq-dump`
- FastQC
- Cutadapt
- MultiQC
- STAR
- SAMtools
- BEDTools, including `flankBed`
- UCSC utilities: `bedToGenePred`, `genePredToGtf`
- Subread `featureCounts`
- `wget`
- `gunzip`

## Usage

Clone the repository:

```bash
git clone https://github.com/rrr259/BDPanalysis.git
cd BDPanalysis
```

Run the Python preparation workflow:

```bash
python3 DATA_PROCESSING.py
```

Then run the R Markdown analysis using `Analysis.Rmd` after `counts.txt` and the required annotation files have been generated.

## Expected input files

Depending on the stage being run, the workflow may require:

- an SRA accession or study identifier
- reference genome FASTA file
- reference annotation GTF file
- promoter-region BED file
- STAR alignment output BAM files
- featureCounts `counts.txt`
- `SraRunTable.txt`
- annotation files such as `overlapping.bed` and `ann_ENS.txt`

## Main output files

The workflow can generate:

- downloaded FASTQ files
- FastQC and MultiQC reports
- STAR genome index files
- sorted BAM alignment files
- modified promoter BED files
- custom GTF files for featureCounts
- `counts.txt`
- positive- and negative-strand bedGraph files
- RPM-normalized count tables
- promoter ranking CSV files
- DESeq2 result tables
- bidirectional promoter catalogue files
- GO functional enrichment result files

## Notes

- The scripts are interactive and ask the user for accessions, file paths, adapter sequences, and output file names.
- Some paths in `Analysis.Rmd` are hard-coded and may need to be changed before running on another system.
- The sample grouping in `Analysis.Rmd` is currently written for OCT4 and BRG1 treated/untreated groups and should be adjusted for other datasets.
- Make sure all external command-line tools are installed and available in your system `PATH`.
