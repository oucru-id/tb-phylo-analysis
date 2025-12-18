# TB FHIR Phylogeny Analysis Pipeline

This pipeline processes FHIR bundle JSON files containing Mycobacterium tuberculosis genomics data to generate SNP distance matrices, phylogenetic trees, and transmission network visualizations. Please refer to our full documentation: https://tb-pipeline-docs.readthedocs.io/en/latest/index.html


## Features

- **Direct FHIR Genomics JSON input**
- **Phylogenetic tree:** Rectangular, circular, and unrooted. 
- **Transmission network visualization:** Graph and statistical plots (histogram, heatmap, violin plot).
- **Clinical metadata integration:** Extracted from Bundle Genomics FHIR.

## Usage

### Requirements

- [Nextflow](https://www.nextflow.io/)
- Python 3.8+
- Python packages: `biopython`, `pandas`, `networkx`, `pyvis`, `matplotlib`, `seaborn`

Install Python dependencies:
```bash
pip install biopython pandas networkx pyvis matplotlib seaborn
```

### Run the Pipeline

```bash
nextflow run main.nf
```

### Input

- FHIR Bundle Genomics JSON files in `data/JSON/`
- Reference genome FASTA in `data/H37Rv.fasta`
