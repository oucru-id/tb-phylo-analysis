# TB FHIR Phylogeny Analysis Pipeline

This pipeline processes FHIR bundle JSON files containing Mycobacterium tuberculosis genomics data to generate SNP distance matrices, phylogenetic trees, transmission network visualizations, and Federated Nextstrain builds. Please refer to our full documentation: https://tb-pipeline-docs.readthedocs.io/en/latest/index.html

## Features

- **FHIR Genomics Data Gateway:** Directly fetch patient FHIR genomic bundles from FHIR servers.
- **Alternative local FHIR Genomics input:** Support for local processing of FHIR Bundles.
- **Phylogenetic Analysis:** Generates SNP distance matrices and phylogenetic trees (Rectangular, Circular, Unrooted).
- **Transmission Network:** Interactive Graph visualization and statistical plots (histogram, heatmap, violin plot).
- **Federated Analysis:** Generates inputs for Nextstrain Augur for phylogenetic federated analytics.
- **Clinical Metadata Integration:** Extracts metadata (location, lineage, patient ID) from FHIR resources.

## Usage

### Requirements

- [Nextflow](https://www.nextflow.io/)
- Python 3.8+
- [Augur](https://docs.nextstrain.org/projects/augur/en/stable/index.html)
- Python packages: `biopython`, `pandas`, `networkx`, `pyvis`, `matplotlib`, `seaborn`, `numpy`, `requests`

Install Python dependencies:
```bash
pip install biopython pandas networkx pyvis matplotlib seaborn numpy requests nextstrain-augur
```

### Run the Pipeline

```bash
nextflow run main.nf
```

### Docker

Requirements:
- Docker and Docker Compose v2

Build and run:
```bash
docker compose build
docker compose up
```

Behavior:
- Code runs inside the container at /opt/pipeline.
- Results are written to /workspace/results in the container, mapped to ./workspace/results on the host.
- Environment variables are loaded from .env automatically by docker-compose.

Override settings:
- To change tree builder or timetree behavior, set in .env:
```
AUGUR_TREE_METHOD=fasttree
AUGUR_ENABLE_TIMETREE=false
```
- Rebuild/run after changes:
```bash
docker compose build && docker compose up
```

### Authentication

This pipeline supports two ways to authenticate to a FHIR server:

- API Key header (X-API-Key)
- OAuth2 Client Credentials (Bearer token from access token URL)

API Key example:
```bash
nextflow run main.nf \
  --fhir_server_url https://your-fhir.example \
  --fhir_server_auth YOUR_API_KEY
```

OAuth2 Client Credentials example:
```bash
nextflow run main.nf \
  --fhir_server_url https://your-fhir.example \
  --fhir_oauth_token_url https://idp.example/oauth2/token \
  --fhir_client_id YOUR_CLIENT_ID \
  --fhir_client_secret YOUR_CLIENT_SECRET \
  --fhir_scope openid
```

Auth precedence:
- If both OAuth2 credentials and API key are provided, OAuth2 (Bearer) is used.
- If only API key is provided, X-API-Key is used.
- If neither is provided, requests are unauthenticated.

#### Using .env

You can provide credentials via environment variables. Create a `.env` file:
```
FHIR_SERVER_URL=https://your-fhir.example
FHIR_API_KEY=YOUR_API_KEY
FHIR_OAUTH_TOKEN_URL=https://idp.example/oauth2/token
FHIR_CLIENT_ID=YOUR_CLIENT_ID
FHIR_CLIENT_SECRET=YOUR_CLIENT_SECRET
FHIR_SCOPE=openid
FHIR_FETCH_SINCE=2025-12-25
```

Load the variables in your shell and run:
```bash
set -a; source .env; set +a
nextflow run main.nf
```

### Input

- FHIR Bundle Genomics JSON files in `data/JSON/` (local) 
- FHIR Bundle Genomics from FHIR server (params.fhir_server_url and params.fhir_server_auth)
- Reference genome FASTA in `data/H37Rv.fasta`
