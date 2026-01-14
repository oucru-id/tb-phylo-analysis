nextflow.enable.dsl = 2

process PREPARE_INPUTS {
    publishDir "${params.results_dir}/federated_export", mode: 'copy'

    input:
    path metadata
    path fasta

    output:
    path "nextstrain_metadata.tsv", emit: meta
    path "nextstrain_sequences.fasta", emit: seqs

    script:
    """
    cp ${fasta} nextstrain_sequences.fasta
    
    python3 -c "
import pandas as pd
import re

df = pd.read_csv('${metadata}', sep='\\t')
df.rename(columns={'sample_id': 'strain'}, inplace=True)
df['date'] = '2023-01-01'

def extract_main_lineage(val):
    match = re.search(r'lineage(\\d+)', str(val), re.IGNORECASE)
    if match:
        return f'Lineage {match.group(1)}'
    return 'Unknown'

df['Lineage'] = df['conclusion'].apply(extract_main_lineage)

def extract_sub_lineage(val):
    match = re.search(r'Lineage lineage([\\d\\.]+)', str(val))
    return match.group(1) if match else 'Unknown'

df['SubLineage'] = df['conclusion'].apply(extract_sub_lineage)

df.to_csv('nextstrain_metadata.tsv', sep='\\t', index=False)
    "
    """
}

process RUN_AUGUR {
    publishDir "${params.results_dir}/nextstrain_build", mode: 'copy'
    
    cpus 6 

    input:
    path sequences
    path metadata
    path reference

    output:
    path "tb_analysis.json"

    script:
    """
    cat ${reference} ${sequences} > aligned.fasta

    augur tree \
      --alignment aligned.fasta \
      --nthreads ${task.cpus} \
      --method ${params.augur_tree_method} \
      --output tree.nwk
    
    EXTRA_FLAGS=""
    if [ "${params.augur_enable_timetree}" = "true" ]; then
      EXTRA_FLAGS="--timetree --coalescent opt --date-confidence --date-inference marginal --clock-rate 1e-7"
    else
      HAS_VALID_DATES=\$(python3 - << 'PY'
import pandas as pd
try:
    df = pd.read_csv('nextstrain_metadata.tsv', sep='\\t')
    if 'date' in df.columns:
        vals = [str(x).strip() for x in df['date'].dropna() if str(x).strip() and str(x).strip().upper() != 'NA']
        uniq = set(vals)
        print("yes" if len(uniq) >= 2 else "no")
    else:
        print("no")
except Exception:
    print("no")
PY
      )
      if [ "\$HAS_VALID_DATES" = "yes" ]; then
        EXTRA_FLAGS="--timetree --coalescent opt --date-confidence --date-inference marginal --clock-rate 1e-7"
      fi
    fi

    augur refine \
      --tree tree.nwk \
      --alignment aligned.fasta \
      --metadata ${metadata} \
      --output-tree refined.nwk \
      --output-node-data branch_lengths.json \
      \$EXTRA_FLAGS

    augur traits \
      --tree refined.nwk \
      --metadata ${metadata} \
      --columns Lineage \
      --confidence \
      --output-node-data traits.json

    # Create Auspice Config
    echo '{
      "title": "TB Federated Analysis",
      "colorings": [
        {
          "key": "Lineage",
          "title": "Main Lineage",
          "type": "categorical"
        },
        {
          "key": "SubLineage",
          "title": "Sub Lineage",
          "type": "categorical"
        },
        {
          "key": "strain",
          "title": "Strain Name",
          "type": "categorical"
        }
      ],
      "display_defaults": {
        "color_by": "Lineage",
        "branch_label": "Lineage"
      }
    }' > auspice_config.json

    echo "NC_000962.3" > exclude.txt
    
    augur filter \
      --metadata ${metadata} \
      --exclude exclude.txt \
      --output-metadata filtered_metadata.tsv
    
    augur export v2 \
      --tree refined.nwk \
      --metadata filtered_metadata.tsv \
      --node-data branch_lengths.json traits.json \
      --auspice-config auspice_config.json \
      --output tb_analysis.json
    """
}

workflow FEDERATED_NEXTSTRAIN {
    take:
    metadata
    fasta
    reference

    main:
    PREPARE_INPUTS(metadata, fasta)
    RUN_AUGUR(PREPARE_INPUTS.out.seqs, PREPARE_INPUTS.out.meta, reference)
}
