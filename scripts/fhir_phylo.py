import json
import argparse
import os
import sys
import csv
import re  
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def load_reference(ref_path):
    record = SeqIO.read(ref_path, "fasta")
    return str(record.seq), record.id

def parse_fhir(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    sample_id = os.path.basename(file_path).replace('.fhir.json', '').replace('.merged', '').replace('.json', '')
    variants = {} 
    
    metadata = {
        "sample_id": sample_id,
        "patient_id": "NA",
        "latitude": "NA",
        "longitude": "NA",
        "conclusion": "NA"
    }
    
    if 'entry' in data:
        for entry in data['entry']:
            res = entry.get('resource', {})
            r_type = res.get('resourceType')

            if r_type == 'Patient':
                metadata["patient_id"] = res.get('id', 'NA')
                for addr in res.get('address', []):
                    for ext in addr.get('extension', []):
                        if ext.get('url') == 'http://hl7.org/fhir/StructureDefinition/geolocation':
                            for sub_ext in ext.get('extension', []):
                                if sub_ext.get('url') == 'latitude':
                                    metadata["latitude"] = sub_ext.get('valueDecimal')
                                elif sub_ext.get('url') == 'longitude':
                                    metadata["longitude"] = sub_ext.get('valueDecimal')

            elif r_type == 'DiagnosticReport':
                conclusions = []
                for cc in res.get('conclusionCode', []):
                    if cc.get('text'): conclusions.append(cc.get('text'))
                if res.get('conclusion'): conclusions.append(res.get('conclusion'))
                if conclusions:
                    seen = set()
                    unique_conclusions = [x for x in conclusions if not (x in seen or seen.add(x))]
                    metadata["conclusion"] = "; ".join(unique_conclusions)

            elif r_type == 'Observation':
                code_coding = res.get('code', {}).get('coding', [])
                is_variant = any(c.get('code') == '69548-6' for c in code_coding)
                
                if is_variant:
                    pos = None
                    alt = None
                    
                    for comp in res.get('component', []):
                        for c in comp.get('code', {}).get('coding', []):
                            code = c.get('code')
                            if code == '81254-5': 
                                if 'valueRange' in comp: pos = comp['valueRange'].get('low', {}).get('value')
                                elif 'valueInteger' in comp: pos = comp.get('valueInteger')

                    if alt is None or pos is None:
                        hgvs_candidates = []
                        vcc = res.get('valueCodeableConcept', {}).get('coding', [])
                        for c in vcc:
                            if 'hgvs' in c.get('system', '') or ':' in c.get('code', ''): hgvs_candidates.append(c.get('code'))
                        for comp in res.get('component', []):
                            vcc = comp.get('valueCodeableConcept', {}).get('coding', [])
                            for c in vcc:
                                if 'hgvs' in c.get('system', '') or ':' in c.get('code', ''): hgvs_candidates.append(c.get('code'))
                        
                        for hgvs_str in hgvs_candidates:
                            if not hgvs_str: continue
                            match = re.search(r'g\.(\d+)[ACGTN]+>([ACGTN]+)', hgvs_str)
                            if match:
                                if pos is None: pos = int(match.group(1))
                                if alt is None: alt = match.group(2)
                                break 

                    if pos is not None and alt:
                        variants[pos] = alt
                        
    return sample_id, variants, metadata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True, help="List of FHIR JSON files")
    parser.add_argument('--reference', required=True, help="Reference genome FASTA")
    parser.add_argument('--anchors', nargs='*', default=[], help="List of Anchor FHIR JSON files")
    args = parser.parse_args()

    ref_seq, ref_name = load_reference(args.reference)
    ref_seq = ref_seq.upper()
    
    all_variant_positions = set()
    sample_variants = {}
    all_metadata = []

    ref_id = "H37Rv" 
    sample_variants[ref_id] = {}
    all_metadata.append({
        "sample_id": ref_id,
        "patient_id": "Reference",
        "latitude": "NA",
        "longitude": "NA",
        "conclusion": "Reference Genome"
    })

    for f in args.anchors:
        sid, vars, meta = parse_fhir(f)
        meta['patient_id'] = "Reference"
        if meta['conclusion'] == "NA":
            meta['conclusion'] = "Anchor"
        else:
            meta['conclusion'] = meta['conclusion'] + " (Anchor)"
            
        sample_variants[sid] = vars
        all_metadata.append(meta)
        all_variant_positions.update(vars.keys())

    for f in args.inputs:
        sid, vars, meta = parse_fhir(f)
        sample_variants[sid] = vars
        all_metadata.append(meta)
        all_variant_positions.update(vars.keys())

    with open("metadata.tsv", "w", newline='') as f:
        fieldnames = ["sample_id", "patient_id", "latitude", "longitude", "conclusion"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(all_metadata)

    samples = list(sample_variants.keys())
    
    sorted_positions = sorted(list(all_variant_positions))
    snp_seqs = []
    
    for sid in samples:
        vars = sample_variants[sid]
        s_chars = []
        for p in sorted_positions:
            if p <= len(ref_seq):
                base = vars.get(p, ref_seq[p-1])
                s_chars.append(base[0])
            else:
                s_chars.append('N')
        snp_seqs.append("".join(s_chars))

    n = len(samples)
    matrix_data = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            val1 = snp_seqs[i]
            val2 = snp_seqs[j]
            dist = 0
            for k in range(len(val1)):
                if val1[k] != 'N' and val2[k] != 'N' and val1[k] != val2[k]:
                    dist += 1
            matrix_data[i][j] = dist
            matrix_data[j][i] = dist
            
    with open("distance_matrix.tsv", "w") as f:
        f.write("snp-dists\t" + "\t".join(samples) + "\n")
        for i, row in enumerate(matrix_data):
            f.write(samples[i] + "\t" + "\t".join(map(str, row)) + "\n")

    if len(snp_seqs) > 0 and len(snp_seqs[0]) > 0:
        aln_records = [SeqRecord(Seq(s), id=i, description="") for s, i in zip(snp_seqs, samples)]
        aln = MultipleSeqAlignment(aln_records)
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        Phylo.write(tree, "phylo_tree.nwk", "newick")
    else:
        with open("phylo_tree.nwk", "w") as f: f.write("();")
    
    full_genome_records = []
    ref_list = list(ref_seq) 
    
    for sid in samples:
        vars = sample_variants[sid]
        sample_seq_list = ref_list[:] 
        
        for pos, alt in vars.items():
            idx = pos - 1
            if 0 <= idx < len(sample_seq_list):
                sample_seq_list[idx] = alt[0] 
        
        full_seq = "".join(sample_seq_list)
        full_genome_records.append(SeqRecord(Seq(full_seq), id=sid, description=""))

    with open("consensus.fasta", "w") as output_handle:
        SeqIO.write(full_genome_records, output_handle, "fasta")

if __name__ == "__main__":
    main()
