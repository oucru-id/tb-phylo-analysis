import requests
import json
import argparse
import os
import sys
from urllib.parse import urljoin, urlparse

def get_headers(api_key):
    return {
        'Accept': 'application/fhir+json',
        'X-API-Key': api_key
    }

def handle_pagination_url(base_url, next_url):
    if not next_url:
        return None
        
    if not next_url.startswith('http'):
        return urljoin(base_url, next_url)
    
    next_parsed = urlparse(next_url)
    base_parsed = urlparse(base_url)
    
    if next_parsed.netloc != base_parsed.netloc:
        return f"{base_url}/Observation?{next_parsed.query}"
    
    return next_url

def fetch_data(base_url, api_key, since=None):
    headers = get_headers(api_key)
    base_url = base_url.rstrip('/')
    
    print(f"Connecting to FHIR Server: {base_url}")
    
    print(f"Searching for Patients with Genetic Data (Code 69548-6)")
    
    patients = set()
    url = f"{base_url}/Observation?code=69548-6&_count=1000"
    
    if since:
        print(f"  > Filtering for data updated after: {since}")
        url += f"&_lastUpdated=gt{since}"
    
    while url:
        try:
            resp = requests.get(url, headers=headers)
            resp.raise_for_status()
            data = resp.json()
            
            if 'entry' in data:
                for entry in data['entry']:
                    res = entry.get('resource', {})
                    subj = res.get('subject', {}).get('reference', '')
                    if subj.startswith('Patient/'):
                        pat_id = subj.split('/')[-1]
                        patients.add(pat_id)
            
            next_link = next((l['url'] for l in data.get('link', []) if l['relation'] == 'next'), None)
            url = handle_pagination_url(base_url, next_link)

        except Exception as e:
            print(f"Error fetching search list: {e}")
            break
            
    print(f"Found {len(patients)} patients with variant data.")
    
    if len(patients) == 0:
        return

    for pid in patients:
        print(f"Fetching full bundle for Patient {pid}")
        patient_resources = []
        
        try:
            p_resp = requests.get(f"{base_url}/Patient/{pid}", headers=headers)
            if p_resp.ok:
                patient_resources.append(p_resp.json())
        except Exception as e:
            print(f"  Error fetching Patient/{pid}: {e}")


        obs_url = f"{base_url}/Observation?patient={pid}&_count=1000"
        
        while obs_url:
            try:
                o_resp = requests.get(obs_url, headers=headers)
                if not o_resp.ok:
                    print(f"  Failed to fetch observations: {o_resp.status_code}")
                    break
                    
                o_data = o_resp.json()
                entries = o_data.get('entry', [])
                if entries:
                    print(f"  - Downloaded {len(entries)} observations...")
                    for e in entries:
                        patient_resources.append(e['resource'])
                
                next_ob_link = next((l['url'] for l in o_data.get('link', []) if l['relation'] == 'next'), None)
                obs_url = handle_pagination_url(base_url, next_ob_link)
                
            except Exception as e:
                print(f"  Error fetching variants for {pid}: {e}")
                break
            
        try:
            d_resp = requests.get(f"{base_url}/DiagnosticReport?patient={pid}", headers=headers)
            if d_resp.ok:
                d_data = d_resp.json()
                for e in d_data.get('entry', []):
                    patient_resources.append(e['resource'])
        except Exception as e:
            print(f"  Error fetching reports for {pid}: {e}")

        bundle = {
            "resourceType": "Bundle",
            "type": "transaction", 
            "entry": [{"resource": r} for r in patient_resources]
        }
        
        fname = f"{pid}.fhir.json"
        with open(fname, 'w') as f:
            json.dump(bundle, f, indent=2)
        print(f"  Saved {len(patient_resources)} resources to {fname}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--url', required=True, help="FHIR Server URL")
    parser.add_argument('--auth', required=True, help="API Key")
    parser.add_argument('--since', help="Filter data updated after YYYY-MM-DD")
    args = parser.parse_args()
    
    fetch_data(args.url, args.auth, args.since)